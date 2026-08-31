#!/usr/bin/env python3
"""calc.py — score an archived submission through all five validator stages, in full detail.

The miner archives every task it builds for under ``data/result/<task created_at>/``
(``submission.json`` and ``last_upload.json``, plus the broadcast ``contract.json`` and
``hbb_reference.json``). That contract
carries ``seed: 0`` — the backend stamps the real seed afterwards — and stages 3 and 4 are keyed on
it, so an archived folder cannot be scored from its own contents alone. This tool pairs the folder
with its task from the public ``/api/v3/tasks`` endpoint, runs the validator's own five stages
against the **real** seed, and writes each stage's working out back into the folder.

A stamped contract carries a comma-joined *list* of round seeds (``"122,321,431"``). Stage 12 is
seed-independent and runs once; stages 3-5 are re-derived per seed and the reported score is their
mean, exactly as ``benchmark_submission`` does it. Per-seed reports and the spread are printed
alongside the mean, and each seed gets its own ``stage{3,4,5}_detail.seed<N>.json``.

    python calc.py                                  # newest folder under data/result/
    python calc.py --folder data/result/2026-08-21T11:42:21
    python calc.py --task-id <uuid>                  # pin the contract, skip fingerprint matching
    python calc.py --submission data/submission.json --out-dir calc   # any loose file

The folder is matched to a task by, in order: the ``task_id`` in the ``last_upload.json`` the
miner archived beside the submission, the folder name against the task's ``created_at`` (how folders
written before that are paired), then the contract fingerprint (active mutations, weights, cell type,
rules). Whatever matches is then *verified* field by field against the archived contract, so a wrong
seed cannot be silently applied to the wrong submission.

Written into the scored folder:

    stage12_detail.json   per-row gate result and structural score, with the off-target k-mer hits
    stage3_detail.json    per-row RNG replay — every draw that decided cut, repair mode and indel
    stage4_detail.json    per-target, per-fold r2 / MAE / residual std / feature importances
    stage5_detail.json    every fidelity term, its share of the geometric mean and what it costs
    validation.json       the headline numbers, the task pairing and the artifact index
    validation.txt        the printed five-stage report
    validation/           the validator's own inputs and outputs, untouched by hand

The archived ``submission.json``, ``contract.json`` and ``hbb_reference.json`` are never written:
the pipeline runs over copies under ``validation/``, because ``stage12.truncate_submission``
rewrites the submission file in place when it caps rows or drops duplicate ids.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import io
import json
import logging
import math
import os
import random
import shutil
import sys
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import niome_subnet.utils.settings as settings  # noqa: E402  (must precede bittensor imports)

# genExp chdir()s to the repo root on import, which is what every relative settings path assumes.
import genExp as G  # noqa: E402

import numpy as np  # noqa: E402
from sklearn.ensemble import RandomForestRegressor  # noqa: E402
from sklearn.metrics import mean_absolute_error, r2_score  # noqa: E402
from sklearn.model_selection import KFold  # noqa: E402

from niome_subnet.genomics.validation import stage12, stage3, stage4, stage5  # noqa: E402


# Fields stage 1 indexes directly. A row missing any of them raises inside the validator rather
# than scoring zero, and run_validation swallows the exception and `continue`s — so the miner gets
# no score at all, not even a zero. Worth refusing to run rather than reproducing the crash.
HARD_FIELDS = ("guideRNA", "target_alignment_start", "mutation", "cas_system")

def parse_seeds(raw) -> list[int]:
    """Split a contract seed field into round seeds.

    The validator stopped scoring at one seed: ``contract["seed"]`` is a comma-joined list
    (``"122,321,431"``), stage 12 runs once because it is seed-independent, and stages 3-5 re-run
    per seed with every breakdown field and ``final_score`` averaged over them
    (``genomics/validation/__init__.py``). Mirrors ``_parse_seeds`` there, and tolerates a bare int
    so archived single-seed contracts still work.
    """
    return [int(part) for part in str(raw).split(",") if part.strip() != ""]


def mean_reports(reports: list[dict]) -> dict:
    """Average the per-seed stage reports the way benchmark_submission averages its breakdown.

    ``stage_report`` nests (``stage4.consistency_factor``, ``stage5.*``), so this recurses and
    averages the numeric leaves. Anything non-numeric — counters, histograms, the seed itself —
    is taken from the first seed so the shape still matches a single-seed report and the existing
    printer keeps working. ``None`` leaves stay ``None``: stage 4 writes them when it cannot fit.
    """
    if not reports:
        return {}
    if len(reports) == 1:
        return reports[0]      # untouched: averaging a single report would float-ify every count

    def merge(values):
        head = values[0]
        if isinstance(head, dict):
            return {k: merge([v.get(k) for v in values if isinstance(v, dict)])
                    for k in head}
        nums = [v for v in values if isinstance(v, (int, float)) and not isinstance(v, bool)]
        if nums and len(nums) == len(values):
            return sum(nums) / len(nums)
        return head

    return merge(reports)


def seed_spread(per_seed: list[tuple], path: tuple[str, ...]):
    """The per-seed values at a dotted path, for showing how much the seed alone moves a number."""
    out = []
    for item in per_seed:
        node = item[3]
        for key in path:
            node = (node or {}).get(key) if isinstance(node, dict) else None
        out.append(node if isinstance(node, (int, float)) else 0.0)
    return out


# Contract fields that must agree between the archived broadcast copy and the task the seed is
# taken from. Everything except the seed itself: the seed is precisely what the archive lacks.
FINGERPRINT = ("active_mutations", "mutation_weights", "mutation_regions", "cell_type", "rules",
               "version")

# Redirected into the staging directory. CHR11_PATH and KMER_CACHE_DIR are deliberately absent: the
# 130 MB FASTA and its k-mer index are read-only and expensive, so they stay shared with data/.
REDIRECTED = (
    "CONTRACT_PATH",
    "HBB_REFERENCE_PATH",
    "MINER_SUBMISSION_PATH",
    "LAST_UPLOAD_PATH",
    "VALID_EXPERIMENTS_PATH",
    "INVALID_EXPERIMENTS_PATH",
    "STAGE3_DATASET",
    "STAGE3_SUMMARY_PATH",
    "FINAL_REWARD_PATH",
    "DISTRIBUTION_FIDELITY_PATH",
)

STAGE4_FEATURES = ("gc", "distance", "gc_score", "dist_score", "consistency", "energy", "mh")
STAGE5_TERMS = ("mutation_coverage_entropy_ratio", "cas_system_coverage_entropy_ratio",
                "strand_coverage_entropy_ratio", "joint_coverage_entropy_ratio",
                "kmer_diversity_entropy_ratio", "distinct_guide_ratio")


def redirect_paths(stage_dir: Path) -> dict[str, str]:
    """Point every writable data path at ``stage_dir``.

    The stage modules bind their paths with ``from ... import NAME`` at import time, so patching
    ``settings`` alone would not reach them — each module's own attribute has to be rebound as well.
    """
    mapping: dict[str, str] = {}
    for name in REDIRECTED:
        original = getattr(settings, name, None)
        if original is None:
            continue
        target = str(stage_dir / Path(original).name)
        setattr(settings, name, target)
        for module in (stage12, stage3, stage4, stage5):
            if hasattr(module, name):
                setattr(module, name, target)
        mapping[name] = target
    return mapping


def read_json(path: str | Path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None


def write_json(path: Path, document) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(document, handle, indent=2)


def captured(report: dict, title: str) -> str:
    """Render a stage report to text and return it, so it can be both printed and archived."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        G.print_stage_report(report, title)
    return buffer.getvalue()


def spread(values) -> dict:
    """min / mean / median / max / std, JSON-safe. Empty input reports n=0 and nothing else."""
    values = [v for v in values if v is not None]
    if not values:
        return {"n": 0}
    array = np.asarray(values, dtype=float)
    return {
        "n": int(array.size),
        "min": float(array.min()),
        "mean": float(array.mean()),
        "median": float(np.median(array)),
        "max": float(array.max()),
        "std": float(array.std()),
    }


# --------------------------------------------------------------------------------------------
# Folder and task resolution
# --------------------------------------------------------------------------------------------

def latest_result_folder(root: Path) -> Path:
    """Newest archive folder holding a submission.

    Folder names are ISO timestamps, so lexicographic order is chronological order. The miner
    writes submission.json first and the contract afterwards, so a folder mid-write is skipped
    rather than scored half-formed.
    """
    if not root.is_dir():
        raise SystemExit(f"{root} does not exist — pass --folder or --submission")
    folders = sorted((p for p in root.iterdir() if p.is_dir()), key=lambda p: p.name, reverse=True)
    if not folders:
        raise SystemExit(f"no archive folders under {root}")
    for folder in folders:
        if (folder / "submission.json").exists():
            if folder is not folders[0]:
                print(f"  ! {folders[0].name} has no submission.json yet; using {folder.name}")
            return folder
    raise SystemExit(f"no folder under {root} contains a submission.json")


def fingerprint(contract: dict) -> dict:
    return {key: contract.get(key) for key in FINGERPRINT}


def resolve_task(task_id: str | None, folder_contract: dict | None, folder_name: str | None,
                 recorded_task_id: str | None, recorded_source: str) -> tuple[dict, str]:
    """The task whose contract this submission should be judged under, and how it was matched.

    ``/api/v3/tasks`` is the only place the real seed is visible: the contract broadcast to miners —
    and archived alongside the submission — carries ``seed: 0``.
    """
    items = G.fetch_all_tasks()      # newest first

    if task_id:
        for item in items:
            if item["id"] == task_id:
                return item, "--task-id"
        raise SystemExit(f"task {task_id} not present in the {len(items)} tasks returned")

    if recorded_task_id:
        for item in items:
            if item["id"] == recorded_task_id:
                return item, f"task id recorded in {recorded_source} ({recorded_task_id})"
        print(f"  ! {recorded_source} records task {recorded_task_id}, which the backend "
              "did not return")

    # The archive folder is named after the task's created_at, truncated to the second.
    if folder_name:
        matches = [item for item in items
                   if str(item.get("created_at", "")).startswith(folder_name)]
        if len(matches) == 1:
            return matches[0], f"folder name against created_at ({folder_name})"
        if len(matches) > 1:
            print(f"  ! {len(matches)} tasks share created_at {folder_name}; "
                  "falling back to the contract fingerprint")

    if folder_contract:
        wanted = fingerprint(folder_contract)
        matches = [item for item in items
                   if fingerprint(item["content"]["contract"]) == wanted]
        stamped = [item for item in matches if item["content"]["contract"].get("seed")]
        if len(stamped) == 1:
            return stamped[0], "contract fingerprint"
        if len(stamped) > 1:
            newest = stamped[0]
            print(f"  ! {len(stamped)} stamped tasks share this contract fingerprint; taking the "
                  f"newest ({newest['id']}, {newest.get('created_at')}). Pass --task-id to pin one.")
            return newest, "contract fingerprint (newest of several)"
        if matches:
            raise SystemExit("the only tasks matching this contract are unstamped (seed 0); "
                             "the backend has not assigned a seed yet")

    for item in items:
        if item["content"]["contract"].get("seed"):
            return item, "newest task with a non-zero seed (no folder match)"

    raise SystemExit("no task with a non-zero contract seed in the returned page")


def verify_pairing(task: dict, folder_contract: dict | None) -> dict:
    """Field-by-field check that the task's contract is the one the submission was built against."""
    if not folder_contract:
        return {"verified": False, "reason": "no archived contract to compare"}
    theirs = fingerprint(task["content"]["contract"])
    ours = fingerprint(folder_contract)
    differences = {key: {"archived": ours[key], "task": theirs[key]}
                   for key in FINGERPRINT if ours[key] != theirs[key]}
    return {"verified": not differences, "differences": differences,
            "archived_seed": folder_contract.get("seed"),
            "task_seed": task["content"]["contract"].get("seed"),
            "compared_against": None}   # filled in by the caller, which knows the file


# --------------------------------------------------------------------------------------------
# Preflight — everything that can be known about the submission before the stages run
# --------------------------------------------------------------------------------------------

def preflight(rows: list, contract: dict, reference: dict) -> dict:
    """Cross-check the submission against the contract it is about to be scored under.

    None of this changes the score — the stages decide that — but it separates "this submission is
    weak" from "this submission was built for a different task", which the stage output alone does
    not distinguish: a contract mismatch shows up only as every row rejected.
    """
    active = list(contract.get("active_mutations", []))
    cas_allowed = list(contract["rules"].get("cas_systems", ["Cas9", "Cas12a"]))
    max_experiments = contract["rules"].get("max_experiments")
    cell_type = contract.get("cell_type")

    malformed = [
        {"index": i, "experiment_id": row.get("experiment_id") if isinstance(row, dict) else None,
         "missing": [f for f in HARD_FIELDS if not isinstance(row, dict) or f not in row]}
        for i, row in enumerate(rows)
        if not isinstance(row, dict) or any(f not in row for f in HARD_FIELDS)
    ]

    ids = [row.get("experiment_id") for row in rows if isinstance(row, dict)]
    usable_ids = [i for i in ids if isinstance(i, str) and i.strip()]
    unusable_ids = len(ids) - len(usable_ids)
    duplicate_ids = len(usable_ids) - len(set(usable_ids))
    over_cap = max(0, len(rows) - max_experiments) if max_experiments else 0

    mutations = Counter(row.get("mutation") for row in rows if isinstance(row, dict))
    cas_systems = Counter(row.get("cas_system") for row in rows if isinstance(row, dict))
    cell_types = Counter(row.get("cell_type") for row in rows if isinstance(row, dict))
    lengths = Counter(len(row["guideRNA"]) for row in rows
                      if isinstance(row, dict) and isinstance(row.get("guideRNA"), str))

    return {
        "rows": len(rows),
        "malformed_rows": malformed,
        "unusable_experiment_ids": unusable_ids,
        "duplicate_experiment_ids": duplicate_ids,
        "max_experiments": max_experiments,
        "rows_over_cap": over_cap,
        "mutation_counts": dict(mutations),
        "mutations_not_in_contract": {m: n for m, n in mutations.items() if m not in active},
        "mutations_unused": [m for m in active if m not in mutations],
        "cas_counts": dict(cas_systems),
        "cas_not_in_rules": {c: n for c, n in cas_systems.items() if c not in cas_allowed},
        "cell_type_counts": dict(cell_types),
        "cell_type_mismatch": {c: n for c, n in cell_types.items() if c != cell_type},
        "guide_lengths": dict(lengths),
        "bad_guide_lengths": {k: v for k, v in lengths.items() if k not in (20, 23)},
        "mutation_map_missing": [m for m in active if m not in reference.get("mutation_map", {})],
    }


def print_preflight(checks: dict, contract: dict) -> None:
    print("\n  PREFLIGHT — submission against this contract")
    print(f"    rows                       {checks['rows']}"
          + (f"  (cap {checks['max_experiments']})" if checks["max_experiments"] else ""))
    print(f"    guide lengths              {checks['guide_lengths']}")
    print(f"    mutations                  {checks['mutation_counts']}")
    print(f"    cas systems                {checks['cas_counts']}")
    print(f"    cell types                 {checks['cell_type_counts']}"
          f"  (contract: {contract.get('cell_type')})")

    warnings = []
    if checks["rows_over_cap"]:
        warnings.append(f"{checks['rows_over_cap']} rows over max_experiments — "
                        "truncate_submission will cut them before scoring")
    if checks["duplicate_experiment_ids"]:
        warnings.append(f"{checks['duplicate_experiment_ids']} duplicate experiment_ids — "
                        "dropped by truncate_submission")
    if checks["unusable_experiment_ids"]:
        warnings.append(f"{checks['unusable_experiment_ids']} rows with a missing or blank "
                        "experiment_id — dropped by truncate_submission")
    if checks["mutations_not_in_contract"]:
        warnings.append(f"mutations not in active_mutations {checks['mutations_not_in_contract']} — "
                        "every one of those rows is rejected as mutation_not_allowed")
    if checks["mutations_unused"]:
        warnings.append(f"active mutations with no rows {checks['mutations_unused']} — "
                        "an empty coverage cell costs stage 5 multiplicatively")
    if checks["cell_type_mismatch"]:
        warnings.append(f"cell_type mismatch {checks['cell_type_mismatch']} — "
                        "rejected as cell_type_mismatch")
    if checks["cas_not_in_rules"]:
        warnings.append(f"cas systems outside rules.cas_systems {checks['cas_not_in_rules']} — "
                        "stage 1 only checks the PAM, but stage 5 scores coverage over the "
                        "rules list, so those rows are dead weight there")
    if checks["bad_guide_lengths"]:
        warnings.append(f"guide lengths outside (20, 23) {checks['bad_guide_lengths']} — "
                        "rejected as invalid_length")
    if checks["mutation_map_missing"]:
        warnings.append(f"active mutations absent from the reference mutation_map "
                        f"{checks['mutation_map_missing']} — contract and hbb_reference disagree")

    if warnings:
        print()
        for warning in warnings:
            print(f"    ! {warning}")
    else:
        print("    ok                         nothing that will cost a row before stage 1")


# --------------------------------------------------------------------------------------------
# Stage 1 + 2 detail
# --------------------------------------------------------------------------------------------

def offtarget_hits(guide: str, cas: str, kmer_index: dict) -> tuple[str, int]:
    """The seed k-mer stage 2 looks up, and how many times it occurs in the indexed window.

    Mirrors ``stage12.offtarget_uniqueness``: Cas9 seeds on the PAM-proximal 12 nt, Cas12a on the
    PAM-distal 12. The hit count is what the 1.0 / 0.7 / 0.4 / 0.1 factor is cut from, so it
    explains a row's off-target penalty in a way the factor alone does not.
    """
    seed = guide[-12:] if cas == "Cas9" else guide[:12]
    return seed, len(kmer_index.get(seed, []))


def stage12_detail(rows: list, valid: list, invalid: list, ctx: G.Context,
                   submitted: int, checks: dict) -> dict:
    """Per-row stage 1 verdict and stage 2 arithmetic. Neither stage reads the seed."""
    total_weighted = sum(entry["stage2"]["weighted_score"] for entry in valid)

    experiments = []
    by_mutation: dict[str, dict] = defaultdict(lambda: {"n": 0, "weighted_score": 0.0})
    by_cell: dict[tuple, dict] = defaultdict(lambda: {"n": 0, "weighted_score": 0.0})

    for entry in valid:
        experiment = entry["experiment"]
        features = entry["features"]
        weighted = entry["stage2"]["weighted_score"]
        seed_kmer, hits = offtarget_hits(experiment["guideRNA"], experiment["cas_system"],
                                         ctx.kmer_index)
        experiments.append({
            "experiment_id": experiment["experiment_id"],
            "design": dict(experiment),
            "stage1": {"passed": True},
            "stage2": {
                "gc": features["gc"],
                "gc_score": features["gc_score"],
                "distance_to_mutation": features["distance_to_mutation"],
                "dist_score": features["dist_score"],
                "consistency": features["consistency"],
                "base_structural_score": (0.625 * features["gc_score"]
                                          + 0.375 * features["dist_score"]),
                "offtarget_seed_kmer": seed_kmer,
                "offtarget_kmer_hits": hits,
                "offtarget_factor": features["offtarget_factor"],
                "structural_score": entry["stage2"]["structural_score"],
                "mutation_weight": features["mutation_weight"],
                "weighted_score": weighted,
                "share_of_total": (weighted / total_weighted) if total_weighted else 0.0,
                "cell_type": features["cell_type"],
                "cell_type_accessibility": features["cell_type_accessibility"],
                "mutation_region": features["mutation_region"],
                "region_energy_offset": features["region_energy_offset"],
            },
        })
        mutation = experiment["mutation"]
        by_mutation[mutation]["n"] += 1
        by_mutation[mutation]["weighted_score"] += weighted
        by_mutation[mutation]["mutation_weight"] = features["mutation_weight"]
        cell = (mutation, experiment["cas_system"], experiment.get("strand"))
        by_cell[cell]["n"] += 1
        by_cell[cell]["weighted_score"] += weighted

    features = [entry["features"] for entry in valid]
    return {
        "stage": "1+2 — structural gate and structural score",
        "seed_dependent": False,
        "formula": "structural_score = (0.625*gc_score + 0.375*dist_score) * offtarget_factor ; "
                   "weighted_score = structural_score * mutation_weight",
        "rows": {
            "submitted": submitted,
            "scored": len(rows),
            "dropped_by_truncation": submitted - len(rows),
            "valid": len(valid),
            "invalid": len(invalid),
        },
        "gate": {
            "reasons": dict(Counter(item["reason"] for item in invalid)),
            "rejected": [{"experiment_id": item["experiment"].get("experiment_id"),
                          "reason": item["reason"], "design": item["experiment"]}
                         for item in invalid],
        },
        "contract_terms": {
            "base_padding": ctx.contract["rules"]["base_padding"],
            "max_mismatches": ctx.contract["rules"].get("max_mismatches"),
            "proximity_gate": ctx.contract["rules"].get("proximity_gate", False),
            "max_experiments": ctx.contract["rules"].get("max_experiments"),
            "mutation_weights": ctx.contract.get("mutation_weights"),
            "cell_type": ctx.contract.get("cell_type"),
            "accessibility": (features[0]["cell_type_accessibility"] if features else None),
        },
        "aggregates": {
            "guide_lengths": dict(Counter(len(row["guideRNA"]) for row in rows)),
            "gc": spread([f["gc"] for f in features]),
            "gc_score": spread([f["gc_score"] for f in features]),
            "distance_to_mutation": spread([f["distance_to_mutation"] for f in features]),
            "dist_score": spread([f["dist_score"] for f in features]),
            "consistency": spread([f["consistency"] for f in features]),
            "offtarget_factor_counts": dict(Counter(f["offtarget_factor"] for f in features)),
            "structural_score": spread([e["stage2"]["structural_score"] for e in valid]),
            "weighted_score": spread([e["stage2"]["weighted_score"] for e in valid]),
            "total_weighted_score": float(total_weighted),
        },
        "by_mutation": {mutation: {**stats, "share_of_total":
                                   (stats["weighted_score"] / total_weighted)
                                   if total_weighted else 0.0}
                        for mutation, stats in by_mutation.items()},
        "by_joint_cell": {" | ".join(str(part) for part in cell): stats
                          for cell, stats in sorted(by_cell.items(), key=lambda kv: str(kv[0]))},
        "preflight": checks,
        "experiments": experiments,
    }


# --------------------------------------------------------------------------------------------
# Stage 3 detail — the RNG stream, replayed
# --------------------------------------------------------------------------------------------

def stage3_replay(entry: dict, round_seed: int) -> dict:
    """Every draw ``stage3.simulate`` makes for one row, in the order it makes them.

    The stream is replayed rather than inferred: ``simulate`` seeds
    ``random.Random(experiment_seed(round_seed, design))`` and then draws the microhomology coin,
    the cut coin, and — only if it cut — the repair-mode uniform, followed by the indel draw. Every
    row's outcome is therefore a threshold comparison this function can show explicitly, which is
    what makes a near miss (``margin_over_cut_p`` a few thousandths) distinguishable from a
    hopeless one.
    """
    experiment = entry["experiment"]
    features = stage3.extract_features(entry)
    seed = stage3.experiment_seed(round_seed, experiment)
    rng = random.Random(seed)

    energy = stage3.sequence_energy(features)
    gc = features["gc"]
    p_mh = min(0.6, (gc * (1 - gc)) * 2.2)
    mh_draw = rng.random()                 # draw 1 — microhomology_trigger
    mh = mh_draw < p_mh
    cas = experiment["cas_system"]
    cut_p = stage3.cut_probability(cas, energy)
    cut_draw = rng.random()                # draw 2 — the cut coin

    detail = {
        "experiment_id": experiment["experiment_id"],
        "design": dict(experiment),
        "experiment_seed": seed,
        "energy": energy,
        "energy_terms": {
            "accessibility": features["cell_type_accessibility"],
            "gc_term": 1.8 * gc,
            "distance_term": 0.6 * math.exp(-features["distance"] / 1500),
            "region_energy_offset": features.get("region_energy_offset", 0.0),
            "clamped": energy in (0.0, 1.0),
        },
        "p_mh": p_mh,
        "mh_draw": mh_draw,
        "mh": mh,
        "cut_p": cut_p,
        "cut_p_ceiling": G.cut_p_ceiling_for(cas),
        "cut_draw": cut_draw,
        "margin_over_cut_p": cut_draw - cut_p,
        "cut": cut_draw <= cut_p,
    }

    if cut_draw > cut_p:
        detail.update({"outcome": "no_cut", "indel_length": 0, "repair": None})
        return detail

    hdr = (0.32 if cas == "Cas9" else 0.24) + 0.35 * energy
    mh_nhej = 0.30 if mh else 0.12
    blunt = 0.35
    total = hdr + mh_nhej + blunt
    mode_draw = rng.random()               # draw 3 — repair_mode
    scaled = mode_draw * total
    if scaled < hdr:
        mode = "HDR"
    elif scaled < hdr + mh_nhej:
        mode = "MH_NHEJ"
    else:
        mode = "BLUNT_NHEJ"
    indel = stage3.sample_indel_length(mode, rng)   # draw 4, except for HDR which needs none

    detail.update({
        "repair": {
            "weights": {"HDR": hdr, "MH_NHEJ": mh_nhej, "BLUNT_NHEJ": blunt, "total": total},
            "probabilities": {"HDR": hdr / total, "MH_NHEJ": mh_nhej / total,
                              "BLUNT_NHEJ": blunt / total},
            "draw": mode_draw,
            "scaled_draw": scaled,
        },
        "outcome": mode,
        "indel_length": indel,
    })
    return detail


def stage3_detail(valid: list, results: list, summary: dict, ctx: G.Context) -> dict:
    """Per-row simulation detail plus the aggregates stage 4 and 5 are then computed from."""
    replays = [stage3_replay(entry, ctx.seed) for entry in valid]

    # The replay must reproduce the pipeline exactly — same seed, same stream. If it ever does not,
    # every number downstream of it in this file is suspect, so it is checked rather than assumed.
    mismatches = [
        {"experiment_id": r["experiment_id"], "replay": [r["outcome"], r["indel_length"]],
         "pipeline": [official["outcome"], official["indel_length"]]}
        for r, official in zip(replays, results)
        if r["outcome"] != official["outcome"] or r["indel_length"] != official["indel_length"]
    ]

    cut_ps = [r["cut_p"] for r in replays]
    outcomes = Counter(r["outcome"] for r in replays)
    by_cas: dict[str, Counter] = defaultdict(Counter)
    for replay in replays:
        by_cas[replay["design"]["cas_system"]][replay["outcome"]] += 1

    return {
        "stage": "3 — biophysical simulation",
        "seed_dependent": True,
        "round_seed": ctx.seed,
        "seeding": "random.Random(sha256(seed|mutation|cas|guide|start|strand) % 2**32)",
        "replay_matches_pipeline": not mismatches,
        "replay_mismatches": mismatches,
        "pipeline_summary": summary,
        "aggregates": {
            "n": len(replays),
            "outcomes": dict(outcomes),
            "cut_rate": (sum(1 for r in replays if r["outcome"] != "no_cut") / len(replays)
                         if replays else 0.0),
            "expected_no_cut": float(sum(1.0 - p for p in cut_ps)),
            "cut_p": spread(cut_ps),
            "cut_p_by_cas": {cas: spread([r["cut_p"] for r in replays
                                          if r["design"]["cas_system"] == cas])
                             for cas in sorted({r["design"]["cas_system"] for r in replays})},
            "energy": spread([r["energy"] for r in replays]),
            "p_mh": spread([r["p_mh"] for r in replays]),
            "mh_true": sum(1 for r in replays if r["mh"]),
            "indel_length": spread([r["indel_length"] for r in replays]),
            "indel_histogram": dict(sorted(Counter(r["indel_length"] for r in replays).items())),
            "outcomes_by_cas": {cas: dict(counts) for cas, counts in sorted(by_cas.items())},
            "outcomes_by_mh": {f"mh={mh}|{outcome}": count for (mh, outcome), count in
                               sorted(Counter((r["mh"], r["outcome"]) for r in replays).items(),
                                      key=lambda kv: str(kv[0]))},
            "mutation_breakdown": stage3.group_by_mutation(results) if results else {},
            "no_cut_rows": [r["experiment_id"] for r in replays if r["outcome"] == "no_cut"],
            "closest_call": min((r["margin_over_cut_p"] for r in replays
                                 if r["outcome"] == "no_cut"), default=None),
        },
        "experiments": replays,
    }


# --------------------------------------------------------------------------------------------
# Stage 4 detail — per target, per fold
# --------------------------------------------------------------------------------------------

def evaluate_folds(X, y, sample_weight, fold_seed: int, n_splits: int = 5) -> dict:
    """``stage4.evaluate``, fold by fold.

    Mirrors it exactly — same ``KFold(shuffle=True, random_state=contract seed)``, same
    ``RandomForestRegressor(n_estimators=200, random_state=42, max_depth=12)``, same
    ``mutation_weight`` sample weights on both fit and metric — and additionally keeps each fold's
    numbers and feature importances. The aggregates it returns are cross-checked against the
    pipeline's own ``model_results`` by the caller.
    """
    n = len(X)
    k = min(n_splits, n) if n > 1 else 1
    kf = KFold(n_splits=max(k, 2), shuffle=True, random_state=fold_seed)

    folds, r2s, maes, residual_stds = [], [], [], []
    importances: dict[str, list[float]] = {column: [] for column in X.columns}

    for index, (train_idx, test_idx) in enumerate(kf.split(X), start=1):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        sw_train = sample_weight.iloc[train_idx] if sample_weight is not None else None
        sw_test = sample_weight.iloc[test_idx] if sample_weight is not None else None

        model = RandomForestRegressor(n_estimators=200, random_state=42, max_depth=12)
        model.fit(X_train, y_train, sample_weight=sw_train)
        pred = model.predict(X_test)

        r2 = float(r2_score(y_test, pred, sample_weight=sw_test))
        mae = float(mean_absolute_error(y_test, pred, sample_weight=sw_test))
        residual_std = float(np.std(y_test - pred))
        r2s.append(r2)
        maes.append(mae)
        residual_stds.append(residual_std)
        for column, importance in zip(X.columns, model.feature_importances_):
            importances[column].append(float(importance))

        folds.append({
            "fold": index,
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "r2": r2,
            "mae": mae,
            "residual_std": residual_std,
            "y_test": {"mean": float(np.mean(y_test)), "std": float(np.std(y_test)),
                       "min": float(np.min(y_test)), "max": float(np.max(y_test))},
            "prediction": {"mean": float(np.mean(pred)), "std": float(np.std(pred)),
                           "min": float(np.min(pred)), "max": float(np.max(pred))},
            "feature_importances": {column: float(importance)
                                    for column, importance in zip(X.columns,
                                                                  model.feature_importances_)},
        })

    return {
        "r2_mean": float(np.mean(r2s)),
        "r2_std": float(np.std(r2s)),
        "mae_mean": float(np.mean(maes)),
        "mae_std": float(np.std(maes)),
        "residual_std_mean": float(np.mean(residual_stds)),
        "n_folds": len(r2s),
        "folds": folds,
        "mean_feature_importances": {column: float(np.mean(values))
                                     for column, values in importances.items()},
    }


def stage4_detail(valid: list, results: list, ctx: G.Context, official: dict,
                  n_folds: int = 5) -> dict:
    """Every quantity stage 4 turns into ``consistency_factor``, per target and per fold.

    ``consistency_score = (0.7*max(avg_r2, 0) + 0.3*(1 - avg_nmae)) * 100`` over the three targets
    ``is_cut``, ``is_hdr`` and ``indel_length``. The two ways it reaches 1.0 are visible here: a
    forest that genuinely predicts the outcomes (positive r2), or targets so degenerate that
    ``r2_score`` hits its zero-variance branch and ``normalized_mae`` short-circuits on
    ``std < 1e-9`` — which is what a single-outcome construction does.
    """
    frame3 = stage4.flatten_stage3(results)
    frame12 = stage4.flatten_stage12(valid)

    detail = {
        "stage": "4 — consistency (RandomForest under KFold)",
        "seed_dependent": True,
        "fold_seed": ctx.seed,
        "formula": "consistency_score = (0.7*max(avg_r2, 0) + 0.3*(1 - avg_nmae)) * 100 ; "
                   "consistency_factor = clamp(consistency_score/100, 0, 1)",
        "model": {"estimator": "RandomForestRegressor", "n_estimators": 200, "max_depth": 12,
                  "random_state": 42, "n_splits": n_folds, "shuffle": True,
                  "sample_weight": "mutation_weight (fit and metrics)"},
        "rows": {"stage12": len(frame12), "stage3": len(frame3)},
        "pipeline_output": official,
    }

    if len(frame12) < 2 or len(frame3) < 2:
        detail["skipped"] = ("fewer than two rows on one side of the join — stage 4 writes a clean "
                             "zero rather than raising")
        detail["targets"] = {}
        return detail

    slim = frame12[["experiment_id", "guideRNA", "start", "stage2_score",
                    "mutation_weight", "weighted_score"]]
    merged = frame3.merge(slim, on="experiment_id", how="inner")
    detail["rows"]["joined"] = len(merged)
    detail["rows"]["join_loss"] = len(frame3) - len(merged)
    if len(merged) == 0:
        detail["skipped"] = "no experiment_id matched across the join — clean zero"
        detail["targets"] = {}
        return detail

    X = stage4.build_X(merged)
    y = stage4.build_y(merged)
    sample_weight = merged["mutation_weight"]

    detail["features"] = {
        "columns": list(X.columns),
        "stats": {column: {**spread(X[column].tolist()),
                           "n_unique": int(X[column].nunique()),
                           "constant": bool(X[column].nunique() <= 1)}
                  for column in X.columns},
        "constant_columns": [column for column in X.columns if X[column].nunique() <= 1],
        "correlation_with_targets": {
            column: {target: (float(np.corrcoef(X[column], y[target])[0, 1])
                              if X[column].nunique() > 1 and y[target].nunique() > 1 else None)
                     for target in y.columns}
            for column in X.columns
        },
    }
    detail["sample_weight"] = {"column": "mutation_weight", **spread(sample_weight.tolist())}

    targets = {}
    for target in y.columns:
        column = y[target]
        scale = float(np.std(column))
        folds = evaluate_folds(X, column, sample_weight=sample_weight,
                               fold_seed=ctx.seed, n_splits=n_folds)
        nmae = stage4.normalized_mae(folds["mae_mean"], column)
        pipeline = (official.get("model_results") or {}).get(target, {})
        targets[target] = {
            "y": {**spread(column.tolist()),
                  "n_unique": int(column.nunique()),
                  "value_counts": {str(k): int(v) for k, v in column.value_counts().items()}
                  if column.nunique() <= 12 else None,
                  "degenerate": scale < 1e-9},
            "r2_mean": folds["r2_mean"],
            "r2_std": folds["r2_std"],
            "mae_mean": folds["mae_mean"],
            "mae_std": folds["mae_std"],
            "residual_std_mean": folds["residual_std_mean"],
            "n_folds": folds["n_folds"],
            "y_std_scale": scale,
            "normalized_mae": nmae,
            "nmae_short_circuited": scale < 1e-9,
            "mean_feature_importances": folds["mean_feature_importances"],
            "folds": folds["folds"],
            "pipeline_model_results": pipeline,
            "delta_vs_pipeline": {
                key: (folds[key] - pipeline[key])
                for key in ("r2_mean", "r2_std", "mae_mean", "mae_std", "residual_std_mean")
                if key in pipeline
            },
        }

    avg_r2 = float(np.mean([t["r2_mean"] for t in targets.values()]))
    avg_nmae = float(np.mean([t["normalized_mae"] for t in targets.values()]))
    consistency_score = (0.7 * max(avg_r2, 0) + 0.3 * (1 - avg_nmae)) * 100
    consistency_factor = 0.0 if math.isnan(consistency_score) else max(
        0.0, min(1.0, consistency_score / 100.0))
    total_weighted = float(merged["weighted_score"].sum())

    detail.update({
        "targets": targets,
        "avg_r2": avg_r2,
        "avg_nmae": avg_nmae,
        "consistency_score": float(consistency_score),
        "consistency_factor": consistency_factor,
        "total_weighted_score": total_weighted,
        "final_reward": total_weighted * consistency_factor,
        "terms": {
            "r2_term": 0.7 * max(avg_r2, 0) * 100,
            "mae_term": 0.3 * (1 - avg_nmae) * 100,
        },
        "delta_vs_pipeline": {
            "consistency_score": float(consistency_score) - official.get("consistency_score", 0.0),
            "consistency_factor": consistency_factor - official.get("consistency_factor", 0.0),
            "total_weighted_score": total_weighted - official.get("total_weighted_score", 0.0),
        },
    })
    return detail


# --------------------------------------------------------------------------------------------
# Stage 5 detail — the six-way geometric mean, term by term
# --------------------------------------------------------------------------------------------

def coverage_detail(counter: Counter, support: list) -> dict:
    """One coverage term: the counts, the entropy they carry and the ratio stage 5 scores."""
    support = list(support)
    counts = [counter.get(category, 0) for category in support]
    h = stage5.shannon_entropy(counts)
    h_max = math.log2(len(support)) if len(support) > 1 else 0.0
    return {
        "support": [list(c) if isinstance(c, tuple) else c for c in support],
        "support_size": len(support),
        "counts": {" | ".join(str(p) for p in c) if isinstance(c, tuple) else str(c):
                   counter.get(c, 0) for c in support},
        "observed_outside_support": {" | ".join(str(p) for p in c) if isinstance(c, tuple)
                                    else str(c): n
                                    for c, n in counter.items() if c not in support},
        "occupied_cells": sum(1 for c in counts if c > 0),
        "empty_cells": [list(c) if isinstance(c, tuple) else c
                        for c, count in zip(support, counts) if count == 0],
        "entropy": h,
        "entropy_max": h_max,
        "ratio": (h / h_max) if h_max > 0 else 1.0,
    }


def stage5_detail(valid: list, results: list, ctx: G.Context, official: dict,
                  total_weighted: float, consistency_factor: float, k: int = 12) -> dict:
    """Every fidelity term, its share of the geometric mean, and what fixing it would be worth.

    The mean is geometric, so the terms multiply: a single empty (mutation x cas x strand) cell
    drags the whole score down by its own sixth root, and a term at the 1e-9 clip takes the
    product to roughly 0.03x. ``cost_if_perfect`` prices each term at the current score.
    """
    contract = ctx.contract
    active = list(contract["active_mutations"])
    cas_systems = list(contract["rules"].get("cas_systems", ["Cas9", "Cas12a"]))
    strands = ["+", "-"]

    experiments = [entry["experiment"] for entry in valid]
    guides = [experiment["guideRNA"] for experiment in experiments]

    mutation = coverage_detail(Counter(e["mutation"] for e in experiments), active)
    cas = coverage_detail(Counter(e["cas_system"] for e in experiments), cas_systems)
    strand = coverage_detail(Counter(e.get("strand") for e in experiments), strands)
    joint = coverage_detail(
        Counter((e["mutation"], e["cas_system"], e.get("strand")) for e in experiments),
        [(m, c, s) for m in active for c in cas_systems for s in strands],
    )

    pool: Counter = Counter()
    total_kmers = 0
    for guide in guides:
        kmers = stage5.extract_kmers(guide, k)
        pool.update(kmers)
        total_kmers += len(kmers)
    kmer_entropy = stage5.shannon_entropy(list(pool.values())) if total_kmers > 1 else 0.0
    kmer_max = math.log2(total_kmers) if total_kmers > 1 else 0.0
    kmer = {
        "k": k,
        "total_kmers": total_kmers,
        "distinct_kmers": len(pool),
        "entropy": kmer_entropy,
        "entropy_max": kmer_max,
        "ratio": (kmer_entropy / kmer_max) if kmer_max > 0 else 0.0,
        "most_repeated": [[kmer_seq, count] for kmer_seq, count in pool.most_common(10)
                          if count > 1],
    }
    duplicates = {guide: count for guide, count in Counter(guides).items() if count > 1}
    guide_term = {
        "guides": len(guides),
        "distinct": len(set(guides)),
        "ratio": (len(set(guides)) / len(guides)) if guides else 0.0,
        "duplicated_guides": duplicates,
    }

    values = {
        "mutation_coverage_entropy_ratio": mutation["ratio"],
        "cas_system_coverage_entropy_ratio": cas["ratio"],
        "strand_coverage_entropy_ratio": strand["ratio"],
        "joint_coverage_entropy_ratio": joint["ratio"],
        "kmer_diversity_entropy_ratio": kmer["ratio"],
        "distinct_guide_ratio": guide_term["ratio"],
    }
    score = G.stage5.geometric_mean([values[name] for name in STAGE5_TERMS])
    factor = max(0.0, min(1.0, score))
    final_score = total_weighted * consistency_factor * factor

    terms = {}
    for name in STAGE5_TERMS:
        value = values[name]
        clipped = max(value, 1e-9)
        # Geometric mean over six terms: each contributes its own sixth root, so replacing one term
        # with 1.0 scales the product by clipped**(-1/6).
        contribution = clipped ** (1.0 / len(STAGE5_TERMS))
        perfect_factor = min(1.0, factor / contribution) if contribution > 0 else factor
        terms[name] = {
            "value": value,
            "clipped_at_1e-9": value < 1e-9,
            "contribution_to_geometric_mean": contribution,
            "log_share": (math.log(clipped) / sum(math.log(max(values[n], 1e-9))
                                                  for n in STAGE5_TERMS))
            if any(values[n] < 1.0 for n in STAGE5_TERMS) else 0.0,
            "factor_if_perfect": perfect_factor,
            "final_score_if_perfect": total_weighted * consistency_factor * perfect_factor,
            "cost_now": total_weighted * consistency_factor * (perfect_factor - factor),
        }

    return {
        "stage": "5 — distribution fidelity",
        "seed_dependent": False,
        "formula": "geometric mean of six ratios, clipped at 1e-9; "
                   "final_score = total_weighted_score * consistency_factor * fidelity_factor",
        "pipeline_output": official,
        "terms": terms,
        "distribution_fidelity_score": score,
        "distribution_fidelity_factor": factor,
        "delta_vs_pipeline": score - (official.get("distribution_fidelity_score") or 0.0),
        "coverage": {
            "mutation": mutation,
            "cas_system": cas,
            "strand": strand,
            "joint": joint,
        },
        "kmer_diversity": kmer,
        "distinct_guides": guide_term,
        "cas_specific_shift_diagnostic": official.get("cas_specific_shift_diagnostic", {}),
        "note": ("the cas_specific_shift diagnostic is reported by the validator but not scored; "
                 "an empty joint cell is scored, and costs the whole product its sixth root"),
    }


# --------------------------------------------------------------------------------------------
# Console output for the two stages that carry the most detail
# --------------------------------------------------------------------------------------------

def print_stage4_detail(detail: dict) -> str:
    lines = ["", "=" * 100,
             f"  STAGE 4 DETAIL — RandomForest under KFold (fold seed {detail['fold_seed']})",
             "=" * 100]
    if not detail.get("targets"):
        lines.append(f"    skipped: {detail.get('skipped')}")
        text = "\n".join(lines)
        print(text)
        return text

    lines.append(f"    rows joined {detail['rows'].get('joined')}   "
                 f"features {', '.join(detail['features']['columns'])}")
    constant = detail["features"]["constant_columns"]
    lines.append(f"    constant feature columns   {constant or 'none'}"
                 + ("   (a constant column carries no signal for the forest)" if constant else ""))

    for target, stats in detail["targets"].items():
        lines.append("")
        degenerate = ("  (degenerate — r2 takes its zero-variance branch and nmae "
                      "short-circuits)" if stats["y"]["degenerate"] else "")
        lines.append(f"    {target}   unique={stats['y']['n_unique']}  "
                     f"y_std={stats['y_std_scale']:.6f}{degenerate}")
        lines.append(f"      r2   mean={stats['r2_mean']:+.6f}  std={stats['r2_std']:.6f}")
        lines.append(f"      mae  mean={stats['mae_mean']:.6f}  std={stats['mae_std']:.6f}  "
                     f"nmae={stats['normalized_mae']:.6f}  "
                     f"residual_std={stats['residual_std_mean']:.6f}")
        lines.append("      fold   n_test     r2          mae      residual_std   pred_mean")
        for fold in stats["folds"]:
            lines.append(f"        {fold['fold']}      {fold['n_test']:>4}   "
                         f"{fold['r2']:+.6f}   {fold['mae']:.6f}   {fold['residual_std']:.6f}     "
                         f"{fold['prediction']['mean']:.6f}")
        importances = sorted(stats["mean_feature_importances"].items(), key=lambda kv: -kv[1])
        lines.append("      importances  " + "  ".join(f"{name}={value:.3f}"
                                                       for name, value in importances))
        delta = stats.get("delta_vs_pipeline", {})
        if delta:
            worst = max(abs(v) for v in delta.values())
            lines.append(f"      matches the validator's own model_results to {worst:.2e}")

    lines.append("")
    lines.append(f"    avg_r2 {detail['avg_r2']:+.6f}   avg_nmae {detail['avg_nmae']:.6f}")
    lines.append(f"    consistency_score {detail['consistency_score']:.6f}  "
                 f"= 0.7*max(avg_r2,0)*100 ({detail['terms']['r2_term']:.4f}) "
                 f"+ 0.3*(1-avg_nmae)*100 ({detail['terms']['mae_term']:.4f})")
    lines.append(f"    CONSISTENCY_FACTOR {detail['consistency_factor']:.6f}")
    text = "\n".join(lines)
    print(text)
    return text


def print_stage5_detail(detail: dict) -> str:
    lines = ["", "=" * 100, "  STAGE 5 DETAIL — distribution fidelity, term by term", "=" * 100,
             "",
             "    term                              value    ^(1/6)    factor if 1.0   "
             "final if 1.0   cost now"]
    for name, term in detail["terms"].items():
        label = name.replace("_entropy_ratio", "").replace("_ratio", "")
        lines.append(f"    {label:<30}{term['value']:>9.6f}{term['contribution_to_geometric_mean']:>10.6f}"
                     f"{term['factor_if_perfect']:>16.6f}{term['final_score_if_perfect']:>15.4f}"
                     f"{term['cost_now']:>11.4f}")
    lines.append("")
    lines.append(f"    DISTRIBUTION_FIDELITY {detail['distribution_fidelity_factor']:.6f}"
                 f"   (validator: {detail['pipeline_output'].get('distribution_fidelity_score')})")

    for name, section in detail["coverage"].items():
        lines.append("")
        lines.append(f"    {name} coverage — {section['occupied_cells']}/{section['support_size']} "
                     f"cells occupied, H={section['entropy']:.6f} / {section['entropy_max']:.6f} "
                     f"= {section['ratio']:.6f}")
        lines.append(f"      counts {section['counts']}")
        if section["empty_cells"]:
            lines.append(f"      ! empty {section['empty_cells']}")
        if section["observed_outside_support"]:
            lines.append(f"      ! outside the scored support (uncounted) "
                         f"{section['observed_outside_support']}")

    kmer = detail["kmer_diversity"]
    lines.append("")
    lines.append(f"    kmer diversity — {kmer['distinct_kmers']} distinct of {kmer['total_kmers']} "
                 f"{kmer['k']}-mers, H={kmer['entropy']:.6f} / {kmer['entropy_max']:.6f} "
                 f"= {kmer['ratio']:.6f}")
    if kmer["most_repeated"]:
        lines.append(f"      most repeated {kmer['most_repeated'][:5]}")
    guides = detail["distinct_guides"]
    lines.append(f"    distinct guides — {guides['distinct']} of {guides['guides']} "
                 f"= {guides['ratio']:.6f}"
                 + (f"   ! duplicated {guides['duplicated_guides']}"
                    if guides["duplicated_guides"] else ""))
    diagnostic = detail["cas_specific_shift_diagnostic"]
    if diagnostic and not diagnostic.get("insufficient_data"):
        lines.append(f"    diagnostic (not scored) — compared {diagnostic.get('compared')}, "
                     f"JSD {diagnostic.get('repair_mode_jensen_shannon_divergence')}, "
                     f"Wasserstein {diagnostic.get('indel_length_wasserstein_distance')}")
    text = "\n".join(lines)
    print(text)
    return text


def print_rejections(invalid: list[dict], submitted: int, examples: int) -> None:
    print(f"\n  REJECTED ROWS — {len(invalid)} of {submitted} scored rows failed stage 1")
    if not invalid:
        print("    none")
        return
    by_reason: dict[str, list] = {}
    for item in invalid:
        by_reason.setdefault(item["reason"], []).append(item["experiment"])
    for reason, group in sorted(by_reason.items(), key=lambda kv: -len(kv[1])):
        print(f"    {reason:<28} {len(group)}")
        for experiment in group[:examples]:
            print(f"      {str(experiment.get('experiment_id')):<12}"
                  f"{str(experiment.get('cas_system')):<8}"
                  f"{str(experiment.get('strand')):<3}"
                  f"start={experiment.get('target_alignment_start')} "
                  f"len={len(experiment.get('guideRNA') or '')} "
                  f"{str(experiment.get('mutation'))[:28]}")
        if len(group) > examples:
            print(f"      ... {len(group) - examples} more")


def print_seed_comparison(scored: dict, alt: dict, scored_seed: int, alt_seed: int) -> str:
    """What the seed alone is worth. Stages 1, 2 and 5 never read it, so only stage 3's outcomes
    and stage 4's fold shuffle can move — the rest of the difference is arithmetic on those."""
    lines = ["", "=" * 100,
             f"  SEED {alt_seed}  ->  SEED {scored_seed}   (same submission, same contract)",
             "=" * 100]
    for label, section, key in (
        ("total_weighted_score", "stage4", "total_weighted_score"),
        ("consistency_factor", "stage4", "consistency_factor"),
        ("distribution_fidelity_factor", "stage5", "distribution_fidelity_factor"),
    ):
        before, after = alt[section][key], scored[section][key]
        moved = "unchanged" if abs(before - after) < 1e-9 else f"{after - before:+.6f}"
        lines.append(f"    {label:<30} {before:>12.6f}  ->  {after:>12.6f}   ({moved})")
    before, after = alt["final_score"], scored["final_score"]
    ratio = after / before if before else 0.0
    lines.append(f"    {'final_score':<30} {before:>12.6f}  ->  {after:>12.6f}   "
                 f"({ratio:.4f}x, {100 * (ratio - 1):+.1f}%)")
    if scored["stage3"]["conforming_rows"] is not None:
        lines.append(f"    {'construction conformance':<30} "
                     f"{alt['stage3']['conforming_rows']:>12} ->  "
                     f"{scored['stage3']['conforming_rows']:>6}")
    lines.append(f"    {'cut_rate':<30} {alt['stage3']['cut_rate']:>12.4f}  ->  "
                 f"{scored['stage3']['cut_rate']:>12.4f}")
    text = "\n".join(lines)
    print(text)
    return text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--result-dir", default="data/result",
                        help="archive root; the newest folder in it is scored by default")
    parser.add_argument("--folder", default=None,
                        help="score this archive folder instead of the newest one")
    parser.add_argument("--submission", default=None,
                        help="score a loose submission file instead of an archive folder")
    parser.add_argument("--out-dir", default=None,
                        help="where the detail files are written (default: the scored folder, "
                             "or calc/ for a loose --submission)")
    parser.add_argument("--task-id", default=None,
                        help="pin the contract to this task instead of matching the folder")
    parser.add_argument("--task", default=None,
                        help="read the task from a local file instead of the backend")
    parser.add_argument("--seed", default=None,
                        help="override the contract seed; accepts a comma-joined list "
                             "(\"122,321,431\") exactly as the stamped contract carries it")
    parser.add_argument("--cell-types", default=None,
                        help="cell-type accessibility JSON; fetched from the backend if omitted")
    parser.add_argument("--construction", default="hdr",
                        choices=sorted(G.CONSTRUCTIONS) + ["none"],
                        help="construction to measure row conformance against (default hdr)")
    parser.add_argument("--compare-seed", type=int, default=None,
                        help="also report at this seed; defaults to the archived contract's seed")
    parser.add_argument("--no-compare", action="store_true",
                        help="skip the second report even if a build seed is available")
    parser.add_argument("--examples", type=int, default=3,
                        help="rejected rows to print per reason (default 3)")
    parser.add_argument("--quiet", action="store_true", help="suppress info logging")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                        format="%(levelname)s %(message)s")

    # ------------------------------------------------------------------ what is being scored
    folder: Path | None = None
    if args.submission:
        submission_path = Path(args.submission).resolve()
        out_dir = Path(args.out_dir or "calc")
    else:
        folder = (Path(args.folder) if args.folder
                  else latest_result_folder(Path(args.result_dir))).resolve()
        submission_path = folder / "submission.json"
        out_dir = Path(args.out_dir) if args.out_dir else folder

    if not submission_path.exists():
        raise SystemExit(f"{submission_path} does not exist")
    rows = read_json(submission_path)
    if not isinstance(rows, list) or not rows:
        raise SystemExit(f"{submission_path} is not a non-empty JSON array of experiments")
    submitted = len(rows)

    folder_contract = read_json(folder / "contract.json") if folder else None
    folder_reference = read_json(folder / "hbb_reference.json") if folder else None
    contract_source = str(folder / "contract.json") if folder else settings.CONTRACT_PATH
    # The miner archives its upload record alongside the submission, and that is the authoritative
    # pairing: folders are named after the archive time, not the task's created_at, so the name
    # match below is only a fallback for folders written before that. Later candidates win.
    folder_upload = read_json(folder / Path(settings.LAST_UPLOAD_PATH).name) if folder else None
    recorded_source = (f"{folder.name}/{Path(settings.LAST_UPLOAD_PATH).name}" if folder
                       else settings.LAST_UPLOAD_PATH)
    recorded_task_id = None
    for candidate in (folder_contract, folder_reference,
                      read_json(folder / "task.json") if folder else None, folder_upload):
        if isinstance(candidate, dict):
            recorded_task_id = (candidate.get("task_id") or candidate.get("id")
                                or recorded_task_id)
    # Read before the redirect, while these still point at data/.
    uploaded = read_json(settings.LAST_UPLOAD_PATH) or {}
    if folder is None:
        # No archive folder: data/contract.json is the broadcast copy and data/last_upload.json
        # records which task the miner built it for.
        folder_contract = read_json(settings.CONTRACT_PATH)
        recorded_task_id = uploaded.get("task_id")

    # The validator's own inputs and outputs live one level down, so the archived submission,
    # contract and reference are never written — truncate_submission rewrites its input in place.
    stage_dir = out_dir / "validation"
    stage_dir.mkdir(parents=True, exist_ok=True)
    mapping = redirect_paths(stage_dir)
    if Path(settings.MINER_SUBMISSION_PATH).resolve() == submission_path:
        raise SystemExit(f"the staging copy would overwrite {submission_path}; pick another "
                         "--out-dir")
    shutil.copyfile(submission_path, settings.MINER_SUBMISSION_PATH)

    # ------------------------------------------------------------------------------- the contract
    if args.task:
        task = read_json(args.task)
        if task is None:
            raise SystemExit(f"could not read {args.task}")
        origin = f"file {args.task}"
    else:
        task, origin = resolve_task(args.task_id, folder_contract,
                                    folder.name if folder else None, recorded_task_id,
                                    recorded_source)

    contract = copy.deepcopy(task["content"]["contract"])
    reference = copy.deepcopy(task["content"]["hbb_reference"])
    real_seed = contract.get("seed")
    if args.seed is not None:
        contract["seed"] = args.seed
    seeds = parse_seeds(contract["seed"])
    seed = seeds[0]                       # the one stages 1/2/5 do not care about anyway
    pairing = verify_pairing(task, folder_contract)
    pairing["compared_against"] = contract_source if folder_contract else None

    print("=" * 100)
    print(f"  submission  {submission_path}  ({submitted} rows)")
    if folder:
        print(f"  folder      {folder}")
    print(f"  task        {task.get('id')}   created {task.get('created_at')}")
    print(f"  matched by  {origin}")
    seed_note = (f"   scored at {contract['seed']} (--seed)" if args.seed is not None else "")
    print(f"  real seed   {real_seed}{seed_note}")
    if len(seeds) > 1:
        print(f"  round seeds {seeds}   ({len(seeds)} of them — stages 3-5 run per seed and the "
              f"final score is their mean)")
    print(f"  cell_type   {contract.get('cell_type')}")
    print(f"  mutations   {contract['active_mutations']}")
    print(f"  weights     {contract.get('mutation_weights')}")
    print(f"  regions     {contract.get('mutation_regions')}")
    print(f"  rules       {contract['rules']}")
    print(f"  detail ->   {out_dir}/    stages -> {stage_dir}/")
    print("=" * 100)

    if pairing.get("differences"):
        print(f"\n  ! {contract_source} and this task's contract disagree:")
        for key, sides in pairing["differences"].items():
            print(f"      {key}\n        archived {sides['archived']}\n        task     {sides['task']}")
        if not args.task_id and not args.task:
            raise SystemExit("refusing to apply this task's seed to a submission built against a "
                             "different contract — pass --task-id to override")
    elif pairing.get("verified"):
        seed_pair = (f"broadcast seed {pairing['archived_seed']}, task seed {pairing['task_seed']}"
                     if pairing["archived_seed"] != pairing["task_seed"]
                     else f"both carry seed {pairing['task_seed']}")
        print(f"\n  contract verified field by field against {contract_source} ({seed_pair})")

    if not real_seed and args.seed is None:
        print("\n  ! this task carries seed 0 — the backend has not stamped it yet, so stages 3 "
              "and 4\n    are being scored against a stream the validator will not use. Re-run "
              "once the\n    task is stamped, or pass --seed to probe a specific one.")

    # ------------------------------------------------------------------------------- preflight
    checks = preflight(rows, contract, reference)
    print_preflight(checks, contract)
    if checks["malformed_rows"]:
        print("\n  FATAL — rows missing a field stage 1 indexes directly:")
        for row in checks["malformed_rows"][:10]:
            print(f"    row {row['index']} ({row['experiment_id']}) missing {row['missing']}")
        print("    The validator raises on these rather than scoring them, and run_validation "
              "swallows\n    the exception — the miner ends up with no score at all. Fix the rows "
              "and re-run.")
        # Not validation.json: that file is the marker for "the five stages ran", and
        # neurons/miner.py skips a folder that already has one.
        write_json(out_dir / "validation_error.json",
                   {"error": "malformed rows", "preflight": checks})
        return 1

    # ------------------------------------------------------------------------------ cell types
    if args.cell_types:
        cell_types = read_json(args.cell_types)
        if cell_types is None:
            raise SystemExit(f"could not read {args.cell_types}")
    else:
        cell_types = G.fetch_cell_types()
    accessibility = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    print(f"\n  accessibility {accessibility} for {contract.get('cell_type')}")

    # -------------------------------------------------------- the validator's own file pipeline
    write_json(Path(settings.CONTRACT_PATH), contract)
    write_json(Path(settings.HBB_REFERENCE_PATH), reference)

    print("\n[1/4] loading chr11 and the k-mer index")
    ctx = G.build_context(contract, reference, cell_types)

    seed_label = f"seed {seed}" if len(seeds) == 1 else f"{len(seeds)} seeds {seeds}"
    print(f"[2/4] running the validator's five stages over {stage_dir}/ at {seed_label}")
    official_score = G.verify_with_validator(ctx)

    scored_rows = read_json(settings.MINER_SUBMISSION_PATH)
    valid = read_json(settings.VALID_EXPERIMENTS_PATH)
    invalid = read_json(settings.INVALID_EXPERIMENTS_PATH)
    results = read_json(settings.STAGE3_DATASET)
    stage3_summary = read_json(settings.STAGE3_SUMMARY_PATH) or {}
    stage4_output = read_json(settings.FINAL_REWARD_PATH) or {}
    stage5_output = read_json(settings.DISTRIBUTION_FIDELITY_PATH) or {}
    dropped = submitted - len(scored_rows)
    if dropped:
        print(f"  truncate_submission cut {dropped} of {submitted} rows "
              f"(cap {checks['max_experiments']}, duplicate or blank experiment_ids)")

    # ------------------------------------------------------ the same numbers, stage by stage
    print("[3/4] measuring every stage quantity from the artifacts it just wrote")
    cfg = None if args.construction == "none" else G.GenConfig(construction=args.construction)

    # Stages 3-5 are re-derived per round seed and averaged, matching benchmark_submission. Stage
    # 12 is seed-independent so the artifacts on disk serve every seed; the per-seed stage-3 results
    # are recomputed in memory rather than read back, because the pipeline leaves only the *last*
    # seed's dataset in data/ and reading that would silently report one seed as if it were the mean.
    per_seed = []
    for sd in seeds:
        sctx = replace(ctx, contract={**contract, "seed": sd})
        sresults = [stage3.simulate(entry, sd) for entry in valid]
        per_seed.append((sd, sctx, sresults, G.stage_report(scored_rows, valid, sresults, sctx, cfg)))

    report = mean_reports([item[3] for item in per_seed])
    if len(seeds) == 1:
        results = per_seed[0][2]
        report_text = captured(report,
                               f"{submission_path.name} — task {task.get('id')} — seed {seed}")
        print(report_text)
    else:
        chunks = []
        for sd, _sctx, _res, rep in per_seed:
            chunks.append(captured(rep, f"{submission_path.name} — task {task.get('id')} — "
                                        f"seed {sd}"))
        header = (f"\n  MEAN OVER {len(seeds)} ROUND SEEDS {seeds} — this is what the validator "
                  f"scores\n" + "  " + "-" * 76 + "\n")
        for label, path in (("total_weighted_score", ("stage4", "total_weighted_score")),
                            ("consistency_factor", ("stage4", "consistency_factor")),
                            ("distribution_fidelity_factor",
                             ("stage5", "distribution_fidelity_factor")),
                            ("final_score", ("final_score",))):
            spread = seed_spread(per_seed, path)
            node = report
            for key in path:
                node = (node or {}).get(key) if isinstance(node, dict) else None
            mean_value = node if isinstance(node, (int, float)) else 0.0
            header += (f"  {label:<32}{mean_value:>12.4f}"
                       f"   per-seed {min(spread):.4f} .. {max(spread):.4f}\n")
        report_text = "\n".join(chunks) + header
        print(report_text)
        results = per_seed[0][2]

    print_rejections(invalid, len(scored_rows), args.examples)

    delta = abs(official_score["final_score"] - report["final_score"])
    mean_note = "" if len(seeds) == 1 else f" (mean of {len(seeds)} seeds)"
    print(f"\n  benchmark_submission final_score = {official_score['final_score']:.6f}{mean_note}")
    print(f"  in-memory replica                = {report['final_score']:.6f}{mean_note}")
    print(f"  delta                            = {delta:.9f}")
    if delta > 1e-6:
        print("  ! the replica and the validator disagree — if this task has several round seeds, "
              "check\n    that both averaged the same set")

    # -------------------------------------------------------------------- per-stage detail files
    print(f"\n[4/4] writing per-stage detail into {out_dir}/")
    detail12 = stage12_detail(scored_rows, valid, invalid, ctx, submitted, checks)

    # Stage 12 is seed-independent, so one file covers every seed. Stages 3-5 are not: with several
    # round seeds each gets its own detail file, and the flat names keep pointing at the first seed
    # so a single-seed run is unchanged.
    extra_detail = {}
    detail3 = detail4 = detail5 = None
    for position, (sd, sctx, sresults, _rep) in enumerate(per_seed):
        d3 = stage3_detail(valid, sresults, stage3_summary if len(seeds) == 1 else {}, sctx)
        d4 = stage4_detail(valid, sresults, sctx, stage4_output if len(seeds) == 1 else {})
        d5 = stage5_detail(valid, sresults, sctx, stage5_output if len(seeds) == 1 else {},
                           total_weighted=d4.get("total_weighted_score", 0.0),
                           consistency_factor=d4.get("consistency_factor", 0.0))
        if position == 0:
            detail3, detail4, detail5 = d3, d4, d5
        if len(seeds) > 1:
            extra_detail[f"stage3_detail.seed{sd}.json"] = d3
            extra_detail[f"stage4_detail.seed{sd}.json"] = d4
            extra_detail[f"stage5_detail.seed{sd}.json"] = d5

    stage4_text = print_stage4_detail(detail4)
    stage5_text = print_stage5_detail(detail5)

    if not detail3["replay_matches_pipeline"]:
        print(f"\n  ! the stage-3 replay disagrees with the pipeline on "
              f"{len(detail3['replay_mismatches'])} rows — see stage3_detail.json")

    # ------------------------------------------------------------ optionally, at a second seed
    compare_seed = args.compare_seed
    if compare_seed is None and not args.no_compare and folder_contract:
        if pairing.get("verified") and folder_contract.get("seed") != seed:
            compare_seed = folder_contract.get("seed")
            print(f"\n  {contract_source} carries seed {compare_seed} — reporting it too, so "
                  "the\n  difference is exactly what the stamped seed was worth "
                  "(--no-compare skips this)")

    comparison_text = ""
    alt_report = None
    if compare_seed is not None and compare_seed != seed:
        alt_ctx = replace(ctx, contract={**contract, "seed": compare_seed})
        alt_results = [stage3.simulate(entry, compare_seed) for entry in valid]
        alt_report = G.stage_report(scored_rows, valid, alt_results, alt_ctx, cfg)
        alt_text = captured(alt_report, f"{submission_path.name} — seed {compare_seed}")
        print(alt_text)
        comparison_text = alt_text + print_seed_comparison(report, alt_report, seed, compare_seed)

    # ------------------------------------------------------------------------------- artifacts
    artifacts = {
        "stage12_detail.json": "per-row stage 1 verdict and stage 2 arithmetic",
        "stage3_detail.json": "per-row RNG replay of the simulation",
        "stage4_detail.json": "per-target, per-fold consistency detail",
        "stage5_detail.json": "distribution fidelity, term by term",
        **({} if len(seeds) == 1 else
           {"stage{3,4,5}_detail.seed<N>.json": f"the same, per round seed ({len(seeds)} of them); "
                                                "the flat files above are the first seed"}),
        "validation.json": "headline numbers and the task pairing",
        "validation.txt": "the printed five-stage report",
        "validation/": "the validator's own inputs and outputs",
    }
    summary = {
        "submission": str(submission_path),
        "folder": str(folder) if folder else None,
        "task_id": task.get("id"),
        "task_created_at": task.get("created_at"),
        "task_matched_by": origin,
        "contract_pairing": pairing,
        "real_seed": real_seed,
        "scored_seed": seed,
        "scored_seeds": seeds,
        "per_seed_final_score": {str(sd): rep.get("final_score")
                                 for sd, _c, _r, rep in per_seed},
        "final_score_is_mean_of": len(seeds),
        "compare_seed": compare_seed if alt_report else None,
        "cell_type": contract.get("cell_type"),
        "accessibility": accessibility,
        "construction": None if cfg is None else cfg.construction,
        "rows": {
            "submitted": submitted,
            "scored": len(scored_rows),
            "dropped_by_truncation": dropped,
            "valid": len(valid),
            "invalid": len(invalid),
            "rejection_reasons": dict(Counter(item["reason"] for item in invalid)),
        },
        "stages": {
            "stage12": {
                "total_weighted_score": detail12["aggregates"]["total_weighted_score"],
                "mean_structural_score": detail12["aggregates"]["structural_score"].get("mean"),
                "offtarget_factor_counts": detail12["aggregates"]["offtarget_factor_counts"],
            },
            "stage3": {
                "outcomes": detail3["aggregates"]["outcomes"],
                "cut_rate": detail3["aggregates"]["cut_rate"],
                "expected_no_cut": detail3["aggregates"]["expected_no_cut"],
                "conforming_rows": report["stage3"]["conforming_rows"],
                "replay_matches_pipeline": detail3["replay_matches_pipeline"],
            },
            "stage4": {
                "avg_r2": detail4.get("avg_r2"),
                "avg_nmae": detail4.get("avg_nmae"),
                "consistency_score": detail4.get("consistency_score"),
                "consistency_factor": detail4.get("consistency_factor"),
                "per_target": {target: {"r2_mean": stats["r2_mean"], "r2_std": stats["r2_std"],
                                        "mae_mean": stats["mae_mean"], "mae_std": stats["mae_std"],
                                        "normalized_mae": stats["normalized_mae"],
                                        "residual_std_mean": stats["residual_std_mean"]}
                               for target, stats in detail4.get("targets", {}).items()},
            },
            "stage5": {
                "distribution_fidelity_factor": detail5["distribution_fidelity_factor"],
                "terms": {name: term["value"] for name, term in detail5["terms"].items()},
                "empty_joint_cells": detail5["coverage"]["joint"]["empty_cells"],
            },
        },
        "score": {
            "total_weighted_score": report["stage4"]["total_weighted_score"],
            "consistency_score": report["stage4"]["consistency_score"],
            "consistency_factor": report["stage4"]["consistency_factor"],
            "distribution_fidelity_score": report["stage5"]["distribution_fidelity_score"],
            "distribution_fidelity_factor": report["stage5"]["distribution_fidelity_factor"],
            "final_score": report["final_score"],
        },
        "benchmark_submission": official_score,
        "replica_delta": delta,
        "at_compare_seed": None if alt_report is None else {
            "seed": compare_seed,
            "total_weighted_score": alt_report["stage4"]["total_weighted_score"],
            "consistency_factor": alt_report["stage4"]["consistency_factor"],
            "distribution_fidelity_factor": alt_report["stage5"]["distribution_fidelity_factor"],
            "final_score": alt_report["final_score"],
        },
        "artifacts": artifacts,
        "redirected_paths": mapping,
    }

    write_json(out_dir / "stage12_detail.json", detail12)
    write_json(out_dir / "stage3_detail.json", detail3)
    write_json(out_dir / "stage4_detail.json", detail4)
    write_json(out_dir / "stage5_detail.json", detail5)
    for name, payload in extra_detail.items():        # one set per round seed, when there are many
        write_json(out_dir / name, payload)
    write_json(out_dir / "validation.json", summary)
    write_json(stage_dir / "task.json", task)
    write_json(stage_dir / "cell_types.json", cell_types)
    write_json(stage_dir / "report.json", report)
    if alt_report is not None:
        write_json(stage_dir / "report_compare_seed.json", alt_report)
    write_json(stage_dir / "miner_score.json", official_score)
    with open(out_dir / "validation.txt", "w") as handle:
        handle.write(report_text + "\n" + stage4_text + "\n" + stage5_text + "\n"
                     + comparison_text + "\n")

    print("\n" + "=" * 100)
    print(f"  {submission_path.name}  ->  final_score {report['final_score']:.6f}  "
          f"=  {report['stage4']['total_weighted_score']:.4f} x "
          f"{report['stage4']['consistency_factor']:.6f} x "
          f"{report['stage5']['distribution_fidelity_factor']:.6f}")
    print(f"  {len(valid)} valid of {len(scored_rows)} scored rows at seed {seed}")
    print("=" * 100)

    print(f"\n  artifacts in {out_dir}/")
    for name in sorted(os.listdir(out_dir)):
        path = out_dir / name
        if path.is_dir():
            print(f"    {name + '/':<40}{'':>10}  {len(os.listdir(path))} files")
            continue
        print(f"    {name:<40}{path.stat().st_size:>10,} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
