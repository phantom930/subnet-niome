#!/usr/bin/env python3
"""test.py — build a submission at one seed, score it at another, report all five stages.

Reads a task file (``task.json`` by default), builds the dataset through **the miner's own
configuration and code path** against a stand-in seed, then runs the validator's five stages against
the task's real seed. That models the unstamped-contract case: the backend broadcasts a task with
``seed: 0`` and assigns a real seed afterwards, so everything the generator engineered against 0 is
judged under a stream it never saw.

Stages 1, 2 and 5 never read the seed and must come out identical either way; stage 3 is keyed on it,
so the difference between the two reports is exactly what the seed was worth.

    python test.py                                 # task.json, build seed 0, into test/
    python test.py --task other.json --out-dir tmp
    python test.py --build-seed 323                # matched seeds; nothing to compare

Every artifact is written under ``--out-dir`` and **nothing under data/ is touched** except the
read-only chr11 FASTA and the shared k-mer cache — so this is safe to run beside a live miner.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import io
import json
import logging
import os
import random
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import niome_subnet.utils.settings as settings  # noqa: E402

# genExp chdir()s to the repo root on import, which is what every relative settings path assumes.
import genExp as G  # noqa: E402

from niome_subnet.genomics.validation import stage12, stage3, stage4, stage5  # noqa: E402


# Redirected into the output directory. CHR11_PATH and KMER_CACHE_DIR are deliberately absent: the
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


def redirect_paths(out_dir: Path) -> dict[str, str]:
    """Point every writable data path at ``out_dir``.

    The stage modules bind their paths with ``from ... import NAME`` at import time, so patching
    ``settings`` alone would not reach them — each module's own attribute has to be rebound as well.
    ``neurons.miner`` reads ``settings.NAME`` at call time, so for that one patching settings is
    enough. Doing both is what keeps a live miner's data/ directory untouched.
    """
    mapping: dict[str, str] = {}
    for name in REDIRECTED:
        original = getattr(settings, name, None)
        if original is None:
            continue
        target = str(out_dir / Path(original).name)
        setattr(settings, name, target)
        for module in (stage12, stage3, stage4, stage5):
            if hasattr(module, name):
                setattr(module, name, target)
        mapping[name] = target
    return mapping


def no_cut_rows(valid: list[dict], results: list[dict], ctx: G.Context) -> list[dict]:
    """Every row stage 3 declined to cut, with the draw that decided it.

    The RNG stream is replayed rather than inferred: ``stage3.simulate`` seeds
    ``random.Random(experiment_seed(round_seed, design))`` and then draws exactly twice before the
    outcome branches — first the microhomology coin, then the cut coin. Replaying those two gives
    the actual uniform that landed above ``cut_p``, so each row reports how far over it fell rather
    than only that it did.

    These rows are *valid*: they pass stage 1, carry a full stage-2 score and still contribute to
    total_weighted_score. What they cost is the construction — is_cut goes to 0 for them, which is
    what stage 4 then cannot recover.
    """
    detail = []
    for entry, result in zip(valid, results):
        if result["outcome"] != "no_cut":
            continue
        experiment = entry["experiment"]
        features = entry["features"]
        seed = stage3.experiment_seed(ctx.seed, experiment)
        rng = random.Random(seed)
        gc = features["gc"]
        p_mh = min(0.6, 2.2 * gc * (1 - gc))
        mh_draw = rng.random()          # draw 1: microhomology coin
        cut_draw = rng.random()         # draw 2: the cut coin, and the reason this row is here
        cut_p = stage3.cut_probability(experiment["cas_system"], result["energy"])
        detail.append({
            "experiment_id": experiment["experiment_id"],
            "design": dict(experiment),
            "stage2": {
                "gc": gc,
                "gc_score": features["gc_score"],
                "distance_to_mutation": features["distance_to_mutation"],
                "dist_score": features["dist_score"],
                "consistency": features["consistency"],
                "offtarget_factor": features["offtarget_factor"],
                "mutation_weight": features["mutation_weight"],
                "cell_type_accessibility": features["cell_type_accessibility"],
                "mutation_region": features["mutation_region"],
                "region_energy_offset": features["region_energy_offset"],
                "structural_score": entry["stage2"]["structural_score"],
                "weighted_score": entry["stage2"]["weighted_score"],
            },
            "stage3": {
                "experiment_seed": seed,
                "energy": result["energy"],
                "p_mh": p_mh,
                "mh": result["mh"],
                "mh_draw": mh_draw,
                "cut_p": cut_p,
                "cut_draw": cut_draw,
                "margin_over_cut_p": cut_draw - cut_p,
                "outcome": result["outcome"],
                "indel_length": result["indel_length"],
            },
        })
    return detail


def no_cut_summary(rows: list[dict], valid: list[dict], results: list[dict], ctx: G.Context,
                   report: dict, detail: list[dict]) -> dict:
    """Aggregate view of the no_cut rows, next to how many were expected."""
    nonconforming = len(rows) - (report["stage3"]["conforming_rows"] or 0)
    return {
        "seed": ctx.seed,
        "rows": len(rows),
        "no_cut": len(detail),
        "expected_no_cut": report["stage3"]["expected_no_cut"],
        "cut_rate": report["stage3"]["cut_rate"],
        "by_cas": dict(Counter(d["design"]["cas_system"] for d in detail)),
        "by_strand": dict(Counter(d["design"]["strand"] for d in detail)),
        "by_mutation": dict(Counter(d["design"]["mutation"] for d in detail)),
        "by_cut_p": dict(Counter(round(d["stage3"]["cut_p"], 4) for d in detail)),
        # no_cut rows keep their stage-2 score; only the construction is lost.
        "weighted_score_carried": sum(d["stage2"]["weighted_score"] for d in detail),
        "weighted_score_total": report["stage2"]["total_weighted_score"],
        "nonconforming_rows": nonconforming,
        "no_cut_share_of_nonconforming": (len(detail) / nonconforming) if nonconforming else 0.0,
        "closest_call": min((d["stage3"]["margin_over_cut_p"] for d in detail), default=None),
        "widest_miss": max((d["stage3"]["margin_over_cut_p"] for d in detail), default=None),
    }


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task", default="task.json",
                        help="task file holding content.contract and content.hbb_reference")
    parser.add_argument("--out-dir", default="test", help="where every artifact is written")
    parser.add_argument("--build-seed", type=int, default=0,
                        help="seed the submission is generated against (default 0)")
    parser.add_argument("--cell-types", default=None,
                        help="cell-type accessibility JSON; fetched from the backend if omitted")
    parser.add_argument("--quiet", action="store_true", help="suppress the miner's own log lines")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                        format="%(levelname)s %(message)s")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mapping = redirect_paths(out_dir)

    # Imported after the redirect so the module picks up the patched settings object.
    from neurons.miner import Miner  # noqa: E402

    task = json.load(open(args.task))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    real_seed = contract.get("seed")

    print("=" * 100)
    print(f"  task {task.get('id')}   created {task.get('created_at')}")
    print(f"  real contract seed {real_seed}   cell_type {contract.get('cell_type')}")
    print(f"  mutations {contract['active_mutations']}")
    print(f"  weights   {contract.get('mutation_weights')}")
    print(f"  rules     {contract['rules']}")
    print(f"  build seed {args.build_seed}   output {out_dir}/")
    print("=" * 100)

    if args.cell_types:
        cell_types = json.load(open(args.cell_types))
    else:
        cell_types = G.fetch_cell_types()
    accessibility = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    print(f"  accessibility {accessibility}")

    # ---------------------------------------------------------------- build, at the stand-in seed
    build_contract = copy.deepcopy(contract)
    build_contract["seed"] = args.build_seed

    miner = Miner.__new__(Miner)          # no wallet, no chain — only _build is exercised
    miner.gen_config = Miner.base_gen_config()
    chosen = Miner.settings_for(contract)
    override = Miner.CELL_TYPE_OVERRIDES.get(contract.get("cell_type"))
    print(f"\n  miner config: " + " ".join(f"{k}={chosen[k]}" for k in Miner.TUNABLE))
    if override:
        print(f"  (cell type {contract.get('cell_type')} overrides {sorted(override)})")

    print(f"\n[1/4] building against seed {args.build_seed} with neurons/miner.py")
    rows = miner._build(build_contract, reference, cell_types)   # writes MINER_SUBMISSION_PATH
    if not rows:
        print("no rows generated", file=sys.stderr)
        return 1

    cfg = Miner.gen_config_for(build_contract)

    # ------------------------------------------------------- report as built, at the stand-in seed
    print(f"\n[2/4] measuring all five stages at the build seed {args.build_seed}")
    build_ctx = G.build_context(build_contract, reference, cell_types)
    _report, build_valid, build_results = G.rescore_under(rows, build_ctx)
    as_built = G.stage_report(rows, build_valid, build_results, build_ctx, cfg)
    built_text = captured(as_built, f"AS BUILT — seed {args.build_seed}")
    print(built_text)

    # -------------------------------------------------------- report as scored, at the real seed
    print(f"[3/4] measuring all five stages at the real contract seed {real_seed}")
    real_ctx = G.build_context(contract, reference, cell_types)
    _report, real_valid, real_results = G.rescore_under(rows, real_ctx)
    as_scored = G.stage_report(rows, real_valid, real_results, real_ctx, cfg)
    scored_text = captured(as_scored, f"AS SCORED — real contract seed {real_seed}")
    print(scored_text)

    # ------------------------------------------- the validator's own file pipeline, for the record
    print(f"[4/4] running the validator's five stages over {out_dir}/ at the real seed")
    write_json(out_dir / "contract.json", contract)             # real seed — what the stages read
    write_json(out_dir / "hbb_reference.json", reference)
    official = G.verify_with_validator(real_ctx)
    print(f"  benchmark_submission final_score = {official['final_score']:.6f}")
    print(f"  in-memory replica                = {as_scored['final_score']:.6f}")
    print(f"  delta                            = "
          f"{abs(official['final_score'] - as_scored['final_score']):.9f}")

    # ------------------------------------------------------------------- the rows that did not cut
    built_no_cut = no_cut_rows(build_valid, build_results, build_ctx)
    scored_no_cut = no_cut_rows(real_valid, real_results, real_ctx)
    no_cut_document = {
        "task_id": task.get("id"),
        "cell_type": contract.get("cell_type"),
        "accessibility": accessibility,
        "build_seed": args.build_seed,
        "real_seed": real_seed,
        "note": ("no_cut rows stay valid and keep their stage-2 weighted_score; what they cost is "
                 "the construction, because is_cut drops to 0 and stage 4 cannot recover it"),
        "summary": {
            "as_built": no_cut_summary(rows, build_valid, build_results, build_ctx,
                                       as_built, built_no_cut),
            "as_scored": no_cut_summary(rows, real_valid, real_results, real_ctx,
                                        as_scored, scored_no_cut),
        },
        "experiments": {"as_built": built_no_cut, "as_scored": scored_no_cut},
    }
    write_json(out_dir / "no_cut.json", no_cut_document)

    print(f"\n  no_cut rows: {len(built_no_cut)} at seed {args.build_seed}, "
          f"{len(scored_no_cut)} at seed {real_seed} "
          f"(expected {as_scored['stage3']['expected_no_cut']:.2f} at either) "
          f"-> {out_dir}/no_cut.json")
    for record in scored_no_cut:
        stage3_detail = record["stage3"]
        print(f"    {record['experiment_id']}  {record['design']['cas_system']:<7}"
              f"{record['design']['strand']}  d={record['stage2']['distance_to_mutation']:>4}  "
              f"cut_p={stage3_detail['cut_p']:.4f}  draw={stage3_detail['cut_draw']:.6f}  "
              f"over by {stage3_detail['margin_over_cut_p']:+.6f}  "
              f"weighted={record['stage2']['weighted_score']:.4f}")

    # ------------------------------------------------------------------------------- the verdict
    lines = ["", "=" * 100, "  WHAT THE SEED COST", "=" * 100]
    for label, section, key in (
        ("total_weighted_score", "stage4", "total_weighted_score"),
        ("consistency_factor", "stage4", "consistency_factor"),
        ("distribution_fidelity_factor", "stage5", "distribution_fidelity_factor"),
    ):
        before, after = as_built[section][key], as_scored[section][key]
        moved = "unchanged" if abs(before - after) < 1e-9 else f"{after - before:+.6f}"
        lines.append(f"    {label:<30} {before:>12.6f}  ->  {after:>12.6f}   ({moved})")
    before, after = as_built["final_score"], as_scored["final_score"]
    ratio = after / before if before else 0.0
    lines.append(f"    {'final_score':<30} {before:>12.6f}  ->  {after:>12.6f}   "
                 f"({ratio:.4f}x, {100 * (ratio - 1):+.1f}%)")
    lines.append(f"    {'construction conformance':<30} "
                 f"{as_built['stage3']['conforming_rows']:>12} / {len(rows)}  ->  "
                 f"{as_scored['stage3']['conforming_rows']:>6} / {len(rows)}")
    verdict = "\n".join(lines)
    print(verdict)

    # --------------------------------------------------------------------------------- artifacts
    write_json(out_dir / "task.json", task)
    write_json(out_dir / "contract_build_seed.json", build_contract)
    write_json(out_dir / "cell_types.json", cell_types)
    write_json(out_dir / "report_as_built.json", as_built)
    write_json(out_dir / "report_as_scored.json", as_scored)
    write_json(out_dir / "miner_score.json", official)
    write_json(out_dir / "config.json", {
        "task_id": task.get("id"),
        "build_seed": args.build_seed,
        "real_seed": real_seed,
        "accessibility": accessibility,
        "miner_constants": {k: (list(v) if isinstance(v, tuple) else v)
                            for k, v in chosen.items()},
        "cell_type_override": override,
        "resolved_config": {k: (list(v) if isinstance(v, tuple) else v)
                            for k, v in vars(cfg).items()},
        "redirected_paths": mapping,
    })
    with open(out_dir / "report.txt", "w") as handle:
        handle.write(built_text + "\n" + scored_text + "\n" + verdict + "\n")

    print(f"\n  artifacts in {out_dir}/")
    for name in sorted(os.listdir(out_dir)):
        size = (out_dir / name).stat().st_size
        print(f"    {name:<40}{size:>10,} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
