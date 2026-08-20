#!/usr/bin/env python3
"""run_task.py — build a submission the way a miner does, score it the way a validator does.

The two halves of a round, in order, against one task file:

    1. the miner builds against the task **with no seed** — the contract it designs on carries
       ``seed: 0``, because that is what a broadcast task hands it
    2. the validator scores those rows **under the seed stamped on the contract**, and that score
       is the one that reaches the leaderboard

Two files come out, beside the task file by default:

    submission.json   the miner's rows — the bare JSON array it PUTs to the presigned S3 URL
    score.json        the validator's score, and nothing the miner claims about itself

``score.json`` is a summary, one screen long. ``score`` is the headline —  the validator's own code
run here, contract and reference persisted to ``data/``, rows dropped where the S3 download would put
them, then stages 1/2 -> 3 -> 4 -> 5 over those files via ``benchmark_submission`` — with
``breakdown`` for the three factors it multiplies, ``stage4_r2`` and ``stage5_ratios`` for where a
weak factor came from, and ``build`` for how the rows were made. Every stage's full output stays in
``data/`` (``final_reward.json``, ``distribution_fidelity_summary.json``,
``invalid_experiments.json``), so nothing summarised here is lost.

``recorded_by_validators`` is the second, independent source: what the real validators actually
published for this task to ``/api/v3/miners/scores`` — readable unsigned, so a score you can go and
look at rather than one this script derived. It carries the field (how many miners, best, median,
top-10 cutoff) and, with ``--uid``, that uid's own recorded score and its delta against the one
computed here. The two can disagree, which is the point of keeping both.

The miner's own predicted score is deliberately absent. It is a self-report, it cannot see
truncation or a stage-1 rejection, and it is computed under whatever seed the miner was handed —
which in production is not the seed the validator scores with.

**Seeds.** Validation always uses the seed in the task file: it is written to ``data/contract.json``,
and both stage 3's outcome draws and stage 4's KFold shuffle read it from there. Generation is
seedless by default — the rows are built against a copy of the contract stamped ``seed: 0``, the
placeholder a miner is actually handed — so the number printed at the end is the honest one, with no
flag to remember.

That gap is the whole design problem. The ``pure`` and ``shaped`` strategies engineer stage-3 outcomes
under the seed they are given and cannot survive it: their collapse *is* the difference between a
miner's predicted score and a recorded one. The default ``robust`` strategy never reads the seed, so
its rows are identical whatever is stamped here — but its score is not, because whether the run lands
on a seed where every row cuts is exactly what the strategy is playing for. **One run of this script
is one sample from that distribution, not the strategy's value.** Run several task files, or
``genExp.py --all-tasks --limit N``, to see the mix.

``--build-with-seed`` builds under the validation seed instead: an oracle run, an upper bound on what
a construction is worth while it holds, and not a score any miner receives. ``--build-seed N``
builds under some other stand-in. Neither changes anything under ``--strategy robust``.

    python scripts/run_task.py                             # build seedless, score under the task's seed
    python scripts/run_task.py --uid 188                   # also pull uid 188's recorded score
    python scripts/run_task.py --task-id e824bae7          # pick one task out of a task list
    python scripts/run_task.py --task-index -1             # the oldest task in that list
    python scripts/run_task.py --build-with-seed           # oracle: build under the scoring seed too
    python scripts/run_task.py --from-submission rows.json # score an existing array, build nothing
    python scripts/run_task.py --offline                   # no backend calls at all

The task file may hold one task or many. One task is either the backend's ``{"id": ...,
"content": {"contract": ..., "hbb_reference": ...}}`` envelope or the bare
``{"contract": ..., "hbb_reference": ...}`` pair. Many is a JSON array of those, or a dict with them
under ``tasks`` (what ``scripts/fetch_tasks.py`` writes) or ``items`` (the backend's own page shape) —
newest first, selected with ``--task-index`` or ``--task-id``. ``testing/submission.json`` and
``testing/score.json`` are overwritten every run whichever task is chosen, so use ``--score-out`` when
comparing tasks.

Like a validator, this **overwrites** ``data/contract.json``, ``data/hbb_reference.json``,
``data/submission.json`` and every stage artifact under ``data/`` — that is how the stages talk to
each other. ``data/`` is gitignored and holds nothing durable.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

# genExp chdir()s to the repo root on import (every settings.py path is relative to it), so resolve
# every CLI path before that happens or a relative --task stops meaning what the caller typed.
REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TASK = REPO_ROOT / "testing" / "task.json"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import requests  # noqa: E402

import niome_subnet.utils.settings as settings  # noqa: E402
import genExp as G  # noqa: E402
import submission as S  # noqa: E402

from niome_subnet.genomics.validation import benchmark_submission  # noqa: E402

# The miner's generation knobs, read off the class rather than re-typed, so the rows scored here are
# the rows the running miner would upload for this contract.
from neurons.miner import Miner  # noqa: E402

TASKS_URL = f"{settings.BASE_URL}/api/v3/tasks"

# The six ratios stage 5 takes the geometric mean of, under names short enough to scan. The stage's
# own diagnostics (per-cas shift, coverage counts) stay in data/distribution_fidelity_summary.json.
STAGE5_RATIOS = (
    ("mutation", "mutation_coverage_entropy_ratio"),
    ("cas_system", "cas_system_coverage_entropy_ratio"),
    ("strand", "strand_coverage_entropy_ratio"),
    ("joint", "joint_coverage_entropy_ratio"),
    ("kmer_diversity", "kmer_diversity_entropy_ratio"),
    ("distinct_guide", "distinct_guide_ratio"),
)


def readable(value: object, places: int = 6) -> object:
    """Round every float in a nested structure, so score.json can be read rather than parsed.

    Six places is finer than any comparison this script makes — the recorded-score delta is judged
    at 1e-6 — so nothing that matters is rounded away.
    """
    if isinstance(value, float):
        return round(value, places)
    if isinstance(value, dict):
        return {key: readable(item, places) for key, item in value.items()}
    if isinstance(value, list):
        return [readable(item, places) for item in value]
    return value


def select_task(path: Path, task_id: str | None, index: int) -> dict:
    """Pick one task out of a task file, which may hold one task or the whole history.

    Four shapes are accepted, because all four are things you end up with: the backend's own
    ``{"items": [...]}`` page, the snapshot ``scripts/fetch_tasks.py`` writes (``{"tasks": [...]}``
    with provenance), a bare JSON array of either, and a single task in envelope or bare form.
    """
    document = json.load(open(path))
    listing = None
    if isinstance(document, list):
        listing = document
    elif isinstance(document, dict):
        for key in ("tasks", "items"):
            if isinstance(document.get(key), list):
                listing = document[key]
                break
    if listing is None:
        return _as_task(document, path, path.stem)

    if not listing:
        raise SystemExit(f"{path}: holds an empty task list")
    if task_id:
        matches = [t for t in listing if str(t.get("id", "")).startswith(task_id)]
        if not matches:
            raise SystemExit(f"{path}: no task id starts with '{task_id}' ({len(listing)} present)")
        if len(matches) > 1:
            raise SystemExit(f"'{task_id}' matches {len(matches)} tasks in {path}; use more of the id")
        chosen = matches[0]
    else:
        if not -len(listing) <= index < len(listing):
            raise SystemExit(f"{path}: --task-index {index} out of range (0..{len(listing) - 1})")
        chosen = listing[index]
    position = listing.index(chosen)
    print(f"  {len(listing)} tasks in {path.name}; using #{position} "
          f"({chosen.get('created_at', '?')[:19]})"
          + ("" if task_id or index else ", the newest — --task-index N or --task-id <id> to pick"))
    return _as_task(chosen, path, f"{path.stem}#{position}")


def _as_task(document: dict, path: Path, fallback_id: str) -> dict:
    """Normalise one task in either the backend envelope or the bare contract/reference shape."""
    if "content" in document:
        content = document["content"]
    elif "contract" in document and "hbb_reference" in document:
        content = document
    else:
        raise SystemExit(
            f"{path}: expected either a 'content' envelope, top-level 'contract' and "
            "'hbb_reference' keys, or a list of tasks under 'tasks'/'items'"
        )
    for key in ("contract", "hbb_reference"):
        if key not in content:
            raise SystemExit(f"{path}: missing '{key}'")
    return {
        "id": document.get("id", fallback_id),
        "created_at": document.get("created_at"),
        "content": content,
    }


def read_artifact(path: str) -> object:
    """Load a stage artifact, tolerating one the run never produced."""
    file = REPO_ROOT / path
    if not file.exists():
        return None
    try:
        return json.load(open(file))
    except json.JSONDecodeError:
        return None


def uid_from_last_upload() -> int | None:
    """The uid the miner last uploaded under, recovered from the presigned key ``niome/<uid>.json``.

    Only a default for ``--uid``: the record is per-machine and says nothing about this task.
    """
    record = read_artifact(settings.LAST_UPLOAD_PATH)
    if not isinstance(record, dict):
        return None
    found = re.search(r"niome/(\d+)\.json", record.get("presigned_url", ""))
    return int(found.group(1)) if found else None


def match_backend_task(contract: dict, reference: dict) -> dict | None:
    """Find the backend task this file came from, so its published scores can be looked up.

    Matched on the contract itself rather than on any id in the file, because a hand-written task
    file may carry no id at all. Falls back to (seed, active_mutations) if the contract has drifted.
    """
    items = requests.get(TASKS_URL, timeout=120).json().get("items", [])
    exact = [i for i in items if i["content"]["contract"] == contract]
    if exact:
        return {"task": exact[0], "match": "contract"}
    loose = [
        i for i in items
        if i["content"]["contract"].get("seed") == contract.get("seed")
        and i["content"]["contract"].get("active_mutations") == contract.get("active_mutations")
    ]
    if loose:
        return {"task": loose[0], "match": "seed+mutations"}
    return None


def fetch_recorded_scores(task_id: str) -> list[dict]:
    """Every score real validators published for one task. Unsigned read, server-side filtered."""
    response = requests.get(f"{settings.MINER_SCORE_URL}?task_id={task_id}", timeout=120)
    response.raise_for_status()
    payload = response.json()
    items = payload["items"] if isinstance(payload, dict) else payload
    return [i for i in items if i.get("task_id") == task_id]


def summarise_recorded(records: list[dict]) -> dict:
    """The shape of the field for one task: what a score was actually worth against real rivals."""
    finals = sorted((r["final_score"] for r in records), reverse=True)
    scored = [f for f in finals if f > 0]
    consistencies = [r["breakdown"]["consistency_factor"] for r in records]
    return {
        "miners_scored": len(records),
        "final_score_max": finals[0] if finals else None,
        "final_score_median": statistics.median(finals) if finals else None,
        "final_score_zero_count": len(finals) - len(scored),
        "top10_cutoff": finals[min(9, len(finals) - 1)] if finals else None,
        "consistency_factor_max": max(consistencies) if consistencies else None,
        "consistency_factor_median": statistics.median(consistencies) if consistencies else None,
        "miners_at_full_consistency": sum(1 for c in consistencies if c >= 0.9999),
        "validator_uids": sorted({r.get("validator_uid") for r in records
                                  if r.get("validator_uid") is not None}),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task", default=str(DEFAULT_TASK))
    parser.add_argument("--task-index", type=int, default=0,
                        help="which task to use when --task holds a list (0 = newest, negative "
                             "indexes from the oldest)")
    parser.add_argument("--task-id", default=None,
                        help="pick a task out of that list by id, or any unique prefix of one")
    parser.add_argument("--submission-out", default=None,
                        help="the miner's rows (default: submission.json beside the task file)")
    parser.add_argument("--score-out", default=None,
                        help="the validator's score (default: score.json beside the task file)")
    parser.add_argument("--from-submission", default=None,
                        help="score this JSON array instead of building one (any miner's upload)")
    # The seed the rows are GENERATED under. Validation is never affected by either of these: it
    # always uses the seed in the task file.
    seeding = parser.add_mutually_exclusive_group()
    seeding.add_argument("--build-seed", type=int, default=0,
                         help="stand-in seed to GENERATE the rows under (default: 0, the "
                              "placeholder a miner is handed — it never sees the real one)")
    seeding.add_argument("--build-with-seed", action="store_true",
                         help="generate under the task file's seed instead of blind: an oracle run "
                              "that no miner gets, useful as an upper bound")
    parser.add_argument("--uid", type=int, default=None,
                        help="miner uid: stamped on the local MinerScore and used to pick your "
                             "record out of the published ones (default: from data/last_upload.json)")
    parser.add_argument("--cell-types", default=None,
                        help="JSON accessibility table to use instead of fetching it")
    parser.add_argument("--offline", action="store_true",
                        help="no backend calls: accessibility defaults to 1.0 and no published "
                             "scores are looked up")
    # Generation knobs; each defaults to the miner's own value. Ignored with --from-submission.
    parser.add_argument("--strategy", choices=("robust", "pure", "shaped"), default=Miner.STRATEGY)
    parser.add_argument("--selection", choices=("packed", "stratified"), default=Miner.SELECTION)
    parser.add_argument("--construction", choices=tuple(G.CONSTRUCTIONS), default=Miner.CONSTRUCTION)
    parser.add_argument("--variants", type=int, default=Miner.VARIANTS)
    parser.add_argument("--flank", type=int, default=Miner.FLANK)
    parser.add_argument("--lengths", default=",".join(str(v) for v in Miner.LENGTHS))
    parser.add_argument("--rows", type=int, default=None, help="override contract max_experiments")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    task_path = Path(args.task).resolve()
    rows_path = (Path(args.submission_out).resolve() if args.submission_out
                 else task_path.with_name("submission.json"))
    score_path = (Path(args.score_out).resolve() if args.score_out
                  else task_path.with_name("score.json"))
    source_path = Path(args.from_submission).resolve() if args.from_submission else None
    cell_types_path = Path(args.cell_types).resolve() if args.cell_types else None
    uid = args.uid if args.uid is not None else uid_from_last_upload()
    started = time.time()

    print(f"[1/5] reading {task_path}")
    task = select_task(task_path, args.task_id, args.task_index)
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    max_experiments = contract.get("rules", {}).get("max_experiments")
    # The seed that scores this submission, and the only one the pipeline will ever see: it is
    # persisted to data/contract.json below, and stage 3 (outcome draws) and stage 4 (KFold shuffle)
    # both read it from there.
    validation_seed = contract.get("seed")
    # The seed the rows are designed against. Blind by default — a broadcast task carries no seed,
    # so the miner builds against 0 and only the validator ever holds the real one.
    build_seed = validation_seed if args.build_with_seed else args.build_seed
    print(f"  task {task['id']} created {task.get('created_at')}")
    print(f"  cell_type={contract.get('cell_type')} "
          f"mutations={len(contract.get('active_mutations', []))} "
          f"max_experiments={max_experiments} uid={uid}")
    print(f"  validation seed={validation_seed} (from the task file, what the score is computed "
          f"under)")
    if not validation_seed:
        print("  ! the task file's seed is 0/absent, so validation runs under seed 0 — a real "
              "validator scores under the stamped seed, which this file does not carry")
    if source_path:
        pass  # nothing is generated, so no build seed applies
    elif build_seed != validation_seed:
        print(f"  build seed={build_seed} — the miner designs blind and is scored under "
              f"{validation_seed}: the outcome construction cannot survive that, and this is the "
              "production case")
    else:
        print(f"  build seed={build_seed} — same as the validation seed. This is an oracle run: no "
              "miner is handed the seed it will be scored under")

    if cell_types_path:
        cell_types = json.load(open(cell_types_path))
        print(f"[2/5] cell types from {cell_types_path} ({len(cell_types)} entries)")
    elif args.offline:
        cell_types = {}
        print("[2/5] offline: accessibility defaults to 1.0")
    else:
        # A validator reads this endpoint signed with its hotkey (api.fetch_cell_types); this reads
        # it unsigned. Same URL and same table today.
        cell_types = G.fetch_cell_types()
        print(f"[2/5] cell types from the backend, unsigned ({len(cell_types)} entries)")
    if contract.get("cell_type") and contract["cell_type"] not in (cell_types or {}):
        print(f"  ! '{contract.get('cell_type')}' is absent from the table: stage 2 falls back to "
              "accessibility 1.0, which a real validator may not do")

    build = None
    if source_path:
        rows = json.load(open(source_path))
        print(f"[3/5] scoring {len(rows)} rows from {source_path} (nothing built)")
    else:
        cfg = G.GenConfig(
            strategy=args.strategy,
            selection=args.selection,
            construction=args.construction,
            variants=args.variants,
            flank=args.flank,
            lengths=tuple(int(v) for v in args.lengths.split(",")),
            rows=args.rows,
            sites_per_cell=Miner.SITES_PER_CELL,
            cas12a_share=Miner.CAS12A_SHARE,
            variant_pool=Miner.VARIANT_POOL,
            refine_passes=Miner.REFINE_PASSES,
        )
        seeding = "under seed " + str(build_seed) if build_seed else "seedless"
        print(f"[3/5] MINER: building {seeding} with strategy={cfg.strategy} "
              f"selection={cfg.selection} construction={cfg.construction} "
              f"variants={cfg.variants} flank={cfg.flank} lengths={list(cfg.lengths)}")
        print("  loading data/chr11.fa and the k-mer index (first call is the slow one)")
        # The task as a miner receives it: same rules, mutations and weights, seed replaced by the
        # stand-in. genExp reads Context.seed off the contract unconditionally, so the field has to
        # exist — 0 is the placeholder value a real broadcast carries. hbb_reference echoes the
        # contract under 'challenge'; nothing reads its seed today, but it is blanked with the
        # contract's so no path can recover the real one from the copy the generator holds.
        # Only this copy is restamped: the task persisted for the pipeline below keeps the task
        # file's seed, which is what the validator scores under.
        build_task = copy.deepcopy(task)
        build_task["content"]["contract"]["seed"] = build_seed
        challenge = build_task["content"]["hbb_reference"].get("challenge")
        if isinstance(challenge, dict) and "seed" in challenge:
            challenge["seed"] = build_seed
        # score=False: the miner's own estimate is not wanted here, and skipping it saves ~5 s.
        meta, rows = S.build_for_task(build_task, cell_types, cfg, score=False)
        build = {
            "blind": build_seed != validation_seed,
            "rows": meta["rows"],
            "weight_skew": meta["weight_skew"],
            "outcome_counts": meta["outcome_counts"],
            # Only the knobs this script can move. The rest of GenConfig is genExp's scoring
            # constants, which are the same on every run and readable there.
            "config": {"strategy": cfg.strategy, "selection": cfg.selection,
                       "construction": cfg.construction, "variants": cfg.variants,
                       "flank": cfg.flank, "lengths": list(cfg.lengths), "rows": cfg.rows},
        }
        if meta["problems"]:
            build["problems"] = meta["problems"]
        print(f"  built {meta['rows']} rows, weight_skew={meta['weight_skew']}, "
              f"outcomes={meta['outcome_counts']}")
        for problem in meta["problems"]:
            print(f"  ! generator invariant violated: {problem}")

    # The miner's artifact: the bare array it would PUT, before any validator touches it.
    rows_path.parent.mkdir(parents=True, exist_ok=True)
    with open(rows_path, "w") as handle:
        json.dump(rows, handle, indent=2)
    print(f"  wrote {len(rows)} rows to {rows_path}")

    # From here on this is the validator's sequence: persist the contract and reference the way
    # api.fetch_task does — carrying the task file's seed, whatever the rows were built under — drop
    # the rows where the S3 download would put them, and benchmark.
    G.persist_task(task)
    with open(settings.MINER_SUBMISSION_PATH, "w") as handle:
        json.dump(rows, handle, indent=2)

    print(f"[4/5] VALIDATOR: scoring {settings.MINER_SUBMISSION_PATH} under seed {validation_seed}")
    score = benchmark_submission(cell_types, uid=uid if uid is not None else 0).model_dump()
    breakdown = score["breakdown"]

    # What the pipeline actually saw, read back from its own artifacts rather than assumed.
    scored_rows = read_artifact(settings.MINER_SUBMISSION_PATH) or []
    invalid = read_artifact(settings.INVALID_EXPERIMENTS_PATH) or []
    stage4 = read_artifact(settings.FINAL_REWARD_PATH) or {}
    stage5 = read_artifact(settings.DISTRIBUTION_FIDELITY_PATH) or {}
    reasons = Counter(item.get("reason") for item in invalid)

    truncated = len(rows) - len(scored_rows)
    if truncated:
        print(f"  ! truncate_submission cut {truncated} of {len(rows)} rows before scoring "
              f"(cap {max_experiments}, unique non-blank experiment_ids)")
    if invalid:
        print(f"  ! stage 1 rejected {len(invalid)} row(s): {dict(reasons)}")
    print(f"  stage 1/2  {breakdown['n_valid_experiments']} valid of {len(scored_rows)} scored, "
          f"total_weighted_score={breakdown['total_weighted_score']:.3f}")
    for target, stats in (stage4.get("model_results") or {}).items():
        print(f"  stage 4    {target:<13} r2={stats['r2_mean']:+.6f} mae={stats['mae_mean']:.3e} "
              f"folds={stats['n_folds']}")
    print(f"  stage 4    consistency_score={breakdown['consistency_score']:.4f} "
          f"-> factor {breakdown['consistency_factor']:.4f}")
    print(f"  stage 5    fidelity_factor={breakdown['distribution_fidelity_factor']:.4f}")
    print(f"  LEADERBOARD SCORE  final_score={score['final_score']:.6f} = "
          f"{breakdown['total_weighted_score']:.3f} x {breakdown['consistency_factor']:.4f} "
          f"x {breakdown['distribution_fidelity_factor']:.4f}")

    recorded: dict | None = None
    if args.offline:
        print("[5/5] offline: skipping the published-score lookup")
    else:
        print(f"[5/5] looking up what real validators published for this contract")
        try:
            matched = match_backend_task(contract, reference)
            if matched is None:
                print("  no backend task carries this contract — nothing was ever scored for it")
            else:
                backend_task = matched["task"]
                records = fetch_recorded_scores(backend_task["id"])
                summary = summarise_recorded(records)
                mine = [r for r in records if r.get("miner_uid") == uid] if uid is not None else []
                recorded = {
                    "task_id": backend_task["id"],
                    "matched_by": matched["match"],
                    "source": f"{settings.MINER_SCORE_URL}?task_id={backend_task['id']}",
                    "miners_scored": summary["miners_scored"],
                    "best": summary["final_score_max"],
                    "median": summary["final_score_median"],
                    "top10_cutoff": summary["top10_cutoff"],
                    "scored_zero": summary["final_score_zero_count"],
                }
                if mine:
                    # One record per validator; there is a single active validator today, and a
                    # second would score the same rows under the same contract.
                    record = mine[0]
                    recorded |= {
                        "uid": uid,
                        "final_score": record["final_score"],
                        "delta_vs_local": score["final_score"] - record["final_score"],
                        "consistency_factor": record["breakdown"]["consistency_factor"],
                        "weight": record.get("weight"),
                        "at": (record.get("created_at") or "")[:19],
                    }
                    if len(mine) > 1:
                        recorded["other_records_for_uid"] = len(mine) - 1
                print(f"  task {backend_task['id']} (created {backend_task.get('created_at')})")
                print(f"  {summary['miners_scored']} miners scored by validators "
                      f"{summary['validator_uids']} | best final={summary['final_score_max']:.3f} "
                      f"median={summary['final_score_median']:.3f} "
                      f"top-10 cutoff={summary['top10_cutoff']:.3f} "
                      f"| {summary['miners_at_full_consistency']} at consistency 1.0")
                for record in mine:
                    b = record["breakdown"]
                    print(f"  REAL VALIDATOR  uid={record['miner_uid']} "
                          f"validator={record.get('validator_hotkey', '?')[:8]} "
                          f"at {record.get('created_at', '')[:19]}")
                    print(f"     final_score={record['final_score']:.6f} = "
                          f"{b['total_weighted_score']:.3f} x {b['consistency_factor']:.4f} "
                          f"x {b['distribution_fidelity_factor']:.4f}  "
                          f"n_valid={b['n_valid_experiments']}  weight={record.get('weight')}")
                    delta = score["final_score"] - record["final_score"]
                    print(f"     local pipeline says {score['final_score']:.6f} "
                          f"(delta {delta:+.6f})")
                    if abs(delta) > 1e-6:
                        seed_note = ""
                        if abs(b["total_weighted_score"] - breakdown["total_weighted_score"]) < 1e-6:
                            # Stages 1, 2 and 5 never read the seed; stage 4 is downstream of it. So
                            # matching weights with a differing consistency isolates the difference
                            # to the seed the rows were built against.
                            seed_note = (" — total_weighted_score matches exactly, so the same rows "
                                         "were scored and only stage 4 differs: the seed the rows "
                                         "were built under is not the seed they were scored under")
                            print(f"     ! consistency {b['consistency_factor']:.4f} recorded vs "
                                  f"{breakdown['consistency_factor']:.4f} here{seed_note}")
                if uid is not None and not mine:
                    print(f"  no published record for uid {uid} on this task")
        except Exception as error:  # noqa: BLE001 - a lookup failure must not lose the local score
            print(f"  ! published-score lookup failed ({type(error).__name__}: {error})")

    # A summary, not an archive: every stage's full output already sits in data/ (final_reward.json,
    # distribution_fidelity_summary.json, invalid_experiments.json), so anything dropped here is one
    # file away. Empty diagnostics are omitted rather than written as zeros, so what remains in the
    # file is what actually happened.
    document = {
        "task": task["id"],
        "cell_type": contract.get("cell_type"),
        "uid": uid,
        "seed": {"build": None if build is None else build_seed, "validation": validation_seed},
        # The one number this script exists to produce: what a validator holding this contract pays
        # for these rows.
        "score": score["final_score"],
        "breakdown": {
            "total_weighted_score": breakdown["total_weighted_score"],
            "consistency_factor": breakdown["consistency_factor"],
            "distribution_fidelity_factor": breakdown["distribution_fidelity_factor"],
        },
        "rows": {"submitted": len(rows), "scored": len(scored_rows),
                 "valid": breakdown["n_valid_experiments"]},
        # Per-target cross-validated R²: the whole of consistency_factor, and the first place a
        # blind build shows up (negative R² where a seed-matched one reads 1.0).
        "stage4_r2": {target: stats["r2_mean"]
                      for target, stats in (stage4.get("model_results") or {}).items()},
        "stage5_ratios": {short: stage5[key] for short, key in STAGE5_RATIOS if key in stage5},
        "recorded_by_validators": recorded,
        "build": build,
        "files": {"task": str(task_path), "submission": str(rows_path)},
        "elapsed_seconds": time.time() - started,
    }
    if truncated:
        document["rows"]["truncated"] = truncated
    if invalid:
        document["rows"]["invalid"] = len(invalid)
        document["rows"]["invalid_reasons"] = dict(reasons)

    score_path.parent.mkdir(parents=True, exist_ok=True)
    with open(score_path, "w") as handle:
        json.dump(readable(document), handle, indent=2)
    print(f"\n  submission -> {rows_path}")
    print(f"  score      -> {score_path}")
    print(f"done in {time.time() - started:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
