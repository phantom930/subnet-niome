#!/usr/bin/env python3
"""run_task.py — build a submission for one task file and score it the way a validator does.

Two files come out, beside the task file by default:

    submission.json   the miner's rows — the bare JSON array it PUTs to the presigned S3 URL
    score.json        the validator's score, and nothing the miner claims about itself

The score has two independent sources, kept apart because they can disagree:

``validator_pipeline`` runs the validator's own code here — contract and reference persisted to
``data/``, rows dropped where the S3 download would put them, then stages 1/2 -> 3 -> 4 -> 5 over
those files via ``benchmark_submission``. It is what a validator holding *this* contract would
compute.

``recorded_by_validators`` is what the real validators actually published for this task to
``/api/v3/miners/scores`` — readable unsigned, so it is a score you can go and look at rather than
one this script derived. If the task file matches a task the backend has issued, every miner's
recorded breakdown for it is summarised, and ``--uid`` picks yours out.

The miner's own predicted score is deliberately absent. It is a self-report, it cannot see
truncation or a stage-1 rejection, and it is computed under whatever seed the miner was handed —
which in production is not the seed the validator scores with.

**Seeds.** Validation always uses the seed in the task file: it is written to ``data/contract.json``,
and both stage 3's outcome draws and stage 4's KFold shuffle read it from there. ``--build-seed``
sets the seed the rows are *generated* under, and defaults to the same value. Passing a different one
reproduces production, where the contract artifact arrives carrying ``seed: 0`` and is scored under
the stamped seed — build under 0, validate under the real seed, and the ``mh`` construction no longer
holds, which is the whole difference between a local score and a recorded one.

    python scripts/run_task.py                             # build + score testing/task.json
    python scripts/run_task.py --uid 188                   # also pull uid 188's recorded score
    python scripts/run_task.py --build-seed 0              # build as the miner really does, score honestly
    python scripts/run_task.py --from-submission rows.json # score an existing array, build nothing
    python scripts/run_task.py --offline                   # no backend calls at all

The task file may be either shape: the backend's ``{"id": ..., "content": {"contract": ...,
"hbb_reference": ...}}`` envelope, or the bare ``{"contract": ..., "hbb_reference": ...}`` pair.

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


def load_task(path: Path) -> dict:
    """Read a task file in either the backend envelope or the bare contract/reference shape."""
    document = json.load(open(path))
    if "content" in document:
        content = document["content"]
    elif "contract" in document and "hbb_reference" in document:
        content = document
    else:
        raise SystemExit(
            f"{path}: expected either a 'content' envelope or top-level "
            "'contract' and 'hbb_reference' keys"
        )
    for key in ("contract", "hbb_reference"):
        if key not in content:
            raise SystemExit(f"{path}: missing '{key}'")
    return {
        "id": document.get("id", path.stem),
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
    parser.add_argument("--submission-out", default=None,
                        help="the miner's rows (default: submission.json beside the task file)")
    parser.add_argument("--score-out", default=None,
                        help="the validator's score (default: score.json beside the task file)")
    parser.add_argument("--from-submission", default=None,
                        help="score this JSON array instead of building one (any miner's upload)")
    parser.add_argument("--build-seed", type=int, default=None,
                        help="seed to GENERATE the rows under (default: the task file's seed). "
                             "Validation always uses the task file's seed, so setting this to 0 "
                             "reproduces production, where the miner is handed an unstamped "
                             "contract and scored under the stamped one")
    parser.add_argument("--uid", type=int, default=None,
                        help="miner uid: stamped on the local MinerScore and used to pick your "
                             "record out of the published ones (default: from data/last_upload.json)")
    parser.add_argument("--cell-types", default=None,
                        help="JSON accessibility table to use instead of fetching it")
    parser.add_argument("--offline", action="store_true",
                        help="no backend calls: accessibility defaults to 1.0 and no published "
                             "scores are looked up")
    # Generation knobs; each defaults to the miner's own value. Ignored with --from-submission.
    parser.add_argument("--strategy", choices=("pure", "shaped"), default=Miner.STRATEGY)
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

    task = load_task(task_path)
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    max_experiments = contract.get("rules", {}).get("max_experiments")
    # The seed that scores this submission, and the only one the pipeline will ever see: it is
    # persisted to data/contract.json below, and stage 3 (outcome draws) and stage 4 (KFold shuffle)
    # both read it from there.
    validation_seed = contract.get("seed")
    build_seed = args.build_seed if args.build_seed is not None else validation_seed
    print(f"[1/5] task {task['id']} from {task_path}")
    print(f"  cell_type={contract.get('cell_type')} "
          f"mutations={len(contract.get('active_mutations', []))} "
          f"max_experiments={max_experiments} uid={uid}")
    print(f"  validation seed={validation_seed} (from the task file)  build seed={build_seed}")
    if not validation_seed:
        print("  ! the task file's seed is 0/absent, so validation runs under seed 0 — a real "
              "validator scores under the stamped seed, which this file does not carry")
    if build_seed != validation_seed:
        print(f"  ! building under {build_seed} and validating under {validation_seed}: the outcome "
              "construction cannot survive that, which is the production case")

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
        )
        print(f"[3/5] building with strategy={cfg.strategy} selection={cfg.selection} "
              f"construction={cfg.construction} variants={cfg.variants} flank={cfg.flank} "
              f"lengths={list(cfg.lengths)}")
        print("  loading data/chr11.fa and the k-mer index (first call is the slow one)")
        # Generate under build_seed. Only the copy handed to the generator is restamped; the task
        # persisted for the pipeline below keeps the task file's seed.
        build_task = copy.deepcopy(task)
        build_task["content"]["contract"]["seed"] = build_seed
        # score=False: the miner's own estimate is not wanted here, and skipping it saves ~5 s.
        meta, rows = S.build_for_task(build_task, cell_types, cfg, score=False)
        build = {
            "seed": build_seed,
            "rows": meta["rows"],
            "weight_skew": meta["weight_skew"],
            "construction": meta["construction"],
            "outcome_counts": meta["outcome_counts"],
            "problems": meta["problems"],
            "config": {k: (list(v) if isinstance(v, tuple) else v) for k, v in vars(cfg).items()},
        }
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

    print(f"[4/5] running the validator pipeline over {settings.MINER_SUBMISSION_PATH}")
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
    print(f"  LOCAL PIPELINE  final_score={score['final_score']:.6f} = "
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
                    "source": f"{settings.MINER_SCORE_URL}?task_id={backend_task['id']}",
                    "matched_by": matched["match"],
                    "task_id": backend_task["id"],
                    "task_created_at": backend_task.get("created_at"),
                    "seed_published_in_task_list": backend_task["content"]["contract"].get("seed"),
                    "field_summary": summary,
                    "for_uid": uid,
                    "records_for_uid": mine,
                }
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

    document = {
        "generated_by": "scripts/run_task.py",
        "task_file": str(task_path),
        "task_id_in_file": task["id"],
        "submission_file": str(rows_path),
        "seed_in_task_file": validation_seed,
        "seeds": {
            "validation": validation_seed,
            "build": build_seed,
            "note": "validation is the seed in the task file, written to data/contract.json and "
                    "read by stage 3 and stage 4; build is the seed the rows were generated under",
        },
        "cell_type": contract.get("cell_type"),
        "uid": uid,
        "elapsed_seconds": time.time() - started,
        # 1. The validator's own code, run here against this contract.
        "validator_pipeline": {
            "scored_by": "niome_subnet.genomics.validation.benchmark_submission "
                         "(stages 1/2->3->4->5)",
            "final_score": score["final_score"],
            "breakdown": breakdown,
            "rows_submitted": len(rows),
            "rows_scored": len(scored_rows),
            "rows_truncated": truncated,
            "rows_invalid": len(invalid),
            "invalid_reasons": dict(reasons),
            "stage4_model_results": stage4.get("model_results"),
            "stage5_summary": stage5,
        },
        # 2. What real validators published for this contract — a score you can go and read.
        "recorded_by_validators": recorded,
        # How the rows were made. Absent with --from-submission.
        "build": build,
    }
    score_path.parent.mkdir(parents=True, exist_ok=True)
    with open(score_path, "w") as handle:
        json.dump(document, handle, indent=2)
    print(f"\n  submission -> {rows_path}")
    print(f"  score      -> {score_path}")
    print(f"done in {time.time() - started:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
