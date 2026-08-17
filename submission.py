#!/usr/bin/env python3
"""submission.py — build the actual experiment datasets for every task.

``genExp.py --all-tasks`` reports what each task *would* score; this writes the rows themselves.
Every entry is a JSON array in the exact submission format stage 1 reads — designs only, no
predicted efficiency, repair mode, indel length or off-target value, because outcome computation
belongs to the validator.

    python submission.py                          # all tasks -> submission.json
    python submission.py --limit 5                # newest 5 only
    python submission.py --task-id <uuid>         # one task
    python submission.py --per-task-dir data/submissions
    python submission.py --no-score               # skip scoring (about 40% faster)

``submission.json`` holds every task keyed by id. That is a convenient archive but **not** an
uploadable artifact: a miner PUTs one bare array per task. ``--per-task-dir`` writes that form,
one ``<task_id>.json`` per task, ready to send as-is.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from dataclasses import replace
from pathlib import Path

import genExp as G


def build_for_task(task: dict, cell_types: dict, cfg: G.GenConfig,
                   score: bool) -> tuple[dict, list[dict]]:
    """Generate one task's submission. Returns (metadata, rows)."""
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, cfg.flank, tuple(sorted(set(cfg.lengths))))

    task_cfg = cfg
    if cfg.strategy == "pure":
        task_cfg = replace(cfg, weight_skew=G.choose_weight_skew(ctx, sites, cfg))

    rows, valid, results = G.generate(ctx, sites, task_cfg)
    rows = G.order_rows(rows, valid)

    meta = {
        "task_id": task["id"],
        "created_at": task.get("created_at"),
        "seed": contract.get("seed"),
        # Unstamped task: stage 3 is keyed on the seed, so the construction stops holding once the
        # backend assigns a real one. The rows remain valid; the expected score does not.
        "seed_provisional": not contract.get("seed"),
        "cell_type": contract.get("cell_type"),
        "active_mutations": contract["active_mutations"],
        "mutation_weights": contract.get("mutation_weights", {}),
        "max_experiments": ctx.max_experiments,
        "weight_skew": task_cfg.weight_skew,
        "rows": len(rows),
        "construction": task_cfg.construction,
        "outcome_counts": dict(Counter(r["outcome"] for r in results)),
        "problems": G.check_invariants(rows, results, task_cfg, valid),
    }

    if score:
        report = G.score_rows(valid, results, ctx) if len(valid) >= 2 else {}
        meta["expected"] = {
            "total_weighted_score": report.get("total_weighted_score", 0.0),
            "consistency_factor": report.get("consistency_factor", 0.0),
            "distribution_fidelity_factor": report.get("distribution_fidelity_factor", 0.0),
            "final_score": report.get("final_score", 0.0),
        }
    return meta, rows


def load_expected_scores(path: str) -> dict[str, float]:
    """Prior per-task scores from a genExp sweep, used to confirm this run reproduces them.

    Generation is deterministic in (contract, config), so any drift means a code or config change
    since the sweep — worth surfacing rather than silently shipping different rows.
    """
    file = Path(path)
    if not file.exists():
        return {}
    document = json.load(open(file))
    records = document["tasks"] if isinstance(document, dict) else document
    return {r["task_id"]: r["final_score"] for r in records if "final_score" in r}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="submission.json")
    parser.add_argument("--per-task-dir", default=None,
                        help="also write one uploadable bare array per task into this directory")
    parser.add_argument("--task-id", default=None, help="build a single task instead of all")
    parser.add_argument("--limit", type=int, default=None, help="only the newest N tasks")
    parser.add_argument("--include-zero-seed", action="store_true",
                        help="keep the placeholder seed==0 tasks")
    parser.add_argument("--strategy", choices=("pure", "shaped"), default="pure")
    parser.add_argument("--construction", choices=tuple(G.CONSTRUCTIONS), default="mh",
                        help="rule every row's outcome must satisfy under the pure strategy")
    parser.add_argument("--flank", type=int, default=G.GenConfig.flank)
    parser.add_argument("--variants", type=int, default=G.GenConfig.variants)
    parser.add_argument("--lengths", default="20,23")
    parser.add_argument("--rows", type=int, default=None, help="override contract max_experiments")
    parser.add_argument("--no-score", action="store_true",
                        help="skip the stage-4/5 scoring pass")
    parser.add_argument("--compare-with", default="result.json",
                        help="prior sweep to cross-check expected scores against")
    parser.add_argument("--indent", type=int, default=2,
                        help="JSON indent for the combined file; 0 writes it compact")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.time()

    print("[1/4] fetching task history")
    tasks = G.fetch_all_tasks()
    if args.task_id:
        tasks = [t for t in tasks if t["id"] == args.task_id]
        if not tasks:
            print(f"task {args.task_id} not found", file=sys.stderr)
            return 1
    else:
        if not args.include_zero_seed:
            kept = [t for t in tasks if t["content"]["contract"].get("seed")]
            print(f"  {len(tasks)} tasks, {len(tasks) - len(kept)} with seed==0 skipped")
            tasks = kept
        if args.limit:
            tasks = tasks[:args.limit]

    cell_types = G.fetch_cell_types()
    lengths = tuple(int(v) for v in args.lengths.split(","))
    cfg = G.GenConfig(strategy=args.strategy,
                      selection="packed" if args.strategy == "pure" else "stratified",
                      flank=args.flank, variants=args.variants, lengths=lengths, rows=args.rows,
                      construction=args.construction)

    print("[2/4] warming reference + site cache")
    warm = G.build_context(tasks[0]["content"]["contract"],
                           tasks[0]["content"]["hbb_reference"], cell_types)
    sites = G.enumerate_sites(warm, cfg.flank, tuple(sorted(set(lengths))))
    print(f"  {len(sites)} PAM sites  {dict(Counter((s.cas, s.strand) for s in sites))}")

    expected = load_expected_scores(args.compare_with) if not args.no_score else {}
    per_task_dir = Path(args.per_task_dir) if args.per_task_dir else None
    if per_task_dir:
        per_task_dir.mkdir(parents=True, exist_ok=True)

    print(f"[3/4] building {len(tasks)} submission(s) with the {cfg.strategy} strategy")
    entries: list[dict] = []
    problems: list[str] = []
    drift: list[str] = []

    for index, task in enumerate(tasks):
        meta, rows = build_for_task(task, cell_types, cfg, score=not args.no_score)
        entry = {**meta, "submission": rows}
        entries.append(entry)

        if meta["problems"]:
            problems.append(f"{task['id'][:8]}: {'; '.join(meta['problems'])}")

        note = ""
        if meta["task_id"] in expected and "expected" in meta:
            delta = meta["expected"]["final_score"] - expected[meta["task_id"]]
            if abs(delta) > 1e-6:
                drift.append(f"{task['id'][:8]}: {delta:+.6f}")
                note = f"  drift={delta:+.4f}"

        if per_task_dir:
            with open(per_task_dir / f"{task['id']}.json", "w") as handle:
                json.dump(rows, handle, indent=2)

        elapsed = time.time() - started
        eta = elapsed / (index + 1) * (len(tasks) - index - 1)
        score_text = (f"final={meta['expected']['final_score']:8.3f}"
                      if "expected" in meta else "not scored")
        print(f"  [{index + 1:>3}/{len(tasks)}] {task['id'][:8]} {meta['cell_type']:<11} "
              f"rows={meta['rows']:>3} {score_text}{note}  eta {eta / 60:.1f}m")

    print(f"[4/4] writing {args.out}")
    document = {
        "generated_by": "submission.py",
        "format": "tasks[].submission is the JSON array a miner PUTs to the presigned S3 URL",
        "elapsed_seconds": time.time() - started,
        "config": {k: (list(v) if isinstance(v, tuple) else v) for k, v in vars(cfg).items()}
        | ({"weight_skew": "fitted per contract — see tasks[].weight_skew"}
           if cfg.strategy == "pure" else {}),
        "tasks": entries,
    }
    with open(args.out, "w") as handle:
        json.dump(document, handle, indent=args.indent or None)

    total_rows = sum(e["rows"] for e in entries)
    size_mb = Path(args.out).stat().st_size / 1e6
    print(f"\n  {len(entries)} task(s), {total_rows} rows, {size_mb:.1f} MB")
    if per_task_dir:
        print(f"  uploadable per-task arrays in {per_task_dir}/")

    if problems:
        print(f"  ! {len(problems)} task(s) violated a submission invariant:")
        for line in problems[:10]:
            print(f"      {line}")
    else:
        print(f"  invariants clean: unique ids, unique designs, '{cfg.construction}' conformance")

    shared = sum(1 for e in entries if e["task_id"] in expected)
    if drift:
        print(f"  ! {len(drift)} of {shared} shared task(s) drifted from {args.compare_with}:")
        for line in drift[:10]:
            print(f"      {line}")
    elif shared:
        print(f"  reproduces {args.compare_with} exactly on all {shared} shared task(s)")
    elif expected:
        # Newer than the sweep — the backend keeps issuing tasks, so this is expected at the head.
        print(f"  no overlap with {args.compare_with}; nothing to cross-check")

    print(f"\ndone in {time.time() - started:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
