#!/usr/bin/env python3
"""robustness.py — option 1: how fragile is each construction to the unknown real seed?

Cas12a rows no_cut under ~4% of seeds and nothing fixes that, so a build made at one seed is scored
under a stream it never saw. This measures the cost distributionally: build a submission at
``--build-seed`` (0, the unstamped case) under each construction, then re-score it across a wide,
out-of-sample set of real seeds and report the spread of consistency_factor and final_score. The
least-fragile construction is the one with the best worst-case and mean across that set — not the one
that happens to peak on any single seed.

    python robustness.py --task test/task.json --constructions mh,hdr --seeds 60

The seeds are spread across 1..~60000 on a prime stride so none coincide with the build seed and the
sample is not clustered in 100..999 (which earlier runs showed a guide can overfit to).
"""

from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import statistics as st
import sys
import time
from collections import Counter
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import niome_subnet.utils.settings as settings  # noqa: E402
import genExp as G  # noqa: E402

STATE: dict = {}


def order_aligned(rows, valid, results):
    """order_rows sorts rows only; realign valid/results so stage 4's KFold sees the same layout."""
    ordered = G.order_rows(rows, valid)
    index = {e["experiment"]["experiment_id"]: i for i, e in enumerate(valid)}
    order = [index[r["experiment_id"]] for r in ordered]
    return ordered, [valid[i] for i in order], [results[i] for i in order]


def build_once(contract, reference, cell_types, construction, build_seed, weight_skew):
    """Build one construction's submission at the build seed — the rows the miner would ship."""
    from neurons.miner import Miner
    bc = copy.deepcopy(contract)
    bc["seed"] = build_seed
    ctx = G.build_context(bc, reference, cell_types)
    cfg = replace(Miner.gen_config_for(bc), construction=construction, weight_skew=weight_skew)
    sites = G.enumerate_sites(ctx, cfg.flank, tuple(sorted(set(cfg.lengths))))
    rows, valid, results = G._generate(ctx, sites, cfg)
    rows, valid, results = order_aligned(rows, valid, results)
    return cfg, rows, valid, results


def worker_init(payload: dict) -> None:
    logging.getLogger().setLevel(logging.ERROR)
    sys.stdout = open(os.devnull, "w")
    task = json.load(open(payload["task"]))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell_types = payload["cell_types"]

    from neurons.miner import Miner
    bc = copy.deepcopy(contract)
    bc["seed"] = payload["build_seed"]
    ctx0 = G.build_context(bc, reference, cell_types)
    sites = G.enumerate_sites(ctx0, Miner.gen_config_for(bc).flank, (20, 23))
    weight_skew = G.choose_weight_skew(ctx0, sites,
                                       Miner.gen_config_for(bc))  # seed-independent, fit once

    builds = {}
    for construction in payload["constructions"]:
        cfg, rows, valid, results = build_once(contract, reference, cell_types,
                                               construction, payload["build_seed"], weight_skew)
        builds[construction] = (rows, valid)
    STATE.update(contract=contract, reference=reference, cell_types=cell_types, builds=builds)


def evaluate(job: tuple[str, int]) -> dict:
    """Re-score one construction's build under one real seed."""
    construction, seed = job
    rows, _valid = STATE["builds"][construction]
    sc = copy.deepcopy(STATE["contract"])
    sc["seed"] = seed
    ctx = G.build_context(sc, STATE["reference"], STATE["cell_types"])
    report, valid, results = G.rescore_under(rows, ctx)
    outcomes = Counter(r["outcome"] for r in results)
    return {
        "construction": construction,
        "seed": seed,
        "final_score": report.get("final_score", 0.0),
        "consistency_factor": report.get("consistency_factor", 0.0),
        "fidelity": report.get("distribution_fidelity_factor", 0.0),
        "total_weighted_score": report.get("total_weighted_score", 0.0),
        "no_cut": outcomes.get("no_cut", 0),
    }


def spread(values: list[float]) -> dict:
    values = sorted(values)
    return {
        "min": values[0], "p10": values[max(0, int(0.10 * len(values)))],
        "median": st.median(values), "mean": st.mean(values),
        "max": values[-1], "std": st.pstdev(values),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task", default="task.json")
    parser.add_argument("--out-dir", default="test/robustness")
    parser.add_argument("--constructions", default="mh,hdr")
    parser.add_argument("--build-seed", type=int, default=0)
    parser.add_argument("--seeds", type=int, default=60, help="number of real seeds to sample")
    parser.add_argument("--stride", type=int, default=997, help="prime stride between sampled seeds")
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--cell-types", default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.ERROR, format="%(levelname)s %(message)s")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    task = json.load(open(args.task))
    contract = task["content"]["contract"]
    cell_types = (json.load(open(args.cell_types)) if args.cell_types else G.fetch_cell_types())
    constructions = args.constructions.split(",")
    # Wide, out-of-sample seeds on a prime stride: never the build seed, never clustered.
    seeds = [args.build_seed + 1 + i * args.stride for i in range(args.seeds)]

    print("=" * 96)
    print(f"  task {task.get('id')}   cell_type {contract.get('cell_type')}   "
          f"real seed (ignored here) {contract.get('seed')}")
    print(f"  build seed {args.build_seed}; scored across {len(seeds)} seeds "
          f"[{seeds[0]}..{seeds[-1]}] stride {args.stride}")
    print(f"  constructions: {constructions}")
    print("=" * 96)

    payload = {"task": str(Path(args.task).resolve()), "cell_types": cell_types,
               "build_seed": args.build_seed, "constructions": constructions}
    jobs = [(c, s) for c in constructions for s in seeds]

    started = time.time()
    records = []
    with Pool(args.jobs, initializer=worker_init, initargs=(payload,)) as pool:
        for i, rec in enumerate(pool.imap_unordered(evaluate, jobs, chunksize=4), 1):
            records.append(rec)
            if i % 40 == 0 or i == len(jobs):
                print(f"    [{i:>4}/{len(jobs)}] {i / (time.time() - started):.0f} scores/s")

    summary = {}
    print(f"\n  {'construction':<14}{'final: min':>11}{'p10':>9}{'median':>9}{'mean':>9}"
          f"{'max':>9}{'std':>8}{'cons.mean':>10}{'nocut.mean':>11}")
    for c in constructions:
        fs = [r["final_score"] for r in records if r["construction"] == c]
        cons = [r["consistency_factor"] for r in records if r["construction"] == c]
        ncut = [r["no_cut"] for r in records if r["construction"] == c]
        s = spread(fs)
        summary[c] = {"final": s, "consistency": spread(cons), "no_cut": spread(ncut),
                      "fidelity_mean": st.mean([r["fidelity"] for r in records
                                                if r["construction"] == c])}
        print(f"  {c:<14}{s['min']:>11.2f}{s['p10']:>9.2f}{s['median']:>9.2f}{s['mean']:>9.2f}"
              f"{s['max']:>9.2f}{s['std']:>8.2f}{st.mean(cons):>10.4f}{st.mean(ncut):>11.2f}")

    # Least fragile: rank by worst-case (min) then mean, since the miner faces one draw of the seed.
    ranked = sorted(constructions, key=lambda c: (summary[c]["final"]["min"],
                                                  summary[c]["final"]["mean"]), reverse=True)
    print(f"\n  least fragile (worst-case final, then mean): {ranked[0]}")
    for c in ranked:
        s = summary[c]["final"]
        print(f"    {c:<8} worst {s['min']:.2f}  mean {s['mean']:.2f}  "
              f"coefficient-of-variation {s['std'] / max(s['mean'], 1e-9):.2f}")

    with open(out_dir / "records.jsonl", "w") as h:
        for r in sorted(records, key=lambda r: (r["construction"], r["seed"])):
            h.write(json.dumps(r) + "\n")
    with open(out_dir / "summary.json", "w") as h:
        json.dump({"task_id": task.get("id"), "cell_type": contract.get("cell_type"),
                   "build_seed": args.build_seed, "seeds": seeds,
                   "constructions": constructions, "least_fragile": ranked[0],
                   "per_construction": summary}, h, indent=2)
    print(f"\n  elapsed {time.time() - started:.0f}s; artifacts in {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
