#!/usr/bin/env python3
"""hedge_task.py — score the hedge builders on one real round, against that round's real field.

The fleet reserves one hotkey for the hedge: the rung that pays a flat score on *every* seed
instead of a spike on ~13 of 900. Two builders can fill it, and CLAUDE.md prices them from
averages over many contracts (all-cut ~+23 over the seed-agnostic hedge on K562). This asks the
narrower question the payout curve actually answers: on one specific round, with its three stamped
seeds and the 248 miners who competed, where would each have *ranked*?

Both are scored through all five stages on the round's own seeds, then dropped into the real field
alongside what our nine all-HDR hotkeys actually scored.
"""
import os

os.environ["NIOME_INSTANCE"] = "hedge"

import json
import sys
import time
import urllib.request
from collections import Counter

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_cut as AC
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

ROWS_CACHE = "hedge_rows.json"     # built rows, so a scoring bug never costs a 205s scan again
TASK_FILE = "task-b9051bc7.json"
SCORES_URL = "https://niome-api.genomes.io/api/v3/miners/scores?limit=40000"
OUT = "hedge_task.json"
# our nine, in miner.sh order, so the round's real all-HDR results can sit in the same table
OURS = [
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW", "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU", "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2", "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",
    "5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3", "5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN",
    "5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2",
]


def score(rows, contract, reference, cell_types, seeds):
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    """Run one row set through all five stages on every stamped seed and average, as a validator
    does. stage12 rewrites data/submission.json in place, so it is written fresh each time."""
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    per = []
    for s in seeds:
        run_stage3(seed=s)
        run_stage4(seed=s)
        per.append(run_stage5())
    keys = ("total_weighted_score", "consistency_factor",
            "distribution_fidelity_factor", "final_score")
    avg = {k: sum(f[k] for f in per) / len(per) for k in keys}
    return avg, [{"seed": s, "consistency_factor": f["consistency_factor"],
                  "final_score": f["final_score"]} for s, f in zip(seeds, per)]


def field_for(task_id):
    """The round's real field: best row per hotkey, ranked by final score."""
    raw = json.load(urllib.request.urlopen(SCORES_URL, timeout=300))
    raw = raw if isinstance(raw, list) else (raw.get("data") or raw.get("items") or [])
    best = {}
    for x in raw:
        if x["task_id"] != task_id:
            continue
        hk = x["miner_hotkey"]
        if hk not in best or x["final_score"] > best[hk]["final_score"]:
            best[hk] = x
    return sorted(best.values(), key=lambda y: -y["final_score"])


def main():
    task = json.load(open(TASK_FILE))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell = contract["cell_type"]
    seeds = [int(x) for x in str(contract["seed"]).split(",")]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    print(f"task {task['id'][:8]}  {cell}  seeds {seeds}  rows "
          f"{contract['rules']['max_experiments']}\n", flush=True)

    built = {}
    cache = {}
    if os.path.exists(ROWS_CACHE):
        blob = json.load(open(ROWS_CACHE))
        if blob.get("task_id") == task["id"]:
            cache = blob.get("rows") or {}
            for name, rows in cache.items():
                built[name] = (rows, blob.get("secs", {}).get(name, 0.0), {"cached": True})
                print(f"{name:<16} reusing {len(rows)} cached rows", flush=True)

    # --- all-cut: pins is_cut across a 900-seed window; consistency ~0.20-0.26 on every seed ---
    t0 = time.monotonic()
    rows, meta = ((None, {"reason": "cached"}) if "all-cut" in built
                  else AC.build_for_cell(contract, reference, cell_types, budget_s=900.0))
    if rows:
        built["all-cut"] = (rows, round(time.monotonic() - t0, 1), meta)
        print(f"all-cut          built {len(rows)} rows in {time.monotonic() - t0:.0f}s "
              f"| {meta.get('clean', '?')} clean", flush=True)
    else:
        print(f"all-cut          DECLINED: {meta.get('reason')}", flush=True)

    # --- the seed-agnostic hedge: the rung below all-cut, strict Cas9 + min-union Cas12a ---
    t0 = time.monotonic()
    try:
        rows, meta = ((None, {"reason": "cached"}) if "seed-agnostic" in built
                      else SA.build_submission(contract, reference, cell_types, task))
    except Exception as exc:
        rows, meta = None, {"reason": str(exc)}
    if rows:
        built["seed-agnostic"] = (rows, round(time.monotonic() - t0, 1), meta)
        print(f"seed-agnostic    built {len(rows)} rows in {time.monotonic() - t0:.0f}s",
              flush=True)
    else:
        print(f"seed-agnostic    DECLINED: {meta.get('reason')}", flush=True)

    json.dump({"task_id": task["id"],
               "rows": {k: v[0] for k, v in built.items()},
               "secs": {k: v[1] for k, v in built.items()}},
              open(ROWS_CACHE, "w"))
    out = {"task_id": task["id"], "cell": cell, "seeds": seeds, "builds": {}}
    print(f"\n{'build':<16} {'rows':>5} {'cells':>6} {'weighted':>9} {'cons':>6} {'fid':>6} "
          f"{'FINAL':>8}  per-seed cons")
    for name, (rows, secs, meta) in built.items():
        avg, per = score(rows, contract, reference, cell_types, seeds)
        cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
        out["builds"][name] = {"rows": len(rows), "cells": cells, "build_s": secs,
                               **avg, "per_seed": per}
        print(f"{name:<16} {len(rows):>5} {cells:>6} {avg['total_weighted_score']:>9.1f} "
              f"{avg['consistency_factor']:>6.4f} {avg['distribution_fidelity_factor']:>6.4f} "
              f"{avg['final_score']:>8.2f}  "
              f"{' '.join(f'{p['consistency_factor']:.4f}' for p in per)}", flush=True)
        json.dump(out, open(OUT, "w"), indent=1)

    # --- drop each build into the round's real field ---
    field = field_for(task["id"])
    out["field_size"] = len(field)
    print(f"\n=== the real field for this round: {len(field)} miners ===")
    print(f"{'rank':>4} {'final':>7} {'weighted':>9} {'cons':>7} {'fid':>6}  who")
    ours = {hk: i for i, hk in enumerate(OURS)}
    for i, x in enumerate(field[:12], 1):
        b = x["breakdown"]
        tag = f"h{ours[x['miner_hotkey']]} (ours, all-HDR)" if x["miner_hotkey"] in ours else ""
        print(f"{i:>4} {x['final_score']:>7.2f} {b['total_weighted_score']:>9.1f} "
              f"{b['consistency_factor']:>7.4f} {b['distribution_fidelity_factor']:>6.4f}  {tag}")
    print("  ...")
    out["our_ranks"] = {f"h{ours[x['miner_hotkey']]}": i for i, x in enumerate(field, 1)
                        if x["miner_hotkey"] in ours}
    print(f"\nour nine all-HDR hotkeys ranked: "
          f"{', '.join(f'h{k[1:]}={v}' for k, v in sorted(out['our_ranks'].items(), key=lambda kv: kv[1]))}")

    finals = [x["final_score"] for x in field]
    print(f"\n=== where each hedge build would have ranked ===")
    print(f"{'build':<16} {'FINAL':>8} {'rank':>6} {'paid?':>6}  vs our best all-HDR this round")
    our_best = min(out["our_ranks"].values())
    our_best_final = field[our_best - 1]["final_score"]
    for name, rec in out["builds"].items():
        f = rec["final_score"]
        rank = sum(1 for v in finals if v > f) + 1
        rec["rank_in_field"] = rank
        paid = "YES" if rank <= 10 else "no"
        print(f"{name:<16} {f:>8.2f} {rank:>6} {paid:>6}  "
              f"{f - our_best_final:+.2f} vs h{[k for k, v in out['our_ranks'].items() if v == our_best][0][1:]} "
              f"({our_best_final:.2f}, rank {our_best})")
    print(f"\nrank-10 cutoff this round: {finals[9]:.2f}   rank-1: {finals[0]:.2f}   "
          f"median: {finals[len(finals) // 2]:.2f}")
    out["cutoff_10"] = finals[9]
    out["rank_1"] = finals[0]
    out["median"] = finals[len(finals) // 2]
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
