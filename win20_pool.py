#!/usr/bin/env python3
"""win20_pool.py — is the width-20 band gain the *window*, or just a bank that fits the load cap?

``load_bank`` truncates to the 60,000 lowest-fail records whatever ``bank_keep`` wrote to disk.
That cap is hard-binding at width 100 (300,000 on disk, a fifth of it used) and nearly non-binding
at width 20 (52,341 of ~149,000 records already fail <= 5 of 20). Since the band bound is
``pool x P(HDR)**B >= group_size``, the pool handed to the greedy is a lever in its own right, and
the two effects are confounded in the width sweep.

Bank build and Cas9 scan are skipped: the band is the complement of the min-union's failed-seed
set, so ``load_bank`` + ``FastGreedy`` is the whole experiment. Banks must already be cached by
win20.py.
"""
import dataclasses
import json
import os
import sys
import time

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics.all_cut import bank_key, load_bank

TASK = "task-k562.json"
JOBS = [((300, 399), 45, [60_000, 150_000, 300_000]),
        ((300, 319), 6, [60_000, 100_000, 150_000])]
OUT = "win20_pool.json"


def main():
    task = json.load(open(TASK))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)
    out = []
    print(f"{'window':>9} {'mf':>3} {'limit':>8} {'pool':>8} {'maxfails':>9} {'band':>5} "
          f"{'frac':>6} {'greedy':>8}")
    for window, max_fail, limits in JOBS:
        cfg = dataclasses.replace(AH.config_for(contract["cell_type"]),
                                  hdr_range=window, main_max_fail=max_fail)
        path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            print(f"{window[0]}-{window[1]:<4} {max_fail:>3}  no cached bank; run win20.py first")
            continue
        width = window[1] - window[0] + 1
        for limit in limits:
            records = load_bank(path, limit=limit)
            started = time.monotonic()
            sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                                per_cell_min=cfg.per_cell_min,
                                caps=AH._group_caps(contract, ctx, cfg))
            index, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
            bad = set()
            for i in index:
                bad.update(int(x) for x in records[i]["fails"])
            band = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
            worst = max(len(records[i]["fails"]) for i in index)
            rec = {"window": list(window), "max_fail": max_fail, "limit": limit,
                   "pool": len(records), "group_worst_fails": worst, "band": band,
                   "clean": len(band), "clean_fraction": round(len(band) / width, 4),
                   "greedy_s": round(time.monotonic() - started, 1)}
            out.append(rec)
            print(f"{window[0]}-{window[1]:<4} {max_fail:>3} {limit:>8} {len(records):>8} "
                  f"{worst:>9} {len(band):>5} {len(band) / width * 100:>5.1f}% "
                  f"{rec['greedy_s']:>7.1f}s", flush=True)
            json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
