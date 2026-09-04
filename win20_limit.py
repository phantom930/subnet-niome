#!/usr/bin/env python3
"""win20_limit.py — does the +1 band seed from a larger load limit survive the Cas9 half?

win20_pool.py shows the width-100 band goes 13 -> 14 when ``load_bank``'s 60,000 cap is lifted, for
half a second of extra greedy. That is only worth having if the build still completes: each seed
added to the band multiplies the conditional Cas9 requirement by P(HDR) ~ 0.57, and the pool at
band 15 was already down to 848 from 1361. Builds and scores the width-100 window at limit 150,000
against the shipped 60,000, on the contract's own seeds.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import dataclasses
import json
import sys
import time
from collections import Counter

import numpy as np

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.all_cut import bank_key, load_bank
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

TASK = "task-k562.json"
WINDOW = (300, 399)
LIMITS = [60_000, 150_000]


def main():
    task = json.load(open(TASK))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    seeds = [int(x) for x in str(contract["seed"]).split(",")]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    cfg = dataclasses.replace(AH.config_for(contract["cell_type"]),
                              hdr_range=WINDOW, main_max_fail=45)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    want = n_rows - cfg.group_size

    print(f"task {task['id'][:8]}  {contract['cell_type']}  window {WINDOW}  seeds {seeds}\n")
    print(f"{'limit':>8} {'band':>5} {'cas9':>6} {'rows':>5} {'cells':>6} {'weighted':>9} "
          f"{'cons':>6} {'fid':>6} {'FINAL':>7} {'build':>7}")
    out = []
    for limit in LIMITS:
        started = time.monotonic()
        records = load_bank(path, limit=limit)
        sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                            per_cell_min=cfg.per_cell_min,
                            caps=AH._group_caps(contract, ctx, cfg))
        index, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
        group = [records[i] for i in index]
        bad = set()
        for r in group:
            bad.update(int(x) for x in r["fails"])
        clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad), dtype=np.int64)
        cas9 = AH.scan_cas9(clean, contract, cell_types, ctx, sites, cfg, want, None)
        MT.free_gpu_memory()
        cells4 = len({(r["mutation"], r["strand"]) for r in cas9})
        if len(cas9) < want or cells4 < 4:
            print(f"{limit:>8} {clean.size:>5} {len(cas9):>6}  DECLINED (short of {want})",
                  flush=True)
            out.append({"limit": limit, "band": [int(x) for x in clean], "cas9_pool": len(cas9),
                        "declined": True})
            continue
        rows = AH.assemble(group, cas9, contract, ctx, cfg, n_rows)
        build_s = time.monotonic() - started
        cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
        json.dump(contract, open(settings.CONTRACT_PATH, "w"))
        json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        per_seed = []
        for seed in seeds:
            run_stage3(seed=seed)
            run_stage4(seed=seed)
            per_seed.append(run_stage5())
        avg = {k: sum(f[k] for f in per_seed) / len(per_seed)
               for k in ("total_weighted_score", "consistency_factor",
                         "distribution_fidelity_factor", "final_score")}
        rec = {"limit": limit, "band": [int(x) for x in clean], "clean": int(clean.size),
               "cas9_pool": len(cas9), "rows": len(rows), "cells": cells,
               "k": len([s for s in seeds if s in set(int(x) for x in clean)]),
               "build_s": round(build_s, 1), **avg}
        out.append(rec)
        print(f"{limit:>8} {clean.size:>5} {len(cas9):>6} {len(rows):>5} {cells:>4}/8 "
              f"{avg['total_weighted_score']:>9.1f} {avg['consistency_factor']:>6.3f} "
              f"{avg['distribution_fidelity_factor']:>6.3f} {avg['final_score']:>7.2f} "
              f"{build_s:>6.1f}s", flush=True)
        print(f"         band {[int(x) for x in clean]}", flush=True)
        json.dump(out, open("win20_limit.json", "w"), indent=1)


if __name__ == "__main__":
    main()
