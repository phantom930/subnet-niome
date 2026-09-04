#!/usr/bin/env python3
"""win20_task.py — score a proposed fleet config on one real round, against the real field.

Builds every hotkey's band for one task under a chosen window width, scores all nine through the
five stages on the round's own stamped seeds, and ranks them in that round's actual field of miner
scores. The point is to see the config as a validator saw the round, not as an average.

Width is the whole question here: a hotkey owns a 100-seed territory, and a narrower band inside it
covers more seeds in total (HEK293: 10.1 of 900 at width 12 against 7.0 at width 100) but confines
them to one end of the territory. A round is decided by where its three seeds actually fall.
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
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

TASK_FILE = "task-29488fed.json"
WIDTH = int(os.getenv("WIN20_W", "12"))
MAX_FAIL = int(os.getenv("WIN20_MF", "6"))
WINDOWS = [(lo, lo + WIDTH - 1) for lo in range(100, 1000, 100)]
OUT = f"win20_task_w{WIDTH}.json"


def build(contract, reference, cell_types, ctx, sites, cfg):
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return None, None, "bank scan produced nothing"
        save_bank(path, bank)
    disk = int(np.load(path, allow_pickle=False)["fails"].shape[0])
    records = load_bank(path)
    if len(records) < cfg.group_size:
        return None, None, f"bank {len(records)} short of group {cfg.group_size}"
    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg))
    index, _ = sel.best(cfg.group_size, restarts=cfg.restarts)
    bad = set()
    for i in index:
        bad.update(int(x) for x in records[i]["fails"])
    band = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    cas9 = AH.scan_cas9(np.array(band, dtype=np.int64), contract, cell_types, ctx, sites, cfg,
                        n_rows - cfg.group_size, None)
    MT.free_gpu_memory()
    if len(cas9) < n_rows - cfg.group_size:
        return band, None, f"Cas9 pool {len(cas9)} short"
    rows = AH.assemble([records[i] for i in index], cas9, contract, ctx, cfg, n_rows)
    return band, (rows, disk, len(cas9)), None


def main():
    task = json.load(open(TASK_FILE))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell = contract["cell_type"]
    seeds = [int(x) for x in str(contract["seed"]).split(",")]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AH.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    print(f"task {task['id'][:8]}  {cell}  seeds {seeds}  width {WIDTH} mf {MAX_FAIL}\n")
    print(f"{'hk':>3} {'window':>9} {'disk':>7} {'band':>5} {'k':>2} {'cas9':>6} {'weighted':>9} "
          f"{'cons':>6} {'fid':>6} {'FINAL':>8} {'per-seed cons':>22} {'build':>7}")
    out = []
    for i, window in enumerate(WINDOWS):
        t0 = time.monotonic()
        cfg = dataclasses.replace(base, hdr_range=window, main_max_fail=MAX_FAIL)
        band, made, err = build(contract, reference, cell_types, ctx, sites, cfg)
        if made is None:
            print(f"h{i:<2} {window[0]}-{window[1]:<4} DECLINED: {err}", flush=True)
            out.append({"hotkey": f"h{i}", "window": list(window), "band": band, "declined": err})
            continue
        rows, disk, pool = made
        json.dump(contract, open(settings.CONTRACT_PATH, "w"))
        json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        per = []
        for s in seeds:
            run_stage3(seed=s)
            run_stage4(seed=s)
            per.append(run_stage5())
        avg = {k: sum(f[k] for f in per) / len(per) for k in
               ("total_weighted_score", "consistency_factor", "distribution_fidelity_factor",
                "final_score")}
        hits = [s for s in seeds if s in set(band)]
        rec = {"hotkey": f"h{i}", "window": list(window), "band": band, "disk": disk,
               "cas9_pool": pool, "k": len(hits), "hits": hits, "rows": len(rows),
               "cells": len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
               **avg, "per_seed": [{"seed": s, "consistency_factor": f["consistency_factor"],
                                    "final_score": f["final_score"]} for s, f in zip(seeds, per)],
               "build_s": round(time.monotonic() - t0, 1)}
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"h{i:<2} {window[0]}-{window[1]:<4} {disk:>7} {len(band):>5} {len(hits):>2} "
              f"{pool:>6} {avg['total_weighted_score']:>9.1f} {avg['consistency_factor']:>6.3f} "
              f"{avg['distribution_fidelity_factor']:>6.3f} {avg['final_score']:>8.2f} "
              f"{' '.join(f'{p['consistency_factor']:.3f}' for p in per):>22} "
              f"{rec['build_s']:>6.1f}s", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
    print("\ndone")


if __name__ == "__main__":
    main()
