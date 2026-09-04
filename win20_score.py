#!/usr/bin/env python3
"""win20_score.py — score each window's all-HDR build through all five stages on the real seeds.

win20.py measures the band; this measures what the band is worth. Both the width-100 baseline and
the five 20-seed windows are rebuilt (banks are cached by then, so this is the greedy + Cas9 scan
+ assemble only) and put through stage12/3/4/5 on the contract's own stamped seeds, per seed rather
than averaged, so the k-of-3 arithmetic is visible instead of inferred.

Runs under ``NIOME_INSTANCE=win20`` so it writes to ``data/inst/win20/`` and cannot disturb a live
hotkey's submission, task artifacts or local scoring.
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

TASK = "task-k562.json"
JOBS = [((300, 399), 45)] + [((lo, lo + 19), 6) for lo in range(300, 400, 20)]
OUT = "win20_scores.json"


def build(contract, reference, cell_types, window, max_fail):
    cfg = dataclasses.replace(AH.config_for(contract["cell_type"]),
                              hdr_range=window, main_max_fail=max_fail)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return None, None, "bank scan produced nothing"
        save_bank(path, bank)
    records = load_bank(path)
    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min,
                        caps=AH._group_caps(contract, ctx, cfg))
    index, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
    group = [records[i] for i in index]
    bad = set()
    for r in group:
        bad.update(int(x) for x in r["fails"])
    clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad), dtype=np.int64)
    if clean.size == 0:
        return None, None, "failures cover the whole window"
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    cas9 = AH.scan_cas9(clean, contract, cell_types, ctx, sites, cfg, n_rows - cfg.group_size, None)
    MT.free_gpu_memory()
    if len(cas9) < n_rows - cfg.group_size:
        return None, [int(x) for x in clean], f"Cas9 pool {len(cas9)} short"
    rows = AH.assemble(group, cas9, contract, ctx, cfg, n_rows)
    return rows, [int(x) for x in clean], None


def score(rows, contract, reference, cell_types, seeds):
    os.makedirs(settings.DATA_DIR, exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)                      # seed-independent, once
    per_seed = []
    for seed in seeds:
        run_stage3(seed=seed)
        run_stage4(seed=seed)
        per_seed.append(run_stage5())
    keys = ("total_weighted_score", "consistency_factor", "distribution_fidelity_factor",
            "final_score", "n_valid_experiments")
    return ({k: sum(f[k] for f in per_seed) / len(per_seed) for k in keys},
            [{"seed": s, "consistency_factor": f["consistency_factor"],
              "final_score": f["final_score"]} for s, f in zip(seeds, per_seed)])


def main():
    task = json.load(open(TASK))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    seeds = [int(x) for x in str(contract["seed"]).split(",")]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    print(f"task {task['id'][:8]}  {contract['cell_type']}  seeds {seeds}\n")
    print(f"{'window':>9} {'mf':>3} {'band':>5} {'k':>2} {'weighted':>9} {'cons':>6} "
          f"{'fid':>6} {'FINAL':>8}   per-seed cons")

    out = {"task_id": task["id"], "seeds": seeds, "results": []}
    for window, max_fail in JOBS:
        started = time.monotonic()
        rows, band, err = build(contract, reference, cell_types, window, max_fail)
        if rows is None:
            print(f"{window[0]}-{window[1]:<4} {max_fail:>3} "
                  f"{len(band or []):>5}  -  DECLINED: {err}", flush=True)
            out["results"].append({"window": list(window), "max_fail": max_fail, "band": band,
                                   "declined": err})
            continue
        avg, per_seed = score(rows, contract, reference, cell_types, seeds)
        hits = [s for s in seeds if s in band]
        rec = {"window": list(window), "max_fail": max_fail, "band": band, "k": len(hits),
               "hits": hits, "rows": len(rows),
               "cas_mix": dict(Counter(r["cas_system"] for r in rows)),
               "wall_s": round(time.monotonic() - started, 1), **avg, "per_seed": per_seed}
        out["results"].append(rec)
        cons_txt = " ".join(f"{p['consistency_factor']:.2f}" for p in per_seed)
        print(f"{window[0]}-{window[1]:<4} {max_fail:>3} {len(band):>5} {len(hits):>2} "
              f"{avg['total_weighted_score']:>9.1f} {avg['consistency_factor']:>6.3f} "
              f"{avg['distribution_fidelity_factor']:>6.3f} {avg['final_score']:>8.2f}   "
              f"{cons_txt}", flush=True)
        json.dump(out, open(OUT, "w"), indent=1)

    json.dump(out, open(OUT, "w"), indent=1)
    ok = [r for r in out["results"] if "final_score" in r]
    base = next((r for r in ok if r["max_fail"] == 45), None)
    narrow = [r for r in ok if r["max_fail"] == 6]
    if base and narrow:
        best = max(narrow, key=lambda r: r["final_score"])
        print(f"\nbest of the five 20-seed windows: {best['window']} k={best['k']} "
              f"final {best['final_score']:.2f}  vs baseline {base['final_score']:.2f}")
        print(f"band seeds covered in 300-399: baseline {len(base['band'])}, "
              f"five narrow windows {sum(len(r['band']) for r in narrow)}")


if __name__ == "__main__":
    main()
