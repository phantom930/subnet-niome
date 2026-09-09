#!/usr/bin/env python3
"""hek_time.py — how long h0's all-cut takes on HEK293, and the band build it competes with.

h0 now runs all-cut on every cell type, and `ALL_CUT_MIN_BUDGET_S` treats HEK293 as the cheap case
(190s against K562's 480s). That gate was set when HEK293 ran its own clustered builder and before
seven siblings screened width-300 band windows on the same GPU, so it is worth timing directly:
HEK293's accessibility is 0.35 against K562's 0.77, which collapses the cut-clean set (~137 of 900
vs ~560) and changes how much work the Cas9 scan does.

Both halves of an HEK293 round are timed cold — h0's all-cut, and one band hotkey's width-300
all-HDR — because the fleet question is whether they coexist, not what either costs alone.

    python hek_time.py [task_id] [band_window]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hektime")

import dataclasses as _dc      # noqa: E402
import json                    # noqa: E402
import logging                 # noqa: E402
import time                    # noqa: E402
from collections import Counter  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np             # noqa: E402

import genExp as G             # noqa: E402
from niome_subnet.genomics import all_cut as AC       # noqa: E402
from niome_subnet.genomics import all_hdr as AH       # noqa: E402
from niome_subnet.genomics import fastgreedy as FG    # noqa: E402
from niome_subnet.genomics import mt19937 as MT       # noqa: E402
from sd_task import task_content                      # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "414dab89-1bf9-4f16-ba6f-93c43186b8cf"
BAND = _ARGV[2] if len(_ARGV) > 2 else "100-399"


def time_all_cut(contract, reference, cell_types, ctx, sites, n_rows):
    cfg = AC.config_for(contract["cell_type"]) or AC.AllCutConfig()
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if os.path.exists(path):
        os.remove(path)                                  # cold, or the timing is meaningless
    t0 = time.monotonic()
    bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
    t_bank = time.monotonic() - t0
    if not bank:
        return {"stage": "bank", "reason": "empty", "bank_s": round(t_bank, 1)}
    AC.save_bank(path, bank)
    t1 = time.monotonic()
    sel = FG.FastGreedy(bank, window_lo=cfg.start_seed, window_hi=cfg.end_seed)
    idx, _u = sel.best(cfg.group_size, restarts=12)
    group = [bank[i] for i in idx]
    bad = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
    t_greedy = time.monotonic() - t1
    t2 = time.monotonic()
    cas9 = AC.scan_cas9(np.array(clean, dtype=np.int64), contract, cell_types, ctx, sites, cfg,
                        n_rows - cfg.group_size)
    t_cas9 = time.monotonic() - t2
    if len(cas9) < n_rows - cfg.group_size:
        return {"stage": "cas9", "reason": f"pool {len(cas9)} short of {n_rows-cfg.group_size}",
                "bank_s": round(t_bank, 1), "greedy_s": round(t_greedy, 1),
                "cas9_s": round(t_cas9, 1), "clean": len(clean)}
    t3 = time.monotonic()
    rows = AC.assemble(group, cas9, contract, ctx, cfg, n_rows)
    t_asm = time.monotonic() - t3
    MT.free_gpu_memory()
    return {"bank": len(bank), "clean": len(clean), "clean_frac": len(clean) / 900.0,
            "cas9": len(cas9), "rows": len(rows),
            "cells": len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
            "bank_s": round(t_bank, 1), "greedy_s": round(t_greedy, 1),
            "cas9_s": round(t_cas9, 1), "assemble_s": round(t_asm, 1),
            "total_s": round(t_bank + t_greedy + t_cas9 + t_asm, 1)}


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    acc = cell_types.get(cell, {}).get("accessibility")
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    print(f"task {TASK[:8]}  {cell}  accessibility {acc}  rows {n_rows}\n")

    print("h0's job — all-cut, whole window 100-999, group 42 (cold bank):")
    ac = time_all_cut(contract, reference, cell_types, ctx, sites, n_rows)
    if ac.get("reason"):
        print(f"  DECLINED at {ac['stage']}: {ac['reason']}   "
              + "  ".join(f"{k}={v}" for k, v in ac.items() if k.endswith('_s')))
    else:
        print(f"  clean {ac['clean']}/900 ({ac['clean_frac']:.1%})  cas9 pool {ac['cas9']}  "
              f"rows {ac['rows']} cells {ac['cells']}/8")
        print(f"  bank {ac['bank_s']}s + greedy {ac['greedy_s']}s + cas9 {ac['cas9_s']}s "
              f"+ assemble {ac['assemble_s']}s = {ac['total_s']}s")

    lo, hi = (int(x) for x in BAND.split("-"))
    print(f"\na band hotkey's job — all-HDR, window {lo}-{hi} (cold bank):")
    cfg = AH.config_for(cell)
    span = hi - lo + 1
    probe = _dc.replace(cfg, hdr_range=(lo, hi),
                        main_max_fail=max(1, round(cfg.main_max_fail * span / 100.0)),
                        variants=min(cfg.variants, AH.WIDE_WINDOW_VARIANTS))
    p = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, probe)}.npz")
    if os.path.exists(p):
        os.remove(p)
    t0 = time.monotonic()
    rows, meta = AH.build_for_cell(contract, reference, cell_types, budget_s=600,
                                   hdr_range=(lo, hi))
    t_hdr = time.monotonic() - t0
    hdr = {"band": meta.get("clean"), "cas9": meta.get("cas9_pool"), "cells": meta.get("cells"),
           "rows": len(rows) if rows else 0, "total_s": round(t_hdr, 1),
           "reason": None if rows else meta.get("reason"),
           "main_max_fail": probe.main_max_fail, "variants": probe.variants}
    if rows:
        print(f"  band {meta['clean']}/{span}  cas9 pool {meta['cas9_pool']}  "
              f"rows {len(rows)} cells {meta['cells']}/8  |  {t_hdr:.0f}s  "
              f"(mf {probe.main_max_fail}, variants {probe.variants})")
    else:
        print(f"  DECLINED: {meta.get('reason')}  ({t_hdr:.0f}s)")

    gate = 190.0
    print(f"\n  ALL_CUT_MIN_BUDGET_S['{cell}'] is {gate:.0f}s.")
    if not ac.get("reason"):
        verdict = ("fits the gate" if ac["total_s"] <= gate
                   else f"EXCEEDS the gate by {ac['total_s']-gate:.0f}s")
        print(f"  h0's all-cut took {ac['total_s']}s single-tenant -> {verdict}.")
        print(f"  Contended, today's K562 round showed all-cut at 338s against 354s idle "
              f"(~1.0x), so treat the gate, not contention, as the risk here.")
    json.dump({"task": TASK, "cell": cell, "accessibility": acc,
               "all_cut": ac, "all_hdr": hdr}, open("hek_time.json", "w"), indent=1)
    print("\nwrote hek_time.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
