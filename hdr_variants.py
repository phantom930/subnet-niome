#!/usr/bin/env python3
"""hdr_variants.py — all-HDR `variants` vs band size and scan time, at a wide band window.

The width-300 fleet layout made the Cas12a HDR bank scan the fleet's binding resource: cost is
roughly `variants x seeds`, so tripling the window triples the scan, and seven hotkeys plus h0's
all-cut can exceed the 900s prefetch budget on one GPU. `variants` is the knob that buys it back —
but it also feeds the band, so the question is what band size survives a cheaper scan.

Timing is only meaningful single-tenant: the same build measured 55s, 185s and 201s depending on how
many siblings shared the GPU. Run this with the fleet idle and read the numbers as a ratio, then
multiply by the concurrency the fleet actually runs at.

`bank_key` folds in `variants`, so each arm builds its own bank; the window defaults to one no
hotkey owns, so every arm measures a COLD scan rather than loading a sibling's cache.

    python hdr_variants.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hdrvar")

import dataclasses as _dc          # noqa: E402
import json                        # noqa: E402
import logging                     # noqa: E402
import time                        # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                 # noqa: E402
from niome_subnet.genomics import all_cut as AC     # noqa: E402
from niome_subnet.genomics import all_hdr as AH     # noqa: E402
from niome_subnet.genomics import mt19937 as MT     # noqa: E402
from sd_task import task_content                    # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "bba85ff0-d2d1-43dc-b47d-a85ac4e2077e"
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "150-449"
VARIANTS = [int(x) for x in (os.getenv("HV_VARIANTS") or "5000,11000,22000,44000").split(",")]


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = (int(x) for x in WINDOW.split("-"))
    span = hi - lo + 1
    base = AH.config_for(cell)
    mf = max(1, round(base.main_max_fail * span / 100.0))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    print(f"task {TASK[:8]}  {cell}  band window {lo}-{hi} (span {span})  "
          f"main_max_fail {mf}  group {base.group_size}  rows {n_rows}\n")
    print(f"  {'variants':>8} {'bank':>7} {'scan s':>7} {'band':>5} {'union':>6} {'cas9':>6} "
          f"{'cells':>5} {'rows':>5} {'rest s':>7} {'total s':>8}  {'s per band seed':>15}")
    out = []
    for v in VARIANTS:
        cfg = _dc.replace(base, hdr_range=(lo, hi), variants=v, main_max_fail=mf)
        path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
        if os.path.exists(path):
            os.remove(path)                      # force a cold scan for an honest timing
        t0 = time.monotonic()
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        t_scan = time.monotonic() - t0
        if not bank:
            print(f"  {v:>8} {'0':>7} {t_scan:>7.0f}   no qualifying guide at mf {mf}")
            out.append({"variants": v, "scan_s": t_scan, "reason": "empty bank"})
            MT.free_gpu_memory()
            continue
        AH.save_bank(path, bank)
        t1 = time.monotonic()
        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=600)
        t_rest = time.monotonic() - t1
        MT.free_gpu_memory()
        rec = {"variants": v, "bank": len(bank), "scan_s": round(t_scan, 1),
               "band": meta.get("clean"), "union": meta.get("union"),
               "cas9": meta.get("cas9_pool"), "cells": meta.get("cells"),
               "rows": len(rows) if rows else 0, "rest_s": round(t_rest, 1),
               "reason": None if rows else meta.get("reason")}
        out.append(rec)
        if not rows:
            print(f"  {v:>8} {len(bank):>7} {t_scan:>7.0f} {str(meta.get('clean')):>5} "
                  f"{str(meta.get('union')):>6} {str(meta.get('cas9_pool')):>6}   "
                  f"declined: {meta.get('reason')}")
            continue
        tot = t_scan + t_rest
        per = tot / meta["clean"] if meta.get("clean") else float("nan")
        print(f"  {v:>8} {len(bank):>7} {t_scan:>7.0f} {meta['clean']:>5} {meta['union']:>6} "
              f"{meta['cas9_pool']:>6} {meta['cells']:>4}/8 {len(rows):>5} {t_rest:>7.0f} "
              f"{tot:>8.0f}  {per:>15.1f}")
    json.dump({"task": TASK, "cell": cell, "window": [lo, hi], "main_max_fail": mf,
               "arms": out}, open("hdr_variants.json", "w"), indent=1)
    good = [r for r in out if r.get("rows")]
    if good:
        print(f"\n  fleet arithmetic — 7 band hotkeys + h0's all-cut (~354s) against a 900s budget:")
        for r in good:
            tot = r["scan_s"] + r["rest_s"]
            print(f"    variants {r['variants']:>6}: 7 x {tot:>5.0f}s + 354s = "
                  f"{7 * tot + 354:>6.0f}s of GPU work   "
                  f"{'FITS' if 7 * tot + 354 <= 900 else 'over budget'}"
                  f"   band {r['band']}")
    print("\nwrote hdr_variants.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
