#!/usr/bin/env python3
"""narrow_width.py — how narrow does the band window have to be to reach the top block's quality?

On HEK293 round a8b9f1bb, ranks 7-10 all scored a SINGLE band hit (cons 0.385-0.386, i.e.
(1 + 2*0.078)/3) and still cleared the rank-10 cutoff of 81.92, because they held weighted ~221 and
fidelity ~0.962. The bar is therefore

    weighted * fidelity >= 81.92 / 0.3854 = 212.6

against our 179.0 after tonight's cas9_cell_target fix -- we need +19% on the product.

The mechanism that could deliver it is the one measured tonight: the Cas9 half must reach HDR on
EVERY band seed, so the candidate pool scales as P(HDR)**band (0.358 per seed on HEK293). A smaller
band means a vastly larger pool, and `assemble` -- which ranks by (distance, |gc-0.50|) and then by
weighted_score -- can take nearly optimal rows. Measured at width 100 (band 7, pool 3986) vs width
300 (band 8, pool 1871): weighted +14%, fidelity +5%.

A width-10 window is what a 90-slot operator can afford, because 90 disjoint width-10 windows tile
the whole 900-seed space. We have 7 band hotkeys, so the same layout would cover ~50 seeds against
our current 56 -- worse. The question this answers is whether the QUALITY half transfers to a width
we can actually run, and where the knee is.

    python narrow_width.py [task_id]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "narrow")

import dataclasses as _dc     # noqa: E402
import json                   # noqa: E402
import logging                # noqa: E402
import statistics as st       # noqa: E402
import time                   # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G            # noqa: E402
from niome_subnet.genomics import all_cut as AC      # noqa: E402
from niome_subnet.genomics import all_hdr as AH      # noqa: E402
from niome_subnet.genomics import mt19937 as MT      # noqa: E402
from niome_subnet.utils import settings              # noqa: E402
from sd_task import score, task_content              # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "a8b9f1bb-f25d-42af-92a1-e0aff16bee3a"
WIDTHS = [int(x) for x in (os.getenv("NW_WIDTHS") or "10,20,30,50,100,300").split(",")]
START = int(os.getenv("NW_START", "300"))
N_HOTKEYS = 7
BAR = 212.6          # weighted * fidelity needed for a k=1 round to clear 81.92
FLOOR_CONS = 0.078   # measured HEK293 off-band floor, from the field's own cons values


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AH.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    cons_k1 = (1 + 2 * FLOOR_CONS) / 3
    print(f"task {TASK[:8]}  {cell}  group {base.group_size}  windows from seed {START}")
    print(f"k=1 consistency {cons_k1:.4f}; bar for 81.92 is weighted*fid >= {BAR:.1f}\n")
    print(f"  {'width':>6} {'mf':>4} {'band':>5} {'cas9':>7} {'weighted':>9} {'fid':>7} "
          f"{'w*fid':>8} {'vs bar':>7} {'k=1':>7} {'cov(7hk)':>9} {'slots for 78%':>14} {'s':>5}")
    out = []
    for w in WIDTHS:
        lo, hi = START, START + w - 1
        mf = AH._scaled_max_fail(cell, base.main_max_fail, w)
        cfg = _dc.replace(base, hdr_range=(lo, hi), main_max_fail=mf,
                          variants=(min(base.variants, AH.WIDE_WINDOW_VARIANTS)
                                    if w > 100 else base.variants))
        p = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
        if os.path.exists(p):
            os.remove(p)
        t0 = time.monotonic()
        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=800)
        dt = time.monotonic() - t0
        MT.free_gpu_memory()
        if not rows:
            print(f"  {w:>6} {mf:>4} {'-':>5} {str(meta.get('cas9_pool') or '-'):>7}"
                  f"   declined: {meta.get('reason')}  ({dt:.0f}s)")
            out.append({"width": w, "mf": mf, "reason": meta.get("reason")})
            continue
        s = score(rows, contract, reference, cell_types, seed=lo)
        band = meta["clean"]
        prod = s["weighted"] * s["fidelity"]
        cov = 1 - (1 - N_HOTKEYS * band / 900.0) ** 3
        # slots needed for a 78% chance that some sibling hits, at this band size
        need = None
        if band:
            need = int(round(900 * (1 - (1 - 0.78) ** (1 / 3)) / band))
        rec = {"width": w, "mf": mf, "band": band, "cas9": meta.get("cas9_pool"),
               "weighted": s["weighted"], "fidelity": s["fidelity"], "w_times_fid": prod,
               "k1_final": s["weighted"] * cons_k1 * s["fidelity"], "fleet_cov_7": cov,
               "slots_for_78pct": need, "build_s": round(dt, 1),
               "cells": meta.get("cells")}
        out.append(rec)
        print(f"  {w:>6} {mf:>4} {band:>5} {meta.get('cas9_pool'):>7} {s['weighted']:>9.2f} "
              f"{s['fidelity']:>7.4f} {prod:>8.1f} {prod/BAR:>6.2f}x "
              f"{rec['k1_final']:>7.2f} {cov:>8.1%} {str(need):>14} {dt:>5.0f}")
    good = [r for r in out if r.get("band")]
    if good:
        best = max(good, key=lambda r: r["w_times_fid"])
        print(f"\n  best quality: width {best['width']}  weighted*fid {best['w_times_fid']:.1f} "
              f"({best['w_times_fid']/BAR:.2f}x the bar)  k=1 {best['k1_final']:.2f}")
        ours = [r for r in good if r["width"] == 300]
        if ours:
            b = ours[0]
            print(f"  our shipped width 300: {b['w_times_fid']:.1f} "
                  f"({b['w_times_fid']/BAR:.2f}x)  k=1 {b['k1_final']:.2f}  cov {b['fleet_cov_7']:.1%}")
        print(f"\n  'slots for 78%' is how many disjoint windows of that width a coldkey needs for")
        print(f"  a 78% chance some sibling hits a band seed -- the 90-slot strategy's arithmetic.")
    json.dump({"task": TASK, "cell": cell, "start": START, "bar": BAR, "arms": out},
              open(os.getenv("NW_OUT", "narrow_width.json"), "w"), indent=1)
    print(f"\nwrote {os.getenv('NW_OUT', 'narrow_width.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
