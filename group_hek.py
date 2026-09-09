#!/usr/bin/env python3
"""group_hek.py — group_size on HEK293, against the quality bar the top block clears.

Ranks 7-10 on HEK293 a8b9f1bb cleared the 81.92 cutoff on a SINGLE band hit, holding weighted ~221
and fidelity ~0.962. At the k=1 consistency of 0.3853 the bar is

    weighted * fidelity >= 81.92 / 0.3853 = 212.6

Window width is now measured and is NOT the lever: weighted*fid is flat at 177-185 across widths
10-300 (narrow_width.py), because narrowing the window grows the band toward the window size rather
than shrinking it. `group_size` is the remaining measured lever with a documented effect on both
terms, and it pulls them in OPPOSITE directions:

  - up:   more Cas12a rows -> higher cas-coverage entropy -> higher fidelity
          (CLAUDE.md's HEK293 sweep: 0.786/0.854/0.886/0.924/0.939 at group 20/40/60/80/100)
  - down: more Cas9 rows, which score better structurally -> higher weighted
          (weighted 326 -> 308 across that same sweep, and K562 258 at g42 vs 250 at g100)

So there is an optimum in the product, and group 80 was chosen on a payout basis that predates
tonight's cas9_cell_target fix -- which changed the weighted term it was balancing against.

`bank_key` excludes group_size, so every arm shares one Cas12a bank and only the first pays for it.

    python group_hek.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "grphek")

import dataclasses as _dc     # noqa: E402
import json                   # noqa: E402
import logging                # noqa: E402
import statistics as st       # noqa: E402
import time                   # noqa: E402
from collections import Counter  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G            # noqa: E402
from niome_subnet.genomics import all_cut as AC      # noqa: E402
from niome_subnet.genomics import all_hdr as AH      # noqa: E402
from niome_subnet.genomics import mt19937 as MT      # noqa: E402
from niome_subnet.utils import settings              # noqa: E402
from sd_task import score, task_content              # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "a8b9f1bb-f25d-42af-92a1-e0aff16bee3a"
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "100-399"
GROUPS = [int(x) for x in (os.getenv("GH_GROUPS") or "42,60,80,100,125,150").split(",")]
N_HOTKEYS = 7
BAR = 212.6
FLOOR_CONS = 0.078


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = (int(x) for x in WINDOW.split("-"))
    span = hi - lo + 1
    base = AH.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    cons_k1 = (1 + 2 * FLOOR_CONS) / 3
    print(f"task {TASK[:8]}  {cell}  window {lo}-{hi} (span {span})  "
          f"mf {AH._scaled_max_fail(cell, base.main_max_fail, span)}")
    print(f"k=1 consistency {cons_k1:.4f}; bar for 81.92 is weighted*fid >= {BAR:.1f}\n")
    print(f"  {'group':>6} {'cas mix':>9} {'target':>7} {'band':>5} {'cas9':>7} {'mean_w':>7} "
          f"{'struct':>7} {'weighted':>9} {'fid':>7} {'w*fid':>8} {'vs bar':>7} {'k=1':>7} "
          f"{'cov':>6} {'s':>5}")
    out = []
    for g in GROUPS:
        cfg = _dc.replace(base, hdr_range=(lo, hi), group_size=g,
                          main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span),
                          variants=(min(base.variants, AH.WIDE_WINDOW_VARIANTS)
                                    if span > 100 else base.variants))
        want = 250 - g
        target = AC.cas9_cell_target(contract, ctx, cfg, want)
        t0 = time.monotonic()
        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=800)
        dt = time.monotonic() - t0
        MT.free_gpu_memory()
        if not rows:
            print(f"  {g:>6} {f'{g}/{want}':>9} {target:>7} {'-':>5} "
                  f"{str(meta.get('cas9_pool') or '-'):>7}   declined: {meta.get('reason')}")
            out.append({"group": g, "reason": meta.get("reason")})
            continue
        s = score(rows, contract, reference, cell_types, seed=lo)
        feats = [d["features"] for d in json.load(open(settings.VALID_EXPERIMENTS_PATH))]
        mw = contract.get("mutation_weights", {})
        band = meta["clean"]
        prod = s["weighted"] * s["fidelity"]
        cov = 1 - (1 - N_HOTKEYS * band / 900.0) ** 3
        rec = {"group": g, "cas_mix": f"{g}/{want}", "cell_target": target, "band": band,
               "cas9": meta.get("cas9_pool"),
               "mean_weight": st.mean(mw.get(r["mutation"], 1.0) for r in rows),
               "structural": st.mean(0.625 * f["gc_score"] + 0.375 * f["dist_score"]
                                     for f in feats),
               "weighted": s["weighted"], "fidelity": s["fidelity"], "w_times_fid": prod,
               "k1_final": s["weighted"] * cons_k1 * s["fidelity"], "fleet_cov_7": cov,
               "cells": meta.get("cells"), "build_s": round(dt, 1)}
        out.append(rec)
        print(f"  {g:>6} {rec['cas_mix']:>9} {target:>7} {band:>5} {meta.get('cas9_pool'):>7} "
              f"{rec['mean_weight']:>7.4f} {rec['structural']:>7.4f} {s['weighted']:>9.2f} "
              f"{s['fidelity']:>7.4f} {prod:>8.1f} {prod/BAR:>6.2f}x {rec['k1_final']:>7.2f} "
              f"{cov:>5.1%} {dt:>5.0f}")
    good = [r for r in out if r.get("band")]
    if good:
        bq = max(good, key=lambda r: r["w_times_fid"])
        bk = max(good, key=lambda r: r["fleet_cov_7"] * r["k1_final"])
        print(f"\n  best quality   : group {bq['group']}  w*fid {bq['w_times_fid']:.1f} "
              f"({bq['w_times_fid']/BAR:.2f}x)  k=1 {bq['k1_final']:.2f}  band {bq['band']}")
        print(f"  best freq*value: group {bk['group']}  "
              f"{bk['fleet_cov_7']*bk['k1_final']:.2f}  (band {bk['band']}, cov "
              f"{bk['fleet_cov_7']:.1%})")
        shipped = [r for r in good if r["group"] == 80]
        if shipped:
            b = shipped[0]
            print(f"  shipped group 80: w*fid {b['w_times_fid']:.1f} "
                  f"({b['w_times_fid']/BAR:.2f}x)  k=1 {b['k1_final']:.2f}  "
                  f"freq*value {b['fleet_cov_7']*b['k1_final']:.2f}")
        print(f"\n  fidelity should rise and weighted fall as group grows; the product is what")
        print(f"  matters, and k=1 must beat 81.92 to place at all on this round.")
    json.dump({"task": TASK, "cell": cell, "window": [lo, hi], "bar": BAR, "arms": out},
              open(os.getenv("GH_OUT", "group_hek.json"), "w"), indent=1)
    print(f"\nwrote {os.getenv('GH_OUT', 'group_hek.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
