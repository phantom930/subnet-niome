#!/usr/bin/env python3
"""conj_grid_report.py — read stage-A grids and show the marginals and the leaders per cell."""
import glob
import json
import statistics as st
import sys

V_CLEAN, V_REST = 0.26, 0.10          # pooled placeholder; stage B measures these per config

for path in sorted(sys.argv[1:] or glob.glob("conj_grid_*.json")):
    d = json.load(open(path))
    ok = [r for r in d if "wxfid" in r]
    if not ok:
        print(f"{path}: no built configs\n"); continue
    cell = ok[0]["cell"]
    bad = [r for r in d if "wxfid" not in r]
    for r in ok:
        other = ((r["clean"] - r["k"]) * V_CLEAN + (900 - r["clean"]) * V_REST) / (900 - r["k"])
        r["k1"] = r["wxfid"] * (1 + 2 * other) / 3
        r["pk1"] = 1 - (1 - r["k"] / 900) ** 3
    print(f"===== {cell}  ({path})  {len(ok)} built / {len(bad)} declined  "
          f"cells 8/8: {all(r['cells'] == 8 for r in ok)}  "
          f"band_cons 1.0: {all(abs(r['band_cons'] - 1) < 1e-9 for r in ok)} =====")
    for key in ("k", "group", "width", "light"):
        vals = sorted({r[key] for r in ok}, key=lambda x: (x is None, x))
        cells = []
        for v in vals:
            rs = [r for r in ok if r[key] == v]
            cells.append(f"{str(v)}:{st.mean(x['wxfid'] for x in rs):.1f}/"
                         f"{st.mean(x['clean'] for x in rs):.0f}")
        print(f"  {key:<6} (w x fid / clean)  " + "  ".join(cells))
    best = max(ok, key=lambda r: r["wxfid"])
    print(f"  BEST w x fid: k={best['k']} group={best['group']} width={best['width']} "
          f"light={best['light']} -> w {best['weighted']:.1f} x fid {best['fidelity']:.4f} "
          f"= {best['wxfid']:.1f}  clean {best['clean']}  cas9 {best['cas9']}  "
          f"modelled k=1 {best['k1']:.1f} at P {best['pk1']:.2%}")
    w = [r["wxfid"] for r in ok]
    print(f"  spread: w x fid {min(w):.1f}-{max(w):.1f} ({max(w)/min(w)-1:.1%})   "
          f"clean {min(r['clean'] for r in ok)}-{max(r['clean'] for r in ok)}\n")
