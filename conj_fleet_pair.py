#!/usr/bin/env python3
"""conj_fleet_pair.py — paired comparison of two conj_fleet runs over the same tasks.

Pairs WITHIN task, which is the whole point: `total_weighted_score` moves 54% with the contract and
the field moves with it, so an unpaired mean of finals across tasks measures contract difficulty
rather than the arm. Every row here is one task's arm-A value against its own arm-B value.
"""
import json
import statistics as st
import sys

A = {r["task"]: r for r in json.load(open(sys.argv[1])) if "reason" not in r}
B = {r["task"]: r for r in json.load(open(sys.argv[2])) if "reason" not in r}
LA = sys.argv[3] if len(sys.argv) > 3 else "A"
LB = sys.argv[4] if len(sys.argv) > 4 else "B"
keys = [k for k in A if k in B]


def best(r):
    built = [h for h in r["hotkeys"] if "final" in h]
    return max(built, key=lambda h: h["final"]) if built else None


print(f"{'task':<10}{'cell':<12}{'arm':<12}{'built':>6}{'|band|':>7}{'bHit':>5}"
      f"{'|clean|':>8}{'cHit':>5}{'ownCln':>7}{'cons':>7}{'wxfid':>8}{'final':>8}"
      f"{'cut10':>7}{'vs':>7}")
rows = []
for k in keys:
    line = []
    for lbl, r in ((LA, A[k]), (LB, B[k])):
        b = best(r)
        f10 = r["field_top10"] or []
        cut = f10[9] if len(f10) >= 10 else float("nan")
        wxf = b["weighted"] * b["fidelity"] if b else 0.0
        fin = b["final"] if b else 0.0
        print(f"{r['task'] if lbl == LA else '':<10}{r['cell'] if lbl == LA else '':<12}"
              f"{lbl + ' g%d k%d' % (r['group'], r['k']):<12}{r['n_built']:>6}"
              f"{len(r['band_total']):>7}{len(r['band_hits']):>5}"
              f"{r['clean_total_n']:>8}{len(r['clean_hits']):>5}"
              f"{b['clean_n'] if b else 0:>7}{b['cons'] if b else 0:>7.3f}{wxf:>8.1f}"
              f"{fin:>8.1f}{cut:>7.1f}{fin / cut if cut == cut and cut else 0:>6.2f}x")
        line.append((fin, cut, r))
    rows.append((k, line[0], line[1]))
    print()

print("=" * 100)
print("PAIRED SUMMARY (within task)")
print("=" * 100)
for lbl, idx in ((LA, 1), (LB, 2)):
    fins = [r[idx][0] for r in rows]
    cuts = [r[idx][1] for r in rows]
    rats = [f / c for f, c in zip(fins, cuts) if c == c and c]
    pl = sum(1 for f, c in zip(fins, cuts) if c == c and f >= c)
    bh = sum(1 for r in rows if len(r[idx][2]["band_hits"]) >= 1)
    bt = st.mean(len(r[idx][2]["band_total"]) for r in rows)
    ct = st.mean(r[idx][2]["clean_total_n"] for r in rows)
    bu = st.mean(r[idx][2]["n_built"] for r in rows)
    print(f"  {lbl:<10} places {pl}/{len(rows)}   median final/cut10 {st.median(rats):.2f}x   "
          f"band>=1 {bh}/{len(rows)}   mean |band| {bt:.0f}  |clean| {ct:.0f}  built {bu:.1f}/10")
wins = sum(1 for r in rows if r[2][0] > r[1][0])
diffs = [r[2][0] - r[1][0] for r in rows]
print(f"\n  {LB} beats {LA} on {wins}/{len(rows)} tasks;  mean paired delta "
      f"{st.mean(diffs):+.1f} final points"
      + (f" (sd {st.stdev(diffs):.1f})" if len(diffs) > 1 else ""))
byc = {}
for k, a, b in rows:
    byc.setdefault(a[2]["cell"], []).append((a, b))
print("\n  by cell type (mean paired delta in final points):")
for cell, sub in sorted(byc.items()):
    d = [b[0] - a[0] for a, b in sub]
    print(f"    {cell:<12} n={len(sub)}  {st.mean(d):+8.1f}   "
          f"{LA} {st.mean(a[0] for a, _b in sub):>7.1f} -> {LB} {st.mean(b[0] for _a, b in sub):>7.1f}")
