#!/usr/bin/env python3
"""conj_fleet_report.py — the detail tables and hit rates from conj_fleet.py."""
import json
import statistics as st
import sys

rows = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "conj_fleet.json"))
rows = [r for r in rows if "reason" not in r]
NHK = max((len(r["hotkeys"]) for r in rows), default=10)

print("=" * 118)
print("PER-TASK SUMMARY")
print("=" * 118)
print(f"{'task':<10}{'cell':<12}{'seeds':<20}{'win':>4}{'g':>5}{'k':>3}"
      f"{'built':>7}{'dist':>5} | {'|band|':>7}{'bHit':>5}{'|clean|':>8}{'cHit':>5} | "
      f"{'best':>7}{'cut10':>7}{'top1':>7}{'rank':>6}")
for r in rows:
    f10 = r["field_top10"] or []
    cut = f10[9] if len(f10) >= 10 else float("nan")
    top = f10[0] if f10 else float("nan")
    bf = r["best_final"] or 0.0
    rank = (sum(1 for x in f10 if x > bf) + 1) if f10 else float("nan")
    print(f"{r['task']:<10}{r['cell']:<12}{str(r['seeds']):<20}"
          f"{r['n_windows']}+{r['padded']:<2}{r['group']:>5}{r['k']:>3}"
          f"{r['n_built']:>7}{r['distinct_builds']:>5} | "
          f"{len(r['band_total']):>7}{len(r['band_hits']):>5}"
          f"{r['clean_total_n']:>8}{len(r['clean_hits']):>5} | "
          f"{bf:>7.1f}{cut:>7.1f}{top:>7.1f}{rank:>6}")

print()
print("=" * 118)
print("HITS PER TASK  (which of the round's three seeds the fleet actually caught)")
print("=" * 118)
print(f"{'task':<10}{'cell':<12}{'seeds':<20} | {'band hits':<22}{'n':>3} | "
      f"{'clean hits':<22}{'n':>3}")
for r in rows:
    print(f"{r['task']:<10}{r['cell']:<12}{str(r['seeds']):<20} | "
          f"{str(r['band_hits']):<22}{len(r['band_hits']):>3} | "
          f"{str(r['clean_hits']):<22}{len(r['clean_hits']):>3}")

print()
print("=" * 118)
print("PER-HOTKEY DETAIL")
print("=" * 118)
for r in rows:
    print(f"-- {r['task']}  {r['cell']}  seeds {r['seeds']}  classes {','.join(r['classes'])}"
          f"  (real windows {r['n_windows']}, padded +{r['padded']})")
    print(f"   {'hk':>3}{'|band|':>7}{'|clean|':>8}{'bHit':>5}{'cHit':>5}"
          f"{'weighted':>10}{'fid':>7}{'cons':>7}{'final':>8}  band")
    for h in r["hotkeys"]:
        if "final" not in h:
            print(f"   {h['hk']:>3}   DECLINED: {h.get('reason')}")
            continue
        print(f"   {h['hk']:>3}{len(h['band']):>7}{h['clean_n']:>8}"
              f"{len(h['band_hits']):>5}{len(h['clean_hits']):>5}"
              f"{h['weighted']:>10.1f}{h['fidelity']:>7.4f}{h['cons']:>7.4f}{h['final']:>8.1f}"
              f"  {h['band']}")
    print()

print("=" * 118)
print("HIT RATES  (fraction of tasks whose 10-hotkey UNION caught >=1 / >=2 / 3 of the 3 seeds)")
print("=" * 118)
n = len(rows)
for label, key in (("band", "band_hits"), ("clean", "clean_hits")):
    s1 = sum(1 for r in rows if len(r[key]) >= 1) / n
    s2 = sum(1 for r in rows if len(r[key]) >= 2) / n
    s3 = sum(1 for r in rows if len(r[key]) >= 3) / n
    print(f"  {label:<6} single {s1:>6.0%}   double {s2:>6.0%}   triple {s3:>6.0%}   (n={n})")
by = {}
for r in rows:
    by.setdefault(r["cell"], []).append(r)
print("\n  by cell type:")
for cell, sub in sorted(by.items()):
    m = len(sub)
    b1 = sum(1 for r in sub if len(r["band_hits"]) >= 1) / m
    b2 = sum(1 for r in sub if len(r["band_hits"]) >= 2) / m
    c1 = sum(1 for r in sub if len(r["clean_hits"]) >= 1) / m
    c2 = sum(1 for r in sub if len(r["clean_hits"]) >= 2) / m
    print(f"    {cell:<12} n={m}  band {b1:>5.0%}/{b2:>5.0%}   clean {c1:>5.0%}/{c2:>5.0%}")

print()
print("=" * 118)
print("SCORE vs THE FIELD  (best of the 10 hotkeys, against the field that played that contract)")
print("=" * 118)
placed = [r for r in rows if r["field_top10"] and len(r["field_top10"]) >= 10
          and (r["best_final"] or 0) >= r["field_top10"][9]]
won = [r for r in rows if r["field_top10"] and (r["best_final"] or 0) >= r["field_top10"][0]]
have = [r for r in rows if r["field_top10"] and len(r["field_top10"]) >= 10]
print(f"  places in the top 10: {len(placed)}/{len(have)}")
print(f"  beats the field's rank 1: {len(won)}/{len(have)}")
if have:
    print(f"  median best/cut10 ratio: "
          f"{st.median((r['best_final'] or 0) / r['field_top10'][9] for r in have):.2f}x")
    print(f"  median best/top1  ratio: "
          f"{st.median((r['best_final'] or 0) / r['field_top10'][0] for r in have):.2f}x")
