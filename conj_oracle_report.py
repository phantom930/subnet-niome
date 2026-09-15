#!/usr/bin/env python3
"""conj_oracle_report.py — hit rates per (group, k) from conj_oracle.py, and what predicts them.

Three things it is careful about, because each would otherwise produce a confident wrong answer:

  * **2-window tasks are a different problem.** When two of the three seeds share a 100-seed class
    the oracle space is 200 seeds rather than 300, so the same band covers 1.5x the fraction. They
    are 30% of tasks (chance, for 9 classes) and are reported split out, never pooled silently.
  * **A config that DECLINES is not a config that misses.** `hotkeys_built` counts how many of the
    10 hotkeys could form a group at all; an arm that built on 0 is excluded from its rate rather
    than scored as a zero, and the count is reported so a thin arm cannot masquerade as a good one.
  * **Per-task argmax is winner's-curse.** With 12 arms, 3 seeds and heavy ties, "the best config
    for this task" is mostly noise. The per-feature section reports the aggregate rate within each
    feature bucket instead, and states the bucket size next to it.
"""
import json
import os
import statistics as st
import sys
from collections import defaultdict

PATH = sys.argv[1] if len(sys.argv) > 1 else "conj_oracle_HEK293.json"
rows = json.load(open(PATH))


def arms_of(r):
    return {(a["group"], a["k"]): a for a in r["arms"]}


def rates(subset, key):
    """-> {(g,k): (n_built, single, double, triple, mean_total, mean_hotkeys)}"""
    acc = defaultdict(list)
    for r in subset:
        for (g, k), a in arms_of(r).items():
            if a["hotkeys_built"] == 0:
                continue
            acc[(g, k)].append((a[f"{key}_hits"], a[f"{key}_total"], a["hotkeys_built"]))
    out = {}
    for gk, v in acc.items():
        n = len(v)
        out[gk] = (n,
                   sum(1 for h, _t, _b in v if h >= 1) / n,
                   sum(1 for h, _t, _b in v if h >= 2) / n,
                   sum(1 for h, _t, _b in v if h >= 3) / n,
                   st.mean(t for _h, t, _b in v),
                   st.mean(b for _h, _t, b in v))
    return out


def table(subset, label):
    print(f"\n===== {label}  (n = {len(subset)} tasks)")
    b, c = rates(subset, "band"), rates(subset, "clean")
    print(f"{'g':>4}{'k':>3}{'arms':>6}{'hk':>5} | {'band>=1':>8}{'band>=2':>8}{'band=3':>8}"
          f"{'|band|':>8} | {'cln>=1':>8}{'cln>=2':>8}{'cln=3':>7}{'|clean|':>9}")
    for gk in sorted(b):
        n, s1, d1, t1, tot, hk = b[gk]
        _n, s2, d2, t2, tot2, _ = c[gk]
        print(f"{gk[0]:>4}{gk[1]:>3}{n:>6}{hk:>5.1f} | {s1:>7.0%}{d1:>8.0%}{t1:>8.0%}{tot:>8.0f}"
              f" | {s2:>7.0%}{d2:>8.0%}{t2:>7.0%}{tot2:>9.0f}")
    if b:
        best = max(b, key=lambda gk: (b[gk][1], b[gk][2], c[gk][1]))
        print(f"  best by band>=1, then band>=2, then clean>=1:  group {best[0]}  k {best[1]}"
              f"   band>=1 {b[best][1]:.0%}  band>=2 {b[best][2]:.0%}  clean>=1 {c[best][1]:.0%}")
    return b, c


ball, call = table(rows, "ALL TASKS")
w3 = [r for r in rows if r["n_windows"] == 3]
w2 = [r for r in rows if r["n_windows"] == 2]
if w3:
    table(w3, "3 distinct seed windows (span 300)")
if w2:
    table(w2, "2 distinct seed windows (span 200)")

# --- chance baselines: is a band of |B| over `span` seeds doing better than |B| random seeds? ---
print("\n===== against chance (hypergeometric: |B| of `span`, 3 seeds drawn)")
print(f"{'g':>4}{'k':>3} | {'|band|':>7}{'span':>6} | {'obs>=1':>8}{'exp>=1':>8} | "
      f"{'obs|clean|':>11}{'exp>=1':>8}{'obs>=1':>8}")
for gk in sorted(ball):
    g, k = gk
    ex_b, ex_c, ob_b, ob_c, spans = [], [], [], [], []
    for r in rows:
        a = arms_of(r).get(gk)
        if not a or a["hotkeys_built"] == 0:
            continue
        n = r["span"]
        for tot, ex, ob, hits in ((a["band_total"], ex_b, ob_b, a["band_hits"]),
                                  (a["clean_total"], ex_c, ob_c, a["clean_hits"])):
            p_none = 1.0
            for i in range(3):
                p_none *= max(0.0, (n - tot - i)) / (n - i)
            ex.append(1 - p_none)
            ob.append(1 if hits >= 1 else 0)
        spans.append(n)
    if ob_b:
        print(f"{g:>4}{k:>3} | {st.mean(a['band_total'] for a in [arms_of(r)[gk] for r in rows if gk in arms_of(r) and arms_of(r)[gk]['hotkeys_built']]):>7.0f}"
              f"{st.mean(spans):>6.0f} | {st.mean(ob_b):>7.0%}{st.mean(ex_b):>8.0%} | "
              f"{st.mean(a['clean_total'] for a in [arms_of(r)[gk] for r in rows if gk in arms_of(r) and arms_of(r)[gk]['hotkeys_built']]):>11.0f}"
              f"{st.mean(ex_c):>8.0%}{st.mean(ob_c):>8.0%}")

# --- band overlap across the 10 hotkeys: how much did the stride actually decorrelate? ---
print("\n===== band overlap across the 10 hotkeys (distinct band seeds vs 10*k slots)")
print(f"{'g':>4}{'k':>3} | {'slots':>6}{'distinct':>9}{'overlap':>9}")
for gk in sorted(ball):
    g, k = gk
    vals = [(a["hotkeys_built"] * k, a["band_total"]) for a in
            (arms_of(r).get(gk) for r in rows) if a and a["hotkeys_built"]]
    if vals:
        sl = st.mean(v[0] for v in vals); ds = st.mean(v[1] for v in vals)
        print(f"{g:>4}{k:>3} | {sl:>6.0f}{ds:>9.0f}{1 - ds / sl:>8.0%}")

# --- features: aggregate rate within each bucket, with bucket size ---
print("\n===== by task feature (aggregate band>=1 within the bucket, not a per-task argmax)")
gks = sorted(ball, key=lambda gk: (-ball[gk][1], -ball[gk][2]))[:4]
feats = {
    "n_windows": lambda r: r["n_windows"],
    "min-gap": lambda r: "<=50" if min(b - a for a, b in zip(r["seeds"], r["seeds"][1:])) <= 50
                         else ">50",
    "bank": lambda r: "300k(capped)" if r["bank"] >= 300000 else "<300k",
}
for fname, fn in feats.items():
    buckets = defaultdict(list)
    for r in rows:
        buckets[fn(r)].append(r)
    print(f"  -- {fname}")
    for bk in sorted(buckets, key=str):
        sub = buckets[bk]
        bb = rates(sub, "band")
        line = "  ".join(f"g{g}k{k} {bb[(g, k)][1]:.0%}" for (g, k) in gks if (g, k) in bb)
        print(f"     {str(bk):<14} n={len(sub):<3} {line}")
