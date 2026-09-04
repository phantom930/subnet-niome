#!/usr/bin/env python3
"""win20_epay3.py — E[pay] across every measured width, per cell type.

Generalises win20_epay2.py (which priced width 20 against the shipped width 100 for the erythroid
types) to arbitrary per-cell width lists, so width 16 can be priced for the erythroid types and
width 12 for HEK293.

Carries forward the four corrections that pricing needed, each of which changed the answer:

* Price per contract, not per cell-type mean -- payout is a step function, so E[pay(mean config)]
  is not mean E[pay].
* Enumerate the seed-draw patterns exactly rather than sampling rounds, and rank our own hotkeys
  against each other so two simultaneous placements cost each other a rank.
* Draw the off-band floor from its **measured** empirical distribution, **pooled per cell type**
  across widths. Welch tests found no difference by width (p = 0.97 / 0.16 / 0.51) and no mechanism
  predicts one -- the floor is a property of the 250-row structure and the seed. Using separate
  per-width floors injected a spurious asymmetry that alone moved the verdict by 6 points.
* Bootstrap with the resample **shared** across widths. configs() orders every width by the same
  sorted contract list, so row i is the same contract; drawing separate indices per width compares
  the configs on different samples and inflates the interval from [-3.9%, +5.9%] to [-48%, +98%].
"""
import json
import math
import os
import statistics as st
import urllib.request
from collections import defaultdict

import numpy as np

from win20_epay import OURS
from win20_epay2 import epay

SAMPLES = 400
PLAN = {"CD34+_HSPC": [100, 20, 16], "HUDEP-2": [100, 20, 16], "K562": [100, 20, 16],
        "HEK293": [100, 20, 12]}
BASELINE = 100
SOURCES = ["win20_cells.json", "win20_narrow.json", "win20_fill.json"]
FLOOR_SOURCES = ["win20_floor.json", "win20_floor2.json"]


def configs():
    recs, seen, uniq = [], set(), []
    for p in SOURCES:
        if os.path.exists(p):
            recs += json.load(open(p))
    for r in recs:
        k = (r["task_id"], r["width"], r["window"][0])
        if k in seen or "declined" in r:
            continue
        seen.add(k)
        uniq.append(r)
    out = {}
    for cell, widths in PLAN.items():
        per = {}
        for w in widths:
            g = [r for r in uniq if r["cell_type"] == cell and r["width"] == w]
            per[w] = {t: [x for x in g if x["task_id"] == t] for t in {x["task_id"] for x in g}}
        common = sorted(set.intersection(*(set(per[w]) for w in widths)))
        for w in widths:
            out[(cell, w)] = [
                {"task_id": t,
                 "band": st.mean(r["band_60k"] for r in per[w][t]),
                 "weighted": st.mean(r["total_weighted_score"] for r in per[w][t]),
                 "fid": st.mean(r["distribution_fidelity_factor"] for r in per[w][t])}
                for t in common]
    return out


def floors():
    """Pooled empirical floor per cell type, reused for every width of that cell type."""
    pool = defaultdict(list)
    for p in FLOOR_SOURCES:
        if os.path.exists(p):
            for r in json.load(open(p)):
                pool[r["cell_type"]] += r["floors"]
    return {(c, w): pool[c] for c, ws in PLAN.items() for w in ws}


def fields():
    tl = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks", timeout=120))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    meta = {}
    for t in tl:
        c = (t.get("content") or {}).get("contract") or {}
        seeds = [s for s in str(c.get("seed") or "").split(",") if s.strip()]
        if len(seeds) == 3 and c.get("cell_type"):
            meta[t["id"]] = c["cell_type"]
    sc = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/miners/scores?limit=40000", timeout=180))
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    best = defaultdict(dict)
    for x in sc:
        if x["task_id"] not in meta or x["miner_hotkey"] in OURS:
            continue
        cur = best[x["task_id"]].get(x["miner_hotkey"])
        if cur is None or x["final_score"] > cur:
            best[x["task_id"]][x["miner_hotkey"]] = x["final_score"]
    out = defaultdict(list)
    for tid, hk in best.items():
        if len(hk) >= 10:
            out[meta[tid]].append(sorted(hk.values(), reverse=True)[:12])
    return out


def main():
    cfg, fld, fl = configs(), fields(), floors()
    rng = np.random.default_rng(20260903)
    draws = rng.random((SAMPLES, 6))

    print("pooled off-band floor (measured, shared across widths within a cell type):")
    for cell in PLAN:
        v = fl[(cell, PLAN[cell][0])]
        print(f"  {cell:<11} mean {st.mean(v):.4f}  sd {st.stdev(v):.4f}  n={len(v)}")

    mat = {}
    for cell, widths in PLAN.items():
        for w in widths:
            rows = cfg[(cell, w)]
            mat[(cell, w)] = ([c["task_id"] for c in rows],
                              np.array([[epay(c["band"], c["weighted"], c["fid"], x,
                                              fl[(cell, w)], draws) for x in fld[cell]]
                                        for c in rows]))

    print(f"\n{'cell':<11} {'w':>4} {'n':>2} {'band':>6} {'weighted':>9} {'fid':>6} "
          f"{'E[pay]':>9} {'vs w100':>9} {'better/worse':>13} {'sign p':>7} {'95% CI':>18}")
    for cell, widths in PLAN.items():
        ids0, m0 = mat[(cell, BASELINE)]
        base = float(m0.mean())
        for w in widths:
            rows = cfg[(cell, w)]
            ids, m = mat[(cell, w)]
            e = float(m.mean())
            if w == BASELINE:
                extra = f"{'--':>9} {'':>13} {'':>7} {'':>18}"
            else:
                order = {t: i for i, t in enumerate(ids)}
                d = [(m[order[t]].mean() / m0[i].mean() - 1) * 100 for i, t in enumerate(ids0)]
                win = sum(1 for x in d if x > 0.5)
                los = sum(1 for x in d if x < -0.5)
                p = (sum(math.comb(win + los, i) for i in range(win, win + los + 1))
                     / 2 ** (win + los)) if win + los else 1.0
                boot = []
                for _ in range(4000):
                    ci = rng.integers(0, m.shape[0], m.shape[0])
                    fi = rng.integers(0, m.shape[1], m.shape[1])
                    boot.append((float(m[np.ix_(ci, fi)].mean())
                                 / float(m0[np.ix_(ci, fi)].mean()) - 1) * 100)
                b = np.sort(np.array(boot))
                extra = (f"{(e / base - 1) * 100:>+8.1f}% {f'{win}/{los}':>13} {p:>7.3f} "
                         f"{f'[{b[100]:+.1f}%, {b[3899]:+.1f}%]':>18}")
            print(f"{cell:<11} {w:>4} {len(rows):>2} {st.mean(c['band'] for c in rows):>6.2f} "
                  f"{st.mean(c['weighted'] for c in rows):>9.1f} "
                  f"{st.mean(c['fid'] for c in rows):>6.3f} {e:>9.5f} {extra}")


if __name__ == "__main__":
    main()
