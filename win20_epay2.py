#!/usr/bin/env python3
"""win20_epay2.py — E[pay] for width 20 vs the shipped width 100, with the floor measured.

Supersedes the flat-floor pricing in win20_epay.py, which assumed 0.10 for a miss and was not
robust to it: sweeping the floor across the range these builds actually occupy moved the verdict
from -4.0% (at 0.086) to +4.1% (0.100) to -1.5% (0.120). The cause is structural — our k=1 score of
86-107 sits in the densest part of the field, against a rank-10 cutoff median of 74-84, so a small
shift crosses several rank steps.

win20_floor.py measured it instead: 38 (contract, width) cells, three off-band seeds each. The
floor is ~0.10-0.12, it is a per-seed *random variable* (sd 0.011-0.032, range 0.073-0.205), and it
is slightly **higher at width 100** on two of three cell types — a small advantage the flat 0.10
was not giving the incumbent.

Three things this pricing does that the flat version did not:

* **Draws the floor from its measured empirical distribution per (cell, width)** rather than using
  a point value, because payout is a step function and E[pay(f)] != pay(E[f]).
* **Pairs the random draws across widths**, so the delta is measured on common random numbers and
  its variance reflects the configs rather than the sampling.
* **Ignores k=0 hotkeys as non-paying**, which is exact here, not an approximation: their best
  possible score is weighted * 0.12 * fid ~ 32 against a rank-10 cutoff of 74-84 in a 247-miner
  field. Only hotkeys that hit at least one seed can place, so at most three scores per round need
  ranking.
"""
import json
import math
import os
import random
import statistics as st
import urllib.request
from collections import defaultdict

import numpy as np

DIST = [0.3, 0.2, 0.2, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
BURNING_RATE = 0.02
HOTKEYS = 9
SEED_SPACE = 900
SAMPLES = 400
CELLS = ["CD34+_HSPC", "HUDEP-2", "K562"]
WIDTHS = [100, 20]

from win20_epay import OURS, configs  # noqa: E402  (same contract pairing and field filtering)


def floors():
    """Empirical off-band consistency draws per (cell type, width)."""
    recs = json.load(open("win20_floor.json"))
    out = defaultdict(list)
    for r in recs:
        out[(r["cell_type"], r["width"])] += r["floors"]
    return out


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


def _pay_vec(scores, cut):
    """Payout share for arrays of our own scores, ranked against one field.

    ``scores`` is a list of equal-length arrays, one per hotkey that holds at least one band seed.
    A hotkey's rank counts the field entries above it plus our own other hotkeys above it, which is
    what makes two simultaneous placements cost each other a rank.
    """
    asc = np.sort(np.asarray(cut))                      # ascending, for searchsorted
    m = len(asc)
    table = np.array(DIST + [0.0]) * (1 - BURNING_RATE)
    total = np.zeros_like(scores[0])
    for i, s in enumerate(scores):
        above = m - np.searchsorted(asc, s, side="right")     # field entries strictly greater
        for j, o in enumerate(scores):
            if j != i:
                above = above + (o > s)
        rank = above + 1
        total = total + table[np.minimum(rank, 11) - 1]
    return total


def epay(band, weighted, fid, field, fl, draws):
    """Expected payout share per round, vectorised over the sampled floor draws.

    Only hotkeys holding at least one band seed can place, so each pattern ranks at most three
    scores: a k=0 hotkey tops out near weighted * 0.12 * fid ~ 32 against a rank-10 cutoff of
    74-84. ``draws`` supplies unit uniforms shared across widths (common random numbers).
    """
    p = band / SEED_SPACE
    a = HOTKEYS * p
    q = 1.0 - a
    fla = np.asarray(fl)
    idx = (draws * len(fla)).astype(np.int64)
    f = fla[idx]                                        # (SAMPLES, 6) floor draws

    def sc(k, cols):
        miss = f[:, cols].sum(axis=1) if cols else 0.0
        return weighted * ((k + miss) / 3) * fid

    patterns = [
        (3 * q * q * a,         lambda: [sc(1, [0, 1])]),
        (3 * q * a * a / 9,     lambda: [sc(2, [0])]),
        (3 * q * a * a * 8 / 9, lambda: [sc(1, [0, 1]), sc(1, [2, 3])]),
        (a ** 3 / 81,           lambda: [sc(3, [])]),
        (a ** 3 * 24 / 81,      lambda: [sc(2, [0]), sc(1, [1, 2])]),
        (a ** 3 * 56 / 81,      lambda: [sc(1, [0, 1]), sc(1, [2, 3]), sc(1, [4, 5])]),
    ]
    total = 0.0
    for weight, build in patterns:
        if weight < 1e-12:
            continue
        cols = build()
        if np.isscalar(cols[0]):
            cols = [np.full(len(draws), c) for c in cols]
        total += weight * float(_pay_vec(cols, field).mean())
    return total


def matrix(cfg, fld, fl, draws):
    """E[pay] for every (cell, width, contract, field) once; the bootstrap then just resamples."""
    out = {}
    for cell in CELLS:
        for w in WIDTHS:
            rows = cfg[(cell, w)]
            out[(cell, w)] = (
                [c["task_id"] for c in rows],
                np.array([[epay(c["band"], c["weighted"], c["fid"], x, fl[(cell, w)], draws)
                           for x in fld[cell]] for c in rows]))
    return out


def main():
    cfg, fld, fl = configs(), fields(), floors()
    rng = np.random.default_rng(20260903)
    draws = rng.random((SAMPLES, 6))

    print("measured off-band floor (empirical draws, not an assumption):")
    for cell in CELLS:
        for w in WIDTHS:
            v = fl[(cell, w)]
            print(f"  {cell:<11} w{w:<4} mean {st.mean(v):.4f}  sd {st.stdev(v):.4f}  n={len(v)}")

    mat = matrix(cfg, fld, fl, draws)
    print(f"\n{'cell':<11} {'w':>4} {'n':>2} {'fields':>7} {'band':>6} {'weighted':>9} "
          f"{'fid':>6} {'E[pay]':>9} {'vs w100':>9}")
    agg = defaultdict(float)
    for cell in CELLS:
        base = None
        for w in WIDTHS:
            rows = cfg[(cell, w)]
            _ids, m = mat[(cell, w)]
            e = float(m.mean())
            agg[w] += e
            rel = "--" if base is None else f"{(e / base - 1) * 100:+8.1f}%"
            if base is None:
                base = e
            print(f"{cell:<11} {w:>4} {len(rows):>2} {m.shape[1]:>7} "
                  f"{st.mean(c['band'] for c in rows):>6.2f} "
                  f"{st.mean(c['weighted'] for c in rows):>9.1f} "
                  f"{st.mean(c['fid'] for c in rows):>6.3f} {e:>9.5f} {rel:>9}")
    delta = (agg[20] / agg[100] - 1) * 100
    print(f"\n{'AGGREGATE':<11} w100 {agg[100]:>9.5f}   w20 {agg[20]:>9.5f}   {delta:+.1f}%")

    wins = losses = flat = 0
    detail = []
    for cell in CELLS:
        ids100, m100 = mat[(cell, 100)]
        ids20, m20 = mat[(cell, 20)]
        order = {t: i for i, t in enumerate(ids20)}
        for i, t in enumerate(ids100):
            a, b = m100[i].mean(), m20[order[t]].mean()
            r = (b / a - 1) * 100 if a else 0.0
            detail.append(r)
            wins += r > 0.5
            losses += r < -0.5
            flat += -0.5 <= r <= 0.5
    k, n = wins, wins + losses
    p = sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n if n else 1.0
    print(f"per contract: {wins} better / {losses} worse / {flat} flat   sign test p = {p:.3f}")
    print(f"  [{', '.join(f'{x:+.0f}%' for x in sorted(detail))}]")

    print("\nbootstrap over contracts and fields (4000 resamples):")
    # The resample must be SHARED across widths: configs() orders both widths by the same sorted
    # contract list, so row i is the same contract in each. Drawing separate indices per width
    # compares the two configs on different samples and inflates the interval to uselessness
    # ([-48%, +98%] rather than the paired spread).
    boot = []
    for _ in range(4000):
        tot = defaultdict(float)
        for cell in CELLS:
            shape = mat[(cell, 100)][1].shape
            ci = rng.integers(0, shape[0], shape[0])
            fi = rng.integers(0, shape[1], shape[1])
            for w in WIDTHS:
                tot[w] += float(mat[(cell, w)][1][np.ix_(ci, fi)].mean())
        boot.append((tot[20] / tot[100] - 1) * 100)
    boot = np.sort(np.array(boot))
    print(f"  mean {boot.mean():+.1f}%   95% CI [{boot[100]:+.1f}%, {boot[3899]:+.1f}%]"
          f"   P(w20 better) = {(boot > 0).mean() * 100:.1f}%")


if __name__ == "__main__":
    main()
