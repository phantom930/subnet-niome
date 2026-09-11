#!/usr/bin/env python3
"""floor_price.py -- what an off-band consistency floor is worth to the fleet, priced against the
fields the field's own 2026-09-08 shift produced.

Why this exists. Until 2026-09-08 our off-band per-seed value (0.101) was the field median, and
CLAUDE.md recorded "nobody has a better floor". Then 205 of 252 hotkeys lifted their median
erythroid `consistency_factor` to 0.1857 on one shared deterministic build (153 of 161 rows tied at
weighted 279.1), paying 2.5% of weighted and 5.0% of fidelity for +72% consistency. Our fleet did
not move, so the floor is now the one term where we are measurably behind the whole field.

The floor pays through the SPIKE, not through itself. A round with k band hits scores
`cons = (k + (3-k)*f) / 3`, so f enters every round we place in. At f = 0.101 a k=1 round is
cons 0.401; at f = 0.237 (all-cut's measured cut-clean seed value) it is 0.491 -- 22% more final
score on exactly the rounds that decide payout. The floor rows themselves never place: the shared
build scores 279.1 * 0.874 * 0.198 = 48.3 against a post-shift rank-10 cutoff median of 95.5.

Three things this prices, because they are not the same question:

  A  floor alone, at our own weighted x fidelity (225.8) -- the value of a mechanism that lifts f
     and changes nothing else. This is the sensitivity curve: it says how much a floor mechanism
     would have to buy to be worth building, BEFORE asking whether one exists.
  B  floor plus the shared build's weighted x fidelity (244.0, +8.1%) -- if whatever lifts f also
     moves row composition the way theirs does.
  C  the shared build itself, no band (f flat, k always 0) -- what 65% of the field now runs, so
     that A and B are priced against the real alternative rather than against nothing.

The floor is CORRELATED across siblings and the band is not. A cut-clean set is a property of the
contract, so every hotkey running the same rule agrees on which seeds are clean: on a round where
no band hits, all eleven hold one score and stack ranks r..r+10. That costs nothing here (floor
rows do not place) but it is modelled rather than assumed, because it is exactly the effect that
inverted the all-cut answer in fleet_price.py.

Fields are post-2026-09-08 erythroid only. Mixing in pre-shift fields prices a field that no
longer exists -- the same error as pooling across contracts, in the time axis. n is small (~23
rounds); the per-cell splits are reported so the reader can see it.
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np

DIST = np.array([0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01])
D = ("/tmp/claude-0/-root-workspace-subnet-niome/"
     "0f021db5-a874-4419-a235-ebcda9baabeb/scratchpad/")
SCORES = os.environ.get("FLP_SCORES", D + "s2.json")
TASKS = os.environ.get("FLP_TASKS", D + "t2.json")
OURS = os.environ.get("FLP_OURS", D + "ours3.json")
SINCE = os.environ.get("FLP_SINCE", "2026-09-08")
TRIALS = int(os.environ.get("FLP_TRIALS", "5000"))
ERY = ("K562", "HUDEP-2", "CD34+_HSPC")

# --- measured, all from the live fleet and the live feed --------------------------------
# PER CELL TYPE, because the contract moves weighted by >50% and pooling applies one cell's
# product to another's rounds -- the same error as pooling fields across contracts, and it is
# decisive here: a k=1 round is 119.6 at HUDEP-2's 298.3 and 88.3 at K562's 220.1, either side
# of most cutoffs. Medians over post-2026-09-08 rows, identified BY HOTKEY ADDRESS (uids are
# recycled: 7 of our 11 uids carried another operator's hotkey during the 09-08 reboot window,
# which contaminated a first pass at these numbers by 57 rows of 228).
OUR_WXF_CELL = {"HUDEP-2": 298.3, "CD34+_HSPC": 233.3, "K562": 220.1}
FIELD_WXF = 244.0     # the shared build: weighted 279.1 x fidelity 0.8739
BAND = 12             # all-HDR band, erythroid, group 80
JOINED = 300          # the joined space window_plan.py confines the fleet to
OUR_F = 0.101         # our measured off-band per-seed value
FIELD_F = 0.198       # the shared build's round consistency, = its per-seed floor (it has no band)
NHK = 11
SEED_LO, SEED_HI = 100, 999
NSEED = SEED_HI - SEED_LO + 1


def fields():
    sc = json.load(open(SCORES))
    tk = json.load(open(TASKS))
    ours = set(json.load(open(OURS)).values())   # ss58 addresses, not uids
    meta = {}
    for t in tk:
        tid = t.get("task_id") or t.get("id")
        c = (t.get("content") or {}).get("contract") or {}
        s = str(c.get("seed", "") or "")
        meta[tid] = (t.get("cell_type") or c.get("cell_type"),
                     [int(x) for x in s.split(",") if x.strip().isdigit()],
                     t.get("created_at") or "")
    best = defaultdict(dict)
    for r in sc:
        cell, seeds, created = meta.get(r["task_id"], (None, [], ""))
        if cell not in ERY or len(seeds) != 3 or created[:10] < SINCE:
            continue
        if r["miner_hotkey"] in ours:
            continue
        hk = r["miner_hotkey"]
        prev = best[r["task_id"]].get(hk)
        if prev is None or (r.get("final_score") or 0) > prev:
            best[r["task_id"]][hk] = r.get("final_score") or 0.0
    out = []
    for tid, v in best.items():
        f = sorted(v.values(), reverse=True)
        if len(f) >= 10 and f[0] > 0:
            out.append((tid, meta[tid][0], meta[tid][2][:16], np.asarray(f)))
    out.sort(key=lambda x: x[2])
    return out


def shares(field, mine):
    """Curve share the fleet captures; ties take consecutive ranks as the validator's sort does."""
    total = 0.0
    for s in mine:
        rank = int(np.sum(field > s)) + 1 + sum(1 for o in mine if o > s)
        if rank <= 10:
            total += DIST[rank - 1]
    return total


def simulate(field, f, wxf, band, rng, trials, nhk=NHK):
    """band=0 -> the flat shared build (correlated, no band). Otherwise nhk decorrelated bands
    over one joined window, with the floor shared across siblings."""
    tot = 0.0
    for _ in range(trials):
        seeds = rng.integers(SEED_LO, SEED_HI + 1, 3)
        if not band:
            tot += shares(field, [f * wxf] * nhk)
            continue
        lo = rng.integers(SEED_LO, SEED_HI - JOINED + 2)
        joined = np.arange(lo, lo + JOINED)
        mine = []
        for _h in range(nhk):
            b = set(rng.choice(joined, band, replace=False).tolist())
            k = sum(1 for s in seeds if int(s) in b)
            mine.append((k + (3 - k) * f) / 3 * wxf)
        tot += shares(field, mine)
    return tot / trials


def run(fs, sweep):
    """Every arm is priced with the ROUND'S OWN cell weighted x fidelity."""
    out = {}
    for f in sweep:
        for tag, wxf_of in (("A", lambda c: OUR_WXF_CELL[c]), ("B", lambda c: FIELD_WXF)):
            vals = []
            for _tid, cell, _d, fld in fs:
                rng = np.random.default_rng(20260910)
                vals.append(simulate(fld, f, wxf_of(cell), BAND, rng, TRIALS))
            out[(tag, f)] = np.asarray(vals)
    vals = []
    for _tid, cell, _d, fld in fs:
        rng = np.random.default_rng(20260910)
        vals.append(simulate(fld, FIELD_F, FIELD_WXF, 0, rng, TRIALS))
    out[("C", None)] = np.asarray(vals)
    return out


def main():
    fs = fields()
    cells = defaultdict(int)
    for _t, c, _d, _f in fs:
        cells[c] += 1
    print("fields: %d post-%s erythroid rounds (%s), our 11 hotkeys excluded by ss58 address"
          % (len(fs), SINCE, ", ".join("%s %d" % kv for kv in sorted(cells.items()))))
    cut = np.array([f[9] for _t, _c, _d, f in fs])
    print("rank-10 cutoff: min %.1f p25 %.1f med %.1f p75 %.1f max %.1f"
          % (cut.min(), np.percentile(cut, 25), np.median(cut),
             np.percentile(cut, 75), cut.max()))
    print("our weighted x fidelity, per cell: %s"
          % "  ".join("%s %.1f" % kv for kv in sorted(OUR_WXF_CELL.items())))
    print("shared build %.1f | band %d in a %d joined space | %d trials per round per arm\n"
          % (FIELD_WXF, BAND, JOINED, TRIALS))

    print("=== round score by k, at each floor, per cell against that cell's OWN fields ===")
    for cell in sorted(OUR_WXF_CELL):
        cc = np.array([f[9] for _t, c, _d, f in fs if c == cell])
        if not len(cc):
            continue
        print("  %s (w x fid %.1f, n=%d fields, cutoff med %.1f)"
              % (cell, OUR_WXF_CELL[cell], len(cc), np.median(cc)))
        for f in (OUR_F, 0.15, FIELD_F, 0.237, 0.30):
            row = []
            for k in (0, 1, 2):
                cons = (k + (3 - k) * f) / 3
                fin = cons * OUR_WXF_CELL[cell]
                row.append("k=%d %6.1f (%3.0f%%)" % (k, fin, 100 * np.mean(fin > cc)))
            print("     f=%.3f  %s" % (f, "  ".join(row)))
    print("  (pct = share of that cell's real post-shift fields the round would have placed in)")

    sweep = [0.101, 0.15, 0.198, 0.237, 0.30]
    res = run(fs, sweep)
    base = res[("A", 0.101)].mean()

    print("\n=== A: floor swept, band 12 kept, each cell at its own weighted x fidelity ===")
    print("  %-8s %-11s %-11s %-7s %-9s" % ("floor", "E[share]", "vs shipped", "ratio", "sign"))
    for f in sweep:
        v = res[("A", f)]
        d = v - res[("A", 0.101)]
        se = d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else 0
        print("  %-8.3f %-11.5f %-11s %-7.2f %d/%d%s"
              % (f, v.mean(), "%+.1f%%" % (100 * (v.mean() / base - 1)), v.mean() / base,
                 int((d > 0).sum()), len(d),
                 "  m/SE %.1f" % (abs(d.mean()) / se) if se else ""))

    print("\n=== B: same, but row composition also moves to the shared build's %.1f ===" % FIELD_WXF)
    print("  %-8s %-11s %-11s %-7s" % ("floor", "E[share]", "vs shipped", "ratio"))
    for f in sweep:
        v = res[("B", f)].mean()
        print("  %-8.3f %-11.5f %-11s %-7.2f"
              % (f, v, "%+.1f%%" % (100 * (v / base - 1)), v / base))

    c = res[("C", None)].mean()
    print("\n=== C: the shared build itself, all 11 hotkeys, no band ===")
    print("  E[share] %.5f  vs shipped %+.1f%%  (round final %.1f flat, cutoff min %.1f med %.1f)"
          % (c, 100 * (c / base - 1) if base else 0, FIELD_F * FIELD_WXF, cut.min(), np.median(cut)))

    print("\n=== per cell, arm A at the shipped floor vs 0.237 ===")
    for cell in sorted(OUR_WXF_CELL):
        idx = [i for i, (_t, cc, _d, _f) in enumerate(fs) if cc == cell]
        if not idx:
            continue
        a = res[("A", 0.101)][idx].mean()
        b = res[("A", 0.237)][idx].mean()
        print("  %-11s n%2d  %.5f -> %.5f  %s"
              % (cell, len(idx), a, b, "%+.1f%%" % (100 * (b / a - 1)) if a else "from zero"))


if __name__ == "__main__":
    main()
