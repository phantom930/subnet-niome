#!/usr/bin/env python3
"""fleet_price.py — correlated all-cut vs decorrelated all-HDR, priced at the FLEET level.

Per single hotkey, all-cut beats all-HDR on P(place) by 3.18x on HUDEP-2 (10.38% vs 3.26%,
measured in conj_test.py's k=0 arm against 28 real fields) and by 3.78x on K562 (cmp_k562.py).
But P(place) per hotkey is the wrong unit for a fleet, because the two constructions differ in
whether siblings can be decorrelated at all:

* all-cut is a deterministic function of the contract -- no window, no band -- so every hotkey
  running it submits IDENTICAL rows. Ten of them hold one score and take ranks r..r+9: they
  collect the TAIL of SCORE_DISTRIBUTION when they place, and all miss together when they don't.
* all-HDR carries a per-hotkey band window, so siblings hold independent 12-seed bands. One hit
  places ALONE, at the top of the curve, and the other nine sit at the floor.

`SCORE_DISTRIBUTION` is [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01], so rank 1
is worth 30x rank 10. That asymmetry is the whole question: does concentrating ten hotkeys on the
higher-P(place) construction beat spreading them across independent lottery tickets?

Method. For each current-regime 3-seed HUDEP-2 round, take the REAL field (every other miner's
final, deduped by hotkey, our own eleven excluded), then Monte-Carlo the round: draw 3 seeds
uniformly from 100-999, draw each all-HDR hotkey's 12-seed band inside a random 300-seed joined
window (which reproduces the measured union of ~104 over ten hotkeys), score every fleet member
from the MEASURED per-seed value ladder, insert them all into the field, and sum the curve shares
they capture. Reports E[curve share] per round for each fleet composition.

Per-seed values are measured, not modelled (CLAUDE.md "there is a third per-seed value"):
band 1.0000 exactly, all-cut cut-clean 0.2404, all-HDR cut-clean 0.162, dirty 0.101-0.106.
"""
import json
import os
import sys
import urllib.request
from collections import defaultdict

import numpy as np

DIST = np.array([0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01])
CELL = os.environ.get("FP_CELL", "HUDEP-2")
TRIALS = int(os.environ.get("FP_TRIALS", "4000"))
# The field lifted its own consistency floor on 2026-09-08 (see CLAUDE.md). Fields before
# that date price a field that no longer exists; "" keeps all current-regime rounds.
SINCE = os.environ.get("FP_SINCE", "")
RNG = int(os.environ.get("FP_RNG", "12345"))
SCORES = os.environ.get("FP_SCORES",
                        "/tmp/claude-0/-root-workspace-subnet-niome/"
                        "0f021db5-a874-4419-a235-ebcda9baabeb/scratchpad/s2.json")
OURS = os.environ.get("FP_OURS",
                      "/tmp/claude-0/-root-workspace-subnet-niome/"
                      "0f021db5-a874-4419-a235-ebcda9baabeb/scratchpad/real.json")

# --- measured construction parameters -------------------------------------------------
# `weighted x fidelity` is PER CELL TYPE. The first version of this script used single constants
# (AC 329.7*0.8686 = 286.4, AH 339.9*0.8931 = 303.6) taken from research builds on HUDEP-2, and
# applied them to K562 as well -- where the live fleet measures 220. That is a 38% inflation of
# BOTH arms on K562, and it matters because all-cut's whole case is that its flat score sometimes
# clears the cutoff: inflate it and it clears fields it cannot actually reach.
#
# These are live medians from the scored feed, identified by ss58 address (uids are recycled).
# all-HDR: every band hotkey, post-2026-09-08. all-cut: h0's own rounds while ALL_CUT_HOTKEYS
# still named it, isolated by the cons signature (all-cut 0.16-0.21 against all-HDR's 0.09-0.10)
# -- n is only 1-3 per cell, so the ratio matters more than the level:
#
#   cell         all-cut   all-HDR   ratio      script previously assumed 0.943
#   HUDEP-2        291.7     298.3   0.978
#   CD34+_HSPC     249.1     233.3   1.068
#   K562           219.4     220.1   0.997      mean 1.014
#
# So all-cut does NOT give up weighted x fidelity, matching CLAUDE.md's "tracks all-cut within 3%".
# AC_WXF is therefore set to the cell's own measured all-cut value where it exists.
AH_WXF_CELL = {"HUDEP-2": 298.3, "CD34+_HSPC": 233.3, "K562": 220.1, "HEK293": 231.3}
AC_WXF_CELL = {"HUDEP-2": 291.7, "CD34+_HSPC": 249.1, "K562": 219.4, "HEK293": 231.3}

# all-cut's clean fraction. 569/900 is the build-time measurement; h0's six live all-cut rounds
# averaged round-cons 0.1774 against the 0.1911 that 569/900 predicts, implying 477/900. Both are
# run -- AC_CLEAN is the shipped assumption, FP_AC_CLEAN=477 the live-implied one.
AC_CLEAN = int(os.environ.get("FP_AC_CLEAN", "569"))
AC_VCLEAN, AC_VDIRTY = 0.2404, 0.1064
AH_BAND, AH_CLEAN, AH_VCLEAN, AH_VDIRTY = 12, 17, 0.162, 0.101
AC_WXF = AC_WXF_CELL[CELL]
AH_WXF = AH_WXF_CELL[CELL]
JOINED = 300          # the joined band space the plan confines all-HDR hotkeys to

# Band size and band SPACE are the two quantities the window layout sets, so both are overridable:
# joined300.py measures the band per cell type at each window width, and the layout question --
# confine every hotkey to a joined 300, or spread them over the whole 900 -- is exactly
# `FP_JOINED=300` against `FP_JOINED=900` at that cell's measured band. `FP_AH_WXF` carries the
# matching `weighted x fidelity`, since a wider window changes it too (measured, not assumed).
AH_BAND = int(os.environ.get("FP_BAND", AH_BAND))
JOINED = int(os.environ.get("FP_JOINED", JOINED))
AH_WXF = float(os.environ.get("FP_AH_WXF", AH_WXF))
SEED_LO, SEED_HI = 100, 999
NSEED = SEED_HI - SEED_LO + 1


def fields():
    sc = json.load(open(SCORES))
    sc = sc if isinstance(sc, list) else sc.get("data") or sc.get("items")
    tk = json.load(urllib.request.urlopen("https://niome-api.genomes.io/api/v3/tasks?limit=500"))
    tk = tk if isinstance(tk, list) else tk.get("data") or tk.get("items")
    meta = {}
    for t in tk:
        tid = t.get("task_id") or t.get("id")
        c = (t.get("content") or {}).get("contract", {})
        s = str(c.get("seed", "") or "")
        meta[tid] = (c.get("cell_type") or t.get("cell_type"),
                     [int(x) for x in s.split(",") if x.strip().isdigit()],
                     (t.get("created_at") or "")[:10])
    ours = set(json.load(open(OURS)).values())
    best = defaultdict(dict)
    for r in sc:
        cell, seeds, created = meta.get(r["task_id"], (None, [], ""))
        if cell != CELL or len(seeds) != 3 or created < SINCE:
            continue
        hk = r["miner_hotkey"]
        if hk in ours:
            continue
        if hk not in best[r["task_id"]] or r["final_score"] > best[r["task_id"]][hk]["final_score"]:
            best[r["task_id"]][hk] = r
    out = []
    for tid, v in best.items():
        f = sorted((x["final_score"] for x in v.values()), reverse=True)
        if len(f) >= 10 and f[0] > 0:
            out.append((tid, np.asarray(f)))
    return out


def shares(field, mine):
    """Curve share our fleet captures. Ties resolve to consecutive ranks, as the validator's sort does.

    Siblings holding the SAME score take consecutive ranks r..r+n, not the same rank n+1 times.
    Counting only strictly-greater siblings gave every tied row the best of the tied ranks, which
    paid four tied rows 1.20 of a curve that sums to 1.00. The overcount scales with how
    CORRELATED a composition is -- identical all-cut submissions, or band hotkeys sharing a seed
    because their windows overlap -- which is exactly the quantity a layout or composition
    comparison is measuring, so it is not a wash between arms.
    """
    total = 0.0
    for i, s in enumerate(sorted(mine, reverse=True)):
        # competitors strictly above, plus every sibling already placed above this one
        rank = int(np.sum(field > s)) + 1 + i
        if rank <= 10:
            total += DIST[rank - 1]
    return total


def simulate(field, comp, rng, trials):
    """comp: (n_allcut, n_allhdr). Returns mean curve share over `trials` seed draws."""
    n_ac, n_ah = comp
    tot = 0.0
    ac_clean_p = AC_CLEAN / NSEED
    for _ in range(trials):
        seeds = rng.integers(SEED_LO, SEED_HI + 1, 3)
        mine = []
        if n_ac:
            # identical rows -> one score, repeated. Which seeds are cut-clean is a property of
            # the build, so draw membership per seed at the measured clean fraction.
            ncl = int(rng.random(3).__lt__(ac_clean_p).sum())
            cons = (ncl * AC_VCLEAN + (3 - ncl) * AC_VDIRTY) / 3
            mine += [cons * AC_WXF] * n_ac
        if n_ah:
            lo = rng.integers(SEED_LO, SEED_HI - JOINED + 2)
            joined = np.arange(lo, lo + JOINED)
            for _h in range(n_ah):
                band = set(rng.choice(joined, AH_BAND, replace=False).tolist())
                k = sum(1 for s in seeds if int(s) in band)
                rest = 3 - k
                ncl = int(rng.random(rest).__lt__(AH_CLEAN / NSEED).sum()) if rest else 0
                cons = (k * 1.0 + ncl * AH_VCLEAN + (rest - ncl) * AH_VDIRTY) / 3
                mine.append(cons * AH_WXF)
        tot += shares(field, mine)
    return tot / trials


def main():
    fs = fields()
    print("%s: %d current-regime 3-seed rounds, real fields (our 11 hotkeys excluded)"
          % (CELL, len(fs)))
    print("all-cut  : clean %d/900, clean-seed %.4f, dirty %.4f, weighted x fid %.1f"
          % (AC_CLEAN, AC_VCLEAN, AC_VDIRTY, AC_WXF))
    print("all-HDR  : band %d, clean %d/900, clean-seed %.3f, dirty %.3f, weighted x fid %.1f"
          % (AH_BAND, AH_CLEAN, AH_VCLEAN, AH_VDIRTY, AH_WXF))
    print("%d trials per round per composition\n" % TRIALS)
    comps = [(0, 11), (1, 10), (2, 9), (4, 7), (6, 5), (8, 3), (11, 0)]
    print("%-14s %-12s %-12s %-10s" % ("all-cut/all-HDR", "E[share]", "vs current", "ratio"))
    res = {}
    for comp in comps:
        rng = np.random.default_rng(RNG)
        vals = [simulate(f, comp, rng, TRIALS) for _tid, f in fs]
        res[comp] = float(np.mean(vals))
    cur = res[(1, 10)]
    for comp in comps:
        v = res[comp]
        print("%-14s %-12.5f %-12s %-10.2f"
              % ("%d / %d" % comp, v, "%+.1f%%" % (100 * (v / cur - 1)) if cur else "-",
                 v / cur if cur else 0))
    print()
    print("current fleet is 1 all-cut (h0) + 10 all-HDR = the '1 / 10' row")


if __name__ == "__main__":
    main()
