#!/usr/bin/env python3
"""gc_tiebreak.py — spend FastGreedy's tie slack on GC, and measure what it costs the band.

HEK293's all-HDR bank is GC-skewed (mean 0.675) because its accessibility of 0.35 leaves `energy`
below the clamp, so P(HDR) keeps rising with GC and the min-union selects the high-GC tail. That is
why HEK293's gc_score is 0.840 against the erythroid types' 0.93 and its weighted is ~166 against
~250. Tightening the GC BOUND was measured and is a net loss (gcdist.py): it shrinks the candidate
pool before the min-union runs, so the band collapses 8 -> 5.

This spends a different resource. The greedy's `argmin` is frequently tied -- measured 71 of 80
picks, median 12 candidates, max 577 -- and ties are currently broken by first-index or uniformly
at random. Preferring the tied candidate closest to GC 0.50 costs nothing in MARGINAL coverage.

**It is not free, though, and that is the point of measuring.** Tied candidates cover DIFFERENT
seeds, so the greedy's path diverges after the first differing pick and the final union can move
either way. Band and weighted are therefore reported together, with the k=1 round score (the value
that actually places) as the decision variable.

    python gc_tiebreak.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "gctb")

import dataclasses as _dc     # noqa: E402
import json                   # noqa: E402
import logging                # noqa: E402
import statistics as st       # noqa: E402
from collections import Counter  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np            # noqa: E402

import genExp as G            # noqa: E402
from niome_subnet.genomics import all_cut as AC      # noqa: E402
from niome_subnet.genomics import all_hdr as AH      # noqa: E402
from niome_subnet.genomics import fastgreedy as FG   # noqa: E402
from niome_subnet.genomics import mt19937 as MT      # noqa: E402
from niome_subnet.utils import settings              # noqa: E402
from sd_task import score, task_content              # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "a8b9f1bb-f25d-42af-92a1-e0aff16bee3a"
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "100-399"
FLOOR = 0.10
RESTARTS = int(os.getenv("GCTB_RESTARTS", "12"))


def gc_of(rec):
    g = rec["guide"]
    return (g.count("G") + g.count("C")) / len(g)


def greedy(recs, cfg, caps, pen=None, rng=None, restarts=1):
    """Min-union greedy, optionally preferring low ``pen`` among argmin ties.

    Mirrors FastGreedy.build's loop (floors, caps, argmin) so the only difference between arms is
    the tie-break. With ``rng`` set, ties are narrowed to the best-``pen`` quartile and then sampled,
    which keeps restarts exploring instead of collapsing to one deterministic path.
    """
    sel = FG.FastGreedy(recs, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=caps)
    fails = [np.asarray(f) for f in sel.fail_lists]
    BIG = 10 ** 9
    n = len(recs)
    hi = max(int(f.max()) for f in fails if f.size) + 1
    best_idx, best_u = None, None
    for r in range(max(1, restarts)):
        rg = np.random.default_rng(r) if (rng is not None and r > 0) else None
        uncovered = np.ones(max(cfg.end_seed + 1, hi), dtype=bool)
        uncovered[:cfg.start_seed] = False
        taken = np.zeros(n, dtype=bool)
        chosen = []
        for _ in range(cfg.group_size):
            cost = np.array([BIG if taken[i] else int(uncovered[fails[i]].sum())
                             for i in range(n)])
            m = int(cost.min())
            if m >= BIG:
                break
            ties = np.flatnonzero(cost == m)
            if pen is not None and len(ties) > 1:
                tp = pen[ties]
                if rg is not None:
                    keep = ties[tp <= np.quantile(tp, 0.25)]
                    pick = int(rg.choice(keep)) if len(keep) else int(ties[0])
                else:
                    pick = int(ties[int(np.argmin(tp))])
            elif rg is not None and len(ties) > 1:
                pick = int(rg.choice(ties))
            else:
                pick = int(ties[0])
            chosen.append(pick)
            taken[pick] = True
            uncovered[fails[pick]] = False
        u = len({int(s) for i in chosen for s in fails[i]})
        if best_u is None or u < best_u:
            best_idx, best_u = chosen, u
    return best_idx


def finish(group, recs, contract, reference, cell_types, ctx, sites, cfg, n_rows, label):
    """scan_cas9 -> assemble -> score, exactly as all_hdr.build_submission does."""
    bad = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad), dtype=np.int64)
    if clean.size == 0:
        return {"arm": label, "reason": "union covers the window"}
    cas9 = AH.scan_cas9(clean, contract, cell_types, ctx, sites, cfg, n_rows - cfg.group_size)
    if len(cas9) < n_rows - cfg.group_size:
        return {"arm": label, "reason": f"cas9 pool {len(cas9)}", "band": int(clean.size)}
    rows = AH.assemble(group, cas9, contract, ctx, cfg, n_rows)
    s = score(rows, contract, reference, cell_types, seed=int(clean[0]))
    detail = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    feats = [d["features"] for d in detail]
    by_cas = {}
    for d, f in zip(detail, feats):
        by_cas.setdefault(d["experiment"]["cas_system"], []).append(f["gc_score"])
    band = int(clean.size)
    k1 = s["weighted"] * ((1 + 2 * FLOOR) / 3) * s["fidelity"]
    return {"arm": label, "band": band, "union": len(bad), "cas9": len(cas9),
            "weighted": s["weighted"], "fidelity": s["fidelity"],
            "gc_score": st.mean(f["gc_score"] for f in feats),
            "dist_score": st.mean(f["dist_score"] for f in feats),
            "gc_mean": st.mean(f["gc"] for f in feats),
            "gc_score_cas12a": st.mean(by_cas.get("Cas12a", [0])),
            "gc_score_cas9": st.mean(by_cas.get("Cas9", [0])),
            "cells": len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
            "k1_final": k1, "fleet_cov": 1 - (1 - 7 * band / 900.0) ** 3}


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = (int(x) for x in WINDOW.split("-"))
    span = hi - lo + 1
    base = AH.config_for(cell)
    cfg = _dc.replace(base, hdr_range=(lo, hi),
                      main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span),
                      variants=min(base.variants, AH.WIDE_WINDOW_VARIANTS))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        if not bank:
            raise SystemExit("bank empty")
        AC.save_bank(path, bank)
    recs = AC.load_bank(path)
    caps = AH._group_caps(contract, ctx, cfg)
    pen = np.array([abs(gc_of(r) - 0.5) for r in recs])
    print(f"task {TASK[:8]}  {cell}  window {lo}-{hi}  mf {cfg.main_max_fail}  "
          f"group {cfg.group_size}  bank {len(recs)} (GC mean {np.mean([gc_of(r) for r in recs]):.3f})\n")

    out = []
    for label, p, rs in (("baseline  (tie = first index)", None, RESTARTS),
                         ("GC tie-break (nearest 0.50)", pen, RESTARTS)):
        idx = greedy(recs, cfg, caps, pen=p, rng=True, restarts=rs)
        got = finish([recs[i] for i in idx], recs, contract, reference, cell_types,
                     ctx, sites, cfg, n_rows, label)
        MT.free_gpu_memory()
        out.append(got)
        if "reason" in got:
            print(f"  {label:<30} declined: {got['reason']}")
            continue
        print(f"  {label:<30} band {got['band']:>3}  weighted {got['weighted']:>7.2f}  "
              f"gc_s {got['gc_score']:.3f} (12a {got['gc_score_cas12a']:.3f} / "
              f"9 {got['gc_score_cas9']:.3f})  fid {got['fidelity']:.4f}  "
              f"k=1 {got['k1_final']:>6.2f}  cov {got['fleet_cov']:.1%}")
    good = [r for r in out if "band" in r and "reason" not in r]
    if len(good) == 2:
        b, g = good
        print(f"\n  delta: band {g['band']-b['band']:+d}   weighted "
              f"{(g['weighted']/b['weighted']-1)*100:+.1f}%   "
              f"gc_score {(g['gc_score']/b['gc_score']-1)*100:+.1f}%   "
              f"k=1 final {(g['k1_final']/b['k1_final']-1)*100:+.1f}%   "
              f"fleet cov {g['fleet_cov']-b['fleet_cov']:+.1%}")
        ev_b = b["fleet_cov"] * b["k1_final"]
        ev_g = g["fleet_cov"] * g["k1_final"]
        print(f"  frequency x value: {ev_b:.2f} -> {ev_g:.2f}  ({(ev_g/ev_b-1)*100:+.1f}%)")
    json.dump({"task": TASK, "cell": cell, "window": [lo, hi], "arms": out},
              open(os.getenv("GCTB_OUT", "gc_tiebreak.json"), "w"), indent=1)
    print(f"\nwrote {os.getenv('GCTB_OUT', 'gc_tiebreak.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
