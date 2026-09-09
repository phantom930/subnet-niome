"""mut_tiebreak.py — spend FastGreedy's tie slack on mutation_weight, not GC.

`total_weighted_score = sum(structural * mutation_weight)`, and our structural is already 0.850 of
a possible 1.0 with `offtarget_factor` at a perfect 1.0. The whole remaining gap to the top miner
is the weight term: measured on HEK293 a8b9f1bb our mean `mutation_weight` is **0.879** against a
reachable ~1.02, and closing it gives 250 * 0.8497 * 1.016 = 215.9 against the rank-1 miner's
216.43.

The leak this targets: the Cas12a group came out **63 light / 17 heavy**, because -- as
FastGreedy's own docstring says -- "the min-union objective is blind to mutation_weight". Ties are
broken by first index or at random, so a light-mutation guide wins as often as a heavy one.

**This is NOT the `light_cell_rows` route**, which all_cut.py records as measured dead: capping the
light mutation drives stage 5's mutation-coverage entropy to its floor and fidelity falls faster
than weighted rises (20 builds, no rank moved). The difference is scale -- shifting only the group's
80 rows leaves ~76 light rows across the Cas9 half, so the entropy term should barely move, where
`light_cell_rows` cuts light to 12. That is the hypothesis under test, and fidelity is the number
that falsifies it.

Ties have equal MARGINAL coverage but cover different seeds, so the band can still move either way
and is reported alongside weighted, with frequency x value as the decision variable.

    python mut_tiebreak.py [task_id] [lo-hi]
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


def weight_pen(recs, contract):
    """Lower is better, so negate the weight: the greedy then prefers the heavier mutation."""
    w = contract.get("mutation_weights", {})
    return np.array([-float(w.get(r["mutation"], 1.0)) for r in recs])


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
    mw = contract.get("mutation_weights", {})
    heavy = max(mw, key=lambda m: mw.get(m, 1.0)) if mw else None
    n_heavy = sum(1 for r in rows if r["mutation"] == heavy)
    grp_heavy = sum(1 for r in group if r["mutation"] == heavy)
    mean_w = st.mean(mw.get(r["mutation"], 1.0) for r in rows)
    return {"arm": label, "band": band, "union": len(bad), "cas9": len(cas9),
            "heavy_rows": n_heavy, "group_heavy": grp_heavy, "mean_weight": mean_w,
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
    wpen = weight_pen(recs, contract)
    for label, p, rs in (("baseline (tie = first index)", None, RESTARTS),
                         ("heavy-mutation tie-break", wpen, RESTARTS)):
        idx = greedy(recs, cfg, caps, pen=p, rng=True, restarts=rs)
        got = finish([recs[i] for i in idx], recs, contract, reference, cell_types,
                     ctx, sites, cfg, n_rows, label)
        MT.free_gpu_memory()
        out.append(got)
        if "reason" in got:
            print(f"  {label:<30} declined: {got['reason']}")
            continue
        print(f"  {label:<28} band {got['band']:>3}  group heavy {got['group_heavy']:>3}/80  "
              f"rows heavy {got['heavy_rows']:>3}/250  mean_w {got['mean_weight']:.4f}  "
              f"weighted {got['weighted']:>7.2f}  fid {got['fidelity']:.4f}  "
              f"k=1 {got['k1_final']:>6.2f}  cov {got['fleet_cov']:.1%}")
    good = [r for r in out if "band" in r and "reason" not in r]
    if len(good) == 2:
        b, g = good
        print(f"\n  delta: band {g['band']-b['band']:+d}   group heavy "
              f"{g['group_heavy']-b['group_heavy']:+d}   mean_w "
              f"{g['mean_weight']-b['mean_weight']:+.4f}   weighted "
              f"{(g['weighted']/b['weighted']-1)*100:+.1f}%   "
              f"fidelity {(g['fidelity']/b['fidelity']-1)*100:+.1f}%   "
              f"k=1 final {(g['k1_final']/b['k1_final']-1)*100:+.1f}%")
        ev_b = b["fleet_cov"] * b["k1_final"]
        ev_g = g["fleet_cov"] * g["k1_final"]
        print(f"  frequency x value: {ev_b:.2f} -> {ev_g:.2f}  ({(ev_g/ev_b-1)*100:+.1f}%)")
    json.dump({"task": TASK, "cell": cell, "window": [lo, hi], "arms": out},
              open(os.getenv("MTB_OUT", "mut_tiebreak.json"), "w"), indent=1)
    print(f"\nwrote {os.getenv('MTB_OUT', 'mut_tiebreak.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
