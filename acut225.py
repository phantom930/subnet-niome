#!/usr/bin/env python3
"""acut225.py — all-cut on four width-225 windows tiling 900, for HEK293 only.

Narrowing an all-cut window normally loses: it raises the clean FRACTION and shrinks the clean SET,
and the set is what pays. That is why width 225/300 lost on K562 (whole-window clean 549 of 900 =
61%, against 0.25 x 0.93 = 23% for a quarter window).

HEK293 inverts it, because its whole-window clean set is nearly empty. Accessibility 0.35 gives
cut_p ~0.875, so cutting on all 900 seeds holds only **137** of them (15.2%, measured live on
634a512c). A quarter window needs the guide to cut on 225 seeds instead of 900, which is a far
weaker requirement, so:

    whole window   P(seed clean) = 137/900                = 0.152
    width 225      P(seed clean) = 0.25 x clean_fraction  = 0.222 at 89% clean

With four hotkeys tiling 100-324 / 325-549 / 550-774 / 775-999, every seed is in exactly one
hotkey's window. Two of a round's three seeds share a window on 1 - (4*3*2)/4**3 = 62.5% of rounds,
and that hotkey then averages two clean seeds against one floor seed — which is the rank 9-10 path
this is for, against HEK293's median rank-10 cutoff of 64.5.

Both halves are measured: the per-window construction (clean set, and consistency in-window-clean /
in-window-dirty / off-window), then the fleet priced on every real HEK293 field.

    python acut225.py [task_id]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else "634a512c-5760-402f-8b57-70d8e1dcc314"
os.environ["NIOME_INSTANCE"] = "acut225"

import dataclasses as _dc  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import random  # noqa: E402
import statistics as st  # noqa: E402
from collections import Counter, defaultdict  # noqa: E402
from itertools import product  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import fastgreedy as FG  # noqa: E402
from niome_subnet.genomics import mt19937 as MT  # noqa: E402
from sd_task import score, task_content  # noqa: E402

WINDOWS = [(100, 324), (325, 549), (550, 774), (775, 999)]
N_SAMPLE = int(os.getenv("A225_SAMPLE", "10"))
OUT = os.getenv("A225_OUT", "acut225.json")
DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]


def measure(lo, hi, task, contract, reference, cell_types, rng):
    """Build all-cut over [lo, hi] and measure consistency in each regime a seed can land in."""
    cell = contract.get("cell_type")
    base = AC.config_for(cell) or AC.AllCutConfig()
    span = hi - lo + 1
    cfg = _dc.replace(base, start_seed=lo, end_seed=hi,
                      cas12a_max_fail=max(1, round(base.cas12a_max_fail * span / 900)))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments

    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        if not bank:
            return None
        AC.save_bank(path, bank)
    records = AC.load_bank(path)

    sel = FG.FastGreedy(records, window_lo=lo, window_hi=hi)
    idx, _u = sel.best(cfg.group_size, restarts=12)
    group = [records[i] for i in idx]
    bad = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = sorted(set(range(lo, hi + 1)) - bad)
    cas9 = AC.scan_cas9(np.array(clean, dtype=np.int64), contract, cell_types, ctx, sites, cfg,
                        n_rows - cfg.group_size)
    if len(cas9) < n_rows - cfg.group_size:
        return {"window": [lo, hi], "reason": f"cas9 pool {len(cas9)}"}
    rows = AC.assemble(group, cas9, contract, ctx, cfg, n_rows)
    cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
    off = [s for s in range(100, 1000) if not (lo <= s <= hi)]

    def sample(pool):
        if not pool:
            return None
        return st.mean(score(rows, contract, reference, cell_types, seed=s)["consistency"]
                       for s in rng.sample(pool, min(N_SAMPLE, len(pool))))

    base_s = score(rows, contract, reference, cell_types, seed=clean[0])
    got = {"window": [lo, hi], "mf": cfg.cas12a_max_fail, "bank": len(records),
           "group": cfg.group_size, "union": len(bad), "clean": len(clean),
           "clean_fraction": len(clean) / span, "cas9_pool": len(cas9), "cells": cells,
           "weighted": base_s["weighted"], "fidelity": base_s["fidelity"],
           "cons_clean": sample(clean), "cons_dirty": sample(sorted(bad)),
           "cons_off": sample(off), "clean_seeds": clean}
    MT.free_gpu_memory()
    return got


def price(per_window, contract, cell_types):
    """Expected curve share of the four-hotkey fleet, on every real HEK293 three-seed field."""
    from sd_task import OURS
    rows = json.load(open("sd_task_scores.json"))
    rows = rows if isinstance(rows, list) else (rows.get("data") or rows.get("items") or [])
    listing = json.load(open("sd_task_listing.json"))
    items = listing if isinstance(listing, list) else (listing.get("items") or [])
    fields, seeds_of = [], []
    for t in items:
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != contract.get("cell_type"):
            continue
        sd = [int(x) for x in str(c.get("seed", "")).split(",") if x.strip().isdigit()]
        if len(sd) != 3:
            continue
        best = {}
        for r in rows:
            if r.get("task_id") != t["id"] or r.get("miner_hotkey") in OURS:
                continue
            best[r["miner_hotkey"]] = max(float(r.get("final_score") or 0),
                                          best.get(r["miner_hotkey"], 0.0))
        if len(best) >= 50:
            fields.append(sorted(best.values(), reverse=True))
            seeds_of.append(sd)

    def share(sc, field):
        rank = sum(1 for f in field if f > sc) + 1
        return DIST[rank - 1] if rank <= len(DIST) else 0.0

    tot_fleet = tot_best = 0.0
    placed = 0
    for field, sd in zip(fields, seeds_of):
        finals = []
        for w in per_window:
            lo, hi = w["window"]
            cleans = set(w["clean_seeds"])
            cons = st.mean([w["cons_clean"] if s in cleans
                            else (w["cons_dirty"] if lo <= s <= hi else w["cons_off"])
                            for s in sd])
            finals.append(w["weighted"] * cons * w["fidelity"])
        # Four hotkeys compete in the same field, so they take consecutive ranks.
        ranks = sorted(finals, reverse=True)
        got = 0.0
        for i, sc in enumerate(ranks):
            rank = sum(1 for f in field if f > sc) + 1 + i
            got += DIST[rank - 1] if rank <= len(DIST) else 0.0
        tot_fleet += got
        tot_best += share(max(finals), field)
        placed += 1 if got > 0 else 0
    n = max(1, len(fields))
    return {"rounds": len(fields), "e_share_fleet": tot_fleet / n,
            "e_share_best_hotkey": tot_best / n, "p_place": placed / n}


def main():
    task, contract, reference = task_content(TASK)
    cell = contract.get("cell_type")
    if cell != "HEK293":
        print(f"note: task is {cell}, not HEK293 — the inversion this tests is HEK293-specific")
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    rng = random.Random(23)
    print(f"task {TASK[:8]}  {cell}  four windows tiling 100-999\n")
    print(f"  {'window':<10} {'mf':>3} {'union':>6} {'clean':>6} {'frac':>6} {'cas9':>6} "
          f"{'cells':>5} {'weighted':>9} {'fid':>7} {'clean cons':>10} {'dirty':>7} {'off':>7}")
    per = []
    for lo, hi in WINDOWS:
        got = measure(lo, hi, task, contract, reference, cell_types, rng)
        if not got or "clean" not in got:
            print(f"  {f'{lo}-{hi}':<10}  declined: {(got or {}).get('reason', 'bank empty')}")
            continue
        per.append(got)
        print(f"  {f'{lo}-{hi}':<10} {got['mf']:>3} {got['union']:>6} {got['clean']:>6} "
              f"{got['clean_fraction']:>5.0%} {got['cas9_pool']:>6} {got['cells']:>5} "
              f"{got['weighted']:>9.1f} {got['fidelity']:>7.4f} {got['cons_clean']:>10.4f} "
              f"{got['cons_dirty']:>7.4f} {got['cons_off']:>7.4f}")
    if not per:
        raise SystemExit("no window produced a submission")

    p = price(per, contract, cell_types)
    print(f"\npriced on {p['rounds']} real {cell} three-seed fields:")
    print(f"  fleet of {len(per)} width-225 all-cut hotkeys : E[share] {p['e_share_fleet']:.4f}"
          f"   places on {p['p_place']:.0%} of rounds")
    print(f"  best single hotkey of the four             : E[share] "
          f"{p['e_share_best_hotkey']:.4f}")
    print(f"\n  references on HEK293 (price_cell.py, same fields):")
    print(f"    all-HDR, per hotkey        E[share] 0.0059   (x4 hotkeys ~ 0.0236)")
    print(f"    all-cut whole window       E[share] 0.0027")
    json.dump({"task": TASK, "cell": cell, "windows": per, "pricing": p},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
