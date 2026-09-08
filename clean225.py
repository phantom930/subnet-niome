#!/usr/bin/env python3
"""clean225.py — Cas12a clean-seed counts for all-cut on a 225-wide seed window.

all-cut today min-unions its Cas12a group over the whole 100-999 window and reaches 557/900 clean
(61.9%) on b9051bc7. Narrowing the window to 225 confines the group's failures to a quarter of the
seed space, so the clean *fraction inside the window* should rise sharply — the question is how far,
and at what screen.

``cas12a_max_fail`` is a screen over the window, so it must scale with the span or the bank is
over-constrained (100 of 900 is 11%; held fixed at 225 it would demand 125 of 225 = 56% clean per
guide). Proportional is 25. Two looser values are measured on one window to check the sensitivity,
because proportional scaling is an assumption, not a result.

Four windows tile 100-999 exactly at width 225.
"""
import os
os.environ["NIOME_INSTANCE"] = "c225"

import dataclasses, json, os.path, sys, time
sys.argv = ["x"]; sys.path.insert(0, "/root/workspace/subnet-niome")
import logging; logging.basicConfig(level=logging.ERROR)

import numpy as np
import genExp as G
from niome_subnet.genomics import all_cut as AC
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank

TASK = os.getenv("C225_TASK", "task-b9051bc7.json")
WIDTH = int(os.getenv("C225_W", "225"))
GROUP = int(os.getenv("C225_GROUP", "42"))
RULE = os.getenv("C225_RULE", "cut")
OUT = os.getenv("C225_OUT", "clean225.json")


def clean_for(contract, reference, cell_types, ctx, sites, base, window, mf):
    """Bank + min-union for one (window, max_fail); returns the clean set inside the window."""
    cfg = dataclasses.replace(base, start_seed=window[0], end_seed=window[1],
                              cas12a_max_fail=mf, group_size=GROUP, rule=RULE)
    path = os.path.join("data/all_cut", f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    t0 = time.monotonic()
    cached = os.path.exists(path)
    if not cached:
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return {"error": "bank scan produced nothing", "build_s": round(time.monotonic()-t0, 1)}
        save_bank(path, bank)
    records = load_bank(path)
    if len(records) < GROUP:
        return {"error": f"bank {len(records)} < group {GROUP}",
                "bank": len(records), "build_s": round(time.monotonic()-t0, 1)}
    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed)
    index, union = sel.best(GROUP, restarts=12)
    bad = set()
    for i in index:
        bad.update(int(x) for x in records[i]["fails"])
    span = window[1] - window[0] + 1
    clean = sorted(set(range(window[0], window[1] + 1)) - bad)
    nf = [len(records[i]["fails"]) for i in index]
    return {"bank": len(records), "cached": cached, "union": len(bad), "clean": len(clean),
            "span": span, "clean_frac": len(clean) / span,
            "group_fails_min": int(min(nf)), "group_fails_max": int(max(nf)),
            "build_s": round(time.monotonic() - t0, 1)}


def main():
    task = json.load(open(TASK))
    contract, reference = task["content"]["contract"], task["content"]["hbb_reference"]
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types(); G.load_sequence()
    base = AC.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    only = os.getenv("C225_WINDOWS")
    windows = ([tuple(int(x) for x in w.split("-")) for w in only.split(",")] if only else
               [(100 + i * WIDTH, 100 + (i + 1) * WIDTH - 1) for i in range(900 // WIDTH)])
    mf_prop = int(os.getenv("C225_MF", str(max(1, round(base.cas12a_max_fail * WIDTH / 900)))))
    print(f"{task['id'][:8]}  {cell}  group {GROUP}  width {WIDTH}  rule {RULE}")
    print(f"baseline (whole 100-999, mf {base.cas12a_max_fail}): clean 557/900 = 61.9%")
    print(f"proportional screen at width {WIDTH}: mf {mf_prop}\n")
    print(f"{'window':>10} {'mf':>4} {'bank':>7} {'union':>6} {'clean':>6} {'/span':>7} "
          f"{'fails':>9} {'build':>8}")
    res = []
    for w in windows:
        r = clean_for(contract, reference, cell_types, ctx, sites, base, w, mf_prop)
        r.update(window=list(w), max_fail=mf_prop)
        res.append(r); json.dump(res, open(OUT, "w"), indent=1)
        if "error" in r:
            print(f"{w[0]:>4}-{w[1]:<5} {mf_prop:>4}  {r['error']}  {r['build_s']}s", flush=True)
        else:
            print(f"{w[0]:>4}-{w[1]:<5} {mf_prop:>4} {r['bank']:>7} {r['union']:>6} "
                  f"{r['clean']:>6} {r['clean_frac']:>6.1%} "
                  f"{r['group_fails_min']:>4}-{r['group_fails_max']:<4} "
                  f"{r['build_s']:>7.1f}s", flush=True)
    # sensitivity to the screen on one window — proportional scaling is an assumption
    print()
    for mf in (int(os.getenv("C225_MF2", "50")), base.cas12a_max_fail):
        if mf == mf_prop:
            continue
        r = clean_for(contract, reference, cell_types, ctx, sites, base, windows[0], mf)
        r.update(window=list(windows[0]), max_fail=mf)
        res.append(r); json.dump(res, open(OUT, "w"), indent=1)
        if "error" in r:
            print(f"{windows[0][0]:>4}-{windows[0][1]:<5} {mf:>4}  {r['error']}", flush=True)
        else:
            print(f"{windows[0][0]:>4}-{windows[0][1]:<5} {mf:>4} {r['bank']:>7} {r['union']:>6} "
                  f"{r['clean']:>6} {r['clean_frac']:>6.1%} "
                  f"{r['group_fails_min']:>4}-{r['group_fails_max']:<4} "
                  f"{r['build_s']:>7.1f}s", flush=True)
    ok = [r for r in res if "clean" in r and r["max_fail"] == mf_prop]
    if ok:
        tot = sum(r["clean"] for r in ok)
        print(f"\nfour windows at mf {mf_prop}: clean {[r['clean'] for r in ok]}  "
              f"total {tot}/900 = {tot/900:.1%}   (whole-window baseline 557/900 = 61.9%)")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
