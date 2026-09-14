#!/usr/bin/env python3
"""band_dupe.py — how much do the eleven hotkeys' clean bands DUPLICATE each other?

The fleet's whole payout case is that eleven bands cover eleven different parts of the seed space:
a coldkey spikes on `1 - (1 - |union|/900)**3`, so what matters is the UNION of the bands, not the
sum. Duplication is the gap between the two, and it has two possible sources:

  * **pigeonhole** — eleven bands of ~13 is 143 seeds, which cannot sit distinctly in a 300-seed
    joined space however independent they are. This is arithmetic and unavoidable at that space.
  * **correlation** — sibling windows overlap heavily (stride 30, width 75-150), so their Cas12a
    candidate pools overlap and `FastGreedy` is deterministic given its bank. If the greedy picks
    the same guides it picks the same band, and eleven hotkeys are worth barely more than one.

Only the second is a defect, and separating them needs the independence model as a control: place
eleven bands of the measured sizes uniformly at random in the same space and compare the unions.
Above the model = correlated. At the model = pigeonhole only, and the fix is a wider space rather
than a different layout.

**What counts as the band.** The clean set of the Cas12a min-union group over the hotkey's window --
`scan_cas9` fills rows onto that band and cannot alter it (band_rebuild.py), so the Cas9 scan and
`assemble` are skipped. That is also what `meta["clean"]` reports and what every band figure in
CLAUDE.md means.

**One bank per (contract, window), deliberately.** The 3x faster shortcut -- one 900-seed bank per
contract, per-window fails derived by filtering -- is measured WRONG in band_rebuild.py:
`main_max_fail` scaled by rate is far stricter over a wide window, the bank collapses ~300k -> ~1k
and every band comes out 4 instead of 12-13.

    BD_LAYOUT=plan|fallback  BD_N=3  BD_CELLS=K562,...  python -u band_dupe.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "banddupe")

import dataclasses as _dc     # noqa: E402
import itertools              # noqa: E402
import json                   # noqa: E402
import logging                # noqa: E402
import random                 # noqa: E402
import statistics as st       # noqa: E402
import time                   # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np            # noqa: E402

import genExp as G            # noqa: E402
import joined_window as JW    # noqa: E402
from niome_subnet.genomics import all_hdr as AH      # noqa: E402
from niome_subnet.genomics import fastgreedy as FG   # noqa: E402
from niome_subnet.genomics import mt19937 as MT      # noqa: E402
from sd_task import task_content                     # noqa: E402

LAYOUT = os.getenv("BD_LAYOUT", "plan")
N_TASK = int(os.getenv("BD_N", "3"))
CELLS = [c for c in os.getenv("BD_CELLS", "K562,HUDEP-2,CD34+_HSPC,HEK293").split(",") if c]
OUT = os.getenv("BD_OUT", "band_dupe.json")
HK = ["niome_hotkey"] + [f"niome_hotkey{i}" for i in range(1, 11)]
RNG = random.Random(int(os.getenv("BD_RNG", "11")))
TRIALS = 2000


def layout_for(cell):
    """{hotkey: [seed, ...]} for the layout under test."""
    if LAYOUT == "fallback":
        got = {h: JW.fallback_for(cell, h) for h in HK}
    else:
        plan = json.load(open("data/window_plan.json"))
        got = {h: sp for h, sp in ((plan.get("assignments") or {}).get(cell) or {}).items()}
    return {h: sorted({s for a, b in sp for s in range(a, b + 1)}) for h, sp in got.items() if sp}


def tasks_for(cell, n):
    """The n most recent three-seed contracts for this cell type."""
    d = json.load(open("sd_task_listing.json"))
    items = d if isinstance(d, list) else (d.get("items") or [])
    got = []
    for t in items:
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != cell:
            continue
        if len([x for x in str(c.get("seed", "")).split(",") if x.strip().isdigit()]) != 3:
            continue
        got.append((t.get("created_at", ""), t["id"]))
    return [tid for _c, tid in sorted(got, reverse=True)[:n]]


def band_for(seeds, contract, reference, cell_types, ctx, sites, cell):
    """The clean band of the min-union Cas12a group over `seeds` — the same cfg build_for_cell makes."""
    base = AH.config_for(cell)
    span = len(seeds)
    cfg = _dc.replace(base, hdr_range=(seeds[0], seeds[-1]), seed_list=tuple(seeds),
                      main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span))
    if span > 100:
        cfg = _dc.replace(cfg, variants=min(cfg.variants, AH.WIDE_WINDOW_VARIANTS))
    cfg = _dc.replace(cfg, light_cell_rows=AH.resolve_light(cfg.light_cell_rows, contract))

    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AH.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        if not bank:
            return None, {"reason": "empty bank"}
        AH.save_bank(path, bank)
    records = AH.load_bank(path)
    if len(records) < cfg.group_size:
        return None, {"reason": f"bank {len(records)} < group {cfg.group_size}"}

    space = cfg.band_seeds()
    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg),
                        seeds=space)
    index, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
    bad = set()
    for rec in (records[i] for i in index):
        bad.update(int(x) for x in rec["fails"])
    MT.free_gpu_memory()
    band = sorted(set(int(x) for x in space) - bad)
    return band, {"bank": len(records), "mf": cfg.main_max_fail, "group": cfg.group_size}


def stats(bands, space):
    """Duplication of a {hotkey: band} mapping against its joined space."""
    sizes = [len(b) for b in bands.values()]
    total = sum(sizes)
    seen = {}
    for b in bands.values():
        for s in b:
            seen[s] = seen.get(s, 0) + 1
    union = len(seen)
    mult = {}
    for c in seen.values():
        mult[c] = mult.get(c, 0) + 1
    pairs = [len(set(a) & set(b)) for a, b in itertools.combinations(bands.values(), 2)]

    # Independence control: the same band sizes placed uniformly at random in the same space.
    sp = list(space)
    sim = []
    for _ in range(TRIALS):
        u = set()
        for k in sizes:
            u.update(RNG.sample(sp, min(k, len(sp))))
        sim.append(len(u))
    exp_u = st.mean(sim) if sim else 0.0

    return {"hotkeys": len(bands), "sizes": sizes, "sum": total, "union": union,
            "space": len(sp),
            "dup_rate": (1 - union / total) if total else 0.0,
            "exp_union_independent": round(exp_u, 1),
            "exp_dup_rate": (1 - exp_u / total) if total else 0.0,
            "multiplicity": {str(k): v for k, v in sorted(mult.items())},
            "mean_pair_overlap": round(st.mean(pairs), 2) if pairs else 0.0,
            "max_pair_overlap": max(pairs) if pairs else 0,
            "p_round_hit": 1 - (1 - union / 900) ** 3,
            "p_round_hit_if_disjoint": 1 - (1 - min(total, 900) / 900) ** 3}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = {"layout": LAYOUT, "n_tasks": N_TASK, "cells": {}}
    print(f"layout: {LAYOUT}   {N_TASK} contracts per cell   band = Cas12a min-union clean set\n")

    for cell in CELLS:
        wins = layout_for(cell)
        if not wins:
            print(f"{cell}: no windows in the {LAYOUT} layout — skipped\n")
            continue
        space = sorted({s for sd in wins.values() for s in sd})
        width = len(next(iter(wins.values())))
        print(f"{cell}  {len(wins)} hotkeys, width {width}, joined space {len(space)} seeds")
        per_task = []
        for tid in tasks_for(cell, N_TASK):
            task, contract, reference = task_content(tid)
            ctx = G.build_context(contract, reference, cell_types)
            sites = G.enumerate_sites(ctx, 3000, (20, 23))
            bands, started = {}, time.monotonic()
            for hk in HK:
                if hk not in wins:
                    continue
                t0 = time.monotonic()
                band, info = band_for(wins[hk], contract, reference, cell_types, ctx, sites, cell)
                if band is None:
                    print(f"    {tid[:8]} {hk:<15} declined: {info.get('reason')}", flush=True)
                    continue
                bands[hk] = band
                print(f"    {tid[:8]} {hk:<15} band {len(band):>3}  bank {info['bank']:>6}"
                      f"  mf {info['mf']:>3}  {time.monotonic() - t0:.0f}s", flush=True)
            if len(bands) < 2:
                continue
            got = stats(bands, space)
            got.update(task=tid, elapsed_s=round(time.monotonic() - started, 1))
            per_task.append(got)
            print(f"  {tid[:8]}  bands {got['sizes']}  sum {got['sum']}  union {got['union']}"
                  f"  dup {got['dup_rate']:.1%}  (independent {got['exp_dup_rate']:.1%})"
                  f"  P(hit) {got['p_round_hit']:.1%}  {got['elapsed_s']:.0f}s")
            out["cells"][cell] = per_task
            json.dump(out, open(OUT, "w"), indent=1)
        if per_task:
            print(f"  -> mean dup {st.mean(t['dup_rate'] for t in per_task):.1%}"
                  f"   independent {st.mean(t['exp_dup_rate'] for t in per_task):.1%}"
                  f"   mean union {st.mean(t['union'] for t in per_task):.1f}"
                  f"   mean P(hit) {st.mean(t['p_round_hit'] for t in per_task):.1%}\n")
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
