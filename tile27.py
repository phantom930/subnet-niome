#!/usr/bin/env python3
"""tile27.py — 11 DISJOINT slices of the joined space, against the shipped overlapping rotation.

Conditioned on the three predicted classes being right, the three seeds are i.i.d. uniform on those
300 seeds -- conditioning i.i.d. draws on a subset leaves them uniform on it -- so WHERE a band sits
is worth nothing and only two things pay:

    P(some hotkey hits >=1)   set by the UNION of the bands
    P(some hotkey hits >=2)   set by the individual band SIZES, as b^2

Both want big bands, and the union also wants them disjoint. The shipped layout is 11 rotated
width-75..150 slices at stride 30, which overlap heavily, so its union is the independence-model
~S(1-(1-b/S)^11) rather than the sum. Eleven windows of 300/11 = 27 tile the space with no overlap,
which makes the bands disjoint BY CONSTRUCTION -- and disjoint bands also make the k>=2 events
mutually exclusive (two hotkeys at k>=2 would need four seeds), so those probabilities simply add.

The open question is whether the band survives a 27-seed window: it is measured flat in width at
group 42 (16 seeds at widths 30/75/100/150 on the erythroid cells) but the Cas9 conditional fill
costs P(HDR) per band seed, so a band that grows toward the window can starve it.

    T27_GROUP=42 T27_WIDTH=27 python -u tile27.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "tile27")

import dataclasses as _dc   # noqa: E402
import json                 # noqa: E402
import logging              # noqa: E402
import math                 # noqa: E402
import time                 # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G          # noqa: E402
from niome_subnet.genomics import all_hdr as AH     # noqa: E402
from niome_subnet.genomics import fastgreedy as FG  # noqa: E402
from niome_subnet.genomics import mt19937 as MT     # noqa: E402
from sd_task import task_content                    # noqa: E402

WIDTH = int(os.getenv("T27_WIDTH", "27"))
GROUP = int(os.getenv("T27_GROUP", "42"))
H, S, NS = 11, 300, 3
OUT = os.getenv("T27_OUT", "tile27.json")
TASKS = {"K562": "a3e302de", "HUDEP-2": "72fb2cd4",
         "CD34+_HSPC": "53671b0f", "HEK293": "b3a077d1"}


def full_id(prefix):
    d = json.load(open("sd_task_listing.json"))
    items = d if isinstance(d, list) else (d.get("items") or [])
    return next(t["id"] for t in items if t["id"].startswith(prefix))


def classes_for(cell):
    plan = json.load(open("data/window_plan.json"))["assignments"][cell]
    seeds = sorted({s for sp in plan.values() for a, b in sp for s in range(a, b + 1)})
    return seeds


def band_for(seeds, contract, reference, cell_types, ctx, sites, cell, group):
    base = AH.config_for(cell)
    span = len(seeds)
    cfg = _dc.replace(base, hdr_range=(seeds[0], seeds[-1]), seed_list=tuple(seeds),
                      group_size=group,
                      main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span,
                                                        AH._native_span(cell)))
    if span > 100:
        cfg = _dc.replace(cfg, variants=min(cfg.variants, AH.WIDE_WINDOW_VARIANTS))
    cfg = _dc.replace(cfg, light_cell_rows=AH.resolve_light(cfg.light_cell_rows, contract))
    p = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AH.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(p):
        bk = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        if not bk:
            return None, {"reason": "empty bank"}
        AH.save_bank(p, bk)
    recs = AH.load_bank(p)
    if len(recs) < cfg.group_size:
        return None, {"reason": f"bank {len(recs)} < group {cfg.group_size}"}
    sp = cfg.band_seeds()
    sel = FG.FastGreedy(recs, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg),
                        seeds=sp)
    idx, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
    bad = set()
    for r in (recs[i] for i in idx):
        bad.update(int(x) for x in r["fails"])
    clean = sorted(set(int(x) for x in sp) - bad)
    # Does the Cas9 conditional fill survive a band this wide? That is what caps the band.
    import numpy as np
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    cas9 = AH.scan_cas9(np.array(clean, dtype=np.int64), contract, cell_types, ctx, sites, cfg,
                        n_rows - cfg.group_size, None)
    MT.free_gpu_memory()
    return clean, {"bank": len(recs), "mf": cfg.main_max_fail, "cas9": len(cas9),
                   "need": n_rows - cfg.group_size, "buildable": len(cas9) >= n_rows - cfg.group_size}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = {"width": WIDTH, "group": GROUP, "cells": {}}
    print(f"11 disjoint width-{WIDTH} slices of the 300-seed joined space, group {GROUP}\n")
    for cell, pref in TASKS.items():
        seeds = classes_for(cell)
        tid = full_id(pref)
        _t, contract, reference = task_content(tid)
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        bands, infos = {}, {}
        for h in range(H):
            sub = seeds[h * WIDTH:(h + 1) * WIDTH]
            if len(sub) < WIDTH:
                break
            t0 = time.monotonic()
            band, info = band_for(sub, contract, reference, cell_types, ctx, sites, cell, GROUP)
            if band is None:
                print(f"  h{h} {sub[0]}-{sub[-1]}  declined: {info.get('reason')}", flush=True)
                continue
            bands[h] = band
            infos[h] = info
            print(f"  h{h:<2} {sub[0]:>3}-{sub[-1]:<3} band {len(band):>3}/{WIDTH}"
                  f"  cas9 {info['cas9']:>5}/{info['need']}"
                  f"  {'OK' if info['buildable'] else 'STARVED'}  {time.monotonic()-t0:.0f}s",
                  flush=True)
        if not bands:
            continue
        ok = [b for h, b in bands.items() if infos[h]["buildable"]]
        union = len({s for b in ok for s in b})
        sizes = [len(b) for b in ok]
        c = math.comb
        p1 = 1 - c(S - union, NS) / c(S, NS) if union < S else 1.0
        p2 = sum((c(b, 2) * (S - b) + c(b, 3)) / c(S, NS) for b in sizes)
        out["cells"][cell] = {"task": tid, "sizes": sizes, "union": union,
                              "p_hit": p1, "p_double": p2,
                              "buildable": len(ok), "of": len(bands)}
        print(f"  -> {len(ok)}/{len(bands)} buildable  sizes {sizes}  union {union}"
              f"  P(>=1) {p1:.1%}  P(>=2) {p2:.2%}\n")
        json.dump(out, open(OUT, "w"), indent=1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
