#!/usr/bin/env python3
"""overlap_test.py — can two hotkeys share a seed window and still cover different seeds?

The fleet tiles disjoint narrow windows because a hotkey's band is ~14 of its 16 seeds — near
saturation. The competing hypothesis, from the leaders' simultaneous-spike spans (79-494, median
242), is that they run WIDE overlapping windows instead.

Overlap only helps if two hotkeys on the same window find DIFFERENT bands. Today they cannot:
FastGreedy.best() takes seed=0, so the build is deterministic and siblings on one window would
submit byte-identical rows. This measures whether the freedom exists at all, by running the
min-union with different selection seeds and comparing the resulting clean sets.

Only the Cas12a half is needed: the band IS the min-union's clean set (Cas9 rows are then required
clean over all of it), so this skips the Cas9 scan and costs one bank per window.

    ARM A  disjoint narrow   6 x width 16, tiled          <- what the fleet runs today
    ARM B  shared wide       6 variants on ONE width-W window
    ARM C  staggered wide    6 x width W, offset by W/2   <- true overlap
"""
import os
os.environ["NIOME_INSTANCE"] = "overlap"

import dataclasses, itertools, json, os.path, sys, time
sys.argv = ["x"]; sys.path.insert(0, "/root/workspace/subnet-niome")
import logging; logging.basicConfig(level=logging.ERROR)

import numpy as np
import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank

TASK = os.getenv("OV_TASK", "task-b9051bc7.json")
WIDE = int(os.getenv("OV_WIDE", "300"))
N = int(os.getenv("OV_N", "6"))
OUT = "overlap_test.json"
ARMS = os.getenv("OV_ARMS", "ABC").upper()


def bands_for(contract, reference, cell_types, ctx, sites, base, window, seeds_list):
    """Clean sets for one window under several selection seeds. One bank, N selections."""
    # main_max_fail is calibrated for a 100-seed window (45 of 100). Held fixed it demands HDR on
    # 255 of 300 seeds at width 300, which no guide meets — the bank comes back empty and the width
    # looks impossible when it is only over-screened. Scale it with the span, as the earlier width
    # sweep had to.
    span = window[1] - window[0] + 1
    cfg = dataclasses.replace(base, hdr_range=window,
                              main_max_fail=max(1, round(base.main_max_fail * span / 100)))
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return None, "bank scan produced nothing"
        save_bank(path, bank)
    records = load_bank(path)
    if len(records) < cfg.group_size:
        return None, f"bank {len(records)} < group {cfg.group_size}"
    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg))
    out = []
    for sd in seeds_list:
        index, _ = sel.best(cfg.group_size, restarts=cfg.restarts, seed=sd)
        bad = set()
        for i in index:
            bad.update(int(x) for x in records[i]["fails"])
        out.append(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad))
    return out, None


def report(name, bands, note=""):
    sets = [set(b) for b in bands]
    union = set().union(*sets) if sets else set()
    total = sum(len(s) for s in sets)
    pair = [len(a & b) for a, b in itertools.combinations(sets, 2)]
    p = 1 - (1 - len(union) / 900) ** 3
    print(f"\n{name}  {note}")
    print(f"  per-hotkey band sizes : {[len(s) for s in sets]}")
    print(f"  sum of bands          : {total}")
    print(f"  UNION (distinct seeds): {len(union)}   duplication {total - len(union)}")
    if pair:
        print(f"  pairwise overlap      : min {min(pair)} median "
              f"{sorted(pair)[len(pair)//2]} max {max(pair)}")
    print(f"  P(>=1 of 3 seeds hits): {p:.1%}")
    return {"bands": [len(s) for s in sets], "sum": total, "union": len(union),
            "pairwise_overlap": pair, "p_hit": p}


def main():
    task = json.load(open(TASK))
    contract, reference = task["content"]["contract"], task["content"]["hbb_reference"]
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types(); G.load_sequence()
    base = AH.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    print(f"{task['id'][:8]}  {cell}  group {base.group_size}  N={N}  wide width {WIDE}")
    res = {}

    t0 = time.monotonic()
    narrow = []
    for i in (range(N) if "A" in ARMS else []):
        lo = 900 + i * 16
        b, err = bands_for(contract, reference, cell_types, ctx, sites, base, (lo, lo + 15), [0])
        if err:
            print(f"  ARM A window {lo}-{lo+15} declined: {err}"); continue
        narrow.append(b[0])
    if narrow:
        res["A_disjoint_narrow"] = report(
            "ARM A  disjoint narrow (width 16, the fleet today)",
            narrow, f"{time.monotonic()-t0:.0f}s")

    t0 = time.monotonic()
    b, err = (bands_for(contract, reference, cell_types, ctx, sites, base,
                        (700, 700 + WIDE - 1), list(range(N)))
              if "B" in ARMS else (None, "skipped"))
    if err:
        print(f"\nARM B declined: {err}")
    else:
        res["B_shared_wide"] = report(
            f"ARM B  shared wide (one {WIDE}-seed window, {N} selection seeds)",
            b, f"{time.monotonic()-t0:.0f}s")

    t0 = time.monotonic()
    stag, step = [], max(1, WIDE // 2)
    for i in (range(N) if "C" in ARMS else []):
        lo = 100 + i * step
        hi = min(999, lo + WIDE - 1)
        if hi - lo + 1 < WIDE:
            break
        bb, err = bands_for(contract, reference, cell_types, ctx, sites, base, (lo, hi), [0])
        if err:
            print(f"  ARM C window {lo}-{hi} declined: {err}"); continue
        stag.append(bb[0])
    if stag:
        res["C_staggered_wide"] = report(
            f"ARM C  staggered wide (width {WIDE}, offset {step})",
            stag, f"{time.monotonic()-t0:.0f}s")

    json.dump(res, open(OUT, "w"), indent=1)
    if len(res) > 1:
        print("\n=== verdict ===")
        for k, v in res.items():
            print(f"  {k:<22} union {v['union']:>4}  P(hit) {v['p_hit']:.1%}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
