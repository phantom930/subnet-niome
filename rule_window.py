#!/usr/bin/env python3
"""rule_window.py — the not_mhnhej rule at a 225-seed window: band, build, and what a hit is worth.

``not_mhnhej`` = {"any": ("HDR", "BLUNT_NHEJ")}: a row complies unless it repairs by MH_NHEJ, and a
no-cut row fails it. Blended per-row compliance is ~0.82 against ``hdr``'s ~0.57, and the band bound
is ``pool * P**B >= group_size``, so the band is far wider -- CLAUDE.md records 45 against hdr's 13.
It was rejected on value: only one of stage 4's three targets is pinned, so a band seed scores
consistency ~0.123, and even three band seeds reach 25.7, under every cutoff.

Two things make it worth re-checking at width 225 rather than 100:

* the width sweep this session showed band size is set by the *pool* clearing the screen, not by
  the window, so a wider window neither helps nor hurts until it clips the band -- at 45 seeds
  ``not_mhnhej`` is nowhere near clipping at 100, and 225 tests whether the band grows when given
  room and a pool sized to match;
* the leaders' windows measure 200-500 wide, so 225 is the width their geometry implies.

``main_max_fail`` is set near the per-row failure mean for each width (rather than carried over),
because the screen only has to admit guides compatible with a wide band and ``load_bank`` keeps the
60,000 lowest-fail regardless.

Reports the band, whether the build completes, and the scored value of a band seed against an
off-band seed -- the number that decides whether a wide band is worth anything.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import dataclasses
import json
import sys
import time
from collections import Counter

import numpy as np

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.all_cut import _params_fn, load_bank, save_bank
from niome_subnet.genomics.validation import stage3
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

TASKS = ["task-k562.json", "task-cd34.json"]
# (window, main_max_fail) -- max_fail near the per-row failure mean for that width
PLAN = [((700, 799), 20), ((700, 924), 45)]
RULE = os.getenv("WIN20_RULE", "not_mhnhej")
BANK_DIR = "data/rule_window"


def bank(contract, reference, cell_types, ctx, sites, cfg, rule):
    cell = contract.get("cell_type")
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    seeds = np.arange(cfg.start_seed, cfg.end_seed + 1, dtype=np.int64)
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == "Cas12a" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    out = []
    for site_index, mutation in jobs:
        site = sites[site_index]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas12a_gc[0], cfg.cas12a_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params = _params_fn(site, distance, acc, offset)
        for guide, fails in MT.screen_guides_rule_gpu(
                guides, seeds, mutation, site.cas, site.start, site.strand, params,
                rule, cfg.main_max_fail).items():
            out.append({"guide": guide, "mutation": mutation, "cas_system": site.cas,
                        "strand": site.strand, "start": site.start, "length": site.length,
                        "fails": fails.astype(np.int16)})
    MT.free_gpu_memory()
    counts = np.asarray([len(b["fails"]) for b in out])
    order = np.argsort(counts)[:cfg.bank_keep]
    return [out[int(i)] for i in order]


def cas9_for(clean, contract, cell_types, ctx, sites, cfg, rule, want):
    cell = contract.get("cell_type")
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    jobs = [(i, m, abs(sites[i].start - ctx.mutation_map[m]))
            for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(sites[i].start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda j: j[2])
    found = []
    for site_index, mutation, distance in jobs:
        site = sites[site_index]
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params = _params_fn(site, distance, acc, offset)
        for guide in MT.screen_guides_rule_gpu(guides, clean, mutation, "Cas9", site.start,
                                               site.strand, params, rule, 0):
            gc, _e, _c = params(guide)
            found.append({"guide": guide, "mutation": mutation, "cas_system": "Cas9",
                          "strand": site.strand, "start": site.start, "length": site.length,
                          "gc": gc, "distance": distance})
        if len(found) >= want * cfg.pool_target and len({(f["mutation"], f["strand"])
                                                         for f in found}) == 4:
            break
    MT.free_gpu_memory()
    return found


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    os.makedirs(BANK_DIR, exist_ok=True)
    print(f"rule = {RULE}\n")
    for tf in TASKS:
        task = json.load(open(tf))
        contract = dict(task["content"]["contract"])
        reference = task["content"]["hbb_reference"]
        cell = contract["cell_type"]
        base = AH.config_for(cell)
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        n_rows = contract["rules"].get("max_experiments") or 250
        print(f"=== {tf}  {cell} ===", flush=True)
        for window, mf in PLAN:
            width = window[1] - window[0] + 1
            cfg = dataclasses.replace(base, hdr_range=window, main_max_fail=mf)
            t0 = time.monotonic()
            path = os.path.join(BANK_DIR,
                                f"{cell}-{RULE}-{window[0]}_{window[1]}-mf{mf}.npz")
            if not os.path.exists(path):
                b = bank(contract, reference, cell_types, ctx, sites, cfg, RULE)
                if not b:
                    print(f"  w{width:<4} mf{mf:<3} empty bank", flush=True)
                    continue
                save_bank(path, b)
            recs = load_bank(path)
            disk = int(np.load(path, allow_pickle=False)["fails"].shape[0])
            if len(recs) < cfg.group_size:
                print(f"  w{width:<4} mf{mf:<3} bank {len(recs)} short of group", flush=True)
                continue
            sel = FG.FastGreedy(recs, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                                per_cell_min=cfg.per_cell_min,
                                caps=AH._group_caps(contract, ctx, cfg))
            idx, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
            badset = set()
            for i in idx:
                badset.update(int(x) for x in recs[i]["fails"])
            clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - badset),
                             dtype=np.int64)
            want = n_rows - cfg.group_size
            cas9 = cas9_for(clean, contract, cell_types, ctx, sites, cfg, RULE, want)
            cells4 = len({(r["mutation"], r["strand"]) for r in cas9})
            line = (f"  w{width:<4} mf{mf:<3} disk {disk:>7} bank {len(recs):>6} "
                    f"band {clean.size:>4}/{width} ({clean.size/width*100:>4.1f}%) "
                    f"cas9 {len(cas9):>5}")
            if len(cas9) < want or cells4 < 4:
                print(line + f"  DECLINED (need {want} over 4 cells, got {cells4})", flush=True)
                continue
            rows = AH.assemble([recs[i] for i in idx], cas9, contract, ctx, cfg, n_rows)
            ncells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
            band_seed = int(clean[0])
            off = next(s for s in range(cfg.start_seed, cfg.end_seed + 1)
                       if s not in set(int(x) for x in clean))
            doc = dict(contract)
            doc["seed"] = f"{band_seed}"
            json.dump(doc, open(settings.CONTRACT_PATH, "w"))
            json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
            json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
            run_stage12(cell_types)
            res = {}
            for label, s in (("band", band_seed), ("off", off)):
                run_stage3(seed=s)
                run_stage4(seed=s)
                res[label] = run_stage5()
            k = lambda lab, key: res[lab][key]
            print(line + f" rows {len(rows)} cells {ncells}/8 {time.monotonic()-t0:>5.0f}s",
                  flush=True)
            print(f"        band seed {band_seed}: cons {k('band','consistency_factor'):.4f} "
                  f"wtd {k('band','total_weighted_score'):.1f} "
                  f"fid {k('band','distribution_fidelity_factor'):.3f} "
                  f"final {k('band','final_score'):.2f}", flush=True)
            print(f"        off  seed {off}: cons {k('off','consistency_factor'):.4f} "
                  f"final {k('off','final_score'):.2f}", flush=True)
            c1, c0 = k('band', 'consistency_factor'), k('off', 'consistency_factor')
            w, f = k('band', 'total_weighted_score'), k('band', 'distribution_fidelity_factor')
            for kk in (1, 2, 3):
                rc = (kk * c1 + (3 - kk) * c0) / 3
                print(f"          k={kk}: round cons {rc:.3f} -> final {w*rc*f:>6.1f}", flush=True)
        print(flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
