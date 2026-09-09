#!/usr/bin/env python3
"""cas12a_union.py — is the all-cut clean set reachable underneath an all-HDR band?

The two builders differ only in what their Cas9 half is required to do:

* all-cut  — Cas9 must be cut-STRICT over the whole Cas12a clean set (559 of 900), so Cas9 rows
  add nothing to the failed-seed union and ``final clean = Cas12a-clean``.
* all-HDR  — Cas9 must be HDR on the 12 band seeds only, so outside the band its rows cut-fail
  freely: 170 rows x ~9 fails = ~1530 events, which is why the measured cut-clean set is 17 of 900.

CLAUDE.md prices what closing that gap is worth: an HDR band sitting on an all-cut clean floor
scores (1 + 2*0.237)/3 = 0.459 at k=1, final 139 on 9ac1d178's field, against the 121.7 (rank 11)
the shipped build got. The falsified table rejects the conjunction on the Cas9 conditional fill
capping the band at 3 seeds -- but that was measured with the scan's early exit in place, so the
pool size it saw is where the scan STOPPED, not what exists.

Pool survival is the whole question. Keeping `want` Cas9 rows out of a pool of N, when each row must
additionally cut on every seed of a clean set C, needs roughly

    N * 0.99**|C| >= want          (Cas9 cut_p is 0.99 at the energy clamp)

so N = 1370 supports |C| ~ 208 and |C| = 559 needs N ~ 46,000. This measures N with the early exit
removed, then computes the exact frontier: for each candidate clean set C, how many pool guides are
genuinely cut-strict over it, seed by seed rather than through that 0.99**|C| approximation.

Reports the achievable |C| at the real per-cell row targets, and what that |C| is worth in round
consistency using the measured per-seed value ladder.
"""
import os

os.environ.setdefault("NIOME_INSTANCE", "cas12a_union")

import dataclasses
import json
import logging
import sys
import time
from collections import Counter

import numpy as np

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
logging.basicConfig(level=logging.ERROR)

from niome_subnet.genomics import all_cut as AC
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics import seed_agnostic as SA
import genExp as G
from niome_subnet.genomics.validation import stage3

ARCH = os.environ.get("CU_ARCH",
                      "data/inst/niome_hotkey3/result/2026-09-09T07:50:28")
GROUPS = [int(x) for x in os.environ.get("CU_GROUPS", "42,80").split(",")]
SEEDS900 = np.arange(100, 1000, dtype=np.int64)
# Per-seed consistency, measured (CLAUDE.md "there is a third per-seed value"):
V_BAND, V_CLEAN_HDR, V_CLEAN_CUT, V_DIRTY = 1.0000, 0.162, 0.237, 0.101


def load_case():
    contract = json.load(open(os.path.join(ARCH, "contract.json")))
    reference = json.load(open(os.path.join(ARCH, "hbb_reference.json")))
    import urllib.request
    cell_types = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/data/cell-types?format=json"))
    return contract, reference, cell_types


def cut_fails_for(records, contract, cell_types, ctx, tag):
    """Cut-failed seeds over 100-999 for each record, grouped by its (mutation, cas, site)."""
    cell = contract.get("cell_type")
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    by_site = {}
    for r in records:
        by_site.setdefault((r["mutation"], r["cas_system"], r["start"], r["strand"],
                            r["length"]), []).append(r)
    out, done, started = {}, 0, time.monotonic()
    for (mutation, cas, start, strand, length), group in by_site.items():
        distance = abs(start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)

        class _S:  # _params_fn only reads these
            pass
        site = _S(); site.cas = cas; site.start = start; site.strand = strand; site.length = length
        params_of = AC._params_fn(site, distance, acc, offset)
        guides = [r["guide"] for r in group]
        res = MT.screen_guides_gpu(guides, SEEDS900, mutation, cas, start, strand,
                                   lambda g: params_of(g)[2], 10 ** 9)
        for g in guides:
            out[(mutation, cas, start, strand, g)] = res.get(g, np.empty(0, dtype=np.int64))
        done += len(guides)
        if done % 5000 < len(guides):
            print("    [%s] cut-fails for %d/%d guides (%.0fs)"
                  % (tag, done, len(records), time.monotonic() - started), flush=True)
    MT.free_gpu_memory()
    return out


def full_cas9_pool(clean_band, contract, cell_types, ctx, sites, cfg):
    """all_hdr.scan_cas9 with the early exit REMOVED — the whole HDR-on-band pool."""
    cell = contract.get("cell_type")
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    jobs = [(i, m, abs(sites[i].start - ctx.mutation_map[m]))
            for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(sites[i].start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda j: j[2])
    found, started = [], time.monotonic()
    for n, (site_index, mutation, distance) in enumerate(jobs, 1):
        site = sites[site_index]
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = AH._params_fn(site, distance, acc, offset)
        for guide in MT.screen_guides_rule_gpu(guides, clean_band, mutation, "Cas9", site.start,
                                               site.strand, params_of, "hdr", 0):
            found.append({"guide": guide, "mutation": mutation, "cas_system": "Cas9",
                          "strand": site.strand, "start": site.start, "length": site.length})
        if n % 100 == 0:
            print("    cas9 scan %d/%d targets | pool %d | %.0fs"
                  % (n, len(jobs), len(found), time.monotonic() - started), flush=True)
    MT.free_gpu_memory()
    return found


def frontier(pool_fails, c12_clean, want, need):
    """Largest clean set C (a subset of the Cas12a clean set) keeping enough strict Cas9 rows.

    Removes from C the seed that the most currently-blocked candidates fail on -- the seed whose
    removal unlocks the most pool -- and records the frontier as C shrinks. Two thresholds:

    * loose      -- `want` strict candidates in total, all four (mutation, strand) cells non-empty
    * production -- every cell also at `need` = cas9_cell_target, which is what scan_cas9's break
                    enforces and what `assemble` actually apportions
    """
    C = set(int(x) for x in c12_clean)
    fails = {k: set(int(x) for x in v) for k, v in pool_fails.items()}
    curve = []
    while True:
        strict = [k for k, f in fails.items() if not (f & C)]
        per_cell = Counter((k[0], k[3]) for k in strict)
        loose = len(strict) >= want and len(per_cell) >= 4
        prod = loose and min(per_cell.values()) >= need
        curve.append((len(C), len(strict), dict(per_cell), loose, prod))
        if prod or not C:
            return curve
        blocked = [f for k, f in fails.items() if (f & C)]
        if not blocked:
            return curve
        tally = Counter()
        for f in blocked:
            for s in (f & C):
                tally[s] += 1
        worst, _ = tally.most_common(1)[0]
        C.discard(worst)


def main():
    contract, reference, cell_types = load_case()
    cell = contract["cell_type"]
    print("contract cell_type=%s mutations=%s" % (cell, list((contract.get('mutation_weights') or {}).keys())))
    base = AH.config_for(cell) or AH.AllHdrConfig()
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments

    for gs in GROUPS:
        cfg = dataclasses.replace(base, group_size=gs)
        print()
        print("=" * 78)
        print("group_size %d   band %d-%d   main_max_fail %d"
              % (gs, cfg.start_seed, cfg.end_seed, cfg.main_max_fail))
        os.makedirs(AH.HDR_BANK_DIR, exist_ok=True)
        path = os.path.join(AH.HDR_BANK_DIR,
                            "cas12a-%s.npz" % AC.bank_key(contract, cell_types, cfg))
        if not os.path.exists(path):
            t = time.monotonic()
            bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
            if not bank:
                print("  bank empty — skipping"); continue
            AC.save_bank(path, bank)
            print("  bank built in %.0fs" % (time.monotonic() - t))
        records = AC.load_bank(path)
        print("  HDR bank: %d guides" % len(records))

        band_space = cfg.band_seeds()
        sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                            per_cell_min=cfg.per_cell_min,
                            caps=AH._group_caps(contract, ctx, cfg))
        idx, _u = sel.best(gs, restarts=cfg.restarts)
        group = [records[i] for i in idx]
        bad = set()
        for r in group:
            bad.update(int(x) for x in r["fails"])
        clean_band = np.array(sorted(set(int(x) for x in band_space) - bad), dtype=np.int64)
        print("  HDR band: %d seeds %s" % (clean_band.size, list(clean_band[:14])))
        if clean_band.size == 0:
            print("  no band — skipping"); continue

        # --- Cas12a half's CUT behaviour: the ceiling on any final clean set -------------
        cf = cut_fails_for(group, contract, cell_types, ctx, "cas12a")
        u12 = set()
        for r in group:
            u12 |= set(int(x) for x in cf[(r["mutation"], r["cas_system"], r["start"],
                                           r["strand"], r["guide"])])
        c12 = sorted(set(int(x) for x in SEEDS900) - u12)
        per = [len(cf[(r["mutation"], r["cas_system"], r["start"], r["strand"], r["guide"])])
               for r in group]
        print("  Cas12a cut-fails: mean %.1f/row  union %d  => Cas12a-clean %d of 900 (%.1f%%)"
              % (np.mean(per), len(u12), len(c12), 100 * len(c12) / 900))
        print("     (all-cut reaches 559 of 900 at group 42 — that is the target)")

        # --- Cas9 pool with the early exit removed --------------------------------------
        want = n_rows - gs
        t = time.monotonic()
        pool = full_cas9_pool(clean_band, contract, cell_types, ctx, sites, cfg)
        print("  FULL Cas9 HDR-on-band pool: %d candidates in %.0fs (the shipped scan stops at "
              "~%d)" % (len(pool), time.monotonic() - t, want * cfg.pool_target))
        if not pool:
            print("  empty pool — skipping"); continue
        need = AC.cas9_cell_target(contract, ctx, cfg, want)
        print("  Cas9 rows needed: %d  (per-cell floor %d)" % (want, need))
        approx = int(np.log(max(want, 1) / len(pool)) / np.log(0.99)) if len(pool) > want else 0
        print("  0.99**|C| approximation says |C| ~ %d is supportable" % approx)

        pf = cut_fails_for(pool, contract, cell_types, ctx, "cas9")
        keyed = {(r["mutation"], r["cas_system"], r["start"], r["strand"], r["guide"]):
                 pf[(r["mutation"], r["cas_system"], r["start"], r["strand"], r["guide"])]
                 for r in pool}
        curve = frontier(keyed, c12, want, need)
        print("  exact frontier (|C| shrinking from the Cas12a-clean set):")
        shown = 0
        for size, strict, percell, loose, prod in curve:
            mark = "<== PRODUCTION-OK" if prod else ("<== loose-OK" if loose else "")
            if shown < 5 or mark or size % 50 == 0:
                print("     |C| = %-4d strict Cas9 = %-6d cells %-2d min/cell %-5d %s"
                      % (size, strict, len(percell),
                         min(percell.values()) if percell else 0, mark))
                shown += 1
            if prod:
                break
        loose_best = next((c for c in curve if c[3]), None)
        prod_best = next((c for c in curve if c[4]), None)
        print("  ACHIEVABLE final clean set: loose %d of 900 | production %d of 900"
              % (loose_best[0] if loose_best else 0, prod_best[0] if prod_best else 0))
        final_clean = prod_best[0] if prod_best else 0

        for label, vclean in (("all-HDR rows (0.162)", V_CLEAN_HDR),
                              ("all-cut rows (0.237)", V_CLEAN_CUT)):
            p = final_clean / 900.0
            off = p * vclean + (1 - p) * V_DIRTY
            for k in (0, 1):
                cons = (k * V_BAND + (3 - k) * off) / 3
                print("     k=%d, %-22s -> off-band seed %.4f  round cons %.4f  final %.1f"
                      % (k, label, off, cons, cons * 303.6))
        print("     (9ac1d178 cutoff was 126.32; the shipped build scored 121.74 at cons 0.401)")


if __name__ == "__main__":
    main()
