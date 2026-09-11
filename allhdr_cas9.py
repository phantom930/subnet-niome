#!/usr/bin/env python3
"""allhdr_cas9.py — the all-HDR band from a Cas9-ONLY group, width 30, K562.

The mirror of `allhdr_cas12a.py`. The shipped construction min-unions **Cas12a** and then requires
HDR of the Cas9 half over the resulting clean band; this asks what a group built the other way
round reaches, with no Cas12a rows at all.

Cas9 is the better half for this rule on both factors:

    cut_p   = min(0.99, 0.86 + 0.18*energy)    vs Cas12a's 0.78 base -> 0.99 against 0.96
    hdr_w   = 0.32 + 0.35*energy               vs Cas12a's 0.24      -> a larger HDR share

so P(HDR) ~ 0.57 against Cas12a's ~0.49 (CLAUDE.md's figure for the conditional Cas9 fill). Since a
bank of N holds about `N * p**k` guides sharing a k-seed band, the same 300,000-guide bank should
support a wider band here -- `ln(250/300000)/ln(0.57)` ~ 13 at group 250 against ~10 for Cas12a.

`all_hdr.build_bank` only banks Cas12a sites, so the scan is mirrored below for Cas9. It is cached
under a `cas9-` filename prefix: `bank_key` does not fold in the cas system (it never needed to),
so a Cas9 bank written under the usual `cas12a-<key>.npz` name would silently collide with the
Cas12a bank of the same config.

**This does not escape the stage-5 problem.** `cas12a_only_score.py` measured a single-cas
submission at fidelity 0.0287-0.0294 -- the cas-coverage term is floored at 1e-9 and costs ~30x, and
joint coverage caps at 4 of 8 cells. A Cas9-only build pays exactly the same penalty; what it may
change is `total_weighted_score`, since Cas9 rows are the structurally better-scoring half (the
Cas12a-only build measured weighted 179.9 against the mixed build's 254.3).

    python allhdr_cas9.py
    AH9_MF=12,14 AH9_GROUPS=80,250 python allhdr_cas9.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "ah9")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA       # noqa: E402
from niome_subnet.genomics.all_cut import _params_fn, bank_key, load_bank, save_bank  # noqa: E402
from niome_subnet.genomics.validation import stage3         # noqa: E402
from allhdr_mf import pick                                  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402

CELL = os.getenv("AH9_CELL", "K562")
WINDOW = tuple(int(x) for x in os.getenv("AH9_WINDOW", "100-129").split("-"))
MFS = [int(x) for x in os.getenv("AH9_MF", "8,10,12,14,16,18,20").split(",")]
GROUPS = [int(x) for x in os.getenv("AH9_GROUPS", "42,80,125,170,250").split(",")]
LOAD = int(os.getenv("AH9_LOAD", "300000"))
OUT = os.getenv("AH9_OUT", "allhdr_cas9.json")


def build_cas9_bank(contract, reference, cell_types, ctx, sites, cfg):
    """`all_hdr.build_bank` for Cas9 sites: HDR on all but `main_max_fail` of the band."""
    cell = contract.get("cell_type")
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    seeds = cfg.band_seeds()
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    bank = []
    for site_index, mutation in jobs:
        site = sites[site_index]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = _params_fn(site, distance, accessibility, offset)
        for guide, fails in MT.screen_guides_rule_gpu(
                guides, seeds, mutation, site.cas, site.start, site.strand, params_of,
                "hdr", cfg.main_max_fail).items():
            bank.append({"guide": guide, "mutation": mutation, "cas_system": site.cas,
                         "strand": site.strand, "start": site.start, "length": site.length,
                         "fails": fails.astype(np.int16)})
    MT.free_gpu_memory()
    counts = np.asarray([len(b["fails"]) for b in bank]) if bank else np.array([])
    order = np.argsort(counts)[:cfg.bank_keep] if bank else []
    return [bank[int(i)] for i in order], len(jobs)


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = WINDOW
    span = hi - lo + 1
    base = AH.config_for(CELL)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_cas9 = len([s for s in sites if s.cas == "Cas9"])
    n_c12 = len([s for s in sites if s.cas == "Cas12a"])
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  window {lo}-{hi} (span {span})"
          f"  rule hdr  CAS9 ONLY  gc {base.cas9_gc}  max_distance {base.max_distance}")
    print(f"sites: {n_cas9} Cas9, {n_c12} Cas12a   "
          f"(P(HDR) ~ 0.57 Cas9 vs ~0.49 Cas12a)\n")
    print(f"{'mf':>4}{'needs':>7}{'raw bank':>10}{'loaded':>8}{'hdr/guide':>12}"
          + "".join(f"{('g'+str(g)):>7}" for g in GROUPS) + f"{'s':>7}")
    out = []
    for mf in MFS:
        t0 = time.monotonic()
        cfg = _dc.replace(base, hdr_range=(lo, hi), main_max_fail=mf)
        rec = {"cell": CELL, "window": [lo, hi], "span": span, "mf": mf, "needs": span - mf,
               "rule": "hdr", "cas": "Cas9", "cas9_gc": list(base.cas9_gc),
               "task": (task.get("task_id") or task["id"])[:8]}
        # `cas9-` prefix: bank_key does not distinguish the cas system, so this must not reuse the
        # Cas12a bank's filename.
        path = os.path.join(AH.HDR_BANK_DIR, f"cas9-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            bank, n_jobs = build_cas9_bank(contract, reference, cell_types, ctx, sites, cfg)
            rec["raw_bank"] = len(bank)
            rec["targets"] = n_jobs
            if not bank:
                rec["empty"] = True
                print(f"{mf:>4}{span-mf:>7}{0:>10}    nothing qualifies{time.monotonic()-t0:>7.0f}")
                out.append(rec); json.dump(out, open(OUT, "w"), indent=1); continue
            save_bank(path, bank)
        records = load_bank(path, limit=LOAD)
        rec["loaded_bank"] = len(records)
        hdrs = np.asarray([span - len(r["fails"]) for r in records])
        rec["hdr_per_guide"] = {"max": int(hdrs.max()), "mean": round(float(hdrs.mean()), 2)}
        line = (f"{mf:>4}{span-mf:>7}{rec.get('raw_bank', 'cached'):>10}{len(records):>8}"
                f"{rec['hdr_per_guide']['mean']:>7.1f}/{rec['hdr_per_guide']['max']:<4}")
        rec["bands"] = {}
        for g in GROUPS:
            if len(records) < g:
                rec["bands"][g] = None
                line += f"{'-':>7}"
                continue
            sel = FG.FastGreedy(records, window_lo=lo, window_hi=hi,
                                per_cell_min=cfg.per_cell_min)
            idx, _u = sel.best(g, restarts=cfg.restarts)
            bad = set()
            for i in idx:
                bad.update(int(x) for x in records[i]["fails"])
            clean = sorted(set(range(lo, hi + 1)) - bad)
            rec["bands"][g] = {"union": len(bad), "clean": len(clean), "seeds": clean}
            line += f"{len(clean):>7}"
        rec["elapsed_s"] = round(time.monotonic() - t0, 1)
        print(line + f"{rec['elapsed_s']:>7.0f}")
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}   (columns under g<N> are the CLEAN BAND for that group size)")


if __name__ == "__main__":
    main()
