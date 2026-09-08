#!/usr/bin/env python3
"""nocut_score.py — what does a no-cut band actually score?

nocut100.py answered the frequency half: a 100-seed window yields a band of **10 seeds** with 250
compliant Cas12a guides (against all-HDR's 12-16), which on band width alone looks competitive.
This answers the half that decides it, and the two pull hard against each other by construction:

    P(no cut) = 1 - min(0.99, base + 0.18*energy),  energy ~ accessibility * 1.8 * gc + ...

so avoiding the cut wants LOW energy, i.e. low GC and large distance — while term 1 wants
``gc_score = max(0, 1 - |gc - 0.5|*2)`` maximal at GC 0.50 and ``dist_score = exp(-d/base_padding)``
maximal at distance 0. The guides that comply are the guides that score worst, and no amount of
searching escapes that: it is one variable driving both.

Reports, per band width k, the best achievable ``sum(weighted_score)`` over the 250 highest-scoring
guides that comply on >= k seeds, against the unconstrained top 250 from the same pool. The ratio is
what a no-cut spike would be worth relative to all-HDR's ~290.

    python nocut_score.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else None
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "100-199"
os.environ["NIOME_INSTANCE"] = "nocutsc"

import json  # noqa: E402
import logging  # noqa: E402
import math  # noqa: E402
import time  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import mt19937 as MT  # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3  # noqa: E402
from nocut100 import pick_task  # noqa: E402

ROWS = 250
OUT = os.getenv("NOCUT_OUT", "nocut_score.json")


def main():
    task, contract, reference = pick_task() if not TASK else (None, None, None)
    if TASK:
        from sd_task import task_content
        task, contract, reference = task_content(TASK)
    lo, hi = (int(x) for x in WINDOW.split("-"))
    seeds = np.arange(lo, hi + 1, dtype=np.int64)
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context({**contract, "seed": 0}, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    cell = contract.get("cell_type")
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    weights = contract.get("mutation_weights") or {}
    regions = contract.get("mutation_regions") or {}
    cfg = AC.AllCutConfig()

    print(f"task {task['id'][:8]}  {cell} (accessibility {acc:.2f})  window {lo}-{hi}  "
          f"base_padding {ctx.base_padding}\n")

    # Screen BOTH cas systems. The Cas9 half is the harder one (base 0.86 against Cas12a's 0.78)
    # and stage 5 needs both present — an empty (mutation x cas x strand) cell costs roughly a
    # 0.03x multiplier on the whole score, so a Cas12a-only band is not a submission.
    rows = []
    for cas in ("Cas12a", "Cas9"):
        jobs = [(i, m) for i, s in enumerate(sites) if s.cas == cas for m in ctx.mutations
                if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
        print(f"screening {len(jobs)} {cas} targets")
        started, seen = time.monotonic(), 0
        for index, (site_index, mutation) in enumerate(jobs, 1):
            site = sites[site_index]
            distance = abs(site.start - ctx.mutation_map[mutation])
            offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
            guides = SA.enumerate_variants(site, ctx, 0.05, 0.95, ctx.max_mismatches, True,
                                           cfg.variants)
            if not guides:
                continue
            params_of = AC._params_fn(site, distance, acc, offset)
            seen += len(guides)
            got = MT.screen_guides_rule_gpu(guides, seeds, mutation, site.cas, site.start,
                                            site.strand, params_of, "nocut", len(seeds))
            dist_score = math.exp(-distance / ctx.base_padding)
            for guide, fails in got.items():
                gc = sum(b in "GC" for b in guide) / site.length
                gc_score = max(0.0, 1.0 - abs(gc - 0.5) * 2.0)
                # offtarget_factor is 1.0: enumerate_variants ran with require_clean=True.
                weighted = ((0.625 * gc_score + 0.375 * dist_score)
                            * weights.get(mutation, 1.0))
                rows.append((len(seeds) - len(fails), weighted, gc, distance, mutation, cas,
                             site.strand))
            if index % 200 == 0:
                MT.free_gpu_memory()
                print(f"  {index}/{len(jobs)} | {seen} guides | {time.monotonic()-started:.0f}s")
        MT.free_gpu_memory()

    clean = np.array([r[0] for r in rows])
    wtd = np.array([r[1] for r in rows])
    cas_arr = np.array([r[5] for r in rows])
    print(f"\n{len(rows)} guides screened "
          f"({int((cas_arr == 'Cas12a').sum())} Cas12a, {int((cas_arr == 'Cas9').sum())} Cas9)")

    # Unconstrained reference: the best 250 by weighted, ignoring the rule entirely.
    top = np.sort(wtd)[::-1][:ROWS]
    ref = float(top.sum())
    print(f"\nunconstrained best {ROWS} rows: total weighted {ref:.1f}"
          f"  (mean {ref / ROWS:.3f}/row)\n")
    print(f"  {'band k':>7} {'compliant':>10} {'Cas9 of them':>13} {'total weighted':>15} "
          f"{'vs uncon.':>10} {'mean gc':>8} {'mean dist':>10}")
    table = {}
    for k in (1, 2, 3, 4, 5, 6, 8, 10, 12):
        m = clean >= k
        n = int(m.sum())
        if n < ROWS:
            print(f"  {k:>7} {n:>10} {'-':>13} {'too few guides':>15}")
            table[k] = {"compliant": n}
            continue
        idx = np.argsort(wtd[m])[::-1][:ROWS]
        sel_w = wtd[m][idx]
        sel_gc = np.array([r[2] for r in rows])[m][idx]
        sel_d = np.array([r[3] for r in rows])[m][idx]
        n9 = int((cas_arr[m][idx] == "Cas9").sum())
        tot = float(sel_w.sum())
        table[k] = {"compliant": n, "total_weighted": tot, "ratio": tot / ref,
                    "cas9_rows": n9, "mean_gc": float(sel_gc.mean()),
                    "mean_distance": float(sel_d.mean())}
        print(f"  {k:>7} {n:>10} {n9:>13} {tot:>15.1f} {tot / ref:>9.0%} "
              f"{sel_gc.mean():>8.3f} {sel_d.mean():>10.0f}")

    print(f"\n  all-HDR on the same window: band 12-16, weighted ~230-370, spike final ~290.")
    print(f"  A no-cut spike is (total weighted) x 1.000 x fidelity, so the 'total weighted'")
    print(f"  column IS its spike score before fidelity (~0.9).")
    json.dump({"task": task["id"], "cell": cell, "accessibility": acc, "window": [lo, hi],
               "guides": len(rows), "unconstrained_top250": ref, "by_band": table},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
