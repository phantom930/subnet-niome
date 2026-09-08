#!/usr/bin/env python3
"""nocut_band.py — how wide a COMMON no-cut band can N Cas12a rows share?

This is the quantity the whole substitution idea rests on, and it is not what nocut100.py measured.
That counted guides with >= k no-cut seeds among *their own best* seeds; the construction needs 80
guides that all skip the SAME k seeds. Those differ enormously, because compliance is an independent
per-(guide, seed) hash draw:

    own best k      ~2.2M guides reach k=8
    a common k      2.2M * p**k, and with p ~ 0.05 that is already < 1 guide by k=4

The greedy choice of which seeds to use is what decides it — put the band on the seeds the most
guides already skip — so this measures the achievable common band directly rather than inferring it.

Compliance over a fixed candidate seed set is stored as a bitmask per guide, so the greedy runs over
integers instead of re-screening.

    python nocut_band.py <task_id> [n_candidate_seeds]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else "83f430e9-1c9a-48b4-b594-e17477a046ea"
N_CAND = int(_ARGV[2]) if len(_ARGV) > 2 else 16
os.environ["NIOME_INSTANCE"] = "ncband"

import json  # noqa: E402
import logging  # noqa: E402
import time  # noqa: E402
from collections import Counter  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import mt19937 as MT  # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3  # noqa: E402
from sd_task import task_content  # noqa: E402

BAND_LO = 100
WANT = [20, 40, 80, 120, 170]        # candidate Cas12a row counts
OUT = os.getenv("NCB_OUT", "nocut_band.json")


def main():
    task, contract, reference = task_content(TASK)
    cell = contract.get("cell_type")
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context({**contract, "seed": 0}, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    cfg = AC.AllCutConfig()
    cand = np.arange(BAND_LO, BAND_LO + N_CAND, dtype=np.int64)
    print(f"task {TASK[:8]}  {cell} (accessibility {acc:.2f})  "
          f"candidate band seeds {cand[0]}-{cand[-1]}\n")

    # One screen per Cas12a target over the candidate seeds, kept as a bitmask of no-cut seeds.
    masks = Counter()
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == "Cas12a" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    started, seen = time.monotonic(), 0
    for index, (site_index, mutation) in enumerate(jobs, 1):
        site = sites[site_index]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas12a_gc[0], cfg.cas12a_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = AC._params_fn(site, distance, acc, offset)
        seen += len(guides)
        got = MT.screen_guides_rule_gpu(guides, cand, mutation, site.cas, site.start,
                                        site.strand, params_of, "nocut", len(cand))
        for guide, fails in got.items():
            bad = set(int(x) for x in fails) if fails is not None else set(int(s) for s in cand)
            bits = 0
            for j, s in enumerate(cand):
                if int(s) not in bad:               # not a fail => it DID no-cut here
                    bits |= 1 << j
            if bits:
                masks[bits] += 1
        if index % 200 == 0:
            MT.free_gpu_memory()
            print(f"  {index}/{len(jobs)} targets | {seen} guides | "
                  f"{time.monotonic() - started:.0f}s")
    MT.free_gpu_memory()

    total = sum(masks.values())
    per_seed = [sum(c for b, c in masks.items() if b >> j & 1) for j in range(N_CAND)]
    print(f"\n{seen} guides screened; {total} skip at least one candidate seed")
    print(f"  single-seed no-cut rate: mean {np.mean(per_seed)/max(1,seen):.3%} "
          f"(best seed {max(per_seed)/max(1,seen):.3%})\n")

    # Greedy: add the seed that keeps the most guides skipping ALL chosen seeds.
    chosen, alive = [], dict(masks)
    print(f"  {'k':>3} {'seed added':>11} {'guides skipping all k':>22}  "
          + "  ".join(f"{w} rows?" for w in WANT))
    table = {}
    for k in range(1, N_CAND + 1):
        best_j, best_keep = None, None
        for j in range(N_CAND):
            if j in chosen:
                continue
            keep = {b: c for b, c in alive.items() if b >> j & 1}
            n = sum(keep.values())
            if best_keep is None or n > sum(best_keep.values()):
                best_j, best_keep = j, keep
        if best_j is None:
            break
        chosen.append(best_j)
        alive = best_keep
        n = sum(alive.values())
        table[k] = {"seed": int(cand[best_j]), "guides": n}
        print(f"  {k:>3} {int(cand[best_j]):>11} {n:>22}  "
              + "  ".join(f"{'yes' if n >= w else 'no':>7}" for w in WANT))
        if n == 0:
            break

    widest = {w: max([k for k, v in table.items() if v["guides"] >= w], default=0) for w in WANT}
    print(f"\n  widest common no-cut band by row count: "
          + ", ".join(f"{w} rows -> k={widest[w]}" for w in WANT))
    print(f"  for reference, all-HDR reaches band 12-16 with 250 rows on this cell type.")
    json.dump({"task": TASK, "cell": cell, "accessibility": acc, "guides": seen,
               "candidate_seeds": [int(x) for x in cand], "greedy": table,
               "widest_by_rows": widest}, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
