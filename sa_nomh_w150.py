#!/usr/bin/env python3
"""sa_nomh_w150.py — not_mhnhej seed-agnostic hedge at a NARROW 150-seed window, with an explicit
per-guide max_fail cutoff during pool selection (mirroring `cas12a_max_fail`), and no GC/mh
theoretical framing this time -- just the empirical min-union numbers.

`sa_nomh_scratch.py` measured the not_mhnhej seed-agnostic hedge over the FULL 900-seed window and
found it collapses catastrophically past group_size ~20 (7.4% clean) regardless of GC/distance
breadth, because its ~20% per-guide failure rate is close to independent across guides -- the
min-union greedy has nothing to exploit there. Two open questions that measurement didn't answer:

1. Does narrowing the WINDOW itself help, independent of any GC/mh mechanism? A smaller seed
   universe means the SAME candidate pool is denser per seed, so pure chance coincidence (two
   guides failing on the same few of 150 seeds, rather than scattered across 900) becomes much more
   likely -- the same mechanism `conjunction.py`'s own width sweeps exploit for its band. This is a
   WINDOW effect, not a design-lever effect, so it can matter even though the GC lever from the last
   script did not.
2. Does an explicit max_fail cutoff during pool selection -- discarding any candidate whose OWN
   not_mhnhej failure count already exceeds a threshold, exactly like the shipped hedge's
   `cas12a_max_fail` -- change anything, given not_mhnhej's population failure rate (~20% of seeds,
   so ~30 of 150 expected per guide) sits far higher than the cut-only hedge's (~2-6%, where
   `cas12a_max_fail=22` keeps nearly everything)?

Enumeration keeps the wide gc_band (0.0-1.0) and no max_distance from the last script, since that
was measured cheap and harmless -- only the window width and the GC/mh theoretical bookkeeping
change here.

    SA_CELL=K562 python sa_nomh_w150.py [task_id]

RESULT (2026-09-21, K562, task 2af71117, window [100,249], 4000-candidate pool): narrowing the
window helps meaningfully; the max_fail cutoff mostly does not, for different reasons -- and
neither closes the gap to plain cut-only.

Per-guide failure rate is unchanged at width 150 (mean 19.8% of 150, matching the full-window rate
-- a per-seed probability, not a window artefact). But clean fraction at group_size 80 rises from
1.6% (900-wide, sa_nomh_scratch.py) to 11-12% (150-wide) -- a ~7x relative improvement, because the
smaller seed universe makes chance coincidence of failures across guides far more likely (denser
pool per seed). It is still nowhere near viable: cut-only on the IDENTICAL pool at the IDENTICAL
width reaches 83-91% clean at group_size 40-150 against not_mhnhej's 9-15%, a ~7x gap that survives
the narrower window too.

max_fail is nearly inert once above the achievable floor: max_fail=10 has ZERO survivors (the best
guide in this pool only reaches 11/150 failures), but 30 through no-filter-at-all give IDENTICAL
min-union results -- the greedy already gravitates to low-failure guides unaided. Set too tight
(max_fail=20, 209 of 4000 survive), it is actively WORSE at larger group sizes (5% clean at
group_size=40 against 15% at max_fail>=30): restricting the pool that hard denies the greedy guides
whose failures happen to OVERLAP an already-chosen guide's, which is what min-union optimises for,
not individual guide quality. See CLAUDE.md's falsified table for the full writeup.
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "sanomhw150")
CELL = os.getenv("SA_CELL", "K562")
os.environ["EH_CELL"] = CELL

import itertools                                          # noqa: E402
import logging                                            # noqa: E402
import time                                                # noqa: E402
from collections import Counter                            # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                          # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA       # noqa: E402
from niome_subnet.genomics.validation import stage3         # noqa: E402
from sd_task import task_content                            # noqa: E402
from energy_hdr_k562 import pick_task                        # noqa: E402

START_SEED, WIDTH = 100, 150
SEEDS = list(range(START_SEED, START_SEED + WIDTH))
FLANK = 3000
LENGTHS = (20, 23)
CAP_PER_SITE = 30
CAP_TOTAL = 4000
MAX_FAIL_SWEEP = (10, 20, 30, 40, 50, 60, 75, 100, WIDTH)   # WIDTH == no filter at all
GROUP_SIZES = (10, 20, 40, 80, 150)


def scan(ctx):
    sites = G.enumerate_sites(ctx, FLANK, LENGTHS)
    jobs = [(s, m) for s in sites for m in ctx.mutations]
    by_cell: dict[tuple, list] = {}
    for job in jobs:
        site, mutation = job
        by_cell.setdefault((mutation, site.cas, site.strand), []).append(job)
    for lst in by_cell.values():
        lst.sort(key=lambda j: abs(j[0].start - ctx.mutation_map[j[1]]))
    interleaved = [job for group in itertools.zip_longest(*by_cell.values())
                   for job in group if job is not None]

    out = []
    scanned_sites = 0
    t0 = time.monotonic()
    for site, mutation in interleaved:
        if len(out) >= CAP_TOTAL:
            break
        scanned_sites += 1
        guides = SA.enumerate_variants(site, ctx, 0.0, 1.0, ctx.max_mismatches, True, CAP_PER_SITE)
        for guide in guides:
            entry = G.build_valid_entry(G.make_experiment(site, guide, mutation, ctx, "cand"), ctx)
            if entry is None:
                continue
            cut_fails, nomh_fails = [], []
            for s in SEEDS:
                rec = stage3.simulate(entry, s)
                if rec["outcome"] == "no_cut":
                    cut_fails.append(s)
                    nomh_fails.append(s)
                elif rec["outcome"] == "MH_NHEJ":
                    nomh_fails.append(s)
            out.append({
                "guide": guide, "mutation": mutation, "cas_system": site.cas,
                "strand": site.strand, "target_alignment_start": site.start, "length": site.length,
                "cut_fails": cut_fails, "fails": nomh_fails, "n_fail": len(nomh_fails),
            })
            if len(out) >= CAP_TOTAL:
                break
    print(f"scanned {scanned_sites} (site,mutation) jobs, banked {len(out)} candidates "
          f"in {time.monotonic() - t0:.1f}s  window [{SEEDS[0]},{SEEDS[-1]}] ({len(SEEDS)} seeds)",
          flush=True)
    return out


def sweep(pool, cfg, label):
    print(f"=== max_fail cutoff x group_size: {label} ===")
    for max_fail in MAX_FAIL_SWEEP:
        survivors = [c for c in pool if c["n_fail"] <= max_fail]
        row = f"  max_fail {max_fail:4d}: {len(survivors):4d}/{len(pool)} survive  "
        if not survivors:
            print(row)
            continue
        for group_size in GROUP_SIZES:
            if len(survivors) < group_size:
                row += f"| g{group_size}: too small  "
                continue
            group = SA.min_union_group(survivors, group_size, cfg, WIDTH)
            union = len({s for c in group for s in c["fails"]})
            row += f"| g{group_size}: {union:3d}/{WIDTH} ({100 * (WIDTH - union) / WIDTH:.0f}% clean)  "
        print(row, flush=True)
    print(flush=True)


def main():
    task_id = _ARGV[1] if len(_ARGV) > 1 else pick_task()
    task, contract, reference = task_content(task_id)
    print(f"task {task_id[:8]}  {CELL}  width-{WIDTH} seed-agnostic window "
          f"[{SEEDS[0]},{SEEDS[-1]}]\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)

    pool = scan(ctx)
    if not pool:
        raise SystemExit("no candidates found")

    n_fails = np.array([c["n_fail"] for c in pool])
    print(f"pool size {len(pool)}  cas mix {dict(Counter(c['cas_system'] for c in pool))}")
    print(f"not_mhnhej failures per guide, of {WIDTH}: mean {n_fails.mean():.1f} "
          f"({100 * n_fails.mean() / WIDTH:.1f}%)  min {n_fails.min()}  max {n_fails.max()}\n",
          flush=True)

    cfg = SA.SeedAgnosticConfig(start_seed=SEEDS[0], end_seed=SEEDS[-1], group_restarts=8)

    sweep(pool, cfg, "not_mhnhej")
    cut_pool = [dict(c, fails=c["cut_fails"], n_fail=len(c["cut_fails"])) for c in pool]
    sweep(cut_pool, cfg, "plain cut-only, same pool, for reference")


if __name__ == "__main__":
    main()
