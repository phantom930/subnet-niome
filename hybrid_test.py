#!/usr/bin/env python3
"""hybrid_test.py — does a no-cut / all-HDR MIXTURE keep consistency, and at what weighted?

The hypothesis under test. A pure pin (every row HDR, or every row no-cut) makes stage 4's three
targets constant, so ``consistency_factor`` is exactly 1.000 — but it also forces all 250 rows
through one rule, which shrinks the guide pool and costs term 1. A *mixture* lets each cas system
do the rule it is relatively best at:

    Cas12a -> no-cut   (P(no cut) 0.22 vs Cas9's 0.14 — base 0.78 against 0.86)
    Cas9   -> HDR      (12.9M candidate guides against Cas12a's 8.7M, at similar P(HDR) ~0.55)

Stage 4 must then PREDICT the outcome rather than see a constant. Whether it can is decided by
``build_X``, which passes exactly:

    ["gc", "distance", "gc_score", "dist_score", "consistency", "energy", "mh"]

``cas_system`` is NOT among them. So a cas-aligned split is only learnable through the fact that gc
is a count over length 20 (Cas9) or 23 (Cas12a) and energy is a function of it — imperfect, which is
a candidate mechanism for the 0.92-0.99 consistencies seen in the field on seed-0 rounds
(19018a0a: 0.9910, 0.9880; e32ec7b3: 0.9890, 0.9220 — all at weighted 325-356). A split aligned to
``dist_score`` instead is directly learnable, so it should score higher still.

Four arms at ONE seed, which is all the consistency question needs:

    pure_hdr     every row HDR                      expect cons 1.000  (the reference)
    pure_nocut   every row no-cut                   expect cons 1.000, fewer Cas9 rows available
    hybrid_cas   Cas12a no-cut, Cas9 HDR            learnable only via gc/energy
    hybrid_dist  split on dist_score, both cas       directly learnable — the interesting arm

    python hybrid_test.py [task_id] [seed]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else None
SEED = int(_ARGV[2]) if len(_ARGV) > 2 else 150
os.environ["NIOME_INSTANCE"] = "hybrid"

import json  # noqa: E402
import logging  # noqa: E402
import math  # noqa: E402
import time  # noqa: E402
from collections import defaultdict  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import mt19937 as MT  # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3  # noqa: E402
from nocut100 import pick_task  # noqa: E402
from sd_task import score, task_content  # noqa: E402

ROWS = 250
CAS12A_ROWS = int(os.getenv("HY_CAS12A", "80"))     # all_hdr's group_size
FLOOR = 4                                           # per (mutation, cas, strand) stage-5 cell
OUT = os.getenv("HY_OUT", "hybrid_test.json")


def candidates(ctx, sites, contract, cell_types, seed, rule, cas):
    """Guides of one cas system that satisfy ``rule`` at ``seed``, with their structural score."""
    acc = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    weights = contract.get("mutation_weights") or {}
    cfg = AC.AllCutConfig()
    seeds = np.array([seed], dtype=np.int64)
    out = []
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == cas for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    for site_index, mutation in jobs:
        site = sites[site_index]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, 0.05, 0.95, ctx.max_mismatches, True,
                                       cfg.variants)
        if not guides:
            continue
        params_of = AC._params_fn(site, distance, acc, offset)
        got = MT.screen_guides_rule_gpu(guides, seeds, mutation, site.cas, site.start,
                                        site.strand, params_of, rule, 0)
        dist_score = math.exp(-distance / ctx.base_padding)
        for guide in got:
            gc = sum(b in "GC" for b in guide) / site.length
            gc_score = max(0.0, 1.0 - abs(gc - 0.5) * 2.0)
            out.append({
                "guide": guide, "mutation": mutation, "cas_system": site.cas,
                "strand": site.strand, "start": site.start, "length": site.length,
                "dist_score": dist_score, "distance": distance,
                "weighted": (0.625 * gc_score + 0.375 * dist_score) * weights.get(mutation, 1.0)})
    return out


def allocate(pool, n_rows, floor=FLOOR):
    """Highest-weighted rows subject to a per-(mutation, cas, strand) floor.

    The floor is not a knob: an empty stage-5 cell puts a zero into a six-way geometric mean that
    clips at 1e-9, which costs ~0.03x on the entire score.
    """
    by_cell = defaultdict(list)
    for rec in pool:
        by_cell[(rec["mutation"], rec["cas_system"], rec["strand"])].append(rec)
    for cell in by_cell:
        by_cell[cell].sort(key=lambda r: -r["weighted"])
    chosen, seen = [], set()

    def take(rec):
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            return False
        seen.add(key)
        chosen.append(rec)
        return True

    for cell, recs in sorted(by_cell.items()):
        got = 0
        for rec in recs:
            if got >= floor:
                break
            got += take(rec)
    rest = sorted((r for r in pool), key=lambda r: -r["weighted"])
    for rec in rest:
        if len(chosen) >= n_rows:
            break
        take(rec)
    return chosen[:n_rows], len(by_cell)


def rows_of(chosen, cell):
    return [{"experiment_id": f"exp-{i:05d}", "guideRNA": r["guide"],
             "target_alignment_start": r["start"],
             "target_alignment_end": r["start"] + r["length"],
             "strand": r["strand"], "mutation": r["mutation"],
             "cas_system": r["cas_system"], "cell_type": cell}
            for i, r in enumerate(chosen)]


def main():
    if TASK:
        task, contract, reference = task_content(TASK)
    else:
        task, contract, reference = pick_task()
    cell = contract.get("cell_type")
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context({**contract, "seed": 0}, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    print(f"task {task['id'][:8]}  {cell}  seed {SEED}  {len(sites)} sites\n")

    t0 = time.monotonic()
    pools = {}
    for rule in ("hdr", "nocut"):
        for cas in ("Cas9", "Cas12a"):
            pools[(rule, cas)] = candidates(ctx, sites, contract, cell_types, SEED, rule, cas)
            print(f"  {rule:<6} {cas:<7} {len(pools[(rule, cas)]):>8} compliant guides "
                  f"({time.monotonic() - t0:.0f}s)")
    print()

    arms = {
        "pure_hdr":   pools[("hdr", "Cas9")] + pools[("hdr", "Cas12a")],
        "pure_nocut": pools[("nocut", "Cas9")] + pools[("nocut", "Cas12a")],
        # The proposal: each cas system on the rule it is relatively better at.
        "hybrid_cas": pools[("hdr", "Cas9")] + pools[("nocut", "Cas12a")],
    }
    # hybrid_dist: the split follows dist_score, which IS a stage-4 feature, so the forest can
    # learn it directly. Near rows take HDR, far rows take no-cut; both cas systems appear in both
    # halves, so the cas-coverage cell stays filled either way.
    every = [(r, "hdr") for r in pools[("hdr", "Cas9")] + pools[("hdr", "Cas12a")]] + \
            [(r, "nocut") for r in pools[("nocut", "Cas9")] + pools[("nocut", "Cas12a")]]
    if every:
        cut = float(np.median([r["dist_score"] for r, _ in every]))
        arms["hybrid_dist"] = [r for r, rule in every
                               if (rule == "hdr") == (r["dist_score"] >= cut)]
        print(f"hybrid_dist splits at dist_score {cut:.4f} "
              f"(>= is HDR, < is no-cut); {len(arms['hybrid_dist'])} eligible\n")

    print(f"  {'arm':<12} {'rows':>5} {'cells':>6} {'cas9':>5} {'weighted':>9} {'cons':>7} "
            f"{'fid':>7} {'final':>8}")
    out = {}
    for name, pool in arms.items():
        if len(pool) < ROWS:
            print(f"  {name:<12} only {len(pool)} compliant guides — cannot fill {ROWS} rows")
            out[name] = {"rows_available": len(pool)}
            continue
        chosen, n_cells = allocate(pool, ROWS)
        rows = rows_of(chosen, cell)
        n9 = sum(1 for r in chosen if r["cas_system"] == "Cas9")
        s = score(rows, contract, reference, cell_types, seed=SEED)
        out[name] = {**{k: s[k] for k in ("weighted", "consistency", "fidelity", "final")},
                     "rows": len(rows), "cells": n_cells, "cas9_rows": n9,
                     "valid": s["valid"]}
        print(f"  {name:<12} {len(rows):>5} {n_cells:>6} {n9:>5} {s['weighted']:>9.1f} "
              f"{s['consistency']:>7.4f} {s['fidelity']:>7.4f} {s['final']:>8.2f}")

    json.dump({"task": task["id"], "cell": cell, "seed": SEED,
               "pools": {f"{r}_{c}": len(v) for (r, c), v in pools.items()}, "arms": out},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
