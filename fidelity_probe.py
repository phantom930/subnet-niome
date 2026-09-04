#!/usr/bin/env python3
"""fidelity_probe.py — can a high-fidelity submission reach consistency 0.45-0.65 with no band?

Across 16,102 current-regime score rows, fidelity gates the elevated consistency regime absolutely:
11.6% of rows above fidelity 0.95 land in consistency 0.45-0.65, against 0 of 6,079 rows at or
below 0.90. It does *not* raise the floor (median 0.103 vs 0.081), and it is a within-miner mode
switch, not a property of who the miner is (43/43 miners with >=4 such rounds show +0.085 fidelity
on exactly those rounds, at -18.7 weighted). Our fleet's median fidelity is 0.899 -- the bucket
where the rate is zero.

This tests whether the regime is *reachable by construction* rather than only correlated. Three
submissions on one contract, differing only in how rows are chosen, each scored on many random
seeds with no band and no seed screening whatsoever:

  A  all-HDR, the shipped build            -- fidelity ~0.89, the incumbent
  B  balanced cells + k-mer greedy         -- maximises stage 5 directly
  C  B, plus gc/distance stratification    -- also widens the stage-4 feature spread

A and B separate "fidelity gates the regime" from "fidelity is a marker of someone else's
construction". B and C separate stage-5 diversity from stage-4 learnability: consistency is
``0.7*max(avg_r2, 0) + 0.3*(1 - avg_nmae)``, so a positive cross-validated R2 needs outcomes that
genuinely vary *with the features*, which k-mer diversity alone does not deliver. Our own strongest
counter-example is h0 on 9ed335da: fidelity 0.973, above the band population's 0.958, and
consistency still 0.397 -- exactly k=1 at our floor.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import dataclasses
import json
import random
import statistics as st
import sys
import time
from collections import Counter, defaultdict

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics import seed_depend as SD
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.genomics.validation import stage5
from niome_subnet.utils import settings

SD.RULES["any"] = lambda rec: True          # no outcome filter: this probe pins nothing
TASKS = ["task-k562.json", "task-cd34.json"]
N_SEEDS = 16
PROBE_SEED = 20260904
OUT = "fidelity_probe.json"


# Bounded per (cell, gc quartile, distance quartile). SD.enumerate_candidates keeps the full
# ``entry`` and ``record`` dict for every candidate, which reached 39 GB and was OOM-killed at
# ~1M candidates -- and this probe needs two such passes. Only four scalars and the 12-mer set
# survive here, and each bucket is capped, so the peak is tens of MB rather than tens of GB.
PER_BUCKET = 300


def candidates(contract, reference, cell_types, gc_band, max_distance, variants):
    pinned = dict(contract)
    pinned["seed"] = 500                     # nothing is screened on it; the rule is "any"
    ctx = G.build_context(pinned, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    jobs = [(s, m) for s in sites for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= max_distance]
    jobs.sort(key=lambda j: abs(j[0].start - ctx.mutation_map[j[1]]))
    buckets = defaultdict(list)
    seen = 0
    for site, mutation in jobs:
        guides = SA.enumerate_variants(site, ctx, gc_band[0], gc_band[1],
                                       ctx.max_mismatches, True, variants)
        for guide in guides:
            entry = G.build_valid_entry(G.make_experiment(site, guide, mutation, ctx, "c"), ctx)
            if entry is None:
                continue
            f = entry["features"]
            gq = min(3, int((f["gc"] - gc_band[0]) / max(1e-9, gc_band[1] - gc_band[0]) * 4))
            dq = min(3, int(f["distance_to_mutation"] / max(1, max_distance) * 4))
            key = ((mutation, site.cas, site.strand), gq, dq)
            if len(buckets[key]) >= PER_BUCKET:
                continue
            seen += 1
            buckets[key].append({
                "guide": guide, "mutation": mutation, "cas_system": site.cas,
                "strand": site.strand, "start": site.start, "length": site.length,
                "weighted_score": entry["stage2"]["weighted_score"],
                "gc": f["gc"], "distance": f["distance_to_mutation"],
                "kmers": frozenset(stage5.extract_kmers(guide, 12))})
    by_cell = defaultdict(list)
    for (cell, _gq, _dq), recs in buckets.items():
        by_cell[cell] += recs
    for recs in by_cell.values():
        recs.sort(key=lambda r: -r["weighted_score"])
    return ctx, by_cell


def balanced(by_cell, n_rows, stratify, kmer_price=1.0, rng=None):
    """n_rows spread evenly over the 8 stage-5 cells, greedy on new 12-mers.

    ``stratify`` additionally forces the picks across gc x distance quartiles inside each cell,
    which is what widens the stage-4 feature spread without touching stage 5.
    """
    cells = sorted(by_cell)
    base, extra = divmod(n_rows, len(cells))
    used, chosen = set(), []
    for i, key in enumerate(cells):
        take = base + (1 if i < extra else 0)
        recs = by_cell[key]
        if stratify:
            gs = sorted({r["gc"] for r in recs})
            ds = sorted({r["distance"] for r in recs})
            def bucket(r):
                return (min(3, int(4 * gs.index(r["gc"]) / max(1, len(gs)))),
                        min(3, int(4 * ds.index(r["distance"]) / max(1, len(ds)))))
            groups = defaultdict(list)
            for r in recs:
                groups[bucket(r)].append(r)
            order = sorted(groups, key=lambda k: -len(groups[k]))
            pool = []
            while len(pool) < take * 12 and any(groups[k] for k in order):
                for k in order:
                    if groups[k]:
                        pool.append(groups[k].pop(0))
            recs = pool
        window = recs[:max(take * 10, 400)]
        taken = [False] * len(window)
        for _ in range(take):
            bestv, besti = None, None
            for j, r in enumerate(window):
                if taken[j]:
                    continue
                v = kmer_price * len(r["kmers"] - used) + r["weighted_score"]
                if bestv is None or v > bestv:
                    bestv, besti = v, j
            if besti is None:
                break
            taken[besti] = True
            chosen.append(window[besti])
            used |= window[besti]["kmers"]
    return chosen


def rows_of(chosen, contract):
    rows, seen = [], set()
    for i, r in enumerate(chosen):
        k = (r["cas_system"], r["start"], r["strand"], r["guide"])
        if k in seen:
            continue
        seen.add(k)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": r["guide"],
                     "target_alignment_start": r["start"],
                     "target_alignment_end": r["start"] + r["length"],
                     "strand": r["strand"], "mutation": r["mutation"],
                     "cas_system": r["cas_system"], "cell_type": contract.get("cell_type")})
    return rows


def score(rows, contract, reference, cell_types, seeds):
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    out = []
    for s in seeds:
        run_stage3(seed=s)
        run_stage4(seed=s)
        out.append(run_stage5())
    return out


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    rng = random.Random(PROBE_SEED)
    seeds = sorted(rng.sample(range(100, 1000), N_SEEDS))
    print(f"probe seeds ({N_SEEDS} random, no band anywhere): {seeds}\n")
    results = []
    for tf in TASKS:
        task = json.load(open(tf))
        contract = dict(task["content"]["contract"])
        reference = task["content"]["hbb_reference"]
        cell = contract["cell_type"]
        contract["seed"] = ",".join(str(s) for s in seeds)
        n_rows = contract["rules"].get("max_experiments") or 250
        print(f"=== {tf}  {cell}  {task['id'][:8]} ===", flush=True)

        builds = {}
        t0 = time.monotonic()
        hdr_rows, meta = AH.build_for_cell(dict(task["content"]["contract"]), reference,
                                           cell_types, budget_s=900.0)
        if hdr_rows:
            builds["A all-HDR (shipped)"] = hdr_rows
        print(f"  A all-HDR: {'built' if hdr_rows else meta.get('reason')} "
              f"({time.monotonic() - t0:.0f}s)", flush=True)

        t0 = time.monotonic()
        _ctx, by = candidates(contract, reference, cell_types, (0.30, 0.70), 600, 1200)
        print(f"  candidates (narrow): {sum(len(v) for v in by.values())} over {len(by)} cells "
              f"({time.monotonic() - t0:.0f}s)", flush=True)
        if len(by) >= 8:
            builds["B balanced + kmer"] = rows_of(balanced(by, n_rows, False), contract)
        by = None                            # released before the second pass is enumerated

        t0 = time.monotonic()
        _ctx, by2 = candidates(contract, reference, cell_types, (0.20, 0.80), 1200, 1200)
        print(f"  candidates (wide):   {sum(len(v) for v in by2.values())} over {len(by2)} cells "
              f"({time.monotonic() - t0:.0f}s)", flush=True)
        if len(by2) >= 8:
            builds["C balanced + kmer + spread"] = rows_of(balanced(by2, n_rows, True), contract)
        by2 = None

        for label, rows in builds.items():
            if len(rows) < n_rows:
                print(f"  {label}: only {len(rows)} rows; skipped")
                continue
            per = score(rows, contract, reference, cell_types, seeds)
            cons = [p["consistency_factor"] for p in per]
            fid = per[0]["distribution_fidelity_factor"]
            wtd = per[0]["total_weighted_score"]
            inband = sum(1 for c in cons if 0.45 <= c < 0.65) / len(cons)
            rec = {"task": task["id"], "cell_type": cell, "variant": label, "rows": len(rows),
                   "fidelity": fid, "weighted": wtd, "seeds": seeds, "consistency": cons,
                   "cons_median": st.median(cons), "cons_max": max(cons),
                   "share_0.45_0.65": inband}
            results.append(rec)
            json.dump(results, open(OUT, "w"), indent=1)
            print(f"  {label:<28} fid {fid:.3f} wtd {wtd:>6.1f} | cons median "
                  f"{st.median(cons):.3f} max {max(cons):.3f} | in 0.45-0.65 {inband*100:>5.1f}%",
                  flush=True)
            print(f"      {' '.join(f'{c:.3f}' for c in cons)}", flush=True)
        print()
    print("done")


if __name__ == "__main__":
    main()
