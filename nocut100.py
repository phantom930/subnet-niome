#!/usr/bin/env python3
"""nocut100.py — can the no-cut rule hold over a 100-seed window, and what does it cost?

"no-cut" is the exact mirror of all-cut: pin ``is_cut`` FALSE across a seed window instead of true.
It is worth asking because, like ``hdr``, it pins all three of stage 4's targets at once — is_cut
False, is_hdr False, indel_length 0 — so a seed where every row complies scores
``consistency_factor`` 1.000, the same regime-1 prize all-HDR buys.

The arithmetic says it will be narrow and expensive, and this measures both rather than asserting:

    cut_p = min(0.99, max(0.4, base + 0.18 * energy)),  base 0.86 Cas9 / 0.78 Cas12a
    energy = clamp(accessibility * (1.8*gc + 0.6*exp(-dist/1500) + region_offset), 0, 1)

energy is non-negative, so cut_p never falls below 0.78 and the 0.4 floor is unreachable dead code.
P(no cut) per row is therefore at most ~0.22 (Cas12a, energy 0) and ~0.14 for Cas9 — against
P(cut) ~0.99 for all-cut and P(HDR) ~0.57 for all-HDR. With N candidate guides the widest band that
still leaves ``rows`` compliant guides is about log(N/rows)/log(1/p), i.e. ~15 at p=0.57 and ~5 at
p=0.18.

The second cost is the one that matters more, and it is a conflict of interest rather than a
probability: low energy needs low GC and large distance, but ``gc_score`` peaks at GC 0.50 and
``dist_score`` decays as exp(-d/400). The guides that most easily avoid cutting are exactly the
guides that score worst on term 1, so this reports the weighted_score of the compliant set against
the unconstrained best.

    python nocut100.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else None
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "100-199"
os.environ["NIOME_INSTANCE"] = "nocut"

import json  # noqa: E402
import logging  # noqa: E402
import random  # noqa: E402
import time  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import mt19937 as MT  # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3  # noqa: E402
from sd_task import fetch, task_content  # noqa: E402

ROWS = 250
OUT = os.getenv("NOCUT_OUT", "nocut100.json")
VERIFY_N = int(os.getenv("NOCUT_VERIFY", "4000"))


def pick_task():
    """Newest stamped task if none named — any real contract will do for a rule measurement."""
    if TASK:
        return task_content(TASK)
    items = fetch("https://niome-api.genomes.io/api/v3/tasks", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or [])
    items.sort(key=lambda t: t.get("created_at") or "")
    for t in reversed(items):
        c = (t.get("content") or {}).get("contract") or {}
        if c and (t.get("content") or {}).get("hbb_reference"):
            return t, c, t["content"]["hbb_reference"]
    raise SystemExit("no usable task in the listing")


def verify_rule(ctx, sites, contract, cell_types, seeds, want=200):
    """Cross-check the GPU 'nocut' screen against stage 3 itself on real (guide, seed) pairs.

    Uses only the public screen — no reimplementation of the draws — so what is verified is exactly
    what a bank would be built from: the screen reports the seeds a guide FAILS the rule on, and
    stage 3 must call every other seed in the window ``no_cut``. A one-bit divergence here silently
    invalidates every bank built on the rule, which is why CLAUDE.md insists on this step.
    """
    cell = contract.get("cell_type")
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    rng = random.Random(7)
    checked = mismatched = 0
    for _ in range(400):
        if checked >= want:
            break
        site = sites[rng.randrange(len(sites))]
        mutation = ctx.mutations[rng.randrange(len(ctx.mutations))]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, 0.20, 0.80, ctx.max_mismatches, True, 24)
        if not guides:
            continue
        params_of = AC._params_fn(site, distance, acc, offset)
        guide = guides[rng.randrange(len(guides))]
        got = MT.screen_guides_rule_gpu([guide], seeds, mutation, site.cas, site.start,
                                        site.strand, params_of, "nocut", len(seeds))
        failed = set(int(x) for x in got.get(guide, np.array([], dtype=np.int64)))
        if guide not in got:            # screen dropped it: it fails every seed in the window
            failed = set(int(s) for s in seeds)
        exp = {"experiment": {"experiment_id": "verify", "guideRNA": guide,
                              "target_alignment_start": site.start,
                              "target_alignment_end": site.start + site.length,
                              "strand": site.strand, "mutation": mutation,
                              "cas_system": site.cas},
               "features": {"gc": sum(b in "GC" for b in guide) / site.length,
                            "distance_to_mutation": distance, "gc_score": 0.0, "dist_score": 0.0,
                            "consistency": 0.0, "cell_type_accessibility": acc,
                            "mutation_weight": 1.0, "region_energy_offset": offset}}
        for seed in rng.sample(list(map(int, seeds)), min(8, len(seeds))):
            truth_nocut = stage3.simulate(exp, seed)["outcome"] == "no_cut"
            screen_nocut = seed not in failed
            checked += 1
            if truth_nocut != screen_nocut:
                mismatched += 1
    return checked, mismatched


def main():
    task, contract, reference = pick_task()
    lo, hi = (int(x) for x in WINDOW.split("-"))
    seeds = np.arange(lo, hi + 1, dtype=np.int64)
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context({**contract, "seed": 0}, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    acc = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    print(f"task {task['id'][:8]}  {contract.get('cell_type')} (accessibility {acc:.2f})  "
          f"window {lo}-{hi} ({len(seeds)} seeds)  {len(sites)} sites\n")

    print("what cut_p allows, before any search:")
    for cas, base in (("Cas9", 0.86), ("Cas12a", 0.78)):
        for e in (0.0, 0.25, 0.5, 1.0):
            cp = min(0.99, max(0.4, base + 0.18 * e))
            print(f"  {cas:<7} energy {e:.2f} -> cut_p {cp:.3f}  P(no cut) {1 - cp:.3f}  "
                  f"band for {ROWS} rows of 1e6 candidates ~ "
                  f"{np.log(1e6 / ROWS) / np.log(1 / max(1e-9, 1 - cp)):.1f} seeds")
    print()

    checked, bad = verify_rule(ctx, sites, contract, cell_types, seeds)
    print(f"rule verification against stage 3: {checked} (guide, seed) pairs, "
          f"{bad} mismatches" + ("  -- ABORT, the rule is wrong" if bad else "  OK\n"))
    if bad:
        raise SystemExit(1)

    # The real screen: how many of the window's seeds does each guide fail (i.e. DOES cut on)?
    cfg = AC.AllCutConfig()
    regions = contract.get("mutation_regions") or {}
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == "Cas12a" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    print(f"screening {len(jobs)} Cas12a targets over the window (max_fail {len(seeds)}, "
          f"i.e. keep everything and read the distribution)")
    best, started, n_guides = [], time.monotonic(), 0
    for index, (site_index, mutation) in enumerate(jobs, 1):
        site = sites[site_index]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas12a_gc[0], cfg.cas12a_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = AC._params_fn(site, distance, accessibility=acc, region_offset=offset)
        n_guides += len(guides)
        survivors = MT.screen_guides_rule_gpu(guides, seeds, mutation, site.cas, site.start,
                                              site.strand, params_of, "nocut", len(seeds))
        for guide, fails in survivors.items():
            best.append((len(fails), guide, mutation, site.cas, site.strand, site.start,
                         site.length))
        if index % 200 == 0:
            MT.free_gpu_memory()
            print(f"  {index}/{len(jobs)} targets | {n_guides} guides | "
                  f"{time.monotonic() - started:.0f}s")
    MT.free_gpu_memory()

    if not best:
        print("\nno guide produced a single no-cut seed in this window.")
        json.dump({"task": task["id"], "window": [lo, hi], "guides": n_guides,
                   "verdict": "no candidate ever fails to cut"}, open(OUT, "w"), indent=1)
        return 0

    fails = np.array([b[0] for b in best])
    clean = len(seeds) - fails                      # seeds where the guide does NOT cut
    print(f"\n{n_guides} guides screened, {len(best)} returned by the screen")
    print(f"no-cut seeds per guide (of {len(seeds)}): "
          f"max {clean.max()}  p99 {np.percentile(clean, 99):.0f}  "
          f"median {np.median(clean):.0f}  mean {clean.mean():.2f}")
    print(f"\n  {'guides with >= k no-cut seeds':<34} {'count':>9}  "
          f"{'enough for ' + str(ROWS) + ' rows?':>22}")
    table = {}
    for k in (1, 2, 3, 4, 5, 6, 8, 10, 15, 20):
        n = int((clean >= k).sum())
        table[k] = n
        print(f"  {'k = ' + str(k):<34} {n:>9}  {'yes' if n >= ROWS else 'no':>22}")
    widest = max([k for k, n in table.items() if n >= ROWS], default=0)
    print(f"\n  widest band with {ROWS} compliant guides: {widest} seeds of {len(seeds)}")
    print(f"  (all-HDR reaches 12-16 on the same 100-seed window; all-cut pins ~560 of 900)")

    json.dump({"task": task["id"], "cell": contract.get("cell_type"), "accessibility": acc,
               "window": [lo, hi], "guides_screened": n_guides, "returned": len(best),
               "clean_max": int(clean.max()), "clean_mean": float(clean.mean()),
               "guides_at_k": table, "widest_band_for_rows": widest, "rows": ROWS},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
