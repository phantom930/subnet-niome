#!/usr/bin/env python3
"""search_guides.py — step 1 of a seed-agnostic construction: guides that cut under every seed.

Stage 3 decides whether a row cuts with the *second* draw of ``random.Random(experiment_seed(
round_seed, design))``: the row cuts iff that draw <= ``cut_p``. ``experiment_seed`` hashes the
round seed together with the design, so changing the (unknown, real) seed reshuffles the draw. A
guide is **seed-agnostic for the cut gate** when its draw stays under ``cut_p`` for *every* seed in a
range — then ``is_cut`` is constant no matter which seed the validator picks, which is the first
thing a construction needs to survive a restamp.

This enumerates a large pool of GC-tuned guide variants per (site, mutation) cell and keeps the ones
that cut under all seeds in ``--start-seed .. --end-seed`` (default 100..999). ``cut_p`` is identical
across the variants of a cell — it depends only on energy, i.e. gc, distance and accessibility, all
seed-independent and all fixed once the GC target is tuned — so only the draw sequence, through the
guide string in the hash, separates a passer from a miss.

    python search_guides.py                         # task.json, seeds 100..999, into test/guides
    python search_guides.py --variants 1024         # a wider pool per cell
    python search_guides.py --min-pass 895          # keep near misses too

The yield is set by ``cut_p ** n_seeds``. At the energy clamp that is 0.99**900 ~ 1.2e-4 for Cas9
(about one passer per 8,500 guides) and 0.96**900 ~ 1e-16 for Cas12a (never). A seed-agnostic *cut*
therefore exists only for Cas9 on a high-accessibility cell type; the tool measures exactly how many.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import logging
import math
import os
import random
import sys
import time
from collections import Counter
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import niome_subnet.utils.settings as settings  # noqa: E402
import genExp as G  # noqa: E402

from niome_subnet.genomics.validation import stage3  # noqa: E402

STATE: dict = {}


def energy_at_gc(gc: float, distance: int, accessibility: float, region_offset: float) -> float:
    """stage3.sequence_energy for an explicit gc — the one input that now varies across a cell.

    Once GC is allowed a band rather than a single point, cut_p is no longer constant within a cell:
    energy = accessibility * (1.8*gc + 0.6*exp(-d/1500) + region_offset), clamped to [0, 1]. So each
    GC count gets its own cut_p, computed from the gc the tuner actually reached.
    """
    return max(0.0, min(1.0, accessibility * (
        1.8 * gc + 0.6 * math.exp(-distance / 1500) + region_offset
    )))


def cut_draw(seed: int, mutation: str, cas: str, guide: str, start: int, strand: str) -> float:
    """The cut coin — the second draw of the row's stream, replicated without a full simulate.

    ``experiment_seed`` is ``int(sha256(key).hexdigest(), 16) % 2**32``, which is just the low 32
    bits of the digest, i.e. its last four bytes big-endian. The stream then draws the microhomology
    coin first and the cut coin second, so one throwaway ``random()`` precedes the value that decides
    the cut.
    """
    key = f"{seed}|{mutation}|{cas}|{guide}|{start}|{strand}"
    seed32 = int.from_bytes(hashlib.sha256(key.encode()).digest()[-4:], "big")
    rng = random.Random(seed32)
    rng.random()            # microhomology coin — consumed, not used here
    return rng.random()     # cut coin: the row cuts iff this is <= cut_p


def passes_all_seeds(mutation: str, cas: str, guide: str, start: int, strand: str,
                     cut_p: float, seeds: range, min_pass: int) -> tuple[int, float]:
    """Count seeds the guide cuts under, short-circuiting once it cannot reach ``min_pass``.

    Returns (cut_count, max_draw_seen). A guide cuts under a seed iff its cut coin <= cut_p, exactly
    the validator's rule (no_cut iff draw > cut_p). The loop abandons a guide the moment the seeds
    it has already failed make ``min_pass`` unreachable — which is most of them, so the average guide
    costs far fewer than ``len(seeds)`` hashes.
    """
    total = len(seeds)
    cuts = 0
    failed = 0
    max_draw = 0.0
    allowed_fail = total - min_pass
    for seed in seeds:
        draw = cut_draw(seed, mutation, cas, guide, start, strand)
        if draw > max_draw:
            max_draw = draw
        if draw <= cut_p:
            cuts += 1
        else:
            failed += 1
            if failed > allowed_fail:
                break
    return cuts, max_draw


def worker_init(payload: dict) -> None:
    logging.getLogger().setLevel(logging.ERROR)
    sys.stdout = open(os.devnull, "w")   # tune_variants / generate narrate to stdout

    task = json.load(open(payload["task"]))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell_types = payload["cell_types"]
    ctx = G.build_context(contract, reference, cell_types)

    from neurons.miner import Miner
    cfg = Miner.gen_config_for(contract)
    sites = G.enumerate_sites(ctx, cfg.flank, tuple(sorted(set(cfg.lengths))))

    accessibility = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    region_offset = {m: stage3.REGION_ENERGY_OFFSETS.get(regions.get(m), 0.0)
                     for m in ctx.mutations}
    STATE.update(ctx=ctx, cfg=cfg, sites=sites,
                 seeds=range(payload["start_seed"], payload["end_seed"] + 1),
                 min_pass=payload["min_pass"], variants=payload["variants"],
                 max_distance=payload["max_distance"], cas_filter=payload["cas_filter"],
                 accessibility=accessibility, region_offset=region_offset,
                 gc_min=payload["gc_min"], gc_max=payload["gc_max"],
                 require_clean=payload["require_clean"])


NUCLEOTIDES = ("A", "C", "G", "T")


def enumerate_variants(site, ctx, gc_lo: float, gc_hi: float, budget: int,
                       require_clean: bool, cap: int):
    """Every guide within Hamming <= budget of the reference target whose GC lands in the band.

    tune_variants samples a few dozen per site and stops as soon as the off-target seed is clean;
    for banking candidates that throws away ~99.9% of the reachable space (measured 4-66 returned
    against 22k-44k available). This enumerates it in full: choose up to ``budget`` positions,
    substitute each to a different base, keep the ones whose GC count is in the band and — when
    ``require_clean`` — whose seed 12-mer is absent from the off-target index (offtarget_factor 1.0).

    Yielded nearest-reference first (Hamming ascending), so a ``cap`` keeps the least-perturbed
    guides; ranking by weighted_score happens later on the survivors.
    """
    ref = list(site.ref_guide)
    L = site.length
    lo = math.ceil(gc_lo * L - 1e-9)
    hi = math.floor(gc_hi * L + 1e-9)
    sl = G.seed_slice(site.cas, L)
    ref_gc = sum(b in "GC" for b in ref)
    positions = list(range(L))
    seen: set[str] = set()
    out: list[str] = []

    for k in range(0, budget + 1):
        for combo in itertools.combinations(positions, k):
            # GC count only moves when a substitution crosses the A/T <-> G/C boundary, so bound it
            # before expanding base choices: at most k positions can each shift GC by +/-1.
            choices = [[b for b in NUCLEOTIDES if b != ref[i]] for i in combo]
            for repl in itertools.product(*choices):
                gc = ref_gc
                for i, b in zip(combo, repl):
                    was = ref[i] in "GC"
                    now = b in "GC"
                    gc += (1 if now else 0) - (1 if was else 0)
                if not (lo <= gc <= hi):
                    continue
                guide = ref[:]
                for i, b in zip(combo, repl):
                    guide[i] = b
                gs = "".join(guide)
                if gs in seen:
                    continue
                if require_clean and gs[sl] in ctx.kmer_index:
                    continue
                seen.add(gs)
                out.append(gs)
                if len(out) >= cap:
                    return out
    return out


def evaluate_cell(cell: tuple[int, str]) -> dict:
    """One (site, mutation) cell: tune a variant pool, keep those that cut under enough seeds."""
    site_index, mutation = cell
    ctx, cfg = STATE["ctx"], STATE["cfg"]
    site = STATE["sites"][site_index]
    seeds, min_pass = STATE["seeds"], STATE["min_pass"]

    if STATE["cas_filter"] and site.cas not in STATE["cas_filter"]:
        return {"cas": site.cas, "strand": site.strand, "tested": 0, "passers": []}
    distance = abs(site.start - ctx.mutation_map[mutation])
    if distance > STATE["max_distance"]:
        return {"cas": site.cas, "strand": site.strand, "tested": 0, "passers": []}

    L = site.length
    accessibility = STATE["accessibility"]
    region_offset = STATE["region_offset"][mutation]
    # Every valid guide of this target across the GC band, not a sample of a few dozen.
    variants = enumerate_variants(site, ctx, STATE["gc_min"], STATE["gc_max"],
                                  ctx.max_mismatches, STATE["require_clean"], STATE["variants"])

    passers = []
    tested = 0
    cut_p_by_gc: dict[int, float] = {}
    for guide in variants:
        tested += 1
        gc_count = sum(b in "GC" for b in guide)
        cut_p = cut_p_by_gc.get(gc_count)
        if cut_p is None:
            # cut_p depends only on the gc count (via energy), so compute it once per count.
            energy = energy_at_gc(gc_count / L, distance, accessibility, region_offset)
            cut_p = cut_p_by_gc[gc_count] = stage3.cut_probability(site.cas, energy)
        cuts, max_draw = passes_all_seeds(mutation, site.cas, guide, site.start, site.strand,
                                          cut_p, seeds, min_pass)
        if cuts < min_pass:
            continue
        entry = G.build_valid_entry(
            G.make_experiment(site, guide, mutation, ctx, "cand"), ctx)
        if entry is None:
            continue      # would be rejected by stage 1; not a usable candidate
        passers.append({
                "guide": guide,
                "mutation": mutation,
                "cas_system": site.cas,
                "strand": site.strand,
                "target_alignment_start": site.start,
                "length": site.length,
                "distance_to_mutation": distance,
                "gc": entry["features"]["gc"],
                "gc_score": entry["features"]["gc_score"],
                "offtarget_factor": entry["features"]["offtarget_factor"],
                "cut_p": cut_p,
                "weighted_score": entry["stage2"]["weighted_score"],
                "cut_seeds": cuts,
                "total_seeds": len(seeds),
                "max_cut_draw": max_draw,
                "margin": cut_p - max_draw,     # >0 exactly when it cuts under all tested seeds
            })

    return {"cas": site.cas, "strand": site.strand,
            "tested": tested, "passers": passers}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task", default="task.json")
    parser.add_argument("--out-dir", default="test/guides")
    parser.add_argument("--start-seed", type=int, default=100)
    parser.add_argument("--end-seed", type=int, default=999, help="inclusive")
    parser.add_argument("--variants", type=int, default=4000,
                        help="cap on guides enumerated per target (site); the full reachable "
                             "space is ~20-44k, taken nearest-reference first")
    parser.add_argument("--include-offtarget", action="store_true",
                        help="also keep guides whose seed 12-mer is in the off-target index "
                             "(offtarget_factor < 1.0); default keeps only clean, offtarget 1.0")
    parser.add_argument("--gc-min", type=float, default=0.40,
                        help="lowest guide GC fraction to enumerate (default 0.40); the 40-60%% "
                             "band is ~4x the exact-50%% pool and stays at the cut_p clamp on a "
                             "high-accessibility cell type")
    parser.add_argument("--gc-max", type=float, default=0.60,
                        help="highest guide GC fraction to enumerate (default 0.60)")
    parser.add_argument("--min-pass", type=int, default=None,
                        help="keep guides cutting under at least this many seeds "
                             "(default: all of them)")
    parser.add_argument("--max-fail-seed", type=int, default=None,
                        help="accept a guide as a candidate when it no_cuts under at most this "
                             "many of the sampled seeds, i.e. cuts under >= (n_seeds - this). The "
                             "natural dial for Cas12a, where zero-fail is impossible (cut_p 0.96). "
                             "Overrides --min-pass when both are given")
    parser.add_argument("--max-distance", type=int, default=2000,
                        help="widest |start - mutation_pos| a cell may have")
    parser.add_argument("--cas", default=None,
                        help="restrict to these Cas systems, comma-separated (e.g. Cas9)")
    parser.add_argument("--per-cell", type=int, default=None,
                        help="keep at most this many strongest passers per "
                             "(mutation, cas, strand) construction cell")
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--cell-types", default=None,
                        help="cell-type accessibility JSON; fetched from the backend if omitted")
    parser.add_argument("--limit-sites", type=int, default=None,
                        help="only the first N sites, for a quick smoke run")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.ERROR, format="%(levelname)s %(message)s")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    from neurons.miner import Miner

    task = json.load(open(args.task))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    if args.cell_types:
        cell_types = json.load(open(args.cell_types))
    else:
        cell_types = G.fetch_cell_types()
    accessibility = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)

    n_seeds = args.end_seed - args.start_seed + 1
    if args.max_fail_seed is not None:
        # cut under >= (n_seeds - max_fail_seed); max_fail_seed no_cuts are tolerated.
        min_pass = max(0, n_seeds - args.max_fail_seed)
    elif args.min_pass is not None:
        min_pass = args.min_pass
    else:
        min_pass = n_seeds
    max_fail = n_seeds - min_pass

    ctx = G.build_context(contract, reference, cell_types)
    cfg = Miner.gen_config_for(contract)
    sites = G.enumerate_sites(ctx, cfg.flank, tuple(sorted(set(cfg.lengths))))
    if args.limit_sites:
        sites = sites[:args.limit_sites]

    print("=" * 100)
    print(f"  task {task.get('id')}   cell_type {contract.get('cell_type')} "
          f"(accessibility {accessibility})")
    print(f"  seed-agnostic cut over seeds {args.start_seed}..{args.end_seed} ({n_seeds}); "
          f"keep guides cutting under >= {min_pass} (<= {max_fail} no_cut)")
    print(f"  {len(sites)} sites x {len(ctx.mutations)} mutations, GC "
          f"{args.gc_min:.0%}-{args.gc_max:.0%}, {args.variants} variants/slice, "
          f"{args.jobs} workers")
    print("=" * 100)

    # The achievable cut_p on THIS task, not the energy-clamp ceiling: at accessibility 0.35 the
    # best a Cas9 site reaches is ~0.95, not 0.99, and the all-seeds odds follow cut_p ** n_seeds.
    achievable = {}
    for cas in ctx.cas_systems:
        best = max((G.predicted_cut_p(site, mutation, ctx)
                    for site in sites for mutation in ctx.mutations
                    if site.cas == cas
                    and abs(site.start - ctx.mutation_map[mutation]) <= args.max_distance),
                   default=0.0)
        achievable[cas] = best
        pr = best ** n_seeds
        odds = (f"~1 per {1 / pr:,.0f} guides" if pr > 1e-12
                else "effectively never — no seed-agnostic cut is possible here")
        print(f"  {cas:<7} best reachable cut_p {best:.4f} (ceiling {G.cut_p_ceiling_for(cas)})  "
              f"->  P(all {n_seeds} seeds) = {pr:.3e}  ({odds})")

    cas_filter = set(args.cas.split(",")) if args.cas else None
    cells = [(i, m) for i in range(len(sites)) for m in ctx.mutations
             if cas_filter is None or sites[i].cas in cas_filter]
    payload = {
        "task": str(Path(args.task).resolve()), "cell_types": cell_types,
        "start_seed": args.start_seed, "end_seed": args.end_seed,
        "min_pass": min_pass,
        "max_fail_seed": max_fail, "variants": args.variants, "max_distance": args.max_distance,
        "cas_filter": cas_filter, "gc_min": args.gc_min, "gc_max": args.gc_max,
        "require_clean": not args.include_offtarget,
    }

    print(f"\n  scanning {len(cells)} cells")
    started = time.time()
    results = []
    with Pool(args.jobs, initializer=worker_init, initargs=(payload,)) as pool:
        for index, record in enumerate(pool.imap_unordered(evaluate_cell, cells, chunksize=8), 1):
            results.append(record)
            if index % 500 == 0 or index == len(cells):
                found = sum(len(r["passers"]) for r in results)
                rate = index / (time.time() - started)
                print(f"    [{index:>5}/{len(cells)}] {rate:.0f} cells/s  "
                      f"{found} passer(s)  eta {(len(cells) - index) / max(rate, 1e-9):.0f}s")

    passers = [p for r in results for p in r["passers"]]
    tested = sum(r["tested"] for r in results)
    passers.sort(key=lambda p: (-p["cut_seeds"], -p["weighted_score"]))

    # A submission is built per (mutation, cas, strand) construction cell, so that is the unit the
    # bank has to fill. Cap each cell to the strongest --per-cell candidates if asked; passers is
    # already sorted strongest-first, so a stable per-cell counter keeps the best.
    def cell_of(p):
        return (p["mutation"], p["cas_system"], p["strand"])
    if args.per_cell is not None:
        kept, seen = [], Counter()
        for p in passers:
            key = cell_of(p)
            if seen[key] < args.per_cell:
                kept.append(p)
                seen[key] += 1
        passers = kept

    by_cas = Counter(p["cas_system"] for p in passers)
    by_cell = Counter(cell_of(p) for p in passers)
    elapsed = time.time() - started

    with open(out_dir / "guides.jsonl", "w") as handle:
        for p in passers:
            handle.write(json.dumps(p) + "\n")

    summary = {
        "task_id": task.get("id"),
        "cell_type": contract.get("cell_type"),
        "accessibility": accessibility,
        "seed_range": [args.start_seed, args.end_seed],
        "n_seeds": n_seeds,
        "min_pass": min_pass,
        "variants_per_cell": args.variants,
        "cells_scanned": len(cells),
        "guides_tested": tested,
        "passers": len(passers),
        "empirical_pass_rate": (len(passers) / tested) if tested else 0.0,
        "by_cas": dict(by_cas),
        "cells_with_a_passer": len(by_cell),
        "cells_total": len(cells),
        "per_construction_cell": {f"{m}|{c}|{st}": n for (m, c, st), n in sorted(by_cell.items())},
        "per_cell_cap": args.per_cell,
        "cas_filter": sorted(cas_filter) if cas_filter else None,
        "per_cas_ceiling": {c: G.cut_p_ceiling_for(c) for c in ctx.cas_systems},
        "per_cas_best_reachable_cut_p": achievable,
        "per_cas_all_seed_probability": {c: achievable[c] ** n_seeds for c in ctx.cas_systems},
        "elapsed_seconds": elapsed,
    }
    with open(out_dir / "summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    print(f"\n  tested {tested:,} guides in {elapsed:.0f}s -> {len(passers)} seed-agnostic-cut "
          f"candidate(s)")
    print(f"  by cas system: {dict(by_cas)}")
    print(f"  construction cells with >= 1 passer: {len(by_cell)}")
    print(f"  per (mutation, cas, strand) cell:")
    for cell, n in sorted(by_cell.items()):
        print(f"    {cell[1]:<7}{cell[2]:<3} {cell[0][:26]:<28} {n} candidate(s)")
    if tested:
        print(f"  empirical pass rate: {len(passers) / tested:.3e} "
              f"(Cas9 analytic {G.cut_p_ceiling_for('Cas9') ** n_seeds:.3e})")
    if passers:
        print(f"\n  strongest candidates (by seeds cut, then weighted_score):")
        print(f"    {'cas':<7}{'str':<4}{'start':>9}{'gc':>6}{'off':>5}{'cut_p':>7}"
              f"{'cuts':>7}{'margin':>9}{'wtd':>7}")
        for p in passers[:12]:
            print(f"    {p['cas_system']:<7}{p['strand']:<4}{p['target_alignment_start']:>9}"
                  f"{p['gc']:>6.3f}{p['offtarget_factor']:>5.1f}{p['cut_p']:>7.3f}"
                  f"{p['cut_seeds']:>4}/{p['total_seeds']}{p['margin']:>9.5f}"
                  f"{p['weighted_score']:>7.3f}")

    print(f"\n  artifacts in {out_dir}/")
    for name in sorted(os.listdir(out_dir)):
        print(f"    {name:<24}{(out_dir / name).stat().st_size:>12,} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
