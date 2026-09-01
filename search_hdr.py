#!/usr/bin/env python
"""search_hdr.py — all-HDR candidate guides for one task, screened over two seed windows.

Research tool, not imported by the neurons. It answers: how many guides survive both

  * **cut** on at least ``900 - max_fail_seed`` of the seeds 100-999, and
  * **HDR** on at least ``100 - max_main_fail_seed`` of the seeds 500-599

for a given cas system. The HDR window is a subset of the cut window, so a guide that reaches HDR
at a seed necessarily cut at it; the wide window is the seed-agnostic cut hedge, the narrow one
pins the repair mode over a band.

**The narrow window is screened first, and that ordering is the whole performance story.** The HDR
tail is ~1e-6 per guide against ~1e-2 for the cut gate, so running it first on 100 seeds kills
almost everything before the 900-seed pass ever runs — about 9x less work than the natural order,
and the second pass sees only a handful of survivors. Both passes reuse
``mt19937.screen_guides_rule_gpu`` rather than adding kernel code: it is already bit-exact against
``random.Random`` and covered by ``verify()``.

**GC is the only design lever on repair mode, and it is worth orders of magnitude here.** Stage 3
draws the microhomology coin before the repair mode, and ``mh`` raises ``mh_nhej`` from 0.12 to
0.30, which pushes probability away from HDR. ``p_mh`` is ``min(0.6, 2.2*gc*(1-gc))``, maximised at
gc 0.50 and falling off either side — so high GC both keeps energy at the clamp *and* suppresses
the competing mode. On Cas12a/CD34+ that moves P(HDR) from 0.4917 at gc 0.50 to 0.5190 at gc 0.90,
which sounds small and is not: over a 100-seed window needing >= 75 hits it is 15x more candidates.
That is why the default band is high rather than centred on 0.50 the way every other builder's is.

Earlier work found repair-mode rules had "no design lever" and topped out at 349-491 failed seeds
of 900. That stands for rules over the *whole* window; it does not transfer to a 100-seed band with
a fail budget, which is a far weaker requirement.

    python search_hdr.py --task task-cd34.json --cas Cas12a            # band 500-599 by default
    python search_hdr.py --task task-hudep.json --cas Cas12a           # band 800-899 by default
    python search_hdr.py --task task-cd34.json --hdr-range 300 399 --groups 42
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import genExp as G                                                      # noqa: E402
from niome_subnet.genomics import mt19937 as MT                         # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA                   # noqa: E402
from niome_subnet.genomics.all_cut import _params_fn                    # noqa: E402
from niome_subnet.genomics.validation import stage3                     # noqa: E402

CUT_WINDOW = (100, 999)          # the seed-agnostic window every other builder uses

# The band the repair mode is pinned over, per cell type. Its *width* is what costs: every seed
# added multiplies the conditional Cas9 requirement by P(HDR) (~0.57), so a 16-seed clean band is
# already near the edge of feasibility and a 31-seed one produced zero Cas9 rows. Its *position*
# is free -- the seeds are independent, so any 100-wide band screens the same and there is no
# reason to expect one to beat another beyond sampling noise.
HDR_RANGE_DEFAULTS: dict[str, tuple[int, int]] = {
    "CD34+_HSPC": (500, 599),
    "HUDEP-2": (800, 899),
    "K562": (700, 799),
    "HEK293": (300, 399),
}
FALLBACK_HDR_RANGE = (500, 599)


def hdr_range_for(cell_type: str, override) -> tuple[int, int]:
    """The pinned band for this cell type, or the explicit override."""
    if override:
        lo, hi = int(override[0]), int(override[1])
        if hi <= lo:
            raise SystemExit(f"--hdr-range needs LO < HI, got {lo} {hi}")
        return lo, hi
    return HDR_RANGE_DEFAULTS.get(cell_type, FALLBACK_HDR_RANGE)


def load_task(path: str) -> tuple[dict, dict]:
    task = json.loads(Path(path).read_text())
    content = task.get("content", task)
    contract = content["contract"]
    reference = content["hbb_reference"]
    return contract, reference


def select_groups(bank: list[dict], sizes: list[int], restarts: int, per_cell_min: int,
                  main_window: tuple[int, int]) -> None:
    """Min-union a group over the HDR failed-seed sets and report the surviving clean band.

    The same selector ``all_cut`` uses on *cut* failures, pointed at HDR failures instead, over the
    100-seed band rather than the 900-seed window. The clean band is the complement of the group's
    failed-seed union: the seeds on which *every* row in the group repairs by HDR. That is the
    prize — with is_cut, is_hdr and indel_length all constant, stage 4's three targets degenerate
    and consistency_factor reaches 1.0 instead of the ~0.20 the cut-only hedge gets.

    The arithmetic is far harsher than the cut hedge's, and it is worth seeing why before reading
    the numbers. There each guide fails 12-22 seeds of 900 (~2%), so unions stay small and a group
    of 42 still leaves a clean majority. Here each guide fails 26-45 of 100 (26-45%), so a union of
    k *random* candidates covers 100*(1 - 0.51**k) seeds — 74% at k=2, 87% at k=3. Only genuine
    coincidence between fail sets keeps a band alive, which is exactly what min-union selects for.
    """
    from niome_subnet.genomics import fastgreedy as FG

    lo, hi = main_window
    print(f"\n### min-union over HDR fail sets, window {lo}-{hi}"
          f" | {len(bank)} candidates | per_cell_min {per_cell_min}\n", flush=True)
    selector = FG.FastGreedy(bank, window_lo=lo, window_hi=hi, per_cell_min=per_cell_min)
    span = hi - lo + 1
    print(f"  {'group':>6}{'union':>8}{'clean':>8}{'clean%':>9}{'cells':>7}"
          f"{'random baseline':>18}", flush=True)
    for k in sizes:
        if k > len(bank):
            print(f"  {k:>6}   only {len(bank)} candidates available", flush=True)
            continue
        index, union = selector.best(k, restarts=restarts)
        chosen = [bank[int(i)] for i in index]
        clean = span - union
        cells = len({(c["mutation"], c["cas_system"], c["strand"]) for c in chosen})
        mean_fail = sum(c["hdr_fails"] for c in chosen) / len(chosen) / span
        expected_union = span * (1 - (1 - mean_fail) ** k)
        print(f"  {k:>6}{union:>8}{clean:>8}{100 * clean / span:>8.0f}%{cells:>5}/4"
              f"{span - expected_union:>15.1f} clean", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", default="task-cd34.json")
    ap.add_argument("--cas", default="Cas12a", choices=("Cas12a", "Cas9"))
    ap.add_argument("--max-fail-seed", type=int, default=100,
                    help="seeds of 900 the guide may fail to cut on (default 100)")
    ap.add_argument("--max-main-fail-seed", type=int, default=35,
                    help="seeds in the HDR band the guide may fail to reach HDR on (default 35)")
    ap.add_argument("--hdr-range", type=int, nargs=2, default=None, metavar=("LO", "HI"),
                    help="the seed band the repair mode is pinned over; defaults per cell type "
                         "(HEK293 300-399, CD34+_HSPC 500-599, K562 700-799, HUDEP-2 800-899)")
    ap.add_argument("--gc", type=float, nargs=2, default=(0.70, 0.95), metavar=("LO", "HI"),
                    help="guide GC band; high suppresses the microhomology coin (default 0.70 0.95)")
    ap.add_argument("--max-distance", type=int, default=400)
    ap.add_argument("--variants", type=int, default=44000)
    ap.add_argument("--flank", type=int, default=3000)
    ap.add_argument("--out", default=None, help="write surviving candidates to this JSON path")
    ap.add_argument("--groups", type=int, nargs="*", default=None, metavar="K",
                    help="after screening, min-union a group of each size over the HDR fail sets")
    ap.add_argument("--restarts", type=int, default=12)
    ap.add_argument("--per-cell-min", type=int, default=2,
                    help="floor per (mutation, cas, strand) cell inside the group")
    args = ap.parse_args()

    contract, reference = load_task(args.task)
    cell = contract.get("cell_type")
    cell_types = G.fetch_cell_types()
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}

    print(f"task {args.task} | {cell} (accessibility {accessibility}) | cas {args.cas}")
    print(f"mutations {contract.get('active_mutations')}")
    main_lo, main_hi = hdr_range_for(cell, args.hdr_range)
    main_span = main_hi - main_lo + 1
    src = "override" if args.hdr_range else f"default for {cell}"
    print(f"cut: <= {args.max_fail_seed} fails of 900 over {CUT_WINDOW[0]}-{CUT_WINDOW[1]}")
    print(f"HDR: <= {args.max_main_fail_seed} fails of {main_span} over {main_lo}-{main_hi} ({src})")
    print(f"GC band {args.gc[0]:.2f}-{args.gc[1]:.2f} | max_distance {args.max_distance} "
          f"| {args.variants} variants/site\n", flush=True)

    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, args.flank, (20, 23))
    cut_seeds = np.arange(CUT_WINDOW[0], CUT_WINDOW[1] + 1, dtype=np.int64)
    main_seeds = np.arange(main_lo, main_hi + 1, dtype=np.int64)

    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == args.cas for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= args.max_distance]
    print(f"{len(jobs)} (site, mutation) targets within {args.max_distance} bp\n", flush=True)

    bank, started, scanned_guides = [], time.monotonic(), 0
    for index, (site_index, mutation) in enumerate(jobs, 1):
        site = sites[site_index]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, args.gc[0], args.gc[1],
                                       ctx.max_mismatches, True, args.variants)
        if not guides:
            continue
        scanned_guides += len(guides)
        params_of = _params_fn(site, distance, accessibility, offset)

        # Pass 1 — the narrow HDR window. Tight, cheap, and it does almost all the killing.
        hdr_ok = MT.screen_guides_rule_gpu(guides, main_seeds, mutation, site.cas, site.start,
                                           site.strand, params_of, "hdr", args.max_main_fail_seed)
        if not hdr_ok:
            continue
        # Pass 2 — the wide cut window, survivors only.
        cut_ok = MT.screen_guides_rule_gpu(list(hdr_ok), cut_seeds, mutation, site.cas, site.start,
                                           site.strand, params_of, "cut", args.max_fail_seed)
        for guide, cut_fails in cut_ok.items():
            gc, energy, cut_p = params_of(guide)
            bank.append({"guide": guide, "mutation": mutation, "cas_system": site.cas,
                         "strand": site.strand, "start": site.start, "length": site.length,
                         "distance": distance, "gc": round(gc, 4), "cut_p": round(cut_p, 4),
                         "cut_fails": int(cut_fails.size),
                         "hdr_fails": int(hdr_ok[guide].size),
                         # the actual failed-seed set is what min-union selects over; keeping only
                         # its size would make the group step impossible
                         "fails": hdr_ok[guide]})
        if index % 25 == 0 or index == len(jobs):
            print(f"  {index}/{len(jobs)} targets | {scanned_guides:,} guides screened | "
                  f"{len(bank)} candidates | {time.monotonic() - started:.0f}s", flush=True)
        MT.free_gpu_memory()

    elapsed = time.monotonic() - started
    print(f"\n### {len(bank)} candidates from {scanned_guides:,} guides in {elapsed:.0f}s")
    if not bank:
        print("### none survived — relax --max-main-fail-seed (the binomial tail is very steep: "
              "35 -> ~5.5e-3, 30 -> ~1.8e-4, 25 -> ~1.8e-6 per guide)")
        return 1
    hf = np.asarray([b["hdr_fails"] for b in bank])
    cf = np.asarray([b["cut_fails"] for b in bank])
    gcs = np.asarray([b["gc"] for b in bank])
    print(f"  hdr_fails  min {hf.min()} median {int(np.median(hf))} max {hf.max()}  "
          f"(of {main_span} seeds in {main_lo}-{main_hi})")
    print(f"  cut_fails  min {cf.min()} median {int(np.median(cf))} max {cf.max()}  "
          f"(of 900 seeds in {CUT_WINDOW[0]}-{CUT_WINDOW[1]})")
    print(f"  gc         min {gcs.min():.3f} median {np.median(gcs):.3f} max {gcs.max():.3f}")
    print(f"  by mutation {dict(Counter(b['mutation'] for b in bank))}")
    print(f"  by strand   {dict(Counter(b['strand'] for b in bank))}")
    best = min(bank, key=lambda b: (b["hdr_fails"], b["cut_fails"]))
    print(f"  best: {best['guide']} hdr_fails={best['hdr_fails']} cut_fails={best['cut_fails']} "
          f"gc={best['gc']:.3f} d={best['distance']}")
    if args.groups:
        select_groups(bank, args.groups, args.restarts, args.per_cell_min,
                      (main_lo, main_hi))
    if args.out:
        Path(args.out).write_text(json.dumps(
            [{k: v for k, v in b.items() if k != "fails"} for b in bank], indent=1))
        print(f"### wrote {len(bank)} candidates to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
