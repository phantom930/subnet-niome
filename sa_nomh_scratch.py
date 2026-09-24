#!/usr/bin/env python3
"""sa_nomh_scratch.py — not_mhnhej as a SEED-AGNOSTIC (window, not single-seed) construction,
enumerated with NO gc_band / max_distance restriction ("all GC, all distance"), from scratch.

Every not_mhnhej test so far in this repo pinned to a KNOWN seed (`seed_depend`'s whole reason to
exist). `seed_agnostic.py` solves the opposite, harder problem -- hedge over a WHOLE seed WINDOW
without knowing which seed will be drawn -- but its machinery (`cut_fail_seeds`, `_scan_site`,
`MT.screen_guides*`) tests ONLY the cut coin; repair mode never enters it, and CLAUDE.md already
records that testing repair-mode rules window-wide tops out at 349-491 failed seeds of 900 and that
"raising max_fail, max_distance or pool size does not move it." not_mhnhej was not named among the
rules that finding covers, and there is a concrete mechanistic reason to check it separately rather
than assume the ceiling transfers.

**Why GC specifically, and why "all" of it.** `stage3.repair_mode` gives
`P(MH_NHEJ | cut) = mh_nhej_w / (hdr_w + mh_nhej_w + blunt_w)`, and BOTH terms in that ratio that can
move are set by GC, not by a coin with no design lever:

    hdr_w      = hdr_base + 0.35 * energy         energy rises with gc (up to the accessibility
                                                    clamp), so hdr_w is largest at high gc
    mh_nhej_w  = 0.30 if mh else 0.12              stage3.microhomology_trigger:
                                                    p(mh=True) = min(0.6, gc*(1-gc)*2.2)
                                                    -- LOWEST at gc near 0 or 1, HIGHEST (0.55) at
                                                    gc = 0.5 -- the CENTRE of the shipped
                                                    SeedAgnosticConfig band (0.40-0.60), not its edge

So the shipped is_cut hedge's own gc_band (chosen to keep weighted_score high, not to help a
repair-mode rule) sits close to the WORST point for not_mhnhej's mh_nhej_w term, and pushing GC
toward either extreme should lower the per-seed failure probability through BOTH channels at once.
Worked example, Cas9, K562-like accessibility (energy already clamped at 1.0 for close guides on
this cell -- see `energy_hdr_k562.py`): at gc=0.5, p_mh=0.55, so
p_fail = 0.55*(0.30/1.32) + 0.45*(0.12/1.14) = 0.55*0.227 + 0.45*0.105 ~ 0.172; at gc=1.0, p_mh=0
exactly, so p_fail = 0.12/1.14 = 0.105 -- a ~39% relative drop, from arithmetic alone, before any
min-union coincidence-picking. Whether real guides near a real PAM site can actually REACH gc near
0/1 within the mismatch budget is an empirical question this script answers rather than assumes --
hence "from scratch": no gc_band, no max_distance, just the raw achievable design space.

**What this measures, not what it ships.** This is a landscape study: enumerate broadly, measure how
the per-guide not_mhnhej failure rate actually varies with the achieved GC (and separately, for
comparison, the plain cut-only failure rate the shipped hedge already solves far better), then run
the same `min_union_group` greedy the shipped hedge uses, unmodified, over this wider pool. It is NOT
running the real seed_agnostic scan (that requires the MT19937 GPU/CPU batched cut-only screen);
`stage3.simulate` is called directly per (candidate, seed) instead, which is only affordable because
not_mhnhej testing is landscape-scale (~1500 candidates), not thousands-of-sites-scale.

    SA_CELL=K562 python sa_nomh_scratch.py [task_id]

RESULT (2026-09-21, K562, task 2af71117, 900-seed window, 1600-candidate pool built with gc_band
0.0-1.0 and no max_distance): dead, and dramatically WORSE than CLAUDE.md's documented repair-mode
ceiling of 349-491 failed seeds -- not inside it, despite the mechanism looking promising on paper.
Achievable GC only reaches 0.150-0.700 (real PAM-adjacent guides can't hit the true 0/1 extremes
within the mismatch budget), and across that whole range the per-guide not_mhnhej failure rate spans
just 0.191-0.222 -- ~3 percentage points, matching the formula's own near-flat theoretical
prediction (0.174-0.181). The aggregate minimum sits INSIDE the shipped 0.40-0.60 gc_band, not at
either extreme, so "all GC" does not even clearly beat the band it was meant to escape.

The decisive result is the group-size curve: `min_union_group` (the shipped hedge's own greedy,
unmodified) collapses explosively as group size grows -- union 122/900 (86.4% clean) at
group_size=1, 479/900 (46.8%) at 5, 833/900 (7.4%) at 20, 886-888/900 (1.3-1.6%) at 40-150, which is
the range any real submission's Cas12a share needs. Plain cut-only on the IDENTICAL pool reaches
532-603/900 (33-41% clean) at the same sizes -- a huge margin at every size. Mechanism: not_mhnhej's
~20% per-guide failure rate is far higher than cut-only's ~2-6% AND nearly uniform across the whole
achievable design space, so unlike whatever gives `hdr` its wider, more exploitable spread, the
min-union greedy has no guide-to-guide heterogeneity here to hunt for -- failures are close to
independent per (guide, seed), so the union tracks a `1-(1-p)^N` independent-failure curve regardless
of pool size or GC/distance breadth. See CLAUDE.md's falsified table and its seed-agnostic ceiling
note for the full writeup.
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "sanomhscratch")
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
from energy_hdr_k562 import pick_task, pearson              # noqa: E402

START_SEED, END_SEED = 100, 999          # matches SeedAgnosticConfig's default window
SEEDS = list(range(START_SEED, END_SEED + 1))
FLANK = 3000                              # "all distance": the whole flank, no max_distance filter
LENGTHS = (20, 23)
CAP_PER_SITE = 30                        # variants kept per (site, mutation), before the seed sweep
CAP_TOTAL = 1600                         # global candidate cap -- bounds the O(candidates x 900) sweep
GROUP_SIZES = (1, 2, 5, 10, 20, 40, 80, 150)


def nomh_p_fail(cas, gc, energy):
    """The formula's own prediction for P(not compliant with not_mhnhej | cut), for comparison."""
    hdr_base = 0.32 if cas == "Cas9" else 0.24
    hdr = hdr_base + 0.35 * energy
    repeat_bias = gc * (1 - gc)
    p_mh = min(0.6, repeat_bias * 2.2)
    mh_nhej = p_mh * 0.30 + (1 - p_mh) * 0.12
    blunt = 0.35
    return mh_nhej / (hdr + mh_nhej + blunt)


def scan(ctx, task, cell_types):
    """Enumerate candidates with NO gc_band / max_distance restriction, interleaved across cells,
    and run the full 900-seed not_mhnhej sweep on each via the real stage3.simulate."""
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
            feat = stage3.extract_features(entry)
            gc, distance = feat["gc"], feat["distance"]
            energy = stage3.sequence_energy(feat)
            cut_p = stage3.cut_probability(site.cas, energy)
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
                "gc": gc, "distance": distance, "energy": energy, "cut_p": cut_p,
                "cut_fails": cut_fails, "fails": nomh_fails, "n_fail": len(nomh_fails),
            })
            if len(out) >= CAP_TOTAL:
                break
    print(f"scanned {scanned_sites} (site,mutation) jobs, banked {len(out)} candidates "
          f"in {time.monotonic() - t0:.1f}s", flush=True)
    return out


def main():
    task_id = _ARGV[1] if len(_ARGV) > 1 else pick_task()
    task, contract, reference = task_content(task_id)
    print(f"task {task_id[:8]}  {CELL}  seed-agnostic window [{START_SEED}, {END_SEED}] "
          f"({len(SEEDS)} seeds), all-GC / all-distance enumeration\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)

    pool = scan(ctx, task, cell_types)
    if not pool:
        raise SystemExit("no candidates found")

    gcs = np.array([c["gc"] for c in pool])
    energies = np.array([c["energy"] for c in pool])
    nomh_rate = np.array([c["n_fail"] / len(SEEDS) for c in pool])
    cut_rate = np.array([len(c["cut_fails"]) / len(SEEDS) for c in pool])
    theo = np.array([nomh_p_fail(c["cas_system"], c["gc"], c["energy"]) for c in pool])

    print(f"achieved GC range: {gcs.min():.3f} - {gcs.max():.3f} (mean {gcs.mean():.3f})")
    print(f"achieved energy range: {energies.min():.3f} - {energies.max():.3f} "
          f"(mean {energies.mean():.3f})")
    print(f"cas mix: {dict(Counter(c['cas_system'] for c in pool))}\n")

    print("=== per-guide not_mhnhej failure rate vs GC (bins away from the shipped 0.40-0.60) ===")
    bins = [(0.0, 0.30), (0.30, 0.40), (0.40, 0.60), (0.60, 0.70), (0.70, 1.01)]
    for lo, hi in bins:
        mask = (gcs >= lo) & (gcs < hi)
        n = int(mask.sum())
        if n == 0:
            print(f"  gc [{lo:.2f},{hi:.2f}): n=0")
            continue
        print(f"  gc [{lo:.2f},{hi:.2f}): n={n:4d}  mean not_mhnhej fail rate "
              f"{nomh_rate[mask].mean():.4f}  mean cut fail rate {cut_rate[mask].mean():.4f}  "
              f"mean theoretical fail rate {theo[mask].mean():.4f}")
    print(f"\n  Pearson(gc, nomh_fail_rate)           = {pearson(gcs, nomh_rate):+.4f}")
    print(f"  Pearson(|gc-0.5|, nomh_fail_rate)      = {pearson(np.abs(gcs - 0.5), nomh_rate):+.4f}"
          "  (predicted negative: further from 0.5 = lower failure)")
    print(f"  Pearson(energy, nomh_fail_rate)        = {pearson(energies, nomh_rate):+.4f}")
    print(f"  Pearson(theoretical, observed) fail rate = {pearson(theo, nomh_rate):+.4f}  "
          "(sanity check on the formula itself)\n", flush=True)

    best = min(pool, key=lambda c: c["n_fail"])
    print(f"best single guide: gc {best['gc']:.3f}  energy {best['energy']:.3f}  "
          f"cas {best['cas_system']}  not_mhnhej fails {best['n_fail']}/{len(SEEDS)}  "
          f"({100 * (1 - best['n_fail'] / len(SEEDS)):.1f}% compliant alone)\n", flush=True)

    cfg = SA.SeedAgnosticConfig(start_seed=START_SEED, end_seed=END_SEED, group_restarts=8)
    print("=== min_union_group: not_mhnhej vs plain cut-only, same pool, same greedy ===")
    for group_size in GROUP_SIZES:
        nomh_group = SA.min_union_group(pool, group_size, cfg, len(SEEDS))
        nomh_union = len({s for c in nomh_group for s in c["fails"]})
        cut_pool = [dict(c, fails=c["cut_fails"]) for c in pool]
        cut_group = SA.min_union_group(cut_pool, group_size, cfg, len(SEEDS))
        cut_union = len({s for c in cut_group for s in c["fails"]})
        print(f"  group_size {group_size:4d}: not_mhnhej union {nomh_union:4d}/900 "
              f"({100 * (900 - nomh_union) / 900:.1f}% clean)   "
              f"cut-only union {cut_union:4d}/900 ({100 * (900 - cut_union) / 900:.1f}% clean)   "
              f"vs CLAUDE.md's documented repair-mode ceiling 349-491: "
              f"{'BEATS it' if nomh_union < 349 else 'inside/above it'}")
    print(flush=True)


if __name__ == "__main__":
    main()
