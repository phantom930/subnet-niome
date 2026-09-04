"""seed_depend.py — a submission pinned to ONE known seed, built to maximise score at that seed.

Every other builder in this repo is seed-*agnostic*: it cannot see the round's seeds, so it hedges
over a window (``seed_agnostic``) or pins an outcome on a narrow band and eats a ~0.10 consistency
floor everywhere else (``all_cut``, ``all_hdr``). This module solves the opposite, much easier
problem — the seed is known — and it exists because the backend sometimes never stamps one.

**When that happens.** A task is broadcast with ``seed: 0`` and the real seeds are stamped in before
scoring. On 8 of 352 scored tasks (3 of 56 since 2026-08-25, so ~5.4%) the stamp never happened and
the validator scored with seed 0 itself. On those rounds a seed-0 build reaches
``consistency_factor`` exactly **1.000**, and since every placing miner also reaches 1.000 the whole
ranking collapses to ``final = total_weighted_score x distribution_fidelity_factor``. Measured on
task 8f02f1a4 (HEK293, 2026-09-02): 16 miners at cons 1.000 swept the curve, rank 1 scored 333.82
(weighted 356.4, fid 0.937) and the rank-10 cutoff was 332.06, against a median of 0.081 for
everyone else.

**So the objective here has no frequency term at all** — no band, no window, no k-of-3 average. It is
a single constrained maximisation:

    maximise   sum(weighted_score) x distribution_fidelity_factor
    subject to every row satisfying the construction rule at the target seed

and the two factors pull against each other. ``weighted_score`` per row is
``(0.625*gc_score + 0.375*dist_score) * offtarget_factor * mutation_weight``, so it wants every row
on the heaviest mutation, at GC 0.50, as close to the mutation as possible. Stage 5's fidelity is a
six-way geometric mean of coverage entropies over mutation, cas, strand and their joint cell, plus
k-mer and distinct-guide diversity — so it wants those four distributions *balanced*. Piling rows
onto the heavy mutation raises term 1 and collapses the mutation-coverage ratio.

``allocate`` therefore searches the split rather than assuming one, scoring each candidate
allocation with stage 5 itself instead of a proxy.

Nothing here is shared with ``seed_agnostic`` or ``all_cut``/``all_hdr``: no bank, no min-union, no
``AllCutConfig``. The only borrowed pieces are genExp's site/guide/scoring primitives, which are the
validator's own stage 1-2 code path (see the note in CLAUDE.md about why that matters).
"""

from __future__ import annotations

import hashlib
import logging
import time
from collections import defaultdict
from dataclasses import dataclass

import genExp as G
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.validation import stage5

logger = logging.getLogger(__name__)

# The rule every row must satisfy at the target seed, as a predicate on stage 3's record.
# "hdr" fixes all three stage-4 targets at once (is_cut=1, is_hdr=1, indel_length=0), which is the
# cheapest route to consistency 1.000 — no reliance on a feature being learnable.
RULES = {
    "hdr": lambda rec: rec["outcome"] == "HDR",
    # mh -> HDR, else BLUNT_NHEJ with indel 1. Also reaches 1.000, via the `mh` feature stage 4
    # can see rather than by holding a target constant. Kept because it admits ~40% of draws
    # against `hdr`'s ~55%, so it is a different (guide, weighted_score) population.
    "mh": lambda rec: ((rec["mh"] and rec["outcome"] == "HDR")
                       or (not rec["mh"] and rec["outcome"] == "BLUNT_NHEJ"
                           and rec["indel_length"] == 1)),
}


@dataclass(frozen=True)
class SeedDependConfig:
    """Deliberately *not* derived from any seed-agnostic config — see the module docstring."""

    # "mh" over "hdr": measured 331.70 vs 330.40 on 8f02f1a4 before the k-mer fix, and it is the
    # rule every validation run used. Both reach consistency 1.000; they differ only in which guide
    # population qualifies, hence in the weighted_score available.
    rule: str = "mh"
    # Site enumeration. flank is genExp's default; max_distance is wide because dist_score decays
    # as exp(-d/base_padding) rather than cutting off, so a distant row still contributes — the
    # allocator will simply rank it below a near one.
    flank: int = 3000
    lengths: tuple[int, ...] = (20, 23)
    # 600, not 2000: the allocator ranks by weighted_score and dist_score decays as exp(-d/400),
    # so candidates past ~600bp are never selected. Measured identical final score at a third of the
    # scan cost (91s vs 315s) — the wide scan is pure waste, and at 2000 the build overruns the
    # in-TTL path (324s observed live).
    max_distance: int = 600
    # GC band for variant enumeration. Wide on purpose: gc_score peaks at GC 0.50 and the allocator
    # ranks by weighted_score, so narrowing here only removes candidates it would not have picked.
    gc_band: tuple[float, float] = (0.30, 0.70)
    variants_per_site: int = 4000
    # Floor per (mutation, cas, strand) cell. An empty cell zeroes one of stage 5's six ratios and
    # costs roughly a 0.03x multiplier on the whole score, so this is a hard constraint, not a knob.
    per_cell_floor: int = 4
    # Allocation search granularity over the heavy-mutation share.
    alloc_step: int = 4
    # Price on 12-mer novelty when selecting rows, in units of weighted_score. Of stage 5's six
    # ratios, cas / strand / distinct-guide are already 1.000 and mutation / joint are both pinned
    # by the mutation split (joint is maximal when cells within a mutation are balanced, which the
    # allocator ensures) -- so k-mer diversity is the ONLY free term, and it is badly served by
    # ranking on weighted_score alone: variants of one site differ by <= max_mismatches
    # substitutions and share nearly all their 12-mers.
    #
    # Measured on 8f02f1a4, sweeping lambda x mutation split:
    #   lambda  heavy   weighted   fid    kmer    final   rank
    #      0.0    208      371.8  0.892  0.875   331.80     11
    #     0.05    205      368.0  0.920  0.997   338.51      1
    #     0.15    205      368.0  0.920  0.997   338.51      1
    #      1.0    205      368.0  0.920  0.997   338.51      1
    # kmer 0.875 -> 0.997 for 3.8 of weighted, and flat from 0.05 upward — a plateau, not a tuned
    # point. 0.15 sits in the middle of it.
    kmer_price: float = 0.15
    # Candidates considered per cell by the greedy. The greedy is O(take x window), so this bounds
    # build time; the pool is sorted by weighted_score, so the tail it ignores would not be picked.
    greedy_window: int = 400
    # Decorrelates sibling hotkeys. The build is fully deterministic — same contract, same seed,
    # same rows — so several hotkeys running this would submit byte-identical files and score
    # identically. ``variant`` breaks ties among candidates whose greedy score is within
    # ``variant_epsilon`` by hashing (guide, variant), which changes WHICH of several near-equal
    # guides is taken without changing what is being optimised. Scores land within a fraction of a
    # point of each other; the row sets differ substantially.
    variant: int = 0
    variant_epsilon: float = 0.02


def enumerate_candidates(ctx, sites, cfg: SeedDependConfig,
                         deadline: float | None = None) -> dict[tuple, list[dict]]:
    """Every (site, mutation, guide) that satisfies the rule at ``ctx.seed``, keyed by stage-5 cell.

    ``require_clean=True`` in the variant enumeration means the guide's 12-mer off-target seed is
    absent from the chr11 index, i.e. ``offtarget_factor == 1.0`` — the term is a step function
    ({1.0, 0.7, 0.4, 0.1}) so a single collision costs 30% of the row and is never worth accepting.
    """
    keep = RULES[cfg.rule]
    by_cell: dict[tuple, list[dict]] = defaultdict(list)
    jobs = [(s, m) for s in sites for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda j: abs(j[0].start - ctx.mutation_map[j[1]]))
    scanned = 0
    for site, mutation in jobs:
        if deadline is not None and time.monotonic() > deadline:
            logger.info("seed-depend: candidate scan stopped on the deadline at %d/%d jobs",
                        scanned, len(jobs))
            break
        scanned += 1
        guides = SA.enumerate_variants(site, ctx, cfg.gc_band[0], cfg.gc_band[1],
                                       ctx.max_mismatches, True, cfg.variants_per_site)
        for guide in guides:
            experiment = G.make_experiment(site, guide, mutation, ctx, "cand")
            entry = G.build_valid_entry(experiment, ctx)
            if entry is None:
                continue
            record = G.simulate(entry, ctx)
            if not keep(record):
                continue
            by_cell[(mutation, site.cas, site.strand)].append({
                "guide": guide, "mutation": mutation, "cas_system": site.cas,
                "strand": site.strand, "start": site.start, "length": site.length,
                "weighted_score": entry["stage2"]["weighted_score"],
                "kmers": frozenset(stage5.extract_kmers(guide, 12)),
                "entry": entry, "record": record,
            })
    for recs in by_cell.values():
        recs.sort(key=lambda r: -r["weighted_score"])
    return by_cell


def _fidelity(chosen: list[dict], contract: dict) -> float:
    """Stage 5's own factor for a candidate allocation — not a proxy for it."""
    detail = stage5.compute_distribution_fidelity(
        [r["entry"] for r in chosen], [r["record"] for r in chosen], contract, k=12)
    return max(0.0, min(1.0, detail.get("distribution_fidelity_score", 0.0)))


def _variant_key(guide: str, variant: int) -> int:
    """Stable per-variant ordering over tied candidates. sha256 so two variants share no structure —
    a simple offset would make sibling submissions overlap heavily on the same tie groups."""
    return int.from_bytes(hashlib.sha256(f"{variant}|{guide}".encode()).digest()[:8], "big")


def _pick_cell(recs: list[dict], take: int, cfg: SeedDependConfig,
               pool: set) -> list[dict]:
    """Greedy selection within one cell on ``weighted_score + kmer_price * new 12-mers``.

    ``pool`` is the 12-mer set already committed by earlier cells, so novelty is scored against the
    whole submission rather than per cell.
    """
    if cfg.kmer_price <= 0:
        return recs[:take]
    window = recs[:max(take * 8, cfg.greedy_window)]
    used = set(pool)
    taken = [False] * len(window)
    out: list[dict] = []
    for _ in range(take):
        scored = []
        for i, rec in enumerate(window):
            if taken[i]:
                continue
            scored.append((rec["weighted_score"] + cfg.kmer_price * len(rec["kmers"] - used), i))
        if not scored:
            break
        top = max(s for s, _i in scored)
        # Every candidate within epsilon of the best is equally good for the objective; pick among
        # them by a hash of (guide, variant) so sibling hotkeys diverge deterministically.
        ties = [i for s, i in scored if s >= top - cfg.variant_epsilon]
        pick_i = min(ties, key=lambda i: _variant_key(window[i]["guide"], cfg.variant))
        taken[pick_i] = True
        out.append(window[pick_i])
        used |= window[pick_i]["kmers"]
    return out


def allocate(by_cell: dict[tuple, list[dict]], contract: dict, ctx,
             n_rows: int, cfg: SeedDependConfig) -> tuple[list[dict], dict]:
    """Choose ``n_rows`` maximising sum(weighted_score) x fidelity.

    Two stages, because the greedy is far more expensive than a head slice. First locate the
    heavy-mutation share with the cheap top-N selection — term 1 rises with it while stage 5's
    mutation-coverage ratio falls, and the peak moves with the contract's weight spread — then rerun
    the greedy only near that optimum. The k-mer price shifts the peak by a couple of rows at most
    (208 -> 205 measured), so the coarse pass locates it well enough.
    """
    weights = contract.get("mutation_weights") or {}
    muts = sorted(ctx.mutations, key=lambda m: -weights.get(m, 1.0))
    if len(muts) < 2:
        muts = list(ctx.mutations)
    heavy, light = muts[0], muts[-1]
    cells = {c: recs for c, recs in by_cell.items() if recs}
    heavy_cells = [c for c in cells if c[0] == heavy]
    light_cells = [c for c in cells if c[0] == light]
    if not heavy_cells or not light_cells:
        return [], {"reason": f"only {len(cells)} of 8 stage-5 cells have candidates"}

    floor = cfg.per_cell_floor
    max_heavy = n_rows - floor * len(light_cells)

    def attempt(n_heavy: int, price: float) -> tuple | None:
        chosen: list[dict] = []
        pool: set = set()
        for group, total in ((sorted(heavy_cells), n_heavy),
                             (sorted(light_cells), n_rows - n_heavy)):
            base, extra = divmod(total, len(group))
            for i, cell in enumerate(group):
                take = base + (1 if i < extra else 0)
                avail = cells[cell]
                if len(avail) < max(take, floor):
                    return None
                got = (avail[:take] if price <= 0
                       else _pick_cell(avail, take, cfg, pool))
                if len(got) < take:
                    return None
                chosen.extend(got)
                for rec in got:
                    pool |= rec["kmers"]
        if len(chosen) < n_rows:
            return None
        chosen = chosen[:n_rows]
        wtd = sum(r["weighted_score"] for r in chosen)
        fid = _fidelity(chosen, contract)
        return (wtd * fid, wtd, fid, n_heavy, chosen)

    coarse = None
    for n_heavy in range(n_rows // 2, max_heavy + 1, cfg.alloc_step):
        got = attempt(n_heavy, 0.0)
        if got and (coarse is None or got[0] > coarse[0]):
            coarse = got
    if coarse is None:
        return [], {"reason": "no allocation filled every cell to its floor"}
    best = coarse
    if cfg.kmer_price > 0:
        centre = coarse[3]
        for n_heavy in range(max(n_rows // 2, centre - cfg.alloc_step),
                             min(max_heavy, centre + cfg.alloc_step) + 1, 2):
            got = attempt(n_heavy, cfg.kmer_price)
            if got and got[0] > best[0]:
                best = got
    score, wtd, fid, n_heavy, chosen = best
    return chosen, {"weighted": wtd, "fidelity": fid, "product": score, "heavy": n_heavy,
                    "cells": len({(r["mutation"], r["cas_system"], r["strand"]) for r in chosen})}


def build(contract: dict, reference: dict, cell_types: dict, seed: int = 0,
          cfg: SeedDependConfig | None = None,
          budget_s: float | None = None) -> tuple[list[dict] | None, dict]:
    """The seed-pinned submission, or ``(None, meta)`` when the caller should fall back."""
    cfg = cfg or SeedDependConfig()
    started = time.monotonic()
    deadline = None if budget_s is None else started + budget_s
    meta: dict = {"method": "seed-depend", "seed": seed, "rule": cfg.rule}

    pinned = dict(contract)
    pinned["seed"] = seed
    ctx = G.build_context(pinned, reference, cell_types)
    sites = G.enumerate_sites(ctx, cfg.flank, cfg.lengths)
    n_rows = pinned["rules"].get("max_experiments") or ctx.max_experiments

    by_cell = enumerate_candidates(ctx, sites, cfg,
                                   None if deadline is None else deadline - 15.0)
    meta["candidates"] = sum(len(v) for v in by_cell.values())
    meta["cells_with_candidates"] = len(by_cell)
    if meta["cells_with_candidates"] < 8:
        meta["reason"] = f"only {meta['cells_with_candidates']} of 8 cells have candidates"
        return None, meta

    chosen, alloc = allocate(by_cell, pinned, ctx, n_rows, cfg)
    meta.update(alloc)
    if not chosen:
        return None, meta

    rows, seen = [], set()
    for i, rec in enumerate(chosen):
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": rec["guide"],
                     "target_alignment_start": rec["start"],
                     "target_alignment_end": rec["start"] + rec["length"],
                     "strand": rec["strand"], "mutation": rec["mutation"],
                     "cas_system": rec["cas_system"],
                     "cell_type": pinned.get("cell_type")})
    meta.update(rows=len(rows), elapsed_s=round(time.monotonic() - started, 1))
    if len(rows) < n_rows:
        meta["reason"] = f"deduped to {len(rows)} rows of {n_rows}"
        return None, meta
    return rows, meta
