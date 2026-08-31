"""all_cut.py — the all-cut construction: every row cuts on a shared seed window.

Named for HEK293, where it was developed, but it is not HEK293-specific and gains far more
elsewhere. All four cell types the backend issues now have a measured config. Against whatever
each shipped before it: +23.45 (t=8.4) on K562, +17.92 (t=4.9, 5/5 contracts) on HUDEP-2 and
+8.56 expected on CD34+_HSPC, all three over the seed-agnostic hedge, against +4 over HEK293's own
clustered builder. See ``CELL_CONFIG`` for what differs per cell type and why.

The construction is conditional, in two stages:

1. **Cas12a** — min-union a group of ``group_size`` guides from a cut-only bank. Their failed-seed
   union defines the seeds this submission cannot save; the complement is the *clean set*.
2. **Cas9** — require strictness only over that clean set, not over all 900 seeds. Strict Cas9
   rows add nothing to the union, so the clean set is preserved exactly:
   ``final clean = Cas12a-clean AND Cas9-clean``.

What that relaxation buys depends on the cell type. On HEK293 (accessibility 0.35) energy never
reaches the clamp, ``0.9564**900`` is 1e-18 and strict Cas9 over the window simply does not exist —
the relaxation makes the hedge *possible*. On K562 it already exists (production finds 350-564
guides), and ``0.99**566`` against ``0.99**900`` is ~36x more candidates — the relaxation makes the
rows *better*, which is worth more.

Measured on the reference task (40 random clean seeds, 40 dirty, scored through all five stages):

    group  clean   Cas9 gc   clean payout   dirty   expected
       20  29.2%     0.724          41.10   18.48      25.09
       38  16.4%     0.500          55.77   22.84      28.26
       42  14.6%     0.494          69.17   21.71      28.61   <- this module
       50  12.7%     0.494          64.56   22.11      27.49
      production        -           28.57   29.36      29.22

Group 42 is chosen for the *tail*, not the mean: expected 28.61 against production's 29.22 is a
statistical tie (t ~ 1.1), but it pays 69.17 on the 14.6% of seeds that come up clean, against a
flat 28.6. Under ``SCORING_SYSTEM = "top"``, which pays only the top 10 miners on a fixed curve, a
heavy right tail should out-earn a flat build of equal mean. That is an argument from the payout
structure, not a measured rank outcome.

Two knobs explain the shape of that table, and neither is obvious:

* **Row GC.** A larger Cas12a group shrinks the clean set, which makes strict Cas9 cheap enough to
  take rows near gc 0.50 where stage 2's ``gc_score`` peaks. GC falls 0.724 -> 0.494 across the
  sweep, lifting *both* legs -- the dirty leg improves too, since better rows score better on every
  seed.
* **Cas mix.** Past group ~34 the GC is already at its optimum and the payout is driven by the cas
  ratio: group 42 is 83/17, inside the 82/18-88/12 plateau the K562 mix sweep found independently.

**The Cas12a bank is the expensive half and it is built per task.** It is cached on disk, keyed on
everything it depends on (cell type, accessibility, mutations, regions, band, distance, thresholds)
and *not* on the seed, which it is independent of -- but every task the backend has issued carries
a distinct mutation set, so across rounds that cache never hits. It earns its keep only within a
session, which is what makes offline sweeps affordable.

Building it per task used to be impossible against the 300 s upload TTL, and that constraint is
what set K562's bank width (see ``CELL_CONFIG``). ``Miner._prefetch_loop`` removed it by building
when the task appears rather than when a validator asks -- a median 30 minutes earlier -- so the
bank now costs ~127 s on HEK293 and 247-412 s on the three erythroid types, inside a 900 s
prepare budget. Only HEK293 still fits the ~225 s in-TTL fallback; see
Miner.ALL_CUT_MIN_BUDGET_S.

The Cas9 stage scans targets nearest the mutation first and stops once it has enough, because the
assembly ranks candidates by ``(distance, |gc - 0.50|)`` anyway -- a full scan produced 71,256
candidates in 248s to fill 208 rows.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import time
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

import genExp as G
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.validation import stage3

logger = logging.getLogger(__name__)

BANK_DIR = "data/all_cut"

# Per cell type. They differ only in cas9_gc, and that one difference is the whole story of how
# cut_p behaves in each: HEK293's accessibility of 0.35 means energy never clamps, so high GC is
# the only route to a usable cut_p at all and the band has to reach 0.80. The erythroid types
# clamp instead -- K562 at 0.77, HUDEP-2 at 0.82, CD34+_HSPC at 0.87, all above ~gc 0.19 -- so
# reaching past 0.60 buys no cut_p and only costs gc_score.
#
# cas12a_max_fail used to be the setting that decided whether the bank fit the upload TTL: on K562
# guides fail only 12-22 times of 900, so a threshold of 100 admits everything and defeats the
# early-out (measured 265-372s against the ~220s budget the in-TTL path had; at 22 it was 31-52s).
# K562 was pinned to 22/d200 for that reason alone, not because the narrower bank scored better.
#
# The round prefetch removed that constraint — Miner._prefetch_loop builds when the task appears,
# a median 30 minutes before a validator asks — so K562 was moved to the same mf100/d400 bank
# HEK293 uses.
#
# All three candidate configs were then scored against each other -- five contracts x 40 seeds,
# paired within (contract, seed), so the same stage-3 stream judges every row set:
#
#   contract    mf22/d200        mf22/d400        mf100/d400
#   8c7174bd     51s  43.64       96s  41.59      407s  39.55
#   a5b89b25          declined   133s  67.09      298s  67.22
#   f35f29b2     43s  41.91       89s  43.42      345s  46.74
#   e0c604a4     31s  40.31       70s  44.24      268s  48.62
#   536fe450     78s  56.95       97s  56.12      374s  59.57
#   -------------------------------------------------------------------------------
#   mf100/d400 vs mf22/d200    +2.36 +- 2.10  t = +1.12  CI [-1.76, +6.47]
#   mf22/d400  vs mf22/d200    +0.64 +- 1.32  t = +0.48  CI [-1.95, +3.23]
#   mf22/d400  vs mf100/d400   -1.85 +- 1.21  t = -1.53  CI [-4.22, +0.52]
#
# **No pair separates.** Point estimates order them mf100/d400 > mf22/d400 > mf22/d200, but every
# interval contains zero, and the wide bank's sign is not even consistent across contracts
# (-4.09, +0.13, +4.82, +8.31, +2.62). mf100/d400 is retained deliberately on that point estimate,
# with the costs below understood -- it is not a measured win, and nothing here should be cited as
# one.
#
# Do not re-derive any of this from clean fraction. It failed as a proxy twice, in both directions:
# the wide bank gains 7-15 clean seeds of 900 over d200 yet scores worse on one contract, and
# mf22/d400 reaches a *higher* clean count than mf100/d400 (558/564/553/574 against
# 557/559/552/566) while scoring 1.85 lower. Term 1 is flat-to-worse on the wide bank
# (weighted 246.8 -> 246.6, 219.5 -> 217.8, 241.0 -> 238.8, 311.4 -> 304.8); the whole delta rides
# on consistency_factor, i.e. which seeds happen to come up clean, which is noisy seed to seed.
#
# What the choice actually costs, since score cannot separate the three:
#
# * The coverage argument for leaving d200 is real but was already solved -- d200 declines 1 in 5
#   on cell coverage, and RETRY_CONFIG's old d400/mf22 recovered exactly that contract. Both d400
#   configs build 5/5, which is why RETRY_CONFIG is now empty.
# * The build is 4-6x slower than mf22/d400 and no longer fits the ~225s in-TTL fallback, so
#   Miner.ALL_CUT_MIN_BUDGET_S had to become per cell type and **K562 all-cut is prefetch-dependent**:
#   a round whose prefetch fails drops to the seed-agnostic hedge, ~64 -> ~41. Granting the wide
#   bank its full +1.85 point estimate over mf22/d400, that trade breaks even at a prefetch failure
#   rate of 8%; above it, mf22/d400 is ahead on expectation.
#
# mf22/d400 is therefore the fallback of record if the prefetch ever proves unreliable: same score
# within noise, same 5/5 coverage, 70-133s, and it fits back under a flat 190s gate.
#
# HUDEP-2 takes K562's config verbatim, and unlike the K562 tuning above it is a *measured* win —
# against the seed-agnostic hedge, which is what every HUDEP-2 round shipped before this entry
# existed. Five contracts x 40 seeds, paired within (contract, seed):
#
#   contract    hedge   all-cut   paired delta
#   a2af3987    40.66    61.51    +20.85 +- 4.29   t = 4.9
#   9503b973    44.17    55.90    +11.73 +- 4.31   t = 2.7
#   beec21a1    31.93    39.15     +7.22 +- 3.04   t = 2.4
#   b25134d3    41.73    68.57    +26.84 +- 4.75   t = 5.6
#   836ed169    41.72    64.66    +22.94 +- 4.10   t = 5.6
#   ------------------------------------------------------------
#   clustered   +17.92 +- 3.65  t = 4.9  95% CI [+10.77, +25.07]; wins 143/200 seed pairs
#
# Positive on 5/5 contracts, each individually significant, and built 5/5. What marks it as real
# rather than the consistency_factor jitter the mf22/mf100 comparison turned up: *all three* score
# terms rise on every contract — weighted 226-342 -> 240-350, consistency 0.14-0.16 -> 0.18-0.23,
# fidelity up throughout.
#
# The band transfers from K562 on the clamp, not on resemblance: at accessibility 0.82 HUDEP-2
# reaches cut_p 0.990 Cas9 / 0.960 Cas12a at gc 0.40 against K562's 0.990 / 0.953, so it is at
# least as clamped and higher GC would only cost gc_score. HEK293's wide band is the wrong
# transfer here. Builds took 247-386s over six cold builds, the same regime as K562, so the same
# prefetch dependency and the same 480s gate apply (the gate only decides whether to start; the
# budget that bounds the build is PREPARE_BUDGET_S).
CELL_CONFIG: dict[str, dict] = {
    "HEK293": {"max_distance": 400, "cas12a_max_fail": 100, "cas9_gc": (0.40, 0.80)},
    "K562": {"max_distance": 400, "cas12a_max_fail": 100, "cas9_gc": (0.40, 0.60)},
    "HUDEP-2": {"max_distance": 400, "cas12a_max_fail": 100, "cas9_gc": (0.40, 0.60)},
    "CD34+_HSPC": {"max_distance": 400, "cas12a_max_fail": 100, "cas9_gc": (0.40, 0.60)},
}

# CD34+_HSPC takes the same config again, and is the weakest of the four -- not because the method
# works less well where it runs, but because it cannot always run. Five contracts x 40 seeds
# against the seed-agnostic hedge:
#
#   contract    hedge   all-cut   delta
#   97fd0b34    30.17    44.03    +13.87   t = 4.3
#   60f5f6a1    33.88    45.03    +11.15   t = 3.0
#   b0fad64f    35.03    52.84    +17.82   t = 5.1
#   d03c7439    25.63   declined (7/8 cells)
#   d2ca8e6e    26.61   declined (7/8 cells)
#
# +14.27 +- 1.93 where it builds, but only 3 of 6 contracts build (a sixth, 470d6d3c, declined on
# cold verification) -> **~+7.1 expected per round**. A decline costs nothing
# (RETRY_CONFIG is empty, so the round ships the hedge it would have shipped anyway); it just
# forgoes the gain.
#
# Shrinking the group recovers coverage and is still not worth it -- measured, not assumed:
#
#   group   mix     builds   gain|built   expected
#      42  83/16      3/5        14.27       8.56   (3/6 incl. verification -> 7.13)
#      38  84/15      3/5        15.28       9.17
#      34  86/13      4/5         9.60       7.68
#      30  88/12      4/5         7.72       6.17
#      25  90/10      5/5         5.59       5.59
#
# Group 25 builds on every contract and still expects less than group 42 building on three of
# five: score falls faster than availability rises. Note this resolves *opposite* to the group
# 75/50 question on HUDEP-2, where the score difference was noise and availability decided it --
# there the spread was +-1, here it is 14.27 -> 5.59. Group 38's +1.00 over 42 is inside noise
# (SE ~2 on each arm), so the shared group_size 42 is kept rather than adding a per-cell override
# the data does not support.
#
# The declines are stage-5 cell coverage (7 of 8 cells), which is contract site geometry -- the
# same contracts fail every time, and accessibility does not predict it. A declining CD34+ round
# is also the slowest prepare in the system: all-cut spends ~287s failing, then the hedge takes
# ~191s, 478s total -- inside PREPARE_BUDGET_S, but past the p10 lead time, so those rounds lean
# on _rows_for_task's wait. CD34+ clamps hardest of
# all four (energy saturates at d200 *and* d400, so cut_p sits at 0.990/0.960 across the whole
# band), which is why the band transfers; the coverage limit is unrelated to that.

# Retried once when the primary config declines, for cell types where a wider config recovers
# something the first one could not.
#
# K562 used to be here as d400/mf22, back when the primary was d200/mf22 and the retry was the
# distance widening. Now that the primary *is* d400/mf100 that retry would be strictly narrower:
# mf22 admits only guides failing <= 22 of 900, a subset of mf100's <= 100, at the same distance,
# and `bank_keep` sorts by fail count ascending so the 300k kept out of ~1.9M candidates contains
# every one of them. min_union_group over a subset cannot beat the same search over its superset,
# so the retry could only ever burn a second ~350s build to arrive somewhere no better. It was
# removed rather than rewritten: the obvious wider retry, d800/mf100, is ~700s on top of a ~350s
# primary and does not fit PREPARE_BUDGET_S.
RETRY_CONFIG: dict[str, dict] = {}


def config_for(cell_type: str, cfg: AllCutConfig | None = None) -> AllCutConfig | None:
    """The tuned config for a cell type, or None if the method has not been measured on it."""
    import dataclasses

    overrides = CELL_CONFIG.get(cell_type)
    if overrides is None:
        return None
    return dataclasses.replace(cfg or AllCutConfig(), **overrides)


def build_for_cell(contract: dict, reference: dict, cell_types: dict,
                   budget_s: float | None = None,
                   reserve_s: float = 0.0) -> tuple[list[dict] | None, dict]:
    """Build for whichever cell type this contract names, retrying once where that helps.

    ``reserve_s`` is budget the caller needs left over for its own fallback. The retry is skipped
    when taking it would spend that reserve: two failed attempts once left only 113 s of a 150 s
    seed-agnostic minimum, dropping the task past the hedge (~41) to the ordinary construction
    (~29). Recovering one contract in nine is not worth risking a 12-point regression on another.

    Returns ``(None, meta)`` for a cell type with no measured config, so the caller falls back
    rather than running an untuned variant.
    """
    import dataclasses

    cell = contract.get("cell_type")
    primary = config_for(cell)
    if primary is None:
        return None, {"reason": f"no measured all-cut config for {cell}"}
    started = time.monotonic()
    rows, meta = build_submission(contract, reference, cell_types, cfg=primary,
                                  budget_s=budget_s)
    if rows:
        return rows, meta
    retry = RETRY_CONFIG.get(cell)
    if retry is None:
        return None, meta
    left = None if budget_s is None else budget_s - (time.monotonic() - started)
    # The retry costs about as much again as the first attempt; only take it if what remains
    # afterwards still covers the caller's fallback.
    if left is not None and left - 120.0 < reserve_s:
        meta["reason"] = (f"{meta.get('reason')}; skipped the retry to keep {reserve_s:.0f}s "
                          f"for the fallback (only {left:.0f}s left)")
        return None, meta
    logger.info("all-cut: %s declined (%s); retrying at d%d",
                cell, meta.get("reason"), retry["max_distance"])
    rows, retry_meta = build_submission(contract, reference, cell_types,
                                        cfg=dataclasses.replace(AllCutConfig(), **retry),
                                        budget_s=left)
    retry_meta["retried"] = True
    retry_meta.setdefault("first_reason", meta.get("reason"))
    return rows, retry_meta


@dataclass(frozen=True)
class AllCutConfig:
    """Every value here was measured; see the module docstring for the sweep it came from."""

    group_size: int = 42               # 83/17 cas mix; peak clean payout (69.17) and expected (28.61)
    start_seed: int = 100
    end_seed: int = 999
    # Cas12a bank: the narrow band keeps the group's own rows near the gc_score peak. A 0-100% band
    # reaches more clean seeds (37.0% vs 29.2% at group 18-20) but pays less on each, and scored the
    # same within noise once the sampling bias was corrected.
    cas12a_gc: tuple[float, float] = (0.40, 0.60)
    cas12a_max_fail: int = 100         # the whole distribution is 62-200; 200 banked 12M unusable rows
    # Cas9: 0.40-0.80 spans everything that can be strict. Below 0.40 strict is ~1e-9 on HEK293, and
    # cut_p clamps at 0.99 for gc >= 0.80, so a wider band adds guides with no better cut_p (it was
    # measured at +2 to +4 candidates -- the excluded targets are low-GC, 117 sites below gc 0.25
    # against 3 above 0.75).
    cas9_gc: tuple[float, float] = (0.40, 0.80)
    # 400, not 1500: every archived contract carries a unique mutation set (40 of 40), and the
    # bank's fail sets are mutation-specific because experiment_seed hashes the mutation string —
    # so a cached bank never hits and it must be built inside the task. Distance is what makes that
    # affordable and it costs nothing: measured group-42 clean was 132 at d400 (120s bank) against
    # 131 at d1500 (~600s). d700 gave 133 in 181s, also no gain.
    max_distance: int = 400
    variants: int = 44000              # the GC band binds first: ~13k/target, well under this cap
    score_cap: int = 2500              # exact stage-2 scoring is only worth it on the ranked head
    bank_keep: int = 300_000
    pool_target: int = 8               # stop the Cas9 scan at this multiple of the rows needed

    @property
    def seeds(self) -> np.ndarray:
        return np.arange(self.start_seed, self.end_seed + 1, dtype=np.int64)


def bank_key(contract: dict, cell_types: dict, cfg: AllCutConfig) -> str:
    """Everything the bank depends on. Deliberately excludes the seed: the bank is seed-independent
    (it is defined *over* the whole window), so one bank serves every task on the same contract
    shape. A changed mutation set or accessibility invalidates it."""
    cell = contract.get("cell_type")
    payload = json.dumps({
        "cell": cell,
        "acc": cell_types.get(cell, {}).get("accessibility"),
        "mutations": sorted(contract.get("active_mutations") or []),
        "regions": contract.get("mutation_regions") or {},
        "gc": cfg.cas12a_gc, "d": cfg.max_distance, "v": cfg.variants,
        "mf": cfg.cas12a_max_fail, "w": [cfg.start_seed, cfg.end_seed],
    }, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _params_fn(site, distance: int, accessibility: float, region_offset: float):
    """(gc, energy, cut_p) per guide, cached on GC count — cut_p depends on the guide only through it."""
    cache: dict[int, tuple[float, float, float]] = {}

    def params_of(guide: str):
        count = sum(b in "GC" for b in guide)
        got = cache.get(count)
        if got is None:
            gc = count / site.length
            energy = max(0.0, min(1.0, accessibility * (
                1.8 * gc + 0.6 * math.exp(-distance / 1500) + region_offset)))
            got = cache[count] = (gc, energy, stage3.cut_probability(site.cas, energy))
        return got

    return params_of


def build_bank(contract: dict, reference: dict, cell_types: dict, ctx, sites,
               cfg: AllCutConfig, deadline: float | None = None) -> list[dict]:
    """Scan the cut-only Cas12a bank over the whole window. ~120s at d400.

    Returns [] if the deadline passes before every target is scanned: a partial bank would silently
    change the construction rather than fail, and the caller can still fall back.
    """
    cell = contract.get("cell_type")
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == "Cas12a" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    logger.info("all-cut: building the Cas12a bank over %d targets (this is the slow half)",
                len(jobs))
    bank, started = [], time.monotonic()
    for index, (site_index, mutation) in enumerate(jobs, 1):
        site = sites[site_index]
        distance = abs(site.start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas12a_gc[0], cfg.cas12a_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = _params_fn(site, distance, accessibility, offset)
        survivors = MT.screen_guides_rule_gpu(guides, cfg.seeds, mutation, site.cas, site.start,
                                              site.strand, params_of, "cut", cfg.cas12a_max_fail)
        for guide, fails in survivors.items():
            bank.append({"guide": guide, "mutation": mutation, "cas_system": site.cas,
                         "strand": site.strand, "start": site.start, "length": site.length,
                         "fails": fails.astype(np.int16)})
        if deadline is not None and time.monotonic() > deadline:
            logger.info("all-cut: bank scan out of budget at %d/%d targets", index, len(jobs))
            return []
        if index % 150 == 0:
            MT.free_gpu_memory()
            logger.info("  %d/%d targets | %d candidates | %.0fs",
                        index, len(jobs), len(bank), time.monotonic() - started)
    MT.free_gpu_memory()
    counts = np.asarray([len(b["fails"]) for b in bank])
    order = np.argsort(counts)[:cfg.bank_keep]
    return [bank[int(i)] for i in order]


def save_bank(path: str, bank: list[dict]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    width = max(len(b["fails"]) for b in bank)
    matrix = np.full((len(bank), width), -1, dtype=np.int16)
    for i, rec in enumerate(bank):
        matrix[i, :len(rec["fails"])] = rec["fails"]
    np.savez_compressed(
        path, fails=matrix,
        guide=np.array([b["guide"] for b in bank]),
        mutation=np.array([b["mutation"] for b in bank]),
        strand=np.array([b["strand"] for b in bank]),
        start=np.array([b["start"] for b in bank], dtype=np.int64),
        length=np.array([b["length"] for b in bank], dtype=np.int16))


def load_bank(path: str, limit: int = 60_000) -> list[dict]:
    """Read a cached bank. Every array is materialised once: ``np.load`` is lazy and each ``d[key]``
    re-decompresses the whole array, so touching them inside the loop is O(rows) decompressions."""
    data = np.load(path, allow_pickle=False)
    fails = np.asarray(data["fails"])
    guide = np.asarray(data["guide"])
    mutation = np.asarray(data["mutation"])
    strand = np.asarray(data["strand"])
    start = np.asarray(data["start"])
    # Banks written before `length` was recorded default to 23: every Cas12a site here is a 23-mer
    # (TTTV PAM), and the field only feeds target_alignment_end.
    length = (np.asarray(data["length"]) if "length" in data.files
              else np.full(fails.shape[0], 23, dtype=np.int16))
    counts = (fails >= 0).sum(axis=1)
    order = np.argsort(counts)[:limit]
    return [{"guide": str(guide[i]), "mutation": str(mutation[i]), "cas_system": "Cas12a",
             "strand": str(strand[i]), "start": int(start[i]), "length": int(length[i]),
             "fails": fails[i][:counts[i]]} for i in order]


def scan_cas9(clean: np.ndarray, contract: dict, cell_types: dict, ctx, sites,
              cfg: AllCutConfig, want: int, deadline: float | None = None) -> list[dict]:
    """Cas9 guides strict over ``clean``, nearest targets first.

    Ordering by distance is what makes the early exit lossless: the assembly ranks candidates by
    ``(distance, |gc - 0.50|)``, so the rows it keeps come from the nearest targets regardless. A
    full scan of all 1,245 targets took 248s and returned 71,256 candidates to fill 208 rows.
    """
    cell = contract.get("cell_type")
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    jobs = [(i, m, abs(sites[i].start - ctx.mutation_map[m]))
            for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda job: job[2])
    found: list[dict] = []
    for site_index, mutation, distance in jobs:
        if deadline is not None and time.monotonic() > deadline:
            logger.info("all-cut: Cas9 scan stopped on the deadline with %d candidates",
                        len(found))
            break
        site = sites[site_index]
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = _params_fn(site, distance, accessibility, offset)
        for guide in MT.screen_guides_rule_gpu(guides, clean, mutation, "Cas9", site.start,
                                               site.strand, params_of, "cut", 0):
            gc, _energy, _cut_p = params_of(guide)
            found.append({"guide": guide, "mutation": mutation, "cas_system": "Cas9",
                          "strand": site.strand, "start": site.start, "length": site.length,
                          "gc": gc, "distance": distance})
        # Enough, and every cell represented: an empty (mutation, strand) cell would zero one of
        # stage 5's six ratios and cost a ~0.03x multiplier on the whole score.
        if len(found) >= want * cfg.pool_target and len({(f["mutation"], f["strand"])
                                                          for f in found}) == 4:
            break
    return found


def assemble(group: list[dict], cas9: list[dict], contract: dict, ctx, cfg: AllCutConfig,
             n_rows: int) -> list[dict]:
    """Cas12a group + mutation-weighted Cas9, mirroring ``seed_agnostic._assemble``'s apportionment."""
    need = n_rows - len(group)
    cas9.sort(key=lambda r: (r["distance"], abs(r["gc"] - 0.50)))
    head = cas9[:cfg.score_cap]
    by_site = {(s.start, s.strand, s.cas, s.length): s for s in
               G.enumerate_sites(ctx, 3000, (20, 23))}
    for rec in head:
        site = by_site.get((rec["start"], rec["strand"], "Cas9", rec["length"]))
        entry = None if site is None else G.build_valid_entry(
            G.make_experiment(site, rec["guide"], rec["mutation"], ctx, "score"), ctx)
        rec["weighted_score"] = entry["stage2"]["weighted_score"] if entry else -1.0
    head = [r for r in head if r["weighted_score"] >= 0]

    by_cell: dict[tuple, list[dict]] = {}
    for rec in head:
        by_cell.setdefault((rec["mutation"], rec["strand"]), []).append(rec)
    for rows in by_cell.values():
        rows.sort(key=lambda r: -r["weighted_score"])
    weights = contract.get("mutation_weights", {})
    shares = {m: max(weights.get(m, 1.0), 1e-9) ** 1.25 for m in ctx.mutations}
    total = sum(shares.values())
    chosen, used = list(group), set()
    for mutation in ctx.mutations:
        per = round(need * shares[mutation] / total)
        for strand, take in (("+", per // 2), ("-", per - per // 2)):
            for rec in by_cell.get((mutation, strand), [])[:take]:
                chosen.append(rec)
                used.add(id(rec))
    # Top up from whatever is left. The per-cell quota can run dry, and a short submission loses
    # term 1 linearly — that silently produced a 137-row build during development.
    if len(chosen) < n_rows:
        spare = sorted((r for r in head if id(r) not in used),
                       key=lambda r: -r["weighted_score"])
        chosen.extend(spare[:n_rows - len(chosen)])

    rows, seen = [], set()
    for i, rec in enumerate(chosen[:n_rows]):
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": rec["guide"],
                     "target_alignment_start": rec["start"],
                     "target_alignment_end": rec["start"] + rec["length"],
                     "strand": rec["strand"], "mutation": rec["mutation"],
                     "cas_system": rec["cas_system"],
                     "cell_type": contract.get("cell_type")})
    return rows


def build_submission(contract: dict, reference: dict, cell_types: dict,
                     cfg: AllCutConfig | None = None, budget_s: float | None = None,
                     allow_bank_build: bool = True) -> tuple[list[dict] | None, dict]:
    """The all-cut submission, or ``(None, meta)`` when the caller should fall back.

    The bank is built inside the task (~120s at the default d400) because contract shapes do not
    repeat — every archived task has a distinct mutation set, so a disk cache never hits. The cache
    is still written and read: it costs nothing and covers a re-broadcast of the same task.

    ``budget_s`` is enforced across both halves. If the bank cannot finish in time the build
    declines rather than overrunning, and the caller falls back.
    """
    cfg = cfg or AllCutConfig()
    started = time.monotonic()
    deadline = None if budget_s is None else started + budget_s
    meta: dict = {"method": "hek293-all-cut", "group_size": cfg.group_size}

    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    path = os.path.join(BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    meta["bank_path"] = path

    if not os.path.exists(path):
        if not allow_bank_build:
            meta["reason"] = "no cached Cas12a bank and bank building is disabled"
            return None, meta
        # Reserve the Cas9 half plus assembly (~10s measured) out of the budget.
        bank_deadline = None if deadline is None else deadline - 20.0
        bank = build_bank(contract, reference, cell_types, ctx, sites, cfg, bank_deadline)
        if not bank:
            meta["reason"] = "Cas12a bank scan ran out of budget"
            return None, meta
        save_bank(path, bank)
    records = load_bank(path)
    meta["bank"] = len(records)

    selector = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed)
    index, union = selector.best(cfg.group_size, restarts=12)
    group = [records[i] for i in index]
    bad = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad), dtype=np.int64)
    window = cfg.end_seed - cfg.start_seed + 1
    meta.update(union=len(bad), clean=len(clean), clean_fraction=len(clean) / window)

    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    cas9 = scan_cas9(clean, contract, cell_types, ctx, sites, cfg, n_rows - cfg.group_size,
                     deadline)
    meta["cas9_pool"] = len(cas9)
    cells = len({(r["mutation"], r["strand"]) for r in cas9})
    if len(cas9) < n_rows - cfg.group_size or cells < 4:
        meta["reason"] = f"Cas9 pool {len(cas9)} over {cells} cells is short of {n_rows - cfg.group_size}"
        MT.free_gpu_memory()
        return None, meta

    rows = assemble(group, cas9, contract, ctx, cfg, n_rows)
    MT.free_gpu_memory()
    meta.update(rows=len(rows), elapsed_s=round(time.monotonic() - started, 1),
                cells=len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
                cas_mix=dict(Counter(r["cas_system"] for r in rows)))
    if len(rows) < n_rows or meta["cells"] < 8:
        meta["reason"] = f"assembled {len(rows)} rows over {meta['cells']} cells"
        return None, meta
    return rows, meta


def main() -> None:
    """Prebuild the Cas12a bank for a task, so the miner only pays the per-task half."""
    import argparse

    parser = argparse.ArgumentParser(description="prebuild an all-cut Cas12a bank")
    parser.add_argument("--task", default="test/task.json")
    parser.add_argument("--cell-types", default="test/cell_types.json")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(message)s")

    task = json.load(open(args.task))
    contract = dict(task["content"]["contract"])
    contract["cell_type"] = "HEK293"
    reference = task["content"]["hbb_reference"]
    cell_types = json.load(open(args.cell_types))
    rows, meta = build_submission(contract, reference, cell_types, allow_bank_build=True)
    print(json.dumps(meta, indent=2, default=str))


if __name__ == "__main__":
    main()
