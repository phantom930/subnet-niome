"""all_hdr.py — the all-HDR construction: every row repairs by HDR on a shared seed band.

The same conditional shape as :mod:`all_cut`, with the pinned outcome moved from *cut* to *HDR*:

1. **Cas12a** — min-union a group of ``group_size`` guides from a bank screened on the ``hdr`` rule
   over a narrow band (``hdr_range``, 100 seeds). Their failed-seed union defines the seeds this
   submission cannot save; the complement is the *clean band*.
2. **Cas9** — require HDR only over that clean band, not over the whole band. Same relaxation
   all-cut makes, one rule along.

On a seed in the clean band every row is HDR, so stage 4's three targets — ``is_cut``, ``is_hdr``
and ``indel_length`` — are all constant and ``consistency_factor`` reaches **exactly 1.0**, against
the ~0.20-0.23 all-cut gets. Verified against ``stage3.simulate`` itself, not just the GPU replica.

**This scores worse than all-cut on the mean, and is shipped deliberately for the tail.** Measured
head to head, five contracts, clean/dirty legs, scored through all five stages:

    cell          band       clean  clean leg  dirty leg  expected  all-cut    delta
    CD34+_HSPC    500-599       16     275.38      29.89     34.26    56.55   -22.29
    K562          700-799       16     254.30      27.22     31.25    61.12   -29.87
    HUDEP-2       800-899       15     191.31      20.07     22.92    41.67   -18.75

all-HDR returns 51-61% of all-cut's expected score. The construction is not at fault — the band is
only ~1.8% of the 100-999 seed space, so the 95% of rounds that miss it pay for the 5% that hit.
The bet being made is that ``SCORING_SYSTEM = "top"`` pays only the top 10 miners on a fixed curve,
so a build worth 111-194 on ~5% of rounds may out-*rank* a flat 56 that never spikes. That is an
argument from the payout structure and is **not** supported by any rank measurement here; if it is
wrong, this costs ~24 points a round across three quarters of task volume. ``Miner.ALL_HDR`` and
``ALL_HDR_CELL_TYPES`` revert it.

Two structural limits, both measured, that no parameter reaches:

* **The band cannot be widened.** Every seed added multiplies the conditional Cas9 requirement by
  P(HDR) ~ 0.57. A 16-seed band yields 484-870 Cas9 candidates against the ~208 needed; a 29-31
  seed band (group 6) yielded **zero** on all three cell types.
* **The band's position is free.** Seeds are independent draws, so 300-399, 500-599, 700-799 and
  800-899 all screen alike — candidate yields 14.0-14.9% and clean bands 15-16 across three cell
  types. Distinct ranges per cell type are a convention here, not a tuning result.

**HEK293 is included, at a larger group.** At accessibility 0.35 energy never clamps, so P(HDR)
falls to ~0.37-0.39 where the enumeration's mass sits: the bank collapses to 1,206 candidates from
2.9M guides (0.041% against ~14.5%) and the clean band is 6-10 seeds against the erythroid types'
15-16. None of that stops the construction — a *narrower* band is easier for the conditional Cas9
half to fill (``0.37**6`` against ``0.37**10``), and every group 20-100 builds 8/8 cells with the
band seed at ``consistency_factor`` exactly 1.0. What HEK293 needs is a bigger group, because the
cas mix is what stage 5 sees:

    group  cas12a/cas9  band  cells  spike cons  fidelity  weighted  spike final
       20       20/230    10    8/8       1.000     0.786       326          256
       40       40/210     8    8/8       1.000     0.854       330          282
       60       60/190     8    8/8       1.000     0.886       326          289
       80       80/170     7    8/8       1.000     0.924       318          294
      100      100/150     6    8/8       1.000     0.939       308          289

Group 80 is shipped: fidelity 0.924 is level with the leaders' 0.947 median on their placing HEK293
rounds (weighted 288, consistency 0.552), and its 7-seed band is hit ~17% more often than group
100's 6 for +0.015 fidelity given up. The weighted score drifts 326 -> 308 across the sweep because
the added Cas12a rows score below the Cas9 rows they displace; that is the whole cost, and 308 still
clears the leaders' 288.

An earlier claim here that HEK293 fidelity was capped near 0.78 was an artefact of only ever
building groups 8-24, where 20 Cas12a rows of 250 drive stage 5's cas-coverage entropy to ~0.65.
Balancing the mix is the same lever that lifted the other three cell types; nothing about
accessibility 0.35 bounds fidelity.
"""

from __future__ import annotations

import dataclasses
import math
import logging
import os
import time
from collections import Counter
from dataclasses import dataclass

import numpy as np

import genExp as G
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.all_cut import (cas9_cell_target, BANK_DIR, _params_fn, assemble, bank_key, load_bank,
                                           save_bank)
from niome_subnet.genomics.validation import stage3

logger = logging.getLogger(__name__)

HDR_BANK_DIR = "data/all_hdr"

# The pinned band per cell type. Position is free (see the module docstring); distinct ranges keep
# the on-disk banks from colliding and make a log line say which cell type it came from.
CELL_CONFIG: dict[str, dict] = {
    # CD34+_HSPC: width 150 at group 100, validated over 5 contracts (spread 1.50-2.71) at
    # E[share] 0.00236 against the shared width-100/group-80 default's 0.00227 (+3.8%). Carries
    # its own `main_max_fail` because `build_for_cell` only rescales the screen for a window
    # passed IN by the caller: 78 is the measured plateau entry at span 150 (the adaptive search
    # had to loosen past the linear 68, unlike K562/HUDEP-2, so this is cell-specific).
    "CD34+_HSPC": {"hdr_range": (500, 649), "main_max_fail": 78, "group_size": 100,
                   "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95)},
    # K562 stays at the shared 100-seed band. A 150-seed default was set earlier in this
    # session on a SINGLE contract (7306c626) and is withdrawn: validated across 5 contracts
    # spanning weight spread 1.50-2.54, width 100 prices at E[share] 0.00204 against width 150's
    # 0.00181 and width 75's 0.00181 — the single-contract grid picked a width that came LAST.
    "K562": {"hdr_range": (700, 799), "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95)},
    "HUDEP-2": {"hdr_range": (800, 899), "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95)},
    # HEK293: width 75 at group 100, the largest validated gain of the four cells — E[share]
    # 0.00228 over 5 contracts (spread 1.50-3.93) against the width-100/group-80 default's
    # 0.00205 (+10.8%). `main_max_fail` 45 is the measured plateau entry at span 75; the z-rule
    # predicted 34 and the adaptive search had to loosen twice, so HEK293 needs a looser screen
    # than its own P(HDR) implies.
    #
    # **group_size 100 is specific to width 75.** At widths 100/150/225 this cell measured group
    # 80 ahead (230.2/228.6/228.2 against 229.2/228.1/224.8), so a build that overrides the window
    # to a wider one should use 80. See the note on CELL_CONFIG below.
    "HEK293": {"hdr_range": (300, 374), "main_max_fail": 45, "group_size": 100,
               "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95)},
}


@dataclass(frozen=True)
class AllHdrConfig:
    """Every value here is measured; see the module docstring for the sweep it came from."""

    # 80/170. This is a *payout* optimum, and it deliberately overrides an earlier *score* optimum
    # — the two disagree, which is the whole reason this comment is long.
    #
    # group_size trades two things against each other, both measured:
    #   band width  (spike frequency)  42:14.6  60:13.7  80:13.0  100:12.1   erythroid means
    #                                  42: 8.3  60: 7.3  80: 7.2  100: 6.5   HEK293
    #     — monotone, 56/56 (contract, window) pairs, no exceptions.
    #   fidelity    (spike score)      42:0.89  60:0.93  80:0.96  100:0.976
    #     — the cas mix: 42/208 drives stage 5's cas-coverage entropy to ~0.65, 100/150 to ~1.0.
    # weighted moves slightly the *other* way (K562 258 at 42 vs 250 at 100), because a smaller
    # group means more Cas9 rows and Cas9 rows score better structurally. So 80 gains a little on
    # band and weighted and gives up only fidelity.
    #
    # Priced as expected payout over 9 disjoint hotkeys — 54 current-regime fields x 25 contracts,
    # each round scored at its k=1 case (one of three seeds in the band, the case that decides
    # placement):
    #
    #   cell            g42      g60      g80     g100    best
    #   K562         0.0079   0.0135   0.0137   0.0116     g80
    #   HUDEP-2      0.0104   0.0122   0.0126   0.0111     g80
    #   CD34+_HSPC   0.0106   0.0126   0.0138   0.0123     g80
    #   HEK293       0.0201   0.0200   0.0197   0.0173     g42 (+2%, inside noise)
    #   AGGREGATE    0.0490   0.0583   0.0598   0.0524     g80  (+14.2% over 100)
    #
    # **Two traps this measurement fell into; do not repeat either.**
    #  1. Scoring against a handful of single fields said group 42 wins by +32%. Those five fields
    #     were all above their cell type's median cutoff. On a hard field the k=1 hit does not
    #     place, so payout rides on the k=2 term, which scales as band**2 and favours wide bands.
    #     Sample many fields, not one.
    #  2. Sampling fields across all backend history says the opposite again, because the subnet
    #     ran at 29-44 miners/task historically against 248 today, and those low cutoffs make every
    #     config look like it places. Only current-regime fields are decision-relevant.
    #
    # The superseded note: group 42 -> 100 was taken on spike-round *score*, 117.1 -> 121.8 (+4.7),
    # which is real and still true. It is the wrong objective — SCORE_DISTRIBUTION is a step
    # function, so a score gain that crosses no rank threshold pays nothing while the narrower band
    # costs frequency on every round.
    group_size: int = 80
    hdr_range: tuple[int, int] = (500, 599)
    # An explicit, possibly NON-CONTIGUOUS band space. When set it replaces `hdr_range` as the seed
    # set the Cas12a group min-unions over, so several disjoint width-100 windows can be joined into
    # one band space. `hdr_range` is still kept in step (min..max) for logging and for anything that
    # only wants a coarse span. `bank_key` folds the set in, so a joined band never shares a bank
    # with the contiguous range that spans it.
    seed_list: tuple[int, ...] | None = None
    # Fails tolerated in the band when banking a Cas12a candidate. 45 of 100 is deliberately loose:
    # the min-union step is what produces the clean band, and a tighter screen shrinks the pool it
    # selects from without improving the group (the binomial tail is far steeper than the gain).
    main_max_fail: int = 45
    max_distance: int = 400
    cas12a_gc: tuple[float, float] = (0.40, 0.95)
    cas9_gc: tuple[float, float] = (0.40, 0.95)
    variants: int = 44000
    score_cap: int = 2500
    bank_keep: int = 300_000
    pool_target: int = 8
    restarts: int = 12
    per_cell_min: int = 2
    # Mutation apportionment. ``light_cell_rows = 6`` is **on**, on a 210-config measurement that
    # supersedes the earlier sweep recorded below; ``light_group_cells`` stays off (group caps were
    # harmful at every setting).
    #   light_group_cells — cap per light (mutation, Cas12a, strand) cell in the min-union group
    #   light_cell_rows   — rows per light (mutation, Cas9, strand) cell in assemble
    #   weight_exponent   — the exponent in assemble's smooth mutation_weight apportionment
    #
    # The motivation was real: the min-union group is blind to mutation_weight and lands near
    # 50/50, so on 9ed335da we shipped 158 of 250 rows on the heavy mutation while the miners
    # taking ranks 8-11 — identical consistency 0.405, *worse* fidelity 0.940 — held ~239, and that
    # 39-point total_weighted_score gap was the whole distance between rank 8 and our rank 12 on
    # the fleet's only spike to date.
    #
    # **The knobs close the gap and it does not help.** Sweep over 4 cell types x 5 settings,
    # ranked on each task's own real field, scoring the k=1 round (one of three seeds in the band —
    # the case that actually decides payout):
    #
    #   cell / task          setting     heavy    fid  weighted  k=1 final   d      rank
    #   K562 9ed335da        None/None     158  0.973     223.9       87.2   -        11
    #                        None/6        194  0.921     242.3       89.3  +2.1      11
    #                        2/4           238  0.764     263.7       80.6  -6.6      11
    #   K562 37737cb7        None/6        185  0.933     214.2       79.9  +1.9      10
    #   CD34+_HSPC           None/6        188  0.932     293.7      109.5  +4.7      11
    #   HUDEP-2              None/6        191  0.928     243.4       90.4  +2.2      10
    #   HEK293               None/6        208  0.871     253.6       88.3  -0.4      11
    #
    # weighted does climb to 263-326 as intended, but stage 5's mutation-coverage entropy falls
    # faster: every backend contract carries exactly 2 mutations, so capping the light one drives
    # that term to its floor. **No rank moved in any of the 20 builds.** Group caps are strictly
    # harmful everywhere; the best setting is light_cell_rows=6 with no group cap, worth +2 to +4.7
    # on three cell types and -0.4 on HEK293 (already skewed to 182/250 at its group 80, because
    # its 170-row Cas9 half is what the exponent governs).
    #
    # **Superseded.** That first sweep held group size and contract sample fixed and priced the k=1
    # final rather than expected payout. Re-measured over 210 configs — 12 contracts across all four
    # cell types x max_distance {400,150,100} x group {42,60,80} x light {None,6} — and priced as
    # E[pay] against each cell type's own current-regime fields, with builds that DECLINE charged as
    # zero (a decline drops to all-cut, which never places):
    #
    #   maxd  grp  light  built  E|built  E[eff]   vs shipped
    #    400   80      6     12   0.0242  0.0242      +15.7%   <- shipped
    #    150   80      6     12   0.0237  0.0237      +13.4%
    #    100   80      6     11   0.0257  0.0235      +12.7%   declines on HEK293/16710fdc
    #    100   80   None     11   0.0239  0.0219       +4.7%
    #    400   80   None     12   0.0209  0.0209          --
    #
    # Two traps in that table. max_distance 100 has the best score *among builds* and is worse
    # overall, because it fails to build on 1 of 12 contracts — always charge declines. And the
    # effect is spread-dependent: paired over 105 configs, light=6 is -0.0008 E[pay] at weight
    # spread <= 2.0 but +0.0026 at 2.0-2.6 and +0.0023 above 2.6. The original rejection was
    # measured on K562 at spread 1.85, the one bucket where it genuinely loses. Gating on spread
    # > 2.0 measured 0.0240, no better than applying it unconditionally, so it is not gated.
    #
    # The cost is real: fidelity falls 0.948 -> 0.889, giving up the one term where this fleet led
    # the field. E[pay] says take it anyway. Solving their row backwards, ranks 8-11 reach
    # weighted 262.8 at fidelity 0.940 with only ~180 heavy rows, which needs base structural
    # ~1.005 against our 0.894 — they have near-perfect gc_score and dist_score, not a heavier
    # skew. dist_score is the lever that does not trade: stage 5 measures mutation/cas/strand/joint
    # coverage plus k-mer and guide diversity, and neither distance nor GC appears in it. With
    # base at 1.0, light_cell_rows=6 then reaches 99.2 on 9ed335da — rank 10 — so these knobs are
    # worth revisiting in that order, and only in that order.
    light_group_cells: int | None = None
    # "auto" resolves per contract from the mutation-weight spread — see `resolve_light`. An int
    # or None pins it, which is what every research script and `all_cut` still pass.
    light_cell_rows: "int | str | None" = "auto"
    weight_exponent: float = 1.25

    @property
    def cas12a_max_fail(self) -> int:
        """Alias ``main_max_fail`` under the name ``all_cut.bank_key`` reads.

        The bank key is borrowed rather than duplicated so the two builders cannot drift on what
        invalidates a cache. It reads ``cfg.cas12a_max_fail``; this config calls the same quantity
        ``main_max_fail`` because that is the parameter's name in search_hdr.py and in the request
        it came from. The alias keeps both names honest instead of renaming one to suit the other.
        """
        return self.main_max_fail

    def band_seeds(self) -> "np.ndarray":
        """The seeds the band is searched over — the explicit set if given, else the range."""
        if self.seed_list:
            return np.asarray(sorted(set(int(x) for x in self.seed_list)), dtype=np.int64)
        return np.arange(self.hdr_range[0], self.hdr_range[1] + 1, dtype=np.int64)

    @property
    def start_seed(self) -> int:
        return self.hdr_range[0]

    @property
    def end_seed(self) -> int:
        return self.hdr_range[1]


# `variants` for a band window wider than the 100 seeds CELL_CONFIG is calibrated at. The cap
# exists because the Cas12a HDR bank scan grows with the span (cost ~`variants x seeds`), so a wide
# window is the one place the tuned per-cell 44000 might not be affordable.
#
# Measured at width 300 on K562 (hdr_variants.py, cold scans, fleet idle; every arm built 250 rows
# at 8/8 cells). Band is far less sensitive to `variants` than to window width -- an 8.8x cut in
# variants costs only 3 band seeds:
#
#   variants |  bank  | scan s | band | build s | fleet coverage (7 hotkeys)
#      5000  | 12619  |    12  |   9  |    27   | 63 seeds  19.6%
#     11000  | 25104  |    21  |  10  |    48   | 70 seeds  21.6%
#     22000  | 43540  |    33  |  11  |    75   | 77 seeds  23.5%
#     44000  | 62144  |    45  |  12  |   105   | 84 seeds  25.5%
#
# **Do not budget this by summing per-build times.** That model said 7 x 105s + h0's 354s all-cut
# = 1089s against the 900s prefetch budget, and it is wrong by ~3x: the builds are substantially
# CPU-bound (variant enumeration, assembly, local scoring) on a 15-core box, so sibling processes
# overlap rather than queue on the GPU. Measured wall clock with all seven concurrent: **~100s** at
# variants 11000 (17:46:45 -> 17:48:25) and **~200s** at 44000 five-way (16:58:34 -> 17:01:55).
# Both fit the budget with room to spare, so the cap is set to the tuned value and buys the band
# back. The real constraint is h0's all-cut, which is heavier and more GPU-bound than any band
# scan; if a round shows it declining alongside seven band builds, step this down to 22000 (band
# 11) before giving up any more coverage.
WIDE_WINDOW_VARIANTS = 44000

# Mean Cas12a fail count over a 100-seed band window, measured per cell with the cap disabled
# (hek_mf.py). Only needed to scale `main_max_fail` to a wider window correctly -- a cell absent
# here falls back to linear scaling, which is what shipped before and is close enough wherever
# `main_max_fail` sits near the mean.
#
# Linear scaling is WRONG for a cell whose P(HDR) is far from the threshold rate, because it scales
# the offset from the mean linearly when that offset should grow as sqrt(span): the binomial
# narrows as the window widens, so the same threshold *rate* is a much deeper tail cut at 300 than
# at 100. On HEK293 (P(HDR) 0.358 measured, against ~0.49 on K562) the tuned mf 48 sits at
# z -3.32 at span 100, and linear scaling to 144 lands at **z -5.66** -- only 53 of 60000 bank
# guides qualified against a group of 80, so every band hotkey declined and the fleet shipped flat
# builds on every HEK293 round. Holding z instead gives 164, which builds: bank 30009, band 8.
#
# Measured at span 300 on HEK293 (414dab89): band is FLAT at 8 for mf 164/184/192/204 while build
# time runs 57/123/190/255s, because `bank_keep` caps the bank and the min-union takes the best
# guides either way. So the z-matched value is both the cheapest that fills the group and the point
# the band saturates -- there is nothing to buy above it.
# K562 is deliberately ABSENT and stays on linear scaling. Its P(HDR) (~0.49) sits close enough to
# the threshold rate that linear works, and measured head to head at span 300 on bba85ff0 the
# z-matched value is slightly WORSE: mf 135 (linear) -> band 12, mf 143 (z-matched) -> band 11,
# both with bank 60000 / 8-of-8 cells. So z-matching is a fix for a cell where linear scaling
# fails outright, not a general improvement -- measure before adding a cell here.
MEAN_FAIL_100 = {"HEK293": 64.2}


def _scaled_max_fail(cell_type: str, mf_100: int, span: int) -> int:
    """``main_max_fail`` for a band of ``span`` seeds, holding the z-score of the tuned value.

    mf = mean*r + z*sd*sqrt(r) with r = span/100, which rearranges to the form below so only the
    mean fail count at span 100 is needed. Falls back to linear when that is unmeasured.
    """
    r = span / 100.0
    mean_100 = MEAN_FAIL_100.get(cell_type)
    if not mean_100:
        return max(1, round(mf_100 * r))
    return max(1, round(mf_100 * r + (mean_100 - mf_100) * (r - math.sqrt(r))))


# `light_cell_rows` is not one number: the right value depends on the CONTRACT, because the knob
# trades `total_weighted_score` against stage 5's mutation-coverage entropy and every backend
# contract carries exactly 2 mutations at a different weight ratio. Measured three ways:
#
#   * 26 distinct K562 contracts, joined-150 band, group 80 (`light_vs_weight.py`): regressing
#     (L12 - L6) on the weight spread gives slope -4.69 per unit spread, r = -0.810, t = -6.76,
#     and a crossover at spread **2.06**. Below it L12 wins 15 of 18; above it L6 wins 5 of 8.
#   * 20 validation contracts across all four cell types: L6 wins every contract at spread >= 2.18,
#     L12 or L25 wins every contract below ~1.9.
#   * The shipped constant 6 is therefore right only for the high-spread half of the distribution;
#     K562's spreads run 1.50-2.54 with a median of 1.85, i.e. mostly BELOW the crossover.
#
# Worth +0.35% on `weighted x fidelity` over a fixed 6 (paired +0.83 +/- 0.22, t = 3.80, winning 16
# of 26 contracts and losing 2). Small, but it is the best-replicated result of the sweep and it
# costs nothing at build time.
#
# Two tiers, not three: L25 beat L12 at spread 1.50 on HUDEP-2/CD34+/HEK293 but LOST to it on
# K562, so the sub-1.6 tier is unresolved and 12 is the safe choice there.
LIGHT_CROSSOVER = 2.06


def resolve_light(value, contract: dict):
    """Turn a ``light_cell_rows`` of ``"auto"`` into a number using the contract's weight spread."""
    if value != "auto":
        return value
    weights = contract.get("mutation_weights") or {}
    if len(weights) != 2 or min(weights.values()) <= 0:
        return 6
    spread = max(weights.values()) / min(weights.values())
    return 12 if spread < LIGHT_CROSSOVER else 6


def config_for(cell_type: str, cfg: AllHdrConfig | None = None) -> AllHdrConfig | None:
    """The tuned config for a cell type, or None where all-HDR has not been measured."""
    overrides = CELL_CONFIG.get(cell_type)
    if overrides is None:
        return None
    return dataclasses.replace(cfg or AllHdrConfig(), **overrides)


def _range_label(seeds) -> str:
    """Compact "100-199,300-399" label for a possibly non-contiguous seed set."""
    xs = sorted(set(int(x) for x in seeds))
    if not xs:
        return "-"
    out, lo, prev = [], xs[0], xs[0]
    for x in xs[1:]:
        if x != prev + 1:
            out.append((lo, prev))
            lo = x
        prev = x
    out.append((lo, prev))
    return ",".join(f"{a}-{b}" for a, b in out)


def build_bank(contract: dict, reference: dict, cell_types: dict, ctx, sites,
               cfg: AllHdrConfig, deadline: float | None = None) -> list[dict]:
    """Cas12a guides reaching HDR on all but ``main_max_fail`` seeds of the band. ~35-40s.

    Far cheaper than all-cut's bank despite the same target count, because the band is 100 seeds
    rather than 900. Returns [] on the deadline: a partial bank would silently change the
    construction rather than fail.
    """
    cell = contract.get("cell_type")
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    seeds = cfg.band_seeds()
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == "Cas12a" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    logger.info("all-hdr: banking Cas12a over %d targets, band %d-%d", len(jobs),
                cfg.start_seed, cfg.end_seed)
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
        for guide, fails in MT.screen_guides_rule_gpu(
                guides, seeds, mutation, site.cas, site.start, site.strand, params_of,
                "hdr", cfg.main_max_fail).items():
            bank.append({"guide": guide, "mutation": mutation, "cas_system": site.cas,
                         "strand": site.strand, "start": site.start, "length": site.length,
                         "fails": fails.astype(np.int16)})
        if deadline is not None and time.monotonic() > deadline:
            logger.info("all-hdr: bank scan out of budget at %d/%d targets", index, len(jobs))
            return []
        if index % 150 == 0:
            MT.free_gpu_memory()
            logger.info("  %d/%d targets | %d candidates | %.0fs", index, len(jobs), len(bank),
                        time.monotonic() - started)
    MT.free_gpu_memory()
    counts = np.asarray([len(b["fails"]) for b in bank])
    order = np.argsort(counts)[:cfg.bank_keep]
    return [bank[int(i)] for i in order]


def scan_cas9(clean: np.ndarray, contract: dict, cell_types: dict, ctx, sites,
              cfg: AllHdrConfig, want: int, deadline: float | None = None) -> list[dict]:
    """Cas9 guides reaching HDR on *every* seed of ``clean``, nearest targets first.

    The relaxation that makes the construction possible: HDR over 15-16 seeds is ~1.1e-4 per guide,
    against ~1e-8 over the full band. Ordering by distance keeps the early exit lossless, since the
    assembly ranks by ``(distance, |gc - 0.50|)`` anyway.
    """
    cell = contract.get("cell_type")
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    jobs = [(i, m, abs(sites[i].start - ctx.mutation_map[m]))
            for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(sites[i].start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda job: job[2])
    found: list[dict] = []
    per_cell: dict[tuple, int] = {}
    cell_target = cas9_cell_target(contract, ctx, cfg, want)
    for site_index, mutation, distance in jobs:
        if deadline is not None and time.monotonic() > deadline:
            logger.info("all-hdr: Cas9 scan stopped on the deadline with %d candidates", len(found))
            break
        site = sites[site_index]
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = _params_fn(site, distance, accessibility, offset)
        for guide in MT.screen_guides_rule_gpu(guides, clean, mutation, "Cas9", site.start,
                                               site.strand, params_of, "hdr", 0):
            gc, _energy, _cut_p = params_of(guide)
            found.append({"guide": guide, "mutation": mutation, "cas_system": "Cas9",
                          "strand": site.strand, "start": site.start, "length": site.length,
                          "gc": gc, "distance": distance})
            key = (mutation, site.strand)
            per_cell[key] = per_cell.get(key, 0) + 1
        # Every cell filled to the quota `assemble` will ask for, not merely non-empty: see
        # all_cut.cas9_cell_target for the HEK293 measurement that made this necessary.
        if (len(found) >= want * cfg.pool_target and len(per_cell) == 4
                and min(per_cell.values()) >= cell_target):
            break
    return found


def _group_caps(contract: dict, ctx, cfg: AllHdrConfig) -> dict | None:
    """Per-cell ceilings for the min-union group: light mutations capped, the heaviest left free.

    Returns None when unset, which is the unconstrained selection the group had before. FastGreedy
    raises any cap below that cell's floor and ignores the whole cap set if it would make the group
    size unreachable, so this can narrow the mutation mix but never empty a stage-5 cell.
    """
    if cfg.light_group_cells is None:
        return None
    weights = contract.get("mutation_weights", {})
    heavy = max(ctx.mutations, key=lambda m: weights.get(m, 1.0))
    return {(m, "Cas12a", strand): cfg.light_group_cells
            for m in ctx.mutations if m != heavy for strand in ("+", "-")}


def build_submission(contract: dict, reference: dict, cell_types: dict,
                     cfg: AllHdrConfig | None = None,
                     budget_s: float | None = None) -> tuple[list[dict] | None, dict]:
    """The all-HDR submission, or ``(None, meta)`` when the caller should fall back."""
    cfg = cfg or AllHdrConfig()
    # Resolve "auto" BEFORE anything reads the field. It is not only `assemble` that does:
    # `scan_cas9` calls `all_cut.cas9_cell_target`, which does arithmetic on light_cell_rows to
    # derive the per-cell Cas9 floor. Resolving late raised TypeError there and would have
    # crash-looped every hotkey on its next restart.
    cfg = dataclasses.replace(cfg, light_cell_rows=resolve_light(cfg.light_cell_rows, contract))
    started = time.monotonic()
    deadline = None if budget_s is None else started + budget_s
    meta: dict = {"method": "all-hdr", "group_size": cfg.group_size,
                  "light_cell_rows": cfg.light_cell_rows,
                  "band": (f"{cfg.start_seed}-{cfg.end_seed}" if not cfg.seed_list
                           else f"joined {len(cfg.seed_list)} seeds "
                                f"{_range_label(cfg.seed_list)}")}

    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    os.makedirs(HDR_BANK_DIR, exist_ok=True)
    # bank_key already folds in the band through start/end_seed, so an all-cut bank and an all-HDR
    # bank for the same contract cannot collide even before the separate directory.
    path = os.path.join(HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    meta["bank_path"] = path

    if not os.path.exists(path):
        bank_deadline = None if deadline is None else deadline - 20.0
        bank = build_bank(contract, reference, cell_types, ctx, sites, cfg, bank_deadline)
        if not bank:
            # An empty bank means either the deadline fired mid-scan or the scan finished and
            # nothing qualified. Those need different fixes — more budget vs a reachable
            # ``main_max_fail`` — so do not report them with one message.
            if bank_deadline is not None and time.monotonic() >= bank_deadline:
                meta["reason"] = "Cas12a HDR bank scan ran out of budget"
            else:
                meta["reason"] = (
                    f"no Cas12a guide reaches HDR on all but {cfg.main_max_fail} of "
                    f"{cfg.band_seeds().size} band seeds")
            return None, meta
        save_bank(path, bank)
    records = load_bank(path)
    meta["bank"] = len(records)
    if len(records) < cfg.group_size:
        meta["reason"] = f"HDR bank {len(records)} short of group {cfg.group_size}"
        return None, meta

    band_space = cfg.band_seeds()
    selector = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                             per_cell_min=cfg.per_cell_min,
                             caps=_group_caps(contract, ctx, cfg),
                             seeds=(band_space if cfg.seed_list else None))
    index, union = selector.best(cfg.group_size, restarts=cfg.restarts)
    group = [records[i] for i in index]
    bad: set[int] = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = np.array(sorted(set(int(x) for x in band_space) - bad), dtype=np.int64)
    span = int(band_space.size)
    meta.update(union=len(bad), clean=len(clean), clean_fraction=len(clean) / span)
    if clean.size == 0:
        meta["reason"] = "the group's HDR failures cover the whole band"
        MT.free_gpu_memory()
        return None, meta

    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    cas9 = scan_cas9(clean, contract, cell_types, ctx, sites, cfg, n_rows - cfg.group_size,
                     deadline)
    meta["cas9_pool"] = len(cas9)
    cells = len({(r["mutation"], r["strand"]) for r in cas9})
    if len(cas9) < n_rows - cfg.group_size or cells < 4:
        meta["reason"] = (f"Cas9 HDR pool {len(cas9)} over {cells} cells is short of "
                          f"{n_rows - cfg.group_size}")
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


def build_for_cell(contract: dict, reference: dict, cell_types: dict,
                   budget_s: float | None = None,
                   hdr_range: tuple[int, int] | None = None,
                   seed_list=None) -> tuple[list[dict] | None, dict]:
    """Build for whichever cell type this contract names, or decline where it is unmeasured.

    ``hdr_range`` overrides the cell type's default band. It is the per-hotkey decorrelation lever:
    the clean band lands inside this window, so running several hotkeys on disjoint windows makes
    their bands disjoint and multiplies the coldkey's round coverage (3 disjoint windows measured
    15.2% against 5.2% for the same window three times). The seeds are independent hashes, so band
    *position* is otherwise free — see the module docstring.
    """
    cell = contract.get("cell_type")
    cfg = config_for(cell)
    if cfg is None:
        return None, {"reason": f"no measured all-HDR config for {cell}"}
    if seed_list:
        # A joined, possibly non-contiguous band space. `hdr_range` is set to its min..max so the
        # span-derived knobs (`_scaled_max_fail`, the wide-window variants cap) see the real size.
        seeds = sorted(set(int(x) for x in seed_list))
        hdr_range = (seeds[0], seeds[-1])
        span = len(seeds)
        cfg = dataclasses.replace(cfg, hdr_range=hdr_range, seed_list=tuple(seeds),
                                  main_max_fail=_scaled_max_fail(cell, cfg.main_max_fail, span))
        if span > 100:
            cfg = dataclasses.replace(cfg, variants=min(cfg.variants, WIDE_WINDOW_VARIANTS))
        return build_submission(contract, reference, cell_types, cfg=cfg, budget_s=budget_s)
    if hdr_range is not None:
        # ``main_max_fail`` is calibrated per cell against a 100-seed band, so it has to scale with
        # the span or a wider window silently demands the impossible. Measured on K562 at 6 targets
        # x 2000 variants: the median guide fails 51 of 100 seeds but 153 of 300, so a flat 45
        # admits 1391 of 10000 guides at width 100 and **0** at width 300 — an empty bank, which
        # build_submission then reports as "ran out of budget". Holding the ratio (45 -> 135) keeps
        # 256 of 10000. Note the qualifying fraction still falls with width (13.9% -> 2.6%) because
        # the binomial concentrates as the window grows; that is the same mechanism that shrinks
        # the band itself from 13 seeds at width 100 to 9 at width 300.
        span = hdr_range[1] - hdr_range[0] + 1
        cfg = dataclasses.replace(cfg, hdr_range=hdr_range)
        if span != 100:
            cfg = dataclasses.replace(
                cfg, main_max_fail=_scaled_max_fail(cell, cfg.main_max_fail, span))
        if span > 100:
            # The scan cost scales with the span, so a wide window has to give some of it back in
            # `variants` or the fleet exceeds its GPU budget. See WIDE_WINDOW_VARIANTS.
            cfg = dataclasses.replace(cfg, variants=min(cfg.variants, WIDE_WINDOW_VARIANTS))
    return build_submission(contract, reference, cell_types, cfg=cfg, budget_s=budget_s)
