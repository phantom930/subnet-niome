"""conjunction.py — all-cut's clean set AND all-HDR's pinned band, on the same rows.

all-cut pins ``is_cut`` over a wide clean set and scores ~0.19-0.24 on each of those seeds.
all-HDR pins all three stage-4 targets over a narrow band and scores exactly 1.0 there, ~0.10
everywhere else. This runs both filters over one row set: a Cas12a group min-unioned on **cut**
over a joined seed space, restricted first to guides that also repair by **HDR** on every seed of
a k-seed band, with the Cas9 half strict on cut over the clean set and on HDR over the band.

That gives three per-seed regimes instead of two:

    band seed        all three targets pinned     cons exactly 1.0        k seeds
    clean, off band  ``is_cut`` pinned            cons 0.19-0.24          ~160-190 of the joined space
    everything else  nothing pinned               cons ~0.10

**CLAUDE.md records this construction as falsified and that entry is stale.** It was measured at
CONTIGUOUS cut windows (100 / 300 / 900) where E[final] fell monotonically in k. Over the JOINED
300-seed space the fleet actually plays, the clean set reaches 160-190 with ``band_cons`` verified
at exactly 1.0000 on every arm — the "spike on an elevated floor" that no earlier construction
produced. The mechanism that killed the contiguous version is still present (the HDR filter shrinks
the Cas12a pool, and a smaller pool min-unions to a larger failed-seed union) and is simply paid
for by the joined window's larger reach.

What it is worth, per cell type, measured over **12 contracts each** and priced against the single
field that played each contract — never a pooled field, which prices field softness instead of the
construction (CLAUDE.md, "Pricing a construction"). Ratios are E[curve share] against a matched
all-HDR build on the same contract and the same joined window, Monte-Carlo'd per seed:

    cell          arm (k/group/width/light)   own-field MC   wins    p      orig 6   new 6
    HEK293        6 / 80 / 300 / 12               2.05x      10/12   0.039   3.00x   1.33x
    CD34+_HSPC    8 / 80 / 100 / 6                1.77x       9/12   0.146   1.65x   1.93x
    K562          8 / 80 / 150 / 6                1.30x      10/12   0.039   1.48x   1.17x
    HUDEP-2       8 / 80 / 225 / 12               0.89x       3/12   0.146   0.90x   0.88x

**HUDEP-2 is deliberately absent from CELL_CONFIG.** It lost on 12 of 12 contracts' worth of
aggregate and 0 of the 6 fresh ones, leave-one-out 0.87-0.91x, so it keeps all-HDR. ``config_for``
returning None for it is the mechanism, and the miner falls through to the all-HDR rung.

**Three honest limits on the table above, all of which should be read before widening this.**

  * **Every effect shrank from n=6 to n=12** (HEK293 3.23 -> 2.05, K562 1.57 -> 1.30, CD34+
    1.92 -> 1.77). These arms were chosen off a 360-config grid per cell, so that is the expected
    regression from selection. The fresh-six halves are the only genuinely out-of-sample numbers:
    5/6, 5/6, 5/6 across the three shipped cells (15/18, p = 0.008). HEK293's halves disagree by
    more than 2x, so treat 2.05x as the least trustworthy of the three point estimates even though
    it is the largest.
  * **There is no single mechanism.** HEK293 wins on term 1 — its ``weighted x fidelity`` comes out
    **1.046x** all-HDR's, 11 of 12 contracts — which is the one cell where all-HDR is documented as
    weak against the leaders' bar. CD34+ and K562 are flat on term 1 (1.007x, 1.006x) and win purely
    on the clean regime. A cut10-median split runs in OPPOSITE directions across the three
    (CD34+ 2.18 soft / 1.46 hard; HEK293 1.72 / 3.36; K562 1.10 / 1.56), so the "soft tail" story
    that motivated this does not generalise and is not why it wins.
  * **This is a PER-HOTKEY measurement, not a fleet one.** Correlated siblings are what sank all-cut
    at fleet level despite winning per hotkey (CLAUDE.md, "Per-hotkey P(place) is the wrong unit for
    a fleet"), and the band here is a deterministic function of (contract, joined space, width), so
    two hotkeys handed the same window build the same band. It is safe as shipped because the fleet
    is down to ONE band hotkey (h0, ``FULL_HK``, the whole joined space) — the per-hotkey number IS
    the fleet number at n=1. **If the fleet regrows, rotate ``band_offset_frac`` per hotkey before
    trusting this**, and re-price it with ``fleet_price.py``.

Build quality is not uniform and the failure is visible in term 1: on 2 of 48 contracts across all
four cells the conjunction's ``weighted x fidelity`` came out far below the matched all-HDR build
(``a8b9f1bb`` HEK293 166.3 vs 180.5; ``33d826fc`` HUDEP-2 204.9 vs 289.9), and both are that cell's
worst contract. Every other contract matched within 3%. Cause unknown; it is not screened for here.

The construction, in order:

  1. bank Cas12a guides on the ``cut`` rule over the joined seed space (an ordinary all-cut bank —
     ``all_cut.bank_key`` folds ``seed_list`` in, so it shares ``data/all_cut`` without colliding)
  2. screen the bank for HDR compliance and greedily take a k-seed band inside a sub-window,
     keeping at each step only guides that repair by HDR on the whole band so far
  3. min-union the survivors on cut -> the clean set
  4. fill the Cas9 half strictly on cut over the clean set AND on HDR over the band
  5. assemble, mutation-weighted, exactly as all-cut does

Declines to ``(None, meta)`` at every step that cannot be met, so the miner falls through to
all-HDR — which is the build the fleet had before this module existed.
"""
from __future__ import annotations

import dataclasses
import math
import logging
import os
import time
from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np

import genExp as G
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.all_cut import (AllCutConfig, BANK_DIR, _params_fn, assemble, bank_key,
                                           bank_slot, cas9_cell_target,
                                           build_bank, config_for as all_cut_config_for, load_bank,
                                           save_bank, scan_cas9)
from niome_subnet.genomics.validation import stage3

logger = logging.getLogger(__name__)

# The tuned arm per cell type, from the 12-contract replication in the module docstring.
# HUDEP-2 is absent on purpose — it measured 0.89x and keeps all-HDR.
#
# **The four values move together and must not be mixed across rows.** `band_width` sets how much
# of the joined space the band may be drawn from, and `band_k` is the depth the Cas9 conditional
# fill can still reach at that width: survivors fall ~0.55x per band seed, so a k that is feasible
# at width 300 is not automatically feasible at width 100. Width 300 means "the whole joined space",
# i.e. no sub-window restriction at all.
#
# **2026-09-15 (later): the four-hotkey fleet arm, set by operator request.** Every hotkey
# min-unions on cut over the same 300-seed joined space and the band is drawn from a width-150
# sub-window whose offset rotates per hotkey (`joined_window.band_offset_frac`, stride 75 -> offsets
# 0/75/150/225). At `band_width` 300 every hotkey would search the identical candidate set and build
# the IDENTICAL band, which is the correlation that sank all-cut at fleet level; at 150 the four
# bands measured 0-3 shared seeds per pair on a live round, a fleet union of 31 of 900 against a
# single hotkey's 8.
#
# **`band_k` is 8 on every cell.** It was 10/8 briefly on 2026-09-15, then 8/6, and HEK293 came
# back to 8 on the arm comparison recorded at the bottom of this block. What is measured, and what
# is not:
#
#   * The only PAIRED fleet measurement is group 80 / k=8 against group 100 / k=10, over the same
#     10 tasks, each hotkey priced in the field that played its contract (`conj_fleet.py`, paired
#     by `conj_fleet_pair.py` -- pairing within task matters because `total_weighted_score` moves
#     54% with the contract):
#
#         arm        places  median final/cut10  band>=1  mean |band|  |clean|  built
#         g80  k8      1/5          0.92x          3/5         56        217    10.0/10
#         g100 k10     0/5          0.00x          0/5         30        119     4.0/10
#
#     g100k10 lost on 0 of 5 tasks, mean paired delta -60.5 final points; per cell CD34+ -17.3
#     (n=1), HUDEP-2 -4.7 (n=1), HEK293 -93.5 (n=3, a BUILD FAILURE -- zero rows on all three, the
#     on-band Cas9 fill cannot reach 150 rows at ~P(HDR)**10 and accessibility 0.35).
#
#   * **`band_k` 8 vs 6, and cut width 300 vs 900, were then measured directly** by rebuilding the
#     HEK293 round `ba815f07` (seeds 580/769/569) on all four band offsets and pricing every arm in
#     that round's own 248-miner field. The two quantities that pay are the fleet band union (spike
#     frequency) and the absolute clean-seed count (the elevated floor):
#
#         arm                band union  P(>=1 hit)  clean seeds of 900  cold build
#         cut300 k=8 (live)      31        9.98%           19-23             51s
#         cut300 k=6             23        7.47%           38-39             41s
#         cut900 k=8             30        9.67%           16-20            180s
#         cut900 k=6             22        7.16%           34-37            180s
#
#     **Two structural results, each consistent across both k values.** (a) The clean-seed count is
#     set by `band_k`, NOT by the width of the cut space -- tripling the cut space to 900 leaves it
#     unchanged or slightly lower, because the min-union's job grows exactly as fast as the space.
#     Band union behaves the same way. So the 900-seed cut space is DOMINATED: equal on both paying
#     terms, 3.5x the cold build. Do not re-run it. (b) `choose_band` is a greedy prefix, so the k=6
#     band is a strict SUBSET of the k=8 band on every hotkey (verified 4/4) -- band coverage is
#     monotone in k and lowering k can only shrink the spike.
#
#     HEK293 therefore sits at k=8: the ~35% more band union is worth more than the ~2x clean set
#     k=6 buys, because HEK293's clean set is only ~38 of 900 even at k=6 -- too small to pay. The
#     erythroid cells hold k=8 for the OPPOSITE reason: their clean sets are 161-167, a floor that
#     does pay, and k=10 collapsed it to ~125 while starving the Cas9 fill.
#
#     Caveat on what this is: one contract, one field. The mechanism half (subset property, clean
#     count independent of cut width) is exact and holds on every contract; the ranking of k=8 over
#     k=6 on HEK293 is a mechanism argument plus a single observation, not a replication.
#
#   * **The arm below is neither of those two.** It pairs group 100 with k=8 on the erythroid cells
#     and group 80 with k=6 on HEK293, so it takes the SHALLOWER band of the winning arm with the
#     LARGER group of the losing one. Nothing has measured that combination. What the k=10 -> 8
#     move is expected to recover is the two things that collapsed at k=10 -- the clean set
#     (119 -> 217 between the two measured arms) and the Cas9 pool, since the on-band filter keeps
#     ~P(HDR)**k -- at the cost of band depth, i.e. spike FREQUENCY. Pre-flight the pool before
#     trusting it: at k=10 CD34+ ran 164-213 candidates against the 150 rows it must fill.
#
# Reverting is one row each. Before the four-hotkey layout the shipped values were HEK293
# k=6/g80/w300, CD34+ and K562 k=8/g80/w300, and no HUDEP-2 entry at all (it kept all-HDR on a
# 0.89x 12-contract replication).
# **2026-09-16: band depth moved to each cell's measured optimum, with `band_cell_aware` ON.**
# Set by operator request, and unlike the two arms above this one IS the measurement — see
# "Band depth, and the per-cell defect that was capping it" in CLAUDE.md. Two things changed:
#
#   * `band_cell_aware` scales `scan_cas9`'s per-cell floor by `P(rule)**k`, because the on-band
#     filter decimates that pool AFTER the scan has stopped. Without it, depth declines: 3 of 12
#     HEK293 contracts at k=8 and 5 of 12 at k=9, and 6 of 12 erythroid at k=11-12. With it, 12/12
#     and 12/12. Paired within contract over 3 erythroid cells x 3 depths it is **20W/2L,
#     sign p = 0.00012**, and it lifts builds that already succeeded (HEK293 `w x fid` +3.25%).
#   * `band_k` is each cell's optimum. Three of the four are the Cas12a BAND-FORMATION WALL -- the
#     depth at which `choose_band`'s survivors fall under `group_size` -- which no Cas9 floor can
#     move and which `group_size` 100 -> 80 does not move either (it is a 0.8x cut against a 0.57x
#     per-seed decay, so it buys the next seed on 1 contract in 4; measured, and it LOSES on
#     E[share]). Only CD34+ has an interior maximum.
#
#     | cell | wall | k here | own-field E[share] | vs the k=8 it replaces |
#     |---|---|---|---|---|
#     | HEK293 | 9 | **9** | 0.000163 | k=8 declined 3/12 |
#     | CD34+_HSPC | 12 | **12**, wall | 0.000145 (k=11 peaks at 0.000153) | 1.43x at its peak |
#     | K562 | 12 | **12**, wall | 0.000220 | 1.36x |
#     | HUDEP-2 | 12 | **12**, wall | 0.000303 | 2.06x |
#
# **CD34+ is set to the wall, not to its peak.** k=11 measures 0.000153 against k=12's 0.000145 --
# a 5% interior maximum at n=4 contracts, which is inside the noise this file's own rules would
# demand more evidence for. Set at 12 by operator request for a uniform erythroid k; revert to 11
# if CD34+ is ever swept properly.
#
# **HUDEP-2 is the load-bearing caveat.** [hud_resolve.py](hud_resolve.py) re-measured the 0.89x
# exclusion on six fresh contracts and **it holds for k=8 (0.94x, 2/6 wins)** -- the arm that was
# shipping. k=12 + the floor is the first HUDEP-2 conjunction arm that does not lose, at **1.12x,
# 4/6, sign p = 0.688**. That is "stops losing", NOT a measured win, and all-HDR remains the
# lower-variance choice on that cell. The mechanism is that HUDEP-2's all-HDR band is **11.5**, so
# any conjunction pinned below ~12 is giving up the dominant term to buy a clean set.
CELL_CONFIG: dict[str, dict] = {
    "HEK293": {"band_k": 9, "group_size": 80, "band_width": 150, "light_cell_rows": 12,
               "band_cell_aware": True},
    "CD34+_HSPC": {"band_k": 12, "group_size": 100, "band_width": 150, "light_cell_rows": 6,
                   "band_cell_aware": True},
    "K562": {"band_k": 12, "group_size": 100, "band_width": 150, "light_cell_rows": 6,
             "band_cell_aware": True},
    "HUDEP-2": {"band_k": 12, "group_size": 100, "band_width": 150, "light_cell_rows": 6,
                "band_cell_aware": True},
}


@dataclass(frozen=True)
class ConjunctionConfig(AllCutConfig):
    """all-cut over an explicit seed space, plus the HDR band that rides on the same rows."""

    # The joined (possibly non-contiguous) cut window. `start_seed`/`end_seed` are kept at its
    # min..max so `bank_key` and anything reading a coarse span stay correct.
    seed_list: tuple[int, ...] = ()
    band_k: int = 8
    band_width: int = 150
    # Let `choose_band` see the Cas9 side before it fixes the band. OFF by default: the band it
    # picks differs from the measured one, so every tuned (k, group, width) row above was measured
    # without it. See `cas9_cell_probe` for what it is for and what it costs.
    band_cell_aware: bool = False
    band_probe_sites: int = 8
    band_probe_guides: int = 2000
    # Ceiling on the band-scaled per-cell Cas9 floor. The scan walks jobs nearest-first, so filling
    # a thin cell means walking most of the distance ordering; this bounds that against the build
    # budget. The deadline still bounds it independently.
    band_fill_cap: int = 50000
    # The rule the band is pinned on. "hdr" is the measured construction; see `hdr_compliance` for
    # what changing it costs. Per-row compliance sets the Cas9 conditional fill, so a rule with
    # lower compliance reaches a smaller feasible `band_k` at the same pool.
    band_rule: str = "hdr"
    # Where the band sub-window starts, as a fraction of the joined space. 2/3 reproduces the
    # measurement, which fixed the start at seed 700 of the joined space [100-199, 400-499,
    # 700-799] — i.e. the start of the third class. The start was held fixed through the whole
    # 360-config grid, so it is NOT a measured parameter; it is carried here only so the shipped
    # slice has the same shape as the measured one. It is also the knob to rotate per hotkey if
    # the fleet regrows past one band hotkey — see the module docstring.
    band_offset_frac: float = 2.0 / 3.0
    # `AllCutConfig` has no `restarts` (only `AllHdrConfig` does) and `all_cut.build_submission`
    # passes its own literal, so the value the research ran at has to be carried here.
    restarts: int = 12
    # `scan_cas9` stops at `pool_target * want` candidates. The band filter then keeps only
    # ~P(HDR)**k of those, so the ordinary 8x leaves ~83 of 21,144 at k=8 against ~170 rows needed
    # and the build declines on an early exit rather than on a real shortage. 500 is effectively a
    # full scan and is the dominant cost of this build.
    pool_target: int = 500

    @property
    def seeds(self) -> np.ndarray:
        if self.seed_list:
            return np.asarray(sorted(set(int(x) for x in self.seed_list)), dtype=np.int64)
        return np.arange(self.start_seed, self.end_seed + 1, dtype=np.int64)


# Fields `ConjunctionConfig` deliberately re-defaults away from `AllCutConfig`'s value. Copying the
# cell's all-cut config wholesale would put them back: it silently reset `pool_target` to 8, the
# Cas9 scan then exited after a handful of sites, the band filter starved the heavy-mutation cell,
# and `total_weighted_score` came out 243 against a reachable 329 — the same starvation mode
# `all_cut.cas9_cell_target` documents, reached by a different route.
_OWN_DEFAULTS = ("pool_target",)


def config_for(cell_type: str, cfg: ConjunctionConfig | None = None) -> ConjunctionConfig | None:
    """The tuned config for a cell type, or None where the conjunction lost or is unmeasured.

    Starts from the cell's all-cut config, so the GC bands, `max_distance` and `cas12a_max_fail`
    are the ones the cut bank was measured with — the conjunction's Cas12a half IS an all-cut bank.
    """
    overrides = CELL_CONFIG.get(cell_type)
    if overrides is None:
        return None
    base = all_cut_config_for(cell_type) or AllCutConfig()
    seeded = ConjunctionConfig(**{f.name: getattr(base, f.name)
                                  for f in dataclasses.fields(AllCutConfig)
                                  if f.name not in _OWN_DEFAULTS})
    return dataclasses.replace(cfg or seeded, **overrides)


def sub_window(seeds: list[int], width: int, offset_frac: float = 2.0 / 3.0) -> list[int]:
    """`width` seeds of the joined space, starting `offset_frac` through it and wrapping."""
    n = len(seeds)
    if width >= n:
        return list(seeds)
    start = int(round(n * offset_frac)) % n
    return [seeds[(start + j) % n] for j in range(width)]


def hdr_compliance(records, contract: dict, cell_types: dict, ctx, candidates,
                   rule: str = "hdr") -> dict:
    """For each banked guide, which of `candidates` it satisfies ``rule`` on.

    ``rule`` defaults to "hdr", the construction as measured. Any key of
    ``mt19937.RULE_SPECS`` works; what changes is how many of stage 4's three targets a band seed
    pins, and therefore what a band seed is WORTH:

        hdr      all 3 pinned                     cons exactly 1.0000
        mh_any   is_cut pinned, is_hdr == mh      cons 0.77-0.85 (measured, K562 and HEK293)

    `mh_any` never makes `is_hdr` constant -- it makes it an exact function of `mh`, which stage 4
    hands the forest as a `build_X` feature, so `r2_score` pays it like a constant. `indel_length`
    is only partly recovered (the BLUNT branch carries an unpinned expovariate), which is the whole
    of the 0.8-vs-1.0 gap.

    Grouped by (site, mutation) because the screen is per target: every guide of one target shares
    the gc/energy/cut_p columns, which is what makes the batched kernel worth using.
    """
    accessibility = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    seeds = np.asarray(sorted(set(int(x) for x in candidates)), dtype=np.int64)
    whole = set(int(s) for s in seeds)
    groups: dict[tuple, list[int]] = defaultdict(list)
    for i, rec in enumerate(records):
        groups[(rec["start"], rec["strand"], rec["cas_system"], rec["length"],
                rec["mutation"])].append(i)
    ok: dict[int, set] = {}
    for (start, strand, cas, length, mutation), idxs in groups.items():
        site = type("S", (), {"start": start, "strand": strand, "cas": cas, "length": length})()
        distance = abs(start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        params_of = _params_fn(site, distance, accessibility, offset)
        guides = [records[i]["guide"] for i in idxs]
        got = MT.screen_guides_rule_gpu(guides, seeds, mutation, cas, start, strand,
                                        params_of, rule, len(seeds))
        for i, guide in zip(idxs, guides):
            fails = got.get(guide)
            bad = set(int(x) for x in fails) if fails is not None else whole
            ok[i] = whole - bad
    return ok


def cas9_cell_probe(contract: dict, cell_types: dict, ctx, sites, cfg, candidates,
                    rule: str) -> dict:
    """Per (mutation, strand) cell, the rule-compliance sets of a SAMPLE of Cas9 guides.

    `choose_band` picks the band from Cas12a guides alone, but the constraint that actually
    declines the build lives on the Cas9 side: stage 5 needs all four cells populated, and the
    on-band filter keeps only ~P(rule)**k of each cell's candidates. A cell that is thin to begin
    with empties first -- on HEK293 the heavy-mutation Cas9 sites sit farther from the mutation
    than the light ones, which is the same starvation `all_cut.cas9_cell_target` documents -- and
    nothing in the Cas12a pool tells `choose_band` it is about to happen. Measured over 12 HEK293
    contracts x depths 6-9, **8 of 9 declines were cell coverage, not pool size**; one had a Cas9
    pool of 1013 against a requirement of 170 and still declined.

    This probes the Cas9 side BEFORE the band is fixed, so the band can be chosen to keep every
    cell alive. It is deliberately OPTIMISTIC: the real pool is further restricted to guides that
    cut on every seed of `clean`, which does not exist until the group does. At `cut_p` 0.99 over a
    ~25-seed clean set that is a ~0.78 factor, so a cell this calls viable is usually viable and a
    cell it calls dead always is. It is a tie-break, not a guarantee.

    Cost is one extra GPU screen of `band_probe_guides` guides over `candidates` -- ~1.2M pairs at
    the defaults, against the 23.9M pairs/s `screen_guides_rule_gpu` sustains.
    """
    cell = contract.get("cell_type")
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    seeds = np.asarray(sorted(set(int(x) for x in candidates)), dtype=np.int64)
    whole = set(int(x) for x in seeds)
    jobs = [(i, m, abs(sites[i].start - ctx.mutation_map[m]))
            for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda job: job[2])
    cap = max(1, cfg.band_probe_guides // max(1, cfg.band_probe_sites))
    seen_sites: dict[tuple, int] = {}
    out: dict[tuple, list] = {}
    for site_index, mutation, distance in jobs:
        site = sites[site_index]
        key = (mutation, site.strand)
        if seen_sites.get(key, 0) >= cfg.band_probe_sites:
            continue
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cap)
        if not guides:
            continue
        seen_sites[key] = seen_sites.get(key, 0) + 1
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        params_of = _params_fn(site, distance, accessibility, offset)
        got = MT.screen_guides_rule_gpu(guides, seeds, mutation, "Cas9", site.start,
                                        site.strand, params_of, rule, len(seeds))
        for guide in guides:
            fails = got.get(guide)
            if fails is None:
                continue
            out.setdefault(key, []).append(whole - set(int(x) for x in fails))
        if len(seen_sites) >= 4 and min(seen_sites.values()) >= cfg.band_probe_sites:
            break
    return out


def choose_band(ok: dict, k: int, candidates, need: int, cell_ok: dict | None = None,
                stats: dict | None = None) -> tuple[list[int], list[int]]:
    """Greedily take a k-seed band, keeping only guides HDR-compliant on all of it.

    One seed per step, taking whichever candidate the most surviving guides comply on, and stopping
    early if the survivors would fall below `need` — a group cannot be formed from fewer guides
    than its own size, so a deeper band than that is worthless. Returns the band reached and the
    indices still alive on it; `len(band) < k` means the fill ran out and the caller should decline.

    `cell_ok` (from `cas9_cell_probe`) makes the step CELL-AWARE: among the seeds that clear
    `need`, take the best one that still leaves every (mutation, strand) cell holding a
    rule-compliant Cas9 probe guide, rather than the best one outright. That is a tie-break on an
    almost-flat surface — the greedy's `argmin`/`argmax` is heavily tied, as `fastgreedy` documents
    for its own picks — so it usually costs nothing in surviving Cas12a pool while avoiding the
    band seeds that empty a thin cell and decline the build at the last step.

    When no seed clears `need` AND keeps every cell, it falls back to the plain argmax rather than
    stopping: a band that might decline still beats no band. `stats` records how often that
    happened, so a build that declines anyway can be told apart from one the probe never steered.
    """
    cand = sorted(set(int(x) for x in candidates))
    col = {s: j for j, s in enumerate(cand)}
    idx = sorted(ok)
    M = np.zeros((len(idx), len(cand)), dtype=bool)
    for r, i in enumerate(idx):
        for s in ok[i]:
            if s in col:
                M[r, col[s]] = True

    # Per-cell probe matrices in the same column space, carried alongside `alive`.
    CM: dict = {}
    calive: dict = {}
    for key, sets in (cell_ok or {}).items():
        if not sets:
            continue
        A = np.zeros((len(sets), len(cand)), dtype=bool)
        for r, ss in enumerate(sets):
            for seed in ss:
                if seed in col:
                    A[r, col[seed]] = True
        CM[key] = A
        calive[key] = np.ones(len(sets), dtype=bool)

    alive = np.ones(len(idx), dtype=bool)
    band: list[int] = []
    taken: set[int] = set()
    steered = fell_back = 0
    for _ in range(k):
        counts = M[alive].sum(axis=0)
        for j in taken:
            counts[j] = -1
        pick = -1
        if CM:
            for cnd in np.argsort(-counts):
                cnd = int(cnd)
                if counts[cnd] < need:
                    break
                if all((calive[key] & CM[key][:, cnd]).any() for key in CM):
                    pick = cnd
                    break
            if pick >= 0 and pick != int(np.argmax(counts)):
                steered += 1
            elif pick < 0:
                fell_back += 1
        if pick < 0:
            pick = int(np.argmax(counts))
            if counts[pick] < need:
                break
        taken.add(pick)
        band.append(cand[pick])
        alive &= M[:, pick]
        for key in CM:
            calive[key] &= CM[key][:, pick]
    if stats is not None:
        stats.update(band_steered=steered, band_fell_back=fell_back,
                     probe_cells=len(CM),
                     probe_alive={f"{m}|{st}": int(v.sum()) for (m, st), v in calive.items()})
    return band, [idx[r] for r in np.flatnonzero(alive)]


def build_submission(contract: dict, reference: dict, cell_types: dict,
                     cfg: ConjunctionConfig | None = None,
                     budget_s: float | None = None) -> tuple[list[dict] | None, dict]:
    """The conjunction submission, or ``(None, meta)`` when the caller should fall back."""
    cfg = cfg or ConjunctionConfig()
    started = time.monotonic()
    deadline = None if budget_s is None else started + budget_s
    joined = [int(x) for x in cfg.seeds]
    meta: dict = {"method": "conjunction", "k": cfg.band_k, "group_size": cfg.group_size,
                  "band_width": cfg.band_width, "band_rule": cfg.band_rule,
                  "band_cell_aware": cfg.band_cell_aware,
                  "light_cell_rows": cfg.light_cell_rows,
                  "window": f"joined {len(joined)} seeds"}

    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    os.makedirs(BANK_DIR, exist_ok=True)
    path = os.path.join(BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    meta["bank_path"] = path

    if not os.path.exists(path):
        bank_deadline = None if deadline is None else deadline - 20.0
        # Every hotkey on this cell folds the same fields into `bank_key`, so without this the
        # four prefetching miners build the identical bank four times and discard three. See
        # `all_cut.bank_slot` for what that costs and what serialising it is measured to buy.
        with bank_slot(path, bank_deadline) as build_it:
            if build_it:
                bank = build_bank(contract, reference, cell_types, ctx, sites, cfg, bank_deadline)
                MT.free_gpu_memory()
                if not bank:
                    # An empty bank is either the deadline firing mid-scan or nothing qualifying.
                    # Those need different fixes, so do not report them with one message.
                    if bank_deadline is not None and time.monotonic() >= bank_deadline:
                        meta["reason"] = "Cas12a cut bank scan ran out of budget"
                    else:
                        meta["reason"] = (f"no Cas12a guide cuts on all but {cfg.cas12a_max_fail} "
                                          f"of {len(joined)} joined seeds")
                    return None, meta
                save_bank(path, bank)
            else:
                meta["bank_shared"] = True
    records = load_bank(path, limit=cfg.bank_keep)
    meta["bank"] = len(records)
    if len(records) < cfg.group_size:
        meta["reason"] = f"cut bank {len(records)} short of group {cfg.group_size}"
        return None, meta

    candidates = sub_window(joined, cfg.band_width, cfg.band_offset_frac)
    ok = hdr_compliance(records, contract, cell_types, ctx, candidates, cfg.band_rule)
    MT.free_gpu_memory()
    cell_ok = None
    p_row = None
    if cfg.band_cell_aware:
        cell_ok = cas9_cell_probe(contract, cell_types, ctx, sites, cfg, candidates,
                                  cfg.band_rule)
        MT.free_gpu_memory()
        meta["probe_guides"] = {f"{m}|{st}": len(v) for (m, st), v in cell_ok.items()}
        # Per-row Cas9 compliance with the band rule, measured rather than assumed. This is what
        # the on-band filter will charge the pool per band seed, so it sets how much bigger than
        # `assemble`'s quota each cell has to be BEFORE the filter.
        frac = [len(x) / max(1, len(candidates)) for v in cell_ok.values() for x in v]
        p_row = float(np.mean(frac)) if frac else None
        meta["probe_p_row"] = p_row
    bstats: dict = {}
    band, alive = choose_band(ok, cfg.band_k, candidates, cfg.group_size, cell_ok, bstats)
    meta.update(bstats)
    meta["band"] = len(band)
    meta["band_seeds"] = list(band)
    if len(band) < cfg.band_k:
        # Fewer band seeds than the arm was measured at is a different construction, not a weaker
        # one — the group/width/light values are tuned at this k. Decline rather than ship it.
        meta["reason"] = (f"band reached {len(band)} of {cfg.band_k} seeds before the surviving "
                          f"pool fell below the group")
        return None, meta
    pool = [records[i] for i in alive]
    meta["pool"] = len(pool)
    if len(pool) < cfg.group_size:
        meta["reason"] = f"HDR-on-band pool {len(pool)} short of group {cfg.group_size}"
        return None, meta

    selector = FG.FastGreedy(pool, window_lo=min(joined), window_hi=max(joined),
                             seeds=np.asarray(joined, dtype=np.int64))
    index, _union = selector.best(cfg.group_size, restarts=cfg.restarts)
    group = [pool[i] for i in index]
    bad: set[int] = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = np.array(sorted(set(joined) - bad), dtype=np.int64)
    meta.update(union=len(bad), clean=int(clean.size),
                clean_fraction=clean.size / max(1, len(joined)))
    if clean.size == 0:
        meta["reason"] = "the group's cut failures cover the whole joined window"
        MT.free_gpu_memory()
        return None, meta

    want = n_rows - cfg.group_size
    cell_floor = None
    if p_row and p_row > 0:
        survive = p_row ** len(band)
        cell_floor = min(cfg.band_fill_cap,
                         int(math.ceil(cas9_cell_target(contract, ctx, cfg, want) / survive)))
        meta["cell_floor"] = cell_floor
    cas9 = scan_cas9(clean, contract, cell_types, ctx, sites, cfg, want, deadline, cell_floor)
    MT.free_gpu_memory()
    meta["cas9_cut"] = len(cas9)
    if cas9:
        ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band, cfg.band_rule)
        cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
        MT.free_gpu_memory()
    meta["cas9_pool"] = len(cas9)
    cells = len({(r["mutation"], r["strand"]) for r in cas9})
    if len(cas9) < n_rows - cfg.group_size or cells < 4:
        meta["reason"] = (f"Cas9 pool {len(cas9)} over {cells} cells is short of "
                          f"{n_rows - cfg.group_size} after the band filter "
                          f"(cut-clean candidates {meta['cas9_cut']})")
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
                   seed_list=None,
                   hdr_range: tuple[int, int] | None = None,
                   band_offset_frac: float | None = None) -> tuple[list[dict] | None, dict]:
    """Build for whichever cell type this contract names, or decline where it is unmeasured.

    Takes the SAME seed space the all-HDR rung would have been given, so the two rungs play one
    window layout and a hotkey's decorrelation carries over unchanged. A contiguous ``hdr_range``
    is accepted and expanded, but the measurement is over a joined space of ~300 seeds and a
    contiguous window is not the operating point any of these arms were tuned at.
    """
    cell = contract.get("cell_type")
    cfg = config_for(cell)
    if cfg is None:
        return None, {"reason": f"no measured conjunction config for {cell}"}
    if seed_list:
        seeds = sorted(set(int(x) for x in seed_list))
    elif hdr_range is not None:
        seeds = list(range(hdr_range[0], hdr_range[1] + 1))
    else:
        return None, {"reason": "the conjunction needs a seed space; none was given"}
    if len(seeds) < cfg.band_width:
        # A window narrower than the band sub-window is not the measured arm: the band would be
        # drawn from the whole space, which is HEK293's configuration and nobody else's.
        return None, {"reason": f"seed space {len(seeds)} narrower than band width "
                                f"{cfg.band_width}"}
    # `cas12a_max_fail` is calibrated against the 900-seed window all-cut banks over, so it has to
    # scale with the span or a 300-seed window silently demands a far stricter guide than the one
    # that was measured. This is the same linear hold the research used (mf 100 -> 33 at span 300).
    mf = max(1, round(cfg.cas12a_max_fail * len(seeds) / 900))
    cfg = dataclasses.replace(cfg, seed_list=tuple(seeds), start_seed=seeds[0], end_seed=seeds[-1],
                              cas12a_max_fail=mf)
    # Per-hotkey band rotation. None keeps `ConjunctionConfig`'s own default, which is the
    # single-hotkey value the 12-contract replication was measured at; the fleet passes
    # `joined_window.band_offset_frac`, so two siblings on one joined window hold different bands.
    if band_offset_frac is not None:
        cfg = dataclasses.replace(cfg, band_offset_frac=float(band_offset_frac))
    return build_submission(contract, reference, cell_types, cfg=cfg, budget_s=budget_s)
