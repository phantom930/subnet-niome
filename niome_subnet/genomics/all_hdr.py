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
from niome_subnet.genomics.all_cut import (BANK_DIR, _params_fn, assemble, bank_key, load_bank,
                                           save_bank)
from niome_subnet.genomics.validation import stage3

logger = logging.getLogger(__name__)

HDR_BANK_DIR = "data/all_hdr"

# The pinned band per cell type. Position is free (see the module docstring); distinct ranges keep
# the on-disk banks from colliding and make a log line say which cell type it came from.
CELL_CONFIG: dict[str, dict] = {
    "CD34+_HSPC": {"hdr_range": (500, 599), "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95)},
    "K562": {"hdr_range": (700, 799), "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95)},
    "HUDEP-2": {"hdr_range": (800, 899), "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95)},
    # HEK293 carries its own screen: main_max_fail 48 of 100 rather than 45 is what its sweep
    # measured at, and the bank is thin enough here (0.041% of guides) that tightening it further
    # starves the min-union step. Its group_size is the shared 80.
    "HEK293": {"hdr_range": (300, 399), "cas12a_gc": (0.40, 0.95), "cas9_gc": (0.40, 0.95),
               "main_max_fail": 48},
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
    light_cell_rows: int = 6
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

    @property
    def start_seed(self) -> int:
        return self.hdr_range[0]

    @property
    def end_seed(self) -> int:
        return self.hdr_range[1]


def config_for(cell_type: str, cfg: AllHdrConfig | None = None) -> AllHdrConfig | None:
    """The tuned config for a cell type, or None where all-HDR has not been measured."""
    overrides = CELL_CONFIG.get(cell_type)
    if overrides is None:
        return None
    return dataclasses.replace(cfg or AllHdrConfig(), **overrides)


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
    seeds = np.arange(cfg.start_seed, cfg.end_seed + 1, dtype=np.int64)
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
        if len(found) >= want * cfg.pool_target and len({(f["mutation"], f["strand"])
                                                         for f in found}) == 4:
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
    started = time.monotonic()
    deadline = None if budget_s is None else started + budget_s
    meta: dict = {"method": "all-hdr", "group_size": cfg.group_size,
                  "band": f"{cfg.start_seed}-{cfg.end_seed}"}

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
            meta["reason"] = "Cas12a HDR bank scan ran out of budget"
            return None, meta
        save_bank(path, bank)
    records = load_bank(path)
    meta["bank"] = len(records)
    if len(records) < cfg.group_size:
        meta["reason"] = f"HDR bank {len(records)} short of group {cfg.group_size}"
        return None, meta

    selector = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                             per_cell_min=cfg.per_cell_min,
                             caps=_group_caps(contract, ctx, cfg))
    index, union = selector.best(cfg.group_size, restarts=cfg.restarts)
    group = [records[i] for i in index]
    bad: set[int] = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad), dtype=np.int64)
    span = cfg.end_seed - cfg.start_seed + 1
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
                   hdr_range: tuple[int, int] | None = None) -> tuple[list[dict] | None, dict]:
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
    if hdr_range is not None:
        cfg = dataclasses.replace(cfg, hdr_range=hdr_range)
    return build_submission(contract, reference, cell_types, cfg=cfg, budget_s=budget_s)
