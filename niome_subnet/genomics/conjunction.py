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
import logging
import os
import time
from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np

import genExp as G
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.all_cut import (AllCutConfig, BANK_DIR, _params_fn, assemble, bank_key,
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
# **2026-09-15: every cell moved to `band_width` 300 by operator request. It is NOT a measured
# improvement and the measurement says the opposite, slightly.** CD34+ was swept at 100/150/225/300
# over six contracts, each priced against the one field that played its contract
# (`conj_replicate.py` with `CR_ARMS`, output `conj_replicate_CD34_HSPC_bw.json`):
#
#     width   own E[share]   vs all-HDR      re-priced   vs width 100
#       100      0.00177        2.34x          0.00125      1.000x     <- was shipped
#       150      0.00122        1.61x          0.00122      0.970x
#       225      0.00105        1.40x          0.00122      0.973x
#       300      0.00101        1.34x          0.00118      0.940x
#
# **Read the RE-PRICED column, not the raw one.** `conj_replicate.py` draws `v_clean`/`v_rest` from
# a shared RNG that advances between arms, so no two arms are scored on the same sampled seeds --
# two arms whose builds were byte-identical (band 8 / clean 183 / wxfid 317.9) priced at 0.00021 and
# 0.00009. Within a contract, across widths, `clean` spreads 2.0% and `wxfid` 0.9% while `v_clean`
# spreads 18%, so the raw 2.34x -> 1.34x ordering is almost entirely that sampling noise.
# `conj_bw_reprice.py` re-prices every arm on a COMMON `v_clean`/`v_rest` per contract, which leaves
# only the deterministic build outputs, and the real effect of 100 -> 300 is **-6%** with 4 of 6
# contracts tying exactly. The paired `v_clean` test behind the raw ordering is t = 1.78, p ~ 0.135.
#
# So widening is ~neutral on score and costs build time (`hdr_compliance` scales with the candidate
# count; `conj_bandwidth.json` measured 52s -> 103s from width 30 -> 300 on K562), which matters
# because these builds are already prefetch-dependent. Band depth is unaffected: k=8 built on 6 of 6
# CD34+ contracts at width 300. Revert by restoring 100 (CD34+) and 150 (K562).
CELL_CONFIG: dict[str, dict] = {
    "HEK293": {"band_k": 6, "group_size": 80, "band_width": 300, "light_cell_rows": 12},
    "CD34+_HSPC": {"band_k": 8, "group_size": 80, "band_width": 300, "light_cell_rows": 6},
    "K562": {"band_k": 8, "group_size": 80, "band_width": 300, "light_cell_rows": 6},
}


@dataclass(frozen=True)
class ConjunctionConfig(AllCutConfig):
    """all-cut over an explicit seed space, plus the HDR band that rides on the same rows."""

    # The joined (possibly non-contiguous) cut window. `start_seed`/`end_seed` are kept at its
    # min..max so `bank_key` and anything reading a coarse span stay correct.
    seed_list: tuple[int, ...] = ()
    band_k: int = 8
    band_width: int = 150
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


def hdr_compliance(records, contract: dict, cell_types: dict, ctx, candidates) -> dict:
    """For each banked guide, which of `candidates` it repairs by HDR on.

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
                                        params_of, "hdr", len(seeds))
        for i, guide in zip(idxs, guides):
            fails = got.get(guide)
            bad = set(int(x) for x in fails) if fails is not None else whole
            ok[i] = whole - bad
    return ok


def choose_band(ok: dict, k: int, candidates, need: int) -> tuple[list[int], list[int]]:
    """Greedily take a k-seed band, keeping only guides HDR-compliant on all of it.

    One seed per step, taking whichever candidate the most surviving guides comply on, and stopping
    early if the survivors would fall below `need` — a group cannot be formed from fewer guides
    than its own size, so a deeper band than that is worthless. Returns the band reached and the
    indices still alive on it; `len(band) < k` means the fill ran out and the caller should decline.
    """
    cand = sorted(set(int(x) for x in candidates))
    col = {s: j for j, s in enumerate(cand)}
    idx = sorted(ok)
    M = np.zeros((len(idx), len(cand)), dtype=bool)
    for r, i in enumerate(idx):
        for s in ok[i]:
            if s in col:
                M[r, col[s]] = True
    alive = np.ones(len(idx), dtype=bool)
    band: list[int] = []
    taken: set[int] = set()
    for _ in range(k):
        counts = M[alive].sum(axis=0)
        for j in taken:
            counts[j] = -1
        j = int(np.argmax(counts))
        if counts[j] < need:
            break
        taken.add(j)
        band.append(cand[j])
        alive &= M[:, j]
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
                  "band_width": cfg.band_width, "light_cell_rows": cfg.light_cell_rows,
                  "window": f"joined {len(joined)} seeds"}

    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    os.makedirs(BANK_DIR, exist_ok=True)
    path = os.path.join(BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    meta["bank_path"] = path

    if not os.path.exists(path):
        bank_deadline = None if deadline is None else deadline - 20.0
        bank = build_bank(contract, reference, cell_types, ctx, sites, cfg, bank_deadline)
        MT.free_gpu_memory()
        if not bank:
            # An empty bank is either the deadline firing mid-scan or nothing qualifying. Those
            # need different fixes, so do not report them with one message.
            if bank_deadline is not None and time.monotonic() >= bank_deadline:
                meta["reason"] = "Cas12a cut bank scan ran out of budget"
            else:
                meta["reason"] = (f"no Cas12a guide cuts on all but {cfg.cas12a_max_fail} of "
                                  f"{len(joined)} joined seeds")
            return None, meta
        save_bank(path, bank)
    records = load_bank(path, limit=cfg.bank_keep)
    meta["bank"] = len(records)
    if len(records) < cfg.group_size:
        meta["reason"] = f"cut bank {len(records)} short of group {cfg.group_size}"
        return None, meta

    candidates = sub_window(joined, cfg.band_width, cfg.band_offset_frac)
    ok = hdr_compliance(records, contract, cell_types, ctx, candidates)
    MT.free_gpu_memory()
    band, alive = choose_band(ok, cfg.band_k, candidates, cfg.group_size)
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

    cas9 = scan_cas9(clean, contract, cell_types, ctx, sites, cfg, n_rows - cfg.group_size,
                     deadline)
    MT.free_gpu_memory()
    meta["cas9_cut"] = len(cas9)
    if cas9:
        ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
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
                   hdr_range: tuple[int, int] | None = None) -> tuple[list[dict] | None, dict]:
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
    return build_submission(contract, reference, cell_types, cfg=cfg, budget_s=budget_s)
