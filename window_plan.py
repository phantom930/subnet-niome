#!/usr/bin/env python3
"""window_plan.py — assign each hotkey its clean-band window for the coming round.

Concentrates several hotkeys onto the 100-seed window the next task is predicted to draw from
(seed_window_model.predict), and leaves the rest spread across their usual territories. The miner
reads the result per build from data/window_plan.json and falls back to its NIOME_HDR_WINDOW pin if
the file is missing, stale or malformed.

**Why this is safe to run before the prediction is proven.** Concentration is exactly EV-neutral
under a uniform generator: six hotkeys at width 16 cover 6 x 14.7 = 88 band seeds, and a uniform
seed hits them with probability 88/900 whether they sit inside one window or across nine. It only
starts costing when the concentrated hotkeys **saturate** the 100-seed window they share -- seven at
width 16 cover 103 band seeds inside 100 and waste the overlap (37.2% against 37.9%) -- or when a
hotkey is spent on something other than a band. So CONCENTRATE is capped below saturation per cell
type, and no hotkey is given up:

    layout                     uniform generator    50% top-1 accuracy
    9 spread (before)                    37.9%              37.9%
    6 concentrated + 3 spread            37.9%              54.5%
    7 concentrated + 2 spread            37.2%              55.9%
    6 + 2 + 1 hedge                      34.3%              51.5%

HEK293 runs width 12 and band 10.1, so eight hotkeys cover 81 seeds and are still below saturation.

Raise CONCENTRATE past the caps here only once the live log in seed_window_model.py clears the ~38%
break-even, because that is where concentration stops being free.

Usage:
    python window_plan.py            # write data/window_plan.json for the next round
    python window_plan.py --dry-run  # print the plan without writing it
"""
import json
import os
import re
import sys
from datetime import datetime, timedelta, timezone

from seed_window_model import (load_log, load_tasks, pending_for, predict, save_log,
                               window as _sw_window,
                               window_label)

PLAN_PATH = "data/window_plan.json"
MINER_SH = "miner.sh"
# Band width per cell type, from the width sweep: band 14.7 of 16 on the erythroid types, 10.1 of
# 12 on HEK293 (168 builds, zero declines).
WIDTH = {"HEK293": 12}
DEFAULT_WIDTH = 16
# Hotkeys concentrated on the window CONC_RANK names. The layout differs by cell type because the
# band width does:
#
#   erythroid (width 16): 6 concentrated + 1 spread on rank 1 + 1 spread on rank 3 + 1 all-cut
#   HEK293    (width 12): 8 concentrated + 1 all-cut, no spread at all
#
# HEK293 takes 8 because its band is 10.1 of 12, so 8 x 12 tiles the whole 100-seed window (81
# band seeds) where 6 would leave 28 seeds uncovered (61 band seeds). That consumes every banded
# hotkey, so HEK293 has no rank-1 or rank-3 exposure -- the concentrated block is the entire bet.
# Set this to 6 to buy that exposure back at the cost of 20 band seeds on the concentrated window.
CONCENTRATE = {"HEK293": 8}
DEFAULT_CONCENTRATE = 6
# Spread hotkeys run the shipped 100-wide window rather than the narrow band width, and sit on the
# NEXT most under-drawn windows rather than their home territories. This deliberately costs
# coverage -- band is 13.00 at width 100 against 14.70/14.60/14.22 at width 16 on
# CD34+/HUDEP-2/K562, and 7.00 against 10.10 at width 12 on HEK293, so each erythroid spread hotkey
# gives up ~9-13% and HEK293's gives up 44% (168-build sweep). Per-hotkey coverage is exactly the
# band count, and P(k=2) = (band/900)**2 whatever the width, so nothing offsets it.
#
# It is kept because it changes one thing at a time: the concentrated block carries the whole
# narrow-window bet while the spread hotkeys stay on the configuration the fleet has always run,
# which is the fallback if the prediction turns out to be worthless.
SPREAD_WIDTH = 100
# Which ranked window the concentrated block tiles, per cell type, as a 1-based rank into
# predict()'s ranked list. 1 = the model's top pick. Any remaining hotkeys take the other ranks in
# order, so with a multi-hotkey fleet at rank 2 the layout is: rank 2 concentrated, rank 1 spread,
# rank 3 spread. With a single hotkey there is no spread — it takes the whole 100-seed window at
# the rank named here and nothing else is covered.
#
# Set per cell type from perrank.py, which measured each rank's *exclusive* hit rate (did that one
# window hold a seed) over 38 scored predictions against the 29.8% chance baseline every individual
# rank shares:
#
#   cell           n      #1      #2      #3      #4
#   CD34+_HSPC    13   46.2%   15.4%   15.4%   46.2%
#   HEK293        13   38.5%   30.8%   15.4%   23.1%
#   HUDEP-2       12   33.3%   50.0%   66.7%    8.3%
#
# So CD34+_HSPC -> 1 and HUDEP-2 -> 3. Read this for what it is: 12-13 rounds per cell, no rank
# significant on its own, and the three cells do not agree on an ordering (HUDEP-2's is inverted).
# It is the best available read, not a proven edge. The aggregate is weakly monotone (39.5% / 31.6%
# / 31.6% / 26.3%), which is why an unmeasured cell type defaults to the top pick.
#
# HEK293 and K562 are not here: they run all-cut instead (miner.sh ALL_CUT_HOTKEYS), because
# HEK293's ranks decline roughly monotonically with no rank worth a lone hotkey and K562 had no
# scored prediction in the window at all.
CELL_RANK = {"CD34+_HSPC": 1, "HUDEP-2": 3}
DEFAULT_RANK = 1
# EXPLICIT MODE, and what actually ships as of 2026-09-08. When this map is non-empty it overrides
# the concentrate/spread layout above entirely: each named hotkey takes the FULL 100-seed window at
# its own ranked position, one hotkey per rank, for every cell type. A hotkey that runs all-cut or
# seed-depend on a given cell is skipped there and its rank simply goes uncovered.
#
#   niome_hotkey1 -> rank 1     niome_hotkey3 -> rank 3
#   niome_hotkey2 -> rank 2     niome_hotkey  -> rank 4  (all-cut on HEK293 instead)
#
# The honest case for this is COVERAGE, not prediction. Per-rank exclusive hit rates over 38 scored
# predictions put rank 1 at 1.33x chance (p=0.130), ranks 2-3 at 1.06x and rank 4 at 0.88x — so
# ranks 2-6 are worth no more than any disjoint window, and that is fine, because six disjoint
# width-100 bands cover 6 x 13 = 78 seeds of 900 and spike on 1-(1-78/900)**3 = 23.8% of rounds
# whether the ordering means anything or not (five windows on K562, where h0 plays all-cut
# instead: 65 seeds, 20.1%). Only rank 1 carries any edge, and h1 holds it.
#
# Width is deliberately the full 100 rather than the 16-seed band width used when several hotkeys
# tile one window. The measured trade is real but small and points the other way: band is 14.70 at
# width 16 against 13.00 at width 100, so four narrow slices would cover 58.8 seeds against 52
# (+13%) — at the cost of leaving 84% of each predicted window uncovered, which throws away
# whatever rank 1 is worth. Switch by giving RANK_WIDTH a value.
RANK_BY_HOTKEY = {"niome_hotkey1": 1, "niome_hotkey2": 2, "niome_hotkey3": 3,
                  "niome_hotkey4": 4, "niome_hotkey5": 5, "niome_hotkey": 6}
RANK_WIDTH = 100

# FIXED MODE — literal per-hotkey windows, overriding both layouts above and the prediction
# entirely. Windows here MAY OVERLAP: the prediction is not consulted, so there is no ranked
# ordering to preserve, and the object that spikes is the ~9-13 seed clean BAND found inside the
# window, not the window itself. Two hotkeys screening overlapping ranges still min-union
# different guide groups and land their bands on different seeds.
#
# What overlap does NOT do is buy coverage. Coverage is `sum(band sizes)` over the fleet, and band
# size falls with window width — measured 13 / 12 / 11 / 9 at widths 100 / 150 / 200 / 300,
# because band size is set by how many rows must agree, not by how wide a range was searched
# (CLAUDE.md, "wider screening window -> wider band", falsified). Measured live at width 300 the
# band is **11-12**, not the 9 on record, so seven width-300 hotkeys cover ~80 band seeds less
# ~2.5 expected collisions in the overlaps = ~77.5, and spike on 1-(1-77.5/900)**3 = 23.7% of
# rounds -- level with the six width-100 disjoint windows' 78 seeds and 23.8%. Hotkey COUNT is
# what buys coverage here, not width: five width-300 hotkeys managed only ~57 seeds (17.8%).
# Set this to {} to return to RANK_BY_HOTKEY.
# CONCENTRATED STRIDE-10 LAYOUT (operator's choice, 2026-09-09). Eleven width-200 windows at
# stride 10 over the 300-seed span 100-399, so h_n covers 100+10n .. 299+10n.
#
# **Measured, so the trade is known.** Bands are effectively INDEPENDENT even at stride 10 —
# adjacent windows share 190 of 200 seeds, yet measured pairwise band overlap is 0.56 seeds against
# an independent-draw expectation of 0.58 (band_overlap.py, K562 1eca4bf1) — because a 10-seed shift
# changes which guides survive the max_fail tail cut, so FastGreedy lands elsewhere. The
# re-correlation failure mode does NOT occur here.
#
# What concentration does cost is room: 11 bands of 12 is 132 seeds, and 132 cannot sit distinctly
# in a 300-seed span, so collisions are forced by pigeonhole rather than by correlation.
#
#   layout                          hotkeys  span  union  P(>=1 of 3 seeds)
#   spread stride 70 over 900            11   900    126        36.4%
#   concentrated stride 10               11   300    104        30.8%
#   current width-300 stride 100          7   900    ~80        24.3%
#   concentrated stride 10                7   260    ~73        22.6%
#
# So this needs ALL ELEVEN slots to beat what it replaces; at seven it is slightly worse. Only h0-h7
# are registered (h8 is keyed but not on the metagraph; h9/h10 have no key at all), and h0 runs
# all-cut, so the seven live band hotkeys are spread across the eleven offsets instead of taking the
# first seven — same layout, best coverage available until the rest register. Set
# `SPREAD_OVER_900 = True` below for the 36.4% variant.
FIXED_WINDOWS = {
    "niome_hotkey1": [100, 299],   # offset 0
    "niome_hotkey2": [120, 319],   # offset 20
    "niome_hotkey3": [130, 329],   # offset 30
    "niome_hotkey4": [150, 349],   # offset 50
    "niome_hotkey5": [170, 369],   # offset 70
    "niome_hotkey6": [180, 379],   # offset 80
    "niome_hotkey7": [200, 399],   # offset 100
    # ---- unfilled slots, assigned as hotkeys register (completes the stride-10 tiling) ----
    "niome_hotkey8": [110, 309],   # offset 10
    "niome_hotkey9": [140, 339],   # offset 40
    "niome_hotkey10": [160, 359],  # offset 60
    "niome_hotkey11": [190, 389],  # offset 90
}
TTL_HOURS = 6                      # survives a missed cron run, expires before it misleads


# ---------------------------------------------------------------------------------------------
# JOINED-WINDOW MODE — the operator's rank-frequency scheme.
#
# Per cell type, over all previous three-seed rounds: count how often each of the nine width-100
# windows was drawn; rank the nine (ties -> most recently updated first, then lower index); record
# the rank of each round's three seed windows; count how often each rank appears; take the three
# most frequent ranks (ties -> most recently updated position, then lower rank); and map those back
# to the windows now holding them. Those three windows are the joined 300-seed band space, and the
# eleven hotkeys take width-200 slices of it at stride 10 in sorted-seed order.
#
# **The prediction has no measured skill** (rank_freq.py, 108 out-of-sample rounds: 102 of 324 seeds
# landed in the predicted window against 108.0 +/- 8.5 expected, z = -0.71; per cell +0.45 / -0.88 /
# +0.00 / -1.07). Seeds are uniform, so counting them cannot forecast them. Band position is free,
# so this is EV-equivalent to any other 300-seed space -- it is kept because it costs nothing and
# would start paying by itself if the generator ever stopped being uniform.
#
# What the layout costs is separate and measured (band_overlap.py): eleven bands of 12 is 132 seeds,
# which cannot sit distinctly in 300, so union is 104 (30.8% of rounds) against 126 (36.4%) for the
# same eleven hotkeys spread over the full 900. Bands themselves stay independent even at stride 10
# (pairwise overlap 0.56 measured against 0.58 expected), so the loss is pigeonhole, not
# correlation. Set JOINED_MODE = False to fall back to FIXED_WINDOWS.
JOINED_MODE = True
# Ten band hotkeys ROTATE around the joined space at stride 30, width 200 — the slice is circular,
# so a hotkey past the end wraps to the front. 10 x 30 = 300 tiles the space exactly once, and every
# seed sits in 200/300 x 10 = 6.67 slices on average. h0 stays on all-cut for every cell type.
JOINED_HK = ["niome_hotkey1", "niome_hotkey2", "niome_hotkey3", "niome_hotkey4",
             "niome_hotkey5", "niome_hotkey6", "niome_hotkey7", "niome_hotkey8",
             "niome_hotkey9", "niome_hotkey10"]
JOINED_SUB_WIDTH = 200
JOINED_STRIDE = 30


def _rank_freq_windows(history, cell):
    """The three windows holding the three most frequent ranks, per the scheme above."""
    counts = [0] * 9
    last_upd = [-1] * 9
    rank_counts = {k: 0 for k in range(10)}
    rank_last = {k: -1 for k in range(10)}
    rows = [t for t in history if t["cell"] == cell]
    for r, t in enumerate(rows):
        order = sorted(range(9), key=lambda w: (-counts[w], -last_upd[w], w))
        rank_of = {w: i + 1 for i, w in enumerate(order)}
        ranks_this = [0, 0, 0] if r == 0 else [rank_of[_sw_window(s)] for s in t["seeds"]]
        for k in ranks_this:
            rank_counts[k] = rank_counts.get(k, 0) + 1
            rank_last[k] = r
        for s in t["seeds"]:
            w = _sw_window(s)
            counts[w] += 1
            last_upd[w] = r
    order = sorted(range(9), key=lambda w: (-counts[w], -last_upd[w], w))
    rank_of = {w: i + 1 for i, w in enumerate(order)}
    top3 = sorted(range(1, 10), key=lambda k: (-rank_counts.get(k, 0), -rank_last.get(k, -1), k))[:3]
    holder = {rk: w for w, rk in rank_of.items()}
    return sorted(holder[k] for k in top3 if k in holder), top3


def _joined_slices(windows):
    """Rotated width-`JOINED_SUB_WIDTH` slices at `JOINED_STRIDE` over the sorted joined seeds.

    The slice is CIRCULAR: hotkey h takes indices (h*stride + k) mod n for k < width, so a slice
    running past the end wraps to the front and every hotkey gets exactly `width` seeds. Returned as
    ascending, non-overlapping [lo, hi] ranges, which is what the plan stores and
    `Miner._window_for` validates.
    """
    seeds = sorted(s for w in windows for s in range(w * 100 + 100, w * 100 + 200))
    n = len(seeds)
    width = min(JOINED_SUB_WIDTH, n)
    out = []
    for h in range(len(JOINED_HK)):
        off = (h * JOINED_STRIDE) % max(1, n)
        chunk = sorted({seeds[(off + k) % n] for k in range(width)})
        spans, lo, prev = [], chunk[0], chunk[0]
        for x in chunk[1:]:
            if x != prev + 1:
                spans.append([lo, prev])
                lo = x
            prev = x
        spans.append([lo, prev])
        out.append(spans)
    return out

def _sh_list(var):
    """A space-separated bash list from miner.sh, so the two files cannot drift."""
    for line in open(MINER_SH):
        m = re.match(rf'\s*{var}=(.*)', line)
        if m:
            return [x for x in m.group(1).strip().strip('"\'').split() if x]
    return []


def hedge_spec():
    """{instance: frozenset(cell types) or None} for the all-cut hotkeys, from ALL_CUT_HOTKEYS.

    ``None`` means every cell type — miner.sh's bare "<hotkey>" form. The "<hotkey>:CELL,CELL" form
    hedges only those cell types and keeps the band elsewhere, which is what a one-hotkey fleet
    runs: the same hotkey is the band bet on the cells where the rank ordering has an edge and the
    flat hedge on the cells where it does not.

    A hotkey has no clean band on the cells it hedges — all-cut pins is_cut across the whole
    900-seed window — so assigning it a window there would waste a window another hotkey could be
    covering, and would make the log read as though a band had existed.
    """
    out = {}
    for token in _sh_list("ALL_CUT_HOTKEYS"):
        name, _, cells = token.partition(":")
        out[name] = frozenset(c for c in cells.split(",") if c) or None
    return out


def hedges_for(cell, spec):
    """The instances running all-cut for this cell type, so they get no window in this round."""
    return {n for n, cells in spec.items() if cells is None or cell in cells}


def seed_depend_hotkeys():
    """Instances pinned to seed 0 by miner.sh's SEED_DEPEND_VARIANTS ("<hotkey>:<variant>").

    Seed-depend is the first rung of the miner's ladder and replaces the construction outright, so
    such a hotkey has no clean band on any cell type. It is excluded for the same reason an all-cut
    hedge is: a window assigned to it is a window nothing covers, and the log would read as though
    a band had been played there.
    """
    return {token.partition(":")[0] for token in _sh_list("SEED_DEPEND_VARIANTS")}


def table():
    """Every (instance, default window) row of miner.sh's HOTKEYS table, unfiltered."""
    rows = []
    for line in open(MINER_SH):
        m = re.match(r'\s*"(\S+)\s+\d+\s+\d+\s+(\d+)-(\d+)"', line)
        if m:
            rows.append((m.group(1), (int(m.group(2)), int(m.group(3)))))
    return rows


def fleet():
    """(instance, default window) in table order, parsed from miner.sh so the two cannot drift.

    Hotkeys hedged on *every* cell type are dropped here: they need no window at all, and leaving
    them in would consume a concentrated slot or a spread window for a build that ignores it. A
    hotkey hedged on only some cell types stays — main() drops it per cell instead.
    """
    skip = ({n for n, cells in hedge_spec().items() if cells is None}
            | set(_sh_list("DEREGISTERED")) | seed_depend_hotkeys())
    rows = [r for r in table() if r[0] not in skip]
    # Concentrated slots are filled from the hotkeys NOT preferred for spread, so a hotkey freed by
    # a deregistration elsewhere in the table lands in the concentrated block rather than pushing
    # a designated spread hotkey out of its role. Spread hotkeys are still used for concentration
    # when the block needs them (HEK293 concentrates 8 and consumes both).
    spread = _sh_list("SPREAD_HOTKEYS")
    rows.sort(key=lambda r: (spread.index(r[0]) + 1) if r[0] in spread else 0)
    return rows


def tile(lo, hi, width, n):
    """n disjoint sub-windows covering [lo, hi] as completely as n slots allow.

    The remainder is spread one seed at a time across the leading slots rather than dropped off the
    end. Packing from ``lo`` at a fixed ``width`` left the tail of the window uncovered -- 8 x 12 =
    96 of 100 on HEK293, 6 x 16 = 96 on the erythroid types -- so a seed landing in the last four
    seeds scored the floor even when the prediction was right. That is not hypothetical: on
    678cf369 the model called 300-399, seed 397 landed in it, and the block tiled only 300-395.

    Widening is free. The width sweep found band flat at 14-15 from width 16 through 50 (168
    builds), so a 17-seed window bands the same as a 16-seed one -- the extra seed is covered at no
    cost to any hotkey's band.

    ``width`` is now a *minimum*: if n slots of that width would overflow the window, the old
    fixed-width packing is kept and the slots that do not fit are dropped.
    """
    span = hi - lo + 1
    if n <= 0 or span <= 0:
        return []
    base, rem = divmod(span, n)
    if base < width:
        # More hotkeys than the window can give `width` each -- keep the intended band width and
        # let the surplus slots fall off, as before.
        return [(lo + i * width, lo + i * width + width - 1)
                for i in range(n) if lo + i * width + width - 1 <= hi]
    out, a = [], lo
    for i in range(n):
        b = a + base + (1 if i < rem else 0) - 1
        out.append((a, min(b, hi)))
        a = b + 1
    return out


def main():
    members = fleet()
    if not table():
        print(f"no HOTKEYS table found in {MINER_SH}", file=sys.stderr)
        return 1
    # An empty ``members`` is not an error: every hotkey may be playing something that has no band
    # (seed-depend, or all-cut on every cell). The plan is still written -- with no assignments --
    # and seed_window_model.py still resolves and emits predictions, so the shadow log keeps
    # accruing the evidence that decides whether a band is worth going back to.
    rows = load_tasks()
    log = load_log()
    now = datetime.now(timezone.utc)
    plan = {"generated_at": now.isoformat(),
            "expires_at": (now + timedelta(hours=TTL_HOURS)).isoformat(),
            "history_through": rows[-1]["at"], "n_tasks": len(rows), "assignments": {}}

    spec = hedge_spec()
    sd = seed_depend_hotkeys()
    hedge_note = ", ".join(f"{n}:{'all cells' if c is None else '/'.join(sorted(c))}"
                           for n, c in sorted(spec.items()))
    print(f"{len(rows)} tasks through {rows[-1]['at'][:16]}; {len(members)} banded hotkeys"
          + (f"; all-cut {hedge_note}" if spec else "")
          + (f"; seed-depend {', '.join(sorted(sd))}" if sd else "") + "\n")
    for cell in sorted({t["cell"] for t in rows}):
        pr = predict(rows, cell)
        hedged = hedges_for(cell, spec)
        # Hotkeys that run all-cut on *this* cell type play no window here, so they are out of the
        # allocation for this cell only and back in it for the next.
        avail = [m for m in members if m[0] not in hedged]
        if JOINED_MODE:
            wins3, top3 = _rank_freq_windows(rows, cell)
            slices = _joined_slices(wins3)
            assign, slots = {}, []
            for name, _default in avail:
                if name in JOINED_HK:
                    assign[name] = slices[JOINED_HK.index(name)]
            top = pr["ranked"][0]
            width = JOINED_SUB_WIDTH
            rank = 0
        elif FIXED_WINDOWS:
            # Fixed mode: literal windows, same for every cell type, overlap permitted.
            assign, slots = {}, []
            for name, _default in avail:
                w = FIXED_WINDOWS.get(name)
                if w:
                    assign[name] = list(w)
            top = pr["ranked"][0]
            width = None
            rank = 0            # no ranked position is consulted in fixed mode
        elif RANK_BY_HOTKEY:
            # Explicit mode: one hotkey per ranked window, same mapping for every cell type.
            assign, slots = {}, []
            for name, _default in avail:
                want = RANK_BY_HOTKEY.get(name)
                if not want or want > len(pr["ranked"]):
                    continue
                a = (pr["ranked"][want - 1]["window"] + 1) * 100
                assign[name] = [a, a + RANK_WIDTH - 1]
            rank = min(RANK_BY_HOTKEY.values()) - 1 if assign else 0
            top = pr["ranked"][rank]
            width = RANK_WIDTH
        else:
            rank = min(CELL_RANK.get(cell, DEFAULT_RANK), len(pr["ranked"])) - 1
            top = pr["ranked"][rank]           # the window the concentrated block tiles
            width = WIDTH.get(cell, DEFAULT_WIDTH)
            conc = min(CONCENTRATE.get(cell, DEFAULT_CONCENTRATE), len(avail))
            lo = (top["window"] + 1) * 100
            slots = tile(lo, lo + 99, width, conc)
            assign = {}
            for (name, _default), w in zip(avail[:len(slots)], slots):
                assign[name] = list(w)
            # The rest take the other ranked windows in order, full width, skipping the one the
            # concentrated block already owns -- so a spread hotkey can never overlap it.
            rest = avail[len(slots):]
            spread_ranks = [r for i, r in enumerate(pr["ranked"]) if i != rank]
            for (name, _default), nxt in zip(rest, spread_ranks[:len(rest)]):
                a = (nxt["window"] + 1) * 100
                assign[name] = [a, a + SPREAD_WIDTH - 1]

        spans = sorted(assign.values())
        if not (FIXED_WINDOWS or JOINED_MODE):
            # Every other layout is built to be disjoint, so an overlap there is a bug that would
            # silently re-correlate two siblings and cost the whole decorrelation. In fixed mode
            # the operator asked for the overlap, so it is reported and kept.
            for (a1, b1), (a2, b2) in zip(spans, spans[1:]):
                if a2 <= b1:
                    print(f"ERROR: {cell} windows {a1}-{b1} and {a2}-{b2} overlap",
                          file=sys.stderr)
                    return 1
        plan["assignments"][cell] = assign
        # Record what was actually applied against this cell's pending prediction, so a resolved
        # round in seed_window_log.json shows predicted window, applied per-hotkey windows, the
        # task id and each hotkey's scored outcome together. Cells where every hotkey hedges still
        # get the block, with no assignments -- the prediction was live and its accuracy is still
        # scored, it just was not played.
        if JOINED_MODE:
            roles = {n: "joined" for n, _d in avail if n in JOINED_HK}
        elif FIXED_WINDOWS:
            roles = {n: "fixed" for n, _d in avail if n in FIXED_WINDOWS}
        elif RANK_BY_HOTKEY:
            roles = {n: f"rank{RANK_BY_HOTKEY[n]}" for n, _d in avail if n in RANK_BY_HOTKEY}
        else:
            roles = {n: ("concentrated" if i < len(slots) else "spread")
                     for i, (n, _d) in enumerate(avail)}
        roles.update({n: "all-cut-hedge" for n in hedged})
        entry = pending_for(log, cell)
        if entry is not None:
            entry["applied"] = {
                "applied_at": now.isoformat(), "predicted_window": top["label"],
                # The width actually assigned, not the configured minimum: tile() widens the slots
                # to cover the window, so a lone hotkey gets the full 100 rather than WIDTH's 16.
                "predicted_p_hit": top["p_hit"],
                "band_width": (slots[0][1] - slots[0][0] + 1) if slots else 0,
                "conc_rank": rank + 1, "rank1_window": pr["ranked"][0]["label"],
                "spread_width": SPREAD_WIDTH, "n_concentrated": len(slots),
                "played": bool(assign),
                "assignments": {n: list(v) for n, v in assign.items()}, "roles": roles}
        head = (f"  {cell:<11} rank1 {pr['ranked'][0]['label']} "
                f"({pr['ranked'][0]['p_hit']:.0%})  | ")
        if not assign:
            why = ", ".join(filter(None, [
                f"all-cut {'/'.join(sorted(hedged))}" if hedged else "",
                f"seed-depend {'/'.join(sorted(sd))}" if sd else ""]))
            print(head + f"no window played ({why or 'no banded hotkey'})")
            continue
        ranks = {(r["window"] + 1) * 100: i for i, r in enumerate(pr["ranked"])}
        if JOINED_MODE:
            wins3, top3 = _rank_freq_windows(rows, cell)
            covered = len({s for spans in assign.values() for a, b in spans
                           for s in range(a, b + 1)})
            print(head + f"beta {pr['beta']:+.2f}  JOINED ranks {top3} -> "
                  + ",".join(f"{w*100+100}-{w*100+199}" for w in wins3)
                  + f"  x{len(assign)} hotkeys, width {JOINED_SUB_WIDTH} stride {JOINED_STRIDE}"
                  f" rotated"
                  + f", span {covered} of 900"
                  + (f"  | all-cut {'/'.join(sorted(hedged))}" if hedged else ""))
            for name in sorted(assign, key=lambda n: JOINED_HK.index(n)):
                print(f"    joined {name:<15} "
                      + ",".join(f"{a}-{b}" for a, b in assign[name]))
            continue
        if FIXED_WINDOWS:
            spans = sorted(assign.values())
            olap = sum(max(0, min(b1, b2) - max(a1, a2) + 1)
                       for (a1, b1), (a2, b2) in zip(spans, spans[1:]))
            covered = len({s for a, b in spans for s in range(a, b + 1)})
            print(head + f"beta {pr['beta']:+.2f}  FIXED windows x{len(assign)}"
                  f"  span {covered} of 900, {olap} seeds overlapped"
                  + (f"  | all-cut {'/'.join(sorted(hedged))}" if hedged else ""))
            for name in sorted(assign, key=lambda n: assign[n][0]):
                lo_, hi_ = assign[name]
                print(f"    fixed  {name:<14} {lo_}-{hi_}  (width {hi_ - lo_ + 1})")
            continue
        if RANK_BY_HOTKEY:
            print(head + f"beta {pr['beta']:+.2f}  width {RANK_WIDTH} x{len(assign)}"
                  + (f"  | all-cut {'/'.join(sorted(hedged))}" if hedged else ""))
            for name in sorted(assign, key=lambda n: RANK_BY_HOTKEY.get(n, 99)):
                lo_, hi_ = assign[name]
                i = ranks.get(lo_)
                print(f"    rank{RANK_BY_HOTKEY[name]}  {name:<14} {lo_}-{hi_}"
                      + (f"  p_hit {pr['ranked'][i]['p_hit']:.0%}" if i is not None else ""))
            continue
        print(head + f"rank{rank + 1} {top['label']} "
              f"({top['p_hit']:.0%}, beta {pr['beta']:+.2f})  "
              f"width {slots[0][1] - slots[0][0] + 1} x{len(slots)}")
        print(f"    banded: "
              + ", ".join(f"{n}:{a}-{b}" for (n, _d), (a, b) in zip(avail, slots)))
        if len(assign) > len(slots):
            print(f"    spread (w{SPREAD_WIDTH}): "
                  + ", ".join(f"{n}:{v[0]}-{v[1]}(rank {ranks.get(v[0], '?')+1},"
                              f" {pr['ranked'][ranks[v[0]]]['p_hit']:.0%})"
                              for n, v in list(assign.items())[len(slots):]))

    if "--dry-run" in sys.argv:
        print("\n--dry-run: not written")
        return 0
    save_log(log)
    os.makedirs(os.path.dirname(PLAN_PATH), exist_ok=True)
    tmp = PLAN_PATH + ".tmp"
    with open(tmp, "w") as handle:            # atomic: miners read this file mid-round
        json.dump(plan, handle, indent=1)
    os.replace(tmp, PLAN_PATH)
    print(f"\nwrote {PLAN_PATH}, valid until {plan['expires_at'][:16]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
