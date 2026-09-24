#!/usr/bin/env python3
"""joined_window.py — the rotated joined-window layout, shared by the plan and the fallback.

[window_plan.py](window_plan.py) writes `data/window_plan.json` once an hour and the miner reads it
per build. When that file is missing, stale or malformed the miner falls back to its own layout —
and the fallback used to be a single contiguous window per hotkey, taken from `NIOME_HDR_WINDOW`.
That was wrong in two ways at once:

  * **The widths collided.** `miner.sh`'s table listed a width per hotkey, but the validated width
    is per CELL TYPE (`all_hdr.CELL_CONFIG`), and `group_size` is tuned AT that width. The pin now
    supplies only the START.
  * **The starts collided.** h0/h1 both started at 100 and h8/h9/h10 all at 900, so eleven hotkeys
    fell back onto EIGHT distinct windows. Three hotkeys drew the identical band and the fleet lost
    their decorrelation on exactly the rounds the plan was unavailable.

So the fallback now runs the SAME construction as the plan — rotated slices at `STRIDE` over a
joined seed space, width from the cell, h0 at a half-stride offset — against a fixed class triple
instead of a predicted one. Both callers import the arithmetic from here rather than restating it,
which is the same reason `window_plan.py` parses `miner.sh` instead of duplicating its lists.

Nothing here imports torch, the seed model or the chain, so the miner can import it at startup.
"""

# The band hotkeys ROTATE around the joined space; the slice is circular, so one running past the
# end wraps to the front and every hotkey gets exactly `width` seeds. `len(ROTATE_HK) * STRIDE`
# should equal the joined space (3 x 100 = 300) so the offsets tile it exactly once: six hotkeys
# at stride 50 do, as ten at stride 30 did before.
#
# **The fleet has shrunk twice on 2026-09-12**, each time with the stride moved to keep the tiling
# exact, verified against the metagraph both times:
#
#   h0-h3 deregistered  -> 7 hotkeys, stride 30 -> 50 (6 rotating x 50 = 300)
#   h4-h7 deregistered  -> 3 hotkeys, stride 50 -> 100 (3 rotating x 100 = 300)
#
# Only h8/h9/h10 remain, holding uids 151/224/10. Leaving the stride behind the hotkey count packs
# the slices into the front of the space and leaves the tail covered only by the circular wrap.
STRIDE = 100
# 2026-09-15: EMPTY. All four registered hotkeys take the FULL joined space (see `FULL_HK`), and
# their decorrelation comes from the BAND sub-window offset below instead of from a joined-window
# slice. `rotated(..., count=0)` returns [], so the plan simply assigns nothing here.
ROTATE_HK: list[str] = []

# One hotkey sits at HALF a stride, i.e. between the first two rotating slices. At offset 0 its
# slice would be byte-identical to the first one's, and two hotkeys on one window draw correlated
# bands — which throws away the whole point of the layout. This is the same collision the fallback
# used to have three ways over.
#
# **This slot holds the SEED-DEPEND hotkey, and that is deliberate.** Whichever hotkey runs
# seed-depend is excluded from the band plan by `window_plan.py`, so parking it here leaves exactly
# `ROTATE_HK` playing bands -- and 6 x STRIDE 50 tiles the 300-seed joined space once, with no gap
# and nothing doubled up. It still gets a half-stride window from the FALLBACK, which is what it
# builds on if seed-depend declines (budget under SEED_DEPEND_MIN_BUDGET_S on a failed prefetch).
# h0 held the slot, then h4, then h10, then it was empty on 2026-09-12. As of 2026-09-13 **h0 holds
# it again**: h0 is the only registered hotkey and it runs seed-depend, so it is excluded from the
# band plan and parking it here costs no rotating slice. It never uses this window on a normal
# round -- seed-depend replaces the construction -- but it builds on it when seed-depend declines
# (budget under SEED_DEPEND_MIN_BUDGET_S of 360s on a failed prefetch), which is the one path that
# would otherwise fall through to miner.sh's last-resort NIOME_HDR_WINDOW column.
#
# ROTATE_HK below is kept pointing at h8/h9/h10 even though all three are deregistered and their
# pm2 apps are deleted: `fallback_for` is keyed on NIOME_INSTANCE, so entries for processes that do
# not exist are inert, and keeping them means the layout is already correct if they re-register.
#
# 2026-09-14: h0 came OFF seed-depend and onto all-HDR, so this slot is now a playing band rather
# than a parking space, and it is the fleet's only one.
# 2026-09-15: h0-h3 are registered (uids 122/41/163/118) and ALL FOUR sit here, i.e. every hotkey
# min-unions on cut over the whole 300-seed joined space. They are NOT correlated by that: the
# conjunction's clean set is a function of the joined space, but its BAND is drawn from a width-150
# sub-window whose offset rotates per hotkey at `BAND_STRIDE` (see `band_offset_frac`). Two hotkeys
# therefore share a clean set and hold different bands, which is where the spike lives.
#
# 2026-09-17: h4 and h5 came back onto the metagraph (uids 190/234, verified against the chain, not
# against miner.sh's stale DEREGISTERED line) and were folded into this layout — six hotkeys now, at
# `BAND_STRIDE` 50 instead of 75 (see below). Set by operator request alongside a second, bigger
# change: for these six hotkeys the CUT computation itself no longer shares the plan's 300-seed
# joined space at all. `conjunction_cut_seeds()` gives them the full 100-999 (900 seeds) regardless
# of what the plan predicts, while the BAND still draws from the plan's 300-seed space exactly as
# before — a hotkey in `BAND_HK` now solves a fundamentally wider cut (more Cas12a pool, likely a
# larger clean set) while keeping the band question ("is this seed inside one of the 3 predicted
# classes") unchanged. This decoupling is new and unmeasured: every previous k/group/light_cell_rows
# row in `conjunction.CELL_CONFIG` was tuned with the cut window EQUAL to the band's 300-seed space,
# and a 900-seed cut bank has a different composition (different `cas12a_max_fail` scaling — it now
# reverts to the cell's native all-cut value instead of being scaled down by seeds/900) that could
# shift where `choose_band`'s formation wall sits. Re-measure before trusting `CELL_CONFIG`'s numbers
# as anything more than the operator's chosen starting point.
FULL_HK: list[str] = ["niome_hotkey", "niome_hotkey1", "niome_hotkey2", "niome_hotkey3",
                      "niome_hotkey4", "niome_hotkey5", "niome_hotkey6", "niome_hotkey7",
                      "niome_hotkey8", "niome_hotkey9"]

# --- the BAND sub-window, which is what decorrelates the fleet now -------------------------------
#
# `conjunction.sub_window` draws the k-seed band from `band_width` seeds of the joined space
# starting `band_offset_frac` through it and wrapping. With every hotkey on the same joined window
# that fraction is the ONLY thing separating two siblings' bands, so it is the fleet layout and
# belongs here rather than in the per-cell config.
#
# `conjunction.py`'s module docstring flags this explicitly: the band is a deterministic function of
# (contract, joined space, width), so two hotkeys handed the same window build the SAME band. That
# was safe while the fleet was one hotkey (n=1 makes the per-hotkey number the fleet number) and
# stops being safe the moment it regrows -- which is what this is.
#
# Four hotkeys at BAND_STRIDE 75 over 300 seeds put the sub-windows at 0 / 75 / 150 / 225. At
# BAND_SUB_WIDTH 150 each covers half the space, so each seed sits in exactly 2 of the 4 -- a
# deliberate 2x overlap, since 4 x 75 = 300 tiles the space once while the width is double the
# stride. Disjoint bands would need width 75, which is not the width the arms were tuned at.
#
# 2026-09-17: six hotkeys (h0-h5) at BAND_STRIDE 50 -- offsets 0/50/100/150/200/250 -- so 6 x 50 =
# 300 still tiles the space exactly once, at BAND_SUB_WIDTH 150 (unchanged) each seed now sits in
# 150/50 = 3 of the 6 slices on average, one more of overlap than the four-hotkey layout.
#
# 2026-09-18 (morning): TEN hotkeys (h0-h9) at BAND_STRIDE 30, all on the predicted 300.
#
# **2026-09-18: the band space is the FULL 900 and the ten hotkeys are SPLIT over it.** Until now
# every band hotkey drew its sub-window from the plan's three predicted classes, so the other 600
# seeds held no band on any hotkey at all -- a round drawing outside the prediction could not spike
# anywhere in the fleet. Two groups now cover the whole space, each tiling its own half exactly
# once:
#
#   BAND_HK  h0-h5  6 x STRIDE  50 = 300   the three PREDICTED classes, 3x overlap at width 150
#   REST_HK  h6-h9  4 x STRIDE 150 = 600   everything the plan did NOT predict, DISJOINT at width 150
#
# `count * stride == span` is the invariant, not a nicety: `band_offset_frac` is
# `(index * stride) % span`, so a count and stride that do not multiply to the span wrap the
# sequence onto itself and two hotkeys draw the identical band from the identical pool. That is the
# collision `BAND_OFFSET_OVERRIDES` was patching by hand at ten hotkeys x stride 50, and it throws
# away the single thing this layout exists to provide. Both groups satisfy it exactly.
#
# **The asymmetry between the two groups is deliberate.** h0-h5 stay concentrated on the predicted
# region at 3x overlap (width 150 / stride 50), which is the layout CLAUDE.md prices and the one the
# 12-contract replication was measured near. h6-h9 spread over twice as much space at ZERO overlap
# (width 150 / stride 150), which is the cheapest possible hedge against the prediction being
# worthless -- and CLAUDE.md measures it as worthless: every strategy in `strategy_rank`'s pool sits
# between 0.93x and 1.06x chance with |z| <= 1.1 over the 160-task cold walk-forward, and SeedFormer
# learned to emit the uniform distribution. Under a uniform generator band position is FREE, so the
# 600-seed half is worth exactly the seeds it holds; nothing is given up by putting four hotkeys
# there, and a round that draws outside the prediction stops being a fleet-wide miss.
#
# What each group costs in sibling correlation, at BAND_SUB_WIDTH 150: inside h0-h5 adjacent
# hotkeys share 100 of 150 candidate seeds (h0 0-149 against h1 50-199), which is LESS correlated
# than the stride-30 layout this replaces (120 of 150). Inside h6-h9 they share NOTHING -- the four
# sub-windows partition the 600. Whether overlapping candidate windows make the resulting BANDS
# correlated is still unmeasured; `widecut_price.py` found bands essentially disjoint (0-1 of 11)
# when the POOLS differed, but within a group here the pool is identical per cell and only the
# window moves. `band_overlap.py` is the tool, and h6-h9 are now the control arm for it.
BAND_SUB_WIDTH = 300

# **2026-09-20: seven hotkeys, width 300 at stride 100, spanning the WHOLE 100-999.** h1/h2/h3 were
# deregistered from netuid 55 overnight (21:11-22:02 on 09-19, verified against the chain at block
# 9106585, not against the logs) and h10 has been gone for longer, so the fleet is h0 plus h4-h9.
# Their pm2 apps are STOPPED rather than deleted.
#
# The layout that replaces the 10 x 30 tiling is not a rotation inside the plan's predicted classes
# at all -- the band space is the full 900 and each hotkey takes a 300-seed CONTIGUOUS window, the
# windows stepping 100 seeds apart:
#
#     h0 100-399   h4 200-499   h5 300-599   h6 400-699   h7 500-799   h8 600-899   h9 700-999
#
# **The invariant is different from every previous layout here and that matters.** The old one was
# `count * stride == span`, which is what a CIRCULAR tiling of non-overlapping slices needs. This
# one is deliberately overlapping and must NOT wrap, so the condition is
# `(count - 1) * stride + width == span` -- 6 x 100 + 300 = 900 exactly. Get it wrong in the other
# direction and the last window runs off the end and `sub_window` wraps it back to seed 100, which
# would silently re-correlate h9 with h0.
#
# Coverage is the full 900 by construction, but it is NOT uniform, and that is inherent to fitting
# seven width-300 windows onto 900 seeds rather than a defect: seeds 100-199 sit in h0's window
# alone and 900-999 in h9's alone, while everything from 300-799 sits in three windows. So the
# edges carry a third of the band mass of the middle. The alternative -- disjoint width-128 windows
# -- would tile evenly but 900/7 is not an integer and the narrower window is a different arm.
# **2026-09-21, operator request: ONE shared window, and the fleet decorrelates on the LOOP AXIS.**
#
# Every previous layout here separated siblings by giving each a different sub-WINDOW -- rotating
# offsets, then contiguous width-300 slices across 900, then a disjoint width-100 tier. All ten
# hotkeys now share the IDENTICAL window and the identical cut, and what separates their bands is
# `ConjunctionConfig.band_loop`: hotkey i plays loop i+1 of the re-banding chain, excluding every
# seed the loops before it banded.
#
#     h0 loop 1   h1 loop 2   h2 loop 3   h3 loop 4   h4 loop 5
#     h5 loop 6   h6 loop 7   h7 loop 8   h8 loop 9   h9 loop 10
#
# **Why this is strictly better than window offsets at separating bands.** Offsets made two bands
# *probably* different -- overlapping candidate windows, measured at 0-1 shared of 11. The loop
# axis makes them *provably* disjoint: loop L bands from candidates with loops 1..L-1's seeds
# removed, so no seed can appear in two hotkeys' bands. Union is exactly `10 * band_k` -- 110 of
# 300 on the erythroid cells, 80 on HEK293 -- which `loop_axis.py` confirmed on 40 of 40 tasks.
#
# **And it is what makes ONE bank possible again.** `bank_key` folds in the seed list, so the
# per-hotkey cut windows of 2026-09-20 meant SEVEN banks per (contract, cell): seven concurrent
# cold Cas12a scans on one GPU at 566-580s each, against a shared bank's 1 x ~410s cold plus warm
# siblings -- and it cost h4 a submission on the first round it ran. One shared window means one
# `bank_key` per (contract, cell), so exactly one hotkey pays the scan and nine load warm. That is
# the whole reason the loop axis can afford the narrow cut the floor measurement wants.
#
# The cost is a chain dependency the offset layout did not have: loop 10 reproduces nine earlier
# `choose_band` calls before its own. They are greedy passes over the already-computed compliance
# matrix -- no extra bank scan, no GPU -- and the measured chains ran all ten in seconds.
# **2026-09-21 (evening): h4 and h5 REMOVED — eight live hotkeys on loops 1-8 contiguously.**
# Both were deregistered from netuid 55 (verified against the chain at block 9117209, not against
# the logs): h4 was crash-looping at 77 restarts on `Wallet ... is not registered`, h5 was a zombie
# -- `check_registered` calls `exit()` but the process survived and kept prefetching for rounds it
# could never be called for, burning GPU. Both pm2 apps are STOPPED.
#
# **Removing them from this list is what reassigns the loops, and leaving them in is not harmless.**
# `band_loop` is `LOOP_HK.index(instance) + 1`, so with the dead pair still present the live fleet
# played loops 1,2,3,4,7,8,9,10 -- loops 5 and 6 were built by nobody, and the band union fell from
# 10*k to 8*k with a HOLE in the middle of the chain. The chain itself stays correct either way
# (each hotkey re-derives its predecessors independently), so this costs coverage, not correctness:
# 110 -> 88 of 300 on the erythroid cells, ~20% of band-hit frequency.
#
#     h0 loop 1   h1 loop 2   h2 loop 3   h3 loop 4
#     h6 loop 5   h7 loop 6   h8 loop 7   h9 loop 8
#
# Note the loop count is now 8, not 10, so the deepest chain re-derives 7 predecessors rather than
# 9 -- marginally cheaper, and the seeds loops 9-10 used to take are simply not banded by anyone.
# Re-add a hotkey here the moment it re-registers; `pm2 start miner-h4 miner-h5` brings the
# processes back, but without a line here they would play no loop at all.
# **2026-09-21 (late): h0 ALSO deregistered — seven live hotkeys, and the loop index is now
# PER CELL TYPE.** h0 lost uid 248 at 19:21 (chain-verified at block 9118707 by ss58, not uid) and
# went the zombie route rather than the crash-loop one: `check_registered` calls `exit()` but the
# process survived and was still prefetching at 20:32, building 250 rows for a round it can never
# be called for. h4/h5 went at ~16:42. Four of ten in nine hours.
# **2026-09-22 00:34/03:50: h6 OUT, h0 BACK IN — a straight swap at the same slot.** h6
# deregistered at 00:34 (which is why no validator called it for round dcbe82b5 — not the
# "isolated miss" it first looked like; its call history and axon were both fine and neither is
# evidence of registration). h0 re-registered with uid **226**, NOT its old 248 — which is exactly
# why this file insists on ss58 over uid.
#
# h0 takes h6's INDEX rather than being prepended, so every other hotkey's loop is untouched on
# both cell groups and only this one slot changes hands. No restart is needed for the other six:
# `band_loop` is `LOOP_ASSIGN[cell][LOOP_HK.index(instance)]` and their indices are unchanged.
# The DEFAULT (erythroid) conjunction list. Per cell, use `conj_hk(cell)` -- HEK293 runs nine.
LOOP_HK: list[str] = ["niome_hotkey", "niome_hotkey1", "niome_hotkey2", "niome_hotkey3",
                      "niome_hotkey4", "niome_hotkey5", "niome_hotkey6"]

# **2026-09-22: the four live hotkeys run the ALL-HDR SPIKE on four windows of a 200-seed space.**
# Set by operator request. `Miner.CONJUNCTION` is False alongside this, so all-HDR is the TOP rung
# rather than the one below the conjunction, and this is the band space it takes.
#
# The space is 200-299 + 400-499 (classes 1 and 3) -- 200 seeds, JOINED, and nothing outside it is
# considered. Four windows of width 100 at stride 50 satisfy the CIRCULAR invariant
# (`count * stride == span`, 4 * 50 = 200), so they wrap the space exactly once:
#
#     h0  200-249     h1  250-299     h2  400-449     h3  450-499
#
# **2026-09-22 (second change): width 100 -> 50 at stride 50.** Operator request. With
# `width == stride` the four windows PARTITION the space rather than overlapping it, so union is
# 200/200 with every seed in exactly ONE window and the four bands are disjoint BY CONSTRUCTION --
# they are drawn from disjoint candidate sets, so no two hotkeys can duplicate a band seat.
#
# The measurement behind it (`hdr_tile.py`, 6 K562 rounds x 10 windows per width, band enumerated
# from the shipped rows): all-HDR's clean band barely moves with window width -- 14-15 at width 20,
# 14 at width 50, 13 at width 100 -- so what the width really sets is band DENSITY inside the
# window: 72% at 20, 28% at 50, 13% at 100. Width 100 at 5x overlap finished LAST of the three arms
# tested (1928.8 summed round final against width 20's 2329.1) with three complete shutouts, because
# five correlated shots at 13% lose to one shot at a high density. Width 50 is the middle point:
# 4 x 14 = 56 disjoint band seats of 200, against width 100's 52 overlapping ones.
#
# Pinned here rather than read from `window_plan.json`, for the same reason `SHARED_CLASSES` is:
# `round_plan.sh` rewrites the plan hourly, and a band space that moves between rounds re-correlates
# the siblings this layout exists to separate.
#
# Width 100 is not a free choice. `all_hdr.CELL_CONFIG` tunes `group_size` AT the cell's own width
# and `_scaled_max_fail` reads its reference span from `hdr_range`, so 100 is the VALIDATED width
# for K562, HUDEP-2 and HEK293 and a rescale for CD34+_HSPC (native 150). Changing it re-prices the
# group.
# 2026-09-23: all-HDR is the rung BELOW the conjunction again, and it plays the SAME seven
# width-15 windows over the same band space, so a round where the conjunction declines does not
# silently re-correlate the siblings on a different tiling.
# **2026-09-23: h7 plays the all-HDR SPIKE, not the conjunction.** Operator request. A hotkey
# listed here skips the conjunction rung entirely (`miner._build` gates on it) and takes all-HDR
# as its top rung, over the WHOLE band space rather than a width-15 slice of it -- all-HDR's band
# is emergent from its own min-union, so a wider space does not mean a wider band: measured live
# on K562 it is 13 of 100 at width 100 and 13 of 15 at width 15.
#
# This is the first MIXED-construction fleet in this file. It is worth stating why that is safe:
# the two rungs read the same `_window_for`, so a hotkey here is decorrelated from the conjunction
# seven by CONSTRUCTION rather than by window, and its band lands wherever its own min-union puts
# it inside 200-299. Overlap with the conjunction bands is therefore possible and unmeasured.
# **2026-09-23: the split between the two constructions is PER CELL TYPE.** Operator request.
#
#   erythroid (K562, HUDEP-2, CD34+_HSPC)   h0-h6 conjunction, 7 x width 15 stride 15
#                                           h7-h9 all-HDR,     3 x width 34 stride 34
#   HEK293                                  h0-h8 conjunction, 9 x width 12 stride 12
#                                           h9    all-HDR,     1 x the whole class
#
# HEK293 gets more conjunction hotkeys and a narrower window because its `band_k` is 8 against the
# erythroid 11, so a width-12 window still holds 8/12 = 67% density where 12 seeds would be far too
# thin for k=11. The arithmetic that has to hold either way is COVERAGE: `conj_hk * stride >= 100`,
# which is 7*15 = 105 and 9*12 = 108.
#
# `ALL_HK` is the fleet in h0..h9 order; the first `conj_hk` play the conjunction and the rest play
# all-HDR. Splitting by position rather than by name means a cell's two groups can never overlap
# and no hotkey can be left out of both.
ALL_HK: list[str] = ["niome_hotkey", "niome_hotkey1", "niome_hotkey2", "niome_hotkey3",
                     "niome_hotkey4", "niome_hotkey5", "niome_hotkey6", "niome_hotkey7",
                     "niome_hotkey8", "niome_hotkey9"]
CELL_LAYOUT: dict[str, dict] = {
    "HEK293":     {"conj_hk": 9, "width": 12, "stride": 12, "hdr_width": 100, "hdr_stride": 100},
    "K562":       {"conj_hk": 7, "width": 15, "stride": 15, "hdr_width": 34, "hdr_stride": 34},
    "HUDEP-2":    {"conj_hk": 7, "width": 15, "stride": 15, "hdr_width": 34, "hdr_stride": 34},
    "CD34+_HSPC": {"conj_hk": 7, "width": 15, "stride": 15, "hdr_width": 34, "hdr_stride": 34},
}
DEFAULT_LAYOUT = CELL_LAYOUT["K562"]


def layout_for(cell) -> dict:
    """The per-cell split and geometry, falling back to the erythroid shape for an unknown cell."""
    return CELL_LAYOUT.get(cell, DEFAULT_LAYOUT)


def conj_hk(cell=None) -> list[str]:
    """The hotkeys that play conjunction+floor on `cell`."""
    return ALL_HK[:layout_for(cell)["conj_hk"]]


def hdr_hk(cell=None) -> list[str]:
    """The hotkeys that play the all-HDR spike on `cell` -- the remainder of the fleet."""
    return ALL_HK[layout_for(cell)["conj_hk"]:]


def is_hdr_only(instance, cell=None) -> bool:
    """True when this hotkey skips the conjunction rung on this cell."""
    return instance in hdr_hk(cell)


# Kept for `window_plan.py`, `fleet_status.py` and the research scripts, which import these by
# name. They describe the ERYTHROID split, which is the default everywhere a cell is not supplied.
HDR_ONLY_HK: list[str] = ALL_HK[DEFAULT_LAYOUT["conj_hk"]:]

# **Which loop each hotkey plays, per cell type. Set by operator request.**
#
#     CD34+_HSPC / HUDEP-2   1, 2, 3, 6, 7, 9, 10
#     K562 / HEK293          1, 3, 5, 6, 7, 8, 10
#
# Positionally aligned with `LOOP_HK`, so entry i is hotkey i's loop on that cell. Both lists are
# NON-CONTIGUOUS, which is a deliberate change from the 1..N packing and has two consequences worth
# stating rather than rediscovering:
#
#   * **The unplayed loops are still computed.** A hotkey at loop 10 reproduces loops 1-9 to build
#     its exclusion set, including loops nobody plays (4, 5, 8 on the erythroid pair; 2, 4, 9 on
#     the other). That is what keeps every played band disjoint from every other, so it is required,
#     not waste — but it does mean the deepest hotkey pays nine `choose_band` passes instead of six,
#     and draws from a pool with 9*k seeds already removed rather than 6*k.
#   * **Coverage is unchanged at 7*k.** The union is the number of hotkeys times k whatever indices
#     they hold; band position is free under a uniform generator, so spreading the indices costs and
#     buys nothing on `P(hit)`. What it does change is the per-hotkey pool margin, which falls with
#     depth — watch the deepest hotkeys for `conjunction declined (band reached N of K)`.
#
# A cell not named here falls back to DEFAULT_LOOPS.
LOOP_ASSIGN: dict[str, tuple[int, ...]] = {
    "CD34+_HSPC": (1, 2, 3, 6, 7, 9, 10),
    "HUDEP-2":    (1, 2, 3, 6, 7, 9, 10),
    "K562":       (1, 3, 5, 6, 7, 8, 10),
    # 2026-09-22: HEK293 gets its own spread, reaching loop 11. The deepest hotkey therefore
    # reproduces TEN predecessor loops to build its exclusion set and bands from candidates with
    # 10*k = 80 of the 300 already removed -- the thinnest pool any hotkey has been asked for.
    # `loop_sweep.py` measured narrow/k8's surviving pool at 2.75x the group (min 2.66x) at these
    # depths over 50 tasks, so the margin is there; watch h9 for `band reached N of K` regardless.
    "HEK293":     (1, 2, 3, 4, 9, 10, 11),
}
DEFAULT_LOOPS: tuple[int, ...] = (1, 2, 3, 4, 5, 6, 7)

# The three width-100 classes the shared window is built from, set explicitly by operator request
# rather than read from `window_plan.json`. Classes are 0-indexed over 100-999, so (0, 1, 3) is
# 100-199, 200-299 and 400-499 -- a JOINED window, not a contiguous 300.
#
# 2026-09-22: the third class moved 4 -> 3 (500-599 -> 400-499) by operator request. Nothing else
# changed: same three-class shape, same 300-seed span, same loop assignment. Note this DOES change
# `bank_key` on every cell whose `CUT_MODE` is `narrow` (HEK293, CD34+_HSPC, HUDEP-2), since the
# cut seeds are the shared window there -- so the first round after the restart pays one cold bank
# scan per (contract, cell) on those three. K562 cuts wide over the full 900 and is unaffected.
#
# Fixing it here rather than following the plan costs nothing measurable and removes a failure
# mode. CLAUDE.md measures every seed-window strategy inside 0.93x-1.06x of chance (|z| <= 1.1 over
# a 160-task cold walk-forward), and under a uniform generator band position is free -- so three
# arbitrary classes are worth exactly what three predicted ones are. What it removes is the churn:
# `round_plan.sh`'s "the plan was rewritten mid-round" cannot happen to a constant.
# **2026-09-23: the band tiles ONE class and the cut spans TWO.** Operator request, on
# `conj_tile.py`: seven width-15 windows at stride 15 over 200-299 caught **7 of 7** reachable
# seeds across six K562 rounds -- the best catch rate of any arm measured this session -- with
# 42/42 builds forming k=11 from only 15 candidates at a pool of 144-205 against group 100.
#
# The two spaces are deliberately different and must not be merged:
#
#   BAND_CLASSES  (1,)    200-299          the 100 seeds the seven windows tile
#   CUT_CLASSES   (1, 3)  200-299+400-499  the 200 seeds the cut min-unions over and clean counts in
#
# The cut is wider than the band because the FLOOR is a property of the cut: a cut restricted to
# 200-299 can never make a 400-499 seed clean, which would leave a third of each round's seeds
# unreachable by either regime. Two classes reaches 15 of 18 drawn seeds across the six rounds.
#
# **Four classes was measured and is WORSE.** The same seven windows at a 400-seed cut
# (200-299+400-499+600-699+700-799) reached 17 of 18 seeds but took band hits 7 -> 3 and summed
# final 1696 -> 1351. `choose_band` picks from the guides that survive the CUT min-union, so a
# harder cut leaves a different, smaller surviving set and the greedy lands on a different band --
# CLAUDE.md puts band overlap between cut arms at 0-1 of 11. Same windows does NOT mean same band.
# Across three spans on identical windows the band hits run 7 (cut 100), 7 (cut 200), 3 (cut 400).
# 2026-09-23 (second move): the whole layout shifts one block down, by operator request.
#   band  200-299 -> 100-199   (class 1 -> 0)
#   cut   200-299+400-499 -> 100-199+500-599   (classes 1,3 -> 0,4)
# Nothing else changes: same seven width-15 windows, same stride, same k, same group. Band
# position is free under a uniform generator (CLAUDE.md measures every seed-window strategy inside
# 0.93x-1.06x of chance), so this is a relocation rather than a re-tuning -- but it DOES change
# `bank_key` on every cell, so the first round after the restart pays one cold bank scan per
# (contract, cell).
BAND_CLASSES: tuple[int, ...] = (3,)
CUT_CLASSES: tuple[int, ...] = (1, 3)

# Compatibility aliases. `window_plan.py`, `fleet_status.py` and several research scripts import
# BAND_HK/BAND_STRIDE/BAND_SPAN by name; keeping them pointed at the single tier means those keep
# working and see one group of ten rather than raising on a missing name.
BAND_HK: list[str] = LOOP_HK
# 2026-09-23: width 15 at stride 15 over a 100-seed band space. `count * stride` is 7 * 15 = 105
# against a span of 100, so this is NOT an exact tiling: the seventh window WRAPS and overlaps the
# first by 5 seeds (w6 = 200-204 + 290-299). Union is still 100/100 and 95 of the 100 seeds sit in
# exactly one window, so the bands are near-disjoint; 5 of the 105 seats can duplicate.
BAND_STRIDE = 15
BAND_SUB_WIDTH = 15
BAND_SPAN = 100

# Group B and group C are both dissolved into the single tier above.
REST_HK: list[str] = []
REST_STRIDE = 100
REST_SPAN = 0
PRED_HK: list[str] = []
PRED_STRIDE = 100
PRED_SPAN = 0
PRED_SUB_WIDTH = 100

# The whole seed space a round can draw from. The band split above and the wide cut window below are
# both expressed against this one list so they cannot drift apart.
FULL_SPACE = list(range(100, 1000))

# Per-hotkey offset overrides. EMPTY and should stay so: the sequence below already places seven
# windows across 900 with no wrap, and an override would pull one off that.
BAND_OFFSET_OVERRIDES: dict[str, int] = {}

# The tiling invariants that guarded the old layouts do not apply: stride 0 means the ten windows
# COINCIDE by design, which every previous layout treated as the one fatal error. What has to hold
# instead is that the shared window is well formed and that the loop indices are distinct -- the
# loop index is now the only thing separating two hotkeys, so a duplicate is exactly the collision
# `BAND_OFFSET_OVERRIDES` used to guard against.
# The offset layout is back, so the invariant that matters again is COVERAGE, not an exact tiling.
# 7 x 15 = 105 over a 100-seed span overshoots by 5, which `sub_window` absorbs by wrapping; what
# must hold is that the seven slices reach every seed of the band space, or a seed exists that no
# hotkey can band. Checked directly rather than inferred from count/stride arithmetic.
assert BAND_STRIDE > 0, "the offset layout needs a nonzero stride; stride 0 is the loop-axis arm"
assert len(BAND_CLASSES) * 100 == BAND_SPAN, \
    "the band span must be exactly the classes the band windows tile"
assert BAND_SUB_WIDTH <= BAND_SPAN, "the sub window is wider than the band space"
assert len(set(BAND_CLASSES)) == len(BAND_CLASSES), "a band class is listed twice"
assert len(set(CUT_CLASSES)) == len(CUT_CLASSES), "a cut class is listed twice"
assert all(0 <= c <= 8 for c in BAND_CLASSES), "a band class is outside 100-999"
assert all(0 <= c <= 8 for c in CUT_CLASSES), "a cut class is outside 100-999"
assert set(BAND_CLASSES) <= set(CUT_CLASSES), \
    "the band space must sit inside the cut space, or `band <= cut` fails"
assert len(set(ALL_HK)) == len(ALL_HK), "a hotkey is listed twice in ALL_HK"
for _cell, _lay in CELL_LAYOUT.items():
    _n = _lay["conj_hk"]
    assert 0 < _n <= len(ALL_HK), f"{_cell}: {_n} conjunction hotkeys of {len(ALL_HK)}"
    assert _n * _lay["stride"] >= BAND_SPAN, \
        (f"{_cell}: {_n} hotkeys x stride {_lay['stride']} = {_n * _lay['stride']} does not reach "
         f"the {BAND_SPAN}-seed band space; some seed sits in no hotkey's window")
    assert (_n - 1) * _lay["stride"] < BAND_SPAN, \
        f"{_cell}: the last conjunction hotkey wraps a full revolution onto another's slice"
    assert _lay["width"] <= BAND_SPAN, f"{_cell}: the band sub window is wider than the space"
    _m = len(ALL_HK) - _n
    assert _m >= 0, f"{_cell}: more conjunction hotkeys than the fleet has"
    if _m:
        assert _m * _lay["hdr_stride"] >= 100, \
            (f"{_cell}: {_m} all-HDR hotkeys x stride {_lay['hdr_stride']} does not reach the "
             f"100-seed class; some seed sits in no all-HDR window")
        assert (_m - 1) * _lay["hdr_stride"] < 100, \
            f"{_cell}: the last all-HDR hotkey wraps a full revolution onto another's slice"
assert BAND_SPAN <= len(FULL_SPACE), "the band span claims more than 100-999"
assert not (set(BAND_HK) & set(REST_HK)), "a hotkey cannot be in both band groups"
assert not (set(PRED_HK) & (set(BAND_HK) | set(REST_HK))), "a hotkey cannot be in two band groups"


def band_window() -> list[int]:
    """The 100 seeds the seven band sub-windows tile (200-299)."""
    return sorted(s for c in BAND_CLASSES for s in range(c * 100 + 100, c * 100 + 200))


def cut_window() -> list[int]:
    """The 200 seeds the cut min-unions over and the clean set is counted in."""
    return sorted(s for c in CUT_CLASSES for s in range(c * 100 + 100, c * 100 + 200))


def shared_window() -> list[int]:
    """The 300 seeds every hotkey bands from and (on a narrow cell) cuts over.

    Built from `SHARED_CLASSES` rather than from the plan, so it is the same list on every hotkey
    and every round -- which is what lets loop L reproduce loops 1..L-1 without coordination.
    """
    return cut_window()


def band_loop(instance, cell=None) -> int | None:
    """This hotkey's position on the re-banding chain for `cell`, 1-based, or None if it plays none.

    Loop L excludes every seed loops 1..L-1 banded, so the played bands are disjoint by
    construction rather than by luck. The mapping is per CELL TYPE (`LOOP_ASSIGN`) since
    2026-09-21; `cell=None` returns the default packing and is only for callers that genuinely
    have no cell in hand (diagnostics), never for a build.
    """
    if instance not in conj_hk(cell):
        return None
    # 2026-09-23: loop 1 for every hotkey. The seven width-15 windows are near-disjoint already
    # (only w6/w0 share 5 seeds), so the loop axis has nothing left to separate -- and stacking it
    # on top would ban seeds a sibling banded in a DIFFERENT window, which is not what the loop
    # exclusion means. `LOOP_ASSIGN` is kept for the layout it belongs to, not read here.
    return 1


def band_layout(instance, cell=None):
    """``(hotkeys, stride, span, on_predicted)`` for a band hotkey, else None.

    ``on_predicted`` is what separates the two groups: True means this hotkey's sub-window is cut
    out of the plan's three predicted classes, False means out of the 600 seeds the plan left over.
    """
    if instance in conj_hk(cell):
        return conj_hk(cell), layout_for(cell)["stride"], BAND_SPAN, True
    if instance in BAND_HK:
        return BAND_HK, BAND_STRIDE, BAND_SPAN, True
    if instance in REST_HK:
        return REST_HK, REST_STRIDE, REST_SPAN, False
    if instance in PRED_HK:
        # `False` in slot 3 means "not the whole-900 space"; `band_space` below distinguishes group
        # C from the (now empty) complement group C is NOT -- it draws from the predicted 300.
        return PRED_HK, PRED_STRIDE, PRED_SPAN, False
    return None


def band_sub_width(instance, cell=None):
    """How wide a sub-window this hotkey draws its band from, or None to keep the cell's own value.

    `conjunction.CELL_CONFIG.band_width` is per CELL, so it cannot separate two hotkeys on different
    tiers -- group A needs 300 (of a 900-seed space) and group C needs 100 (of a 300-seed one). The
    per-cell value stays the default for anything outside a band group.
    """
    if instance in PRED_HK:
        return PRED_SUB_WIDTH
    if instance in conj_hk(cell):
        return layout_for(cell)["width"]
    if instance in REST_HK:
        return BAND_SUB_WIDTH
    return None


def band_offset(instance, cell=None):
    """This hotkey's band sub-window offset, in seeds of ITS OWN space, or None."""
    layout = band_layout(instance, cell)
    if layout is None:
        return None
    hotkeys, stride, span, _ = layout
    if instance in conj_hk(cell):
        hotkeys, stride = conj_hk(cell), layout_for(cell)["stride"]
    off = BAND_OFFSET_OVERRIDES.get(instance, hotkeys.index(instance) * stride)
    # 2026-09-23: wrapping is INTENDED again. The 09-20 layout PLACED windows across 100-999 with
    # no wrap, so an offset past `span - width` meant a broken layout and this raised. The seven
    # width-15 slices rotate inside a 100-seed space and 7 * 15 = 105 overshoots by 5 on purpose --
    # `sub_window` wraps the last window to 200-204 + 290-299, which is exactly the arm
    # `conj_tile.py` measured at 7/7 seeds caught. What must not happen is a FULL revolution, which
    # would land two hotkeys on the identical slice; that is what the bound below checks.
    width = band_sub_width(instance, cell) or BAND_SUB_WIDTH
    if off >= span:
        raise ValueError(f"band offset {off} reaches or passes the {span}-seed band space for "
                         f"{instance}; it would wrap a full revolution onto another hotkey's slice")
    return off


def band_offset_frac(instance, cell=None):
    """This hotkey's band sub-window offset as a fraction of its own space, or None.

    None means "not a band hotkey here" and the caller keeps `ConjunctionConfig`'s own default,
    which is the single-hotkey value the 12-contract replication was measured at.

    **The fraction is against the hotkey's OWN group span (300 for h0-h5, 600 for h6-h9), not
    against the 900-seed whole.** `conjunction.sub_window` turns it back into an index with
    `int(round(n * frac))` on whatever list it is handed, so the only correct pairing is this
    fraction with the list `band_space()` returns for the same hotkey. Pair it with the other
    group's space and the offset lands in the wrong place — which is why `band_space` exists rather
    than leaving each caller to slice the plan window itself.
    """
    layout = band_layout(instance, cell)
    if layout is None:
        return None
    return band_offset(instance, cell) / float(layout[2])


def band_space(instance, predicted, cell=None):
    """The seeds this hotkey cuts its band sub-window out of, or None if it plays no band.

    ``predicted`` is the plan's joined space — the three predicted width-100 classes, 300 seeds, as
    a sorted seed list. A group A hotkey draws from exactly that; a group B hotkey draws from its
    complement in 100-999, the 600 seeds no group A hotkey can reach.

    Returning the complement here rather than at the call site is what keeps the two groups from
    drifting: the 600 is DEFINED as "whatever the plan did not predict", so a plan that moves its
    classes moves both groups together and the fleet still covers 100-999 exactly once per group.
    """
    layout = band_layout(instance, cell)
    if layout is None:
        return None
    # 2026-09-20: the band space is the FULL 100-999 for every band hotkey. `predicted` is accepted
    # and ignored -- the plan's three classes no longer steer the band at all, only the all-HDR rung
    # below the conjunction still reads them. Kept in the signature so callers do not have to know
    # which layout is live.
    if layout[3]:
        # 2026-09-23: the band space is BAND_CLASSES (200-299, 100 seeds) and each hotkey takes a
        # width-15 slice of it at stride 15. `predicted` is accepted and ignored -- the band space
        # is fixed here, not by the plan. Separation is the OFFSET again, not the loop axis: seven
        # width-15 slices of 100 seeds are already near-disjoint (w6 wraps and overlaps w0 by 5).
        return band_window()
    if instance in PRED_HK:
        # Group C is the one tier the plan still steers. With `JOINED_SOURCE` on "uniform" that is a
        # deterministic 100-399 every round (stable-sort tie-break over nine tied classes), so these
        # three are effectively fixed windows -- which costs nothing, because band position is free
        # under a uniform generator.
        pred = sorted({int(s) for s in (predicted or [])})
        return pred or list(FULL_SPACE[:PRED_SPAN])
    pred = sorted({int(s) for s in (predicted or [])})
    rest = [s for s in FULL_SPACE if s not in set(pred)]
    return rest or list(FULL_SPACE)


# 2026-09-17: the CUT window for band hotkeys, decoupled from the band space above. Every
# previously-tuned `conjunction.CELL_CONFIG` row assumed cut == band space (300 seeds); this widens
# the cut to the full 100-999 for exactly the hotkeys that also get a rotating band offset, i.e.
# every hotkey `band_layout()` names -- BOTH groups since 2026-09-18, which is why the gate tests
# `band_layout` rather than `BAND_HK` (group B would otherwise silently lose the wide cut).
# Which CELL TYPES take the wide cut is the second gate below. Both must pass, so this is a per-cell
# veto on a per-hotkey layout.
#
# 2026-09-17, measured (`conj_wall.py`, 4 contracts per cell, wide arm against a narrow control):
# the wide cut LOWERS HEK293's band-formation wall from 9 to 8 and leaves all three erythroid cells
# at 12. Unanimous, 4/4 within every cell, and the narrow control returns CLAUDE.md's recorded
# 9/12/12/12 on the same code path -- so the probe reproduces the known answer before being trusted
# on the new one.
#
# The mechanism is the `bank_keep` cap, which is why one cell moved and three did not. The wall sits
# where `bank * P(rule)**k` falls under `group_size`; min-unioning the cut across 900 seeds instead
# of 300 is a strictly harder constraint, so fewer guides qualify. HEK293's bank falls to 140k-182k,
# BELOW the 300,000 cap, and the wall drops a seed. The erythroid banks are still pinned AT the cap
# even over 900 seeds, so nothing changes for them. (Implied per-seed decay backs this up: 0.38-0.40
# on HEK293 against its documented P(HDR) ~ 0.37, 0.513 on the erythroid cells.)
#
# So the wide cut costs HEK293 specifically, twice: it puts `CELL_CONFIG`'s `band_k` 8 exactly AT
# the wall rather than one step below it -- the build then succeeds only while the wall holds at 8,
# with no slack for a contract whose bank comes in thinner -- and it makes k=9 unreachable, which is
# the arm CLAUDE.md prices as HEK293's best (0.000163 against k=8's 0.000136, -17%).
#
# HEK293 was therefore excluded and kept the 300-seed joined cut. The erythroid cells keep the wide
# cut by operator request: their wall is 12 either way, so the measurement gives no reason to move
# them, and their gates in `Miner.CONJUNCTION_MIN_BUDGET_S` are already calibrated for it.
#
# **2026-09-18: HEK293 is back IN, by operator request, so all four cells now take the 900-seed cut
# and this frozenset is inert.** It is kept as a set of all four rather than deleted so the veto can
# be re-narrowed to one cell without rebuilding the mechanism. The measurement above is unchanged
# and is the evidence AGAINST this: `band_k` 8 now sits exactly AT the wide-cut wall of 8 with no
# step of slack, and k=9 -- the arm CLAUDE.md prices as HEK293's best, 0.000163 against k=8's
# 0.000136 -- is unreachable. What changed on the other side of the ledger is the 2026-09-18 band
# split: with `REST_HK` drawing from the 600 seeds the plan did NOT predict, a NARROW HEK293 cut put
# h6-h9's band entirely outside their own cut window, and widening the cut restores band ⊆ cut for
# all ten hotkeys. Measured margin at the wide cut, so the "no slack" is a number rather than a
# worry: see the live check recorded in CLAUDE.md.
#
# **What excluding a cell here does to the BAND, which is not nothing and was once recorded wrong.**
# The band CANDIDATE set is identical either way -- `band_space()` picks it, and that function does
# not know about the cut at all -- but the chosen band is not: `choose_band` picks from the guides
# that survive the cut min-union, and those pools differ between a 300-seed and a 900-seed cut.
# Measured overlap is 0 or 1 of 11 on 8 of 8 contracts, so the two arms build essentially DISJOINT
# bands of the same size. (An earlier form of this comment claimed they "produce the identical band";
# that is wrong and CLAUDE.md records the correction.) Under a uniform generator band position is
# free, so same-size-different-place costs nothing -- which is what `widecut_price.py` then measured
# directly, at +1.1% with sign p = 0.688.
#
# **One consequence of the 2026-09-18 split, specific to HEK293.** A group B hotkey (h6-h9) draws
# its band from the 600 seeds the plan did NOT predict, while HEK293 -- excluded here -- keeps the
# predicted 300 as its CUT window. So on HEK293 those four hotkeys hold a band that lies entirely
# OUTSIDE their own cut space. That is mechanically sound rather than a bug: the band is pinned by
# `hdr_compliance` over `band_candidates`, which is computed per seed against the bank and never
# consults `cfg.seeds`, and a guide satisfying the `hdr` rule on a seed necessarily cut on it -- so
# every band row still pins all three of stage 4's targets there. What it does mean is that the
# band and the off-band clean floor now come from disjoint seed sets on that one cell, and that the
# wall (`bank * P(rule)**k` against `group_size`) is being asked about seeds the bank was not
# selected on. Unmeasured; HEK293's live `band_k` is 8 against a narrow-cut wall of 9, so there is
# one step of slack for it to spend.
# **2026-09-18 (evening): the cut span is now a per-cell MODE, not a per-cell boolean.**
# `hek_cutspan.py` priced all three spans on HEK293 and found E[own-field share] cannot separate
# them at k=8 (±4%, nothing above 3 of 4 paired) -- so the decision fell to the non-score terms, and
# there the 900-seed cut is the one arm that is dominated: it walls at 8 where the 300- and
# 450-seed spans wall at 9, and it holds a fifth less off-band floor (cut-clean 19.2 against 24.0 of
# 900, enumerated from the shipped rows over all 900, not from `meta["clean"]`). `hek_k9.py` then
# measured the step the wall controls: k=8 -> k=9 is worth **+10.4% at cut 450** and +12.2% at cut
# 300, 4W/0L on both, building 4/4. So HEK293 moves to `union` + k=9, worth ~1.13x the shipped
# wide/k=8 arm (chained across two runs at an anchor that reproduces to 0.0%).
#
#   "wide"   the full 100-999. Erythroid cells: their wall is 12 either way, so the wide cut costs
#            them no depth, and `widecut_price.py` priced it a wash (+1.1%, p = 0.688) on score.
#   "union"  the predicted space PLUS this hotkey's own band sub-window. For `BAND_HK` that IS the
#            predicted space (their band already sits inside it); for `REST_HK` it is 300 + 150 =
#            450 seeds, which is the narrowest span that still satisfies band ⊆ cut.
#
# What `union` costs: TWO bank keys per (contract, cell) instead of one, since h0-h5 min-union over
# 300 seeds and h6-h9 over 450. `bank_slot` means one hotkey pays each scan and the rest load warm,
# but it is two cold scans per round rather than one.
# **2026-09-20: every cell is on "union" — the erythroid cells drop the 900-seed cut.** Set by
# operator request, and it is the change this file has recommended since 2026-09-17 rather than a
# new direction. With the band space now the full 900 and each hotkey holding a 300-seed contiguous
# window, "union" resolves to that window, so the cut and the band coincide on all four cells.
#
# What the wide cut was measured to be worth on erythroid, and why dropping it is not a loss:
#   * `widecut_price.py`, 4 contracts x 2 cells paired within contract: **+1.1%, 4W/2L/2T,
#     sign p = 0.688** — a wash on own-field E[share].
#   * `coldbuild.py`: the wide cut roughly DOUBLES the cold build (K562 199s -> 410s, CD34+
#     182s -> 386s), and essentially all of it is the Cas12a bank scan (2.48x).
#   * The prefetch LEAD is the binding constraint, not the gate: at a p10 lead of 395s a 410s cold
#     build misses the window on about the bottom decile of rounds — and because `bank_slot` shares
#     one bank per (contract, cell), when it misses EVERY erythroid hotkey falls through to all-HDR
#     together.
# So this trades ~1% of score, inside noise, for halving the build and removing a fleet-wide
# correlated failure mode.
#
# The one real cost is the floor, and it points the other way from HEK293: the wide cut WIDENS the
# erythroid cut-clean set +30-40% (86-104 off-band clean seeds narrow against 115-133 wide), where
# on HEK293 it NARROWS it. That widening is exactly what `floor_price.py` arm A prices at ~+3% and
# `widecut_price.py` measured at +1%, so it is already counted above — do not subtract it twice.
#
# The wall does not move: `conj_wall.py` measured the erythroid band-formation wall at **12 on both
# the narrow and the wide cut**, unanimous 4/4 per cell, so `band_k` 11 stays one full step below it
# either way and the narrow arm's surviving-pool margin is the healthier of the two.
# **2026-09-20 (reverted same day): the erythroid cells go back to the 900-seed shared cut.** The
# all-union config above was live for one round and the first one measured its real cost, which the
# prior evidence could not have shown:
#
#   * With a PER-HOTKEY cut window there is no shared bank. `bank_key` folds in the seed list, so
#     seven different 300-seed cuts are **seven different banks** — confirmed, 7 bank files written
#     in one round — and every build is therefore COLD. There is no warm sibling any more.
#   * Seven concurrent Cas12a scans on one GPU (100% utilisation) took **566-580s each**, against
#     the shared wide bank's 1 x ~410s cold plus 6 x ~56s warm.
#   * It cost a submission on the very first round: h4 was called 196s into its build, waited, ran
#     out of window and shipped the fallback; its conjunction only finished 3 minutes later.
#
# **The change was meant to remove a fleet-wide correlated failure and GPU contention re-created
# it.** Not sharing a bank does not decouple the hotkeys when they all contend for the same device —
# they still miss together, now at 578s against a p10 prefetch lead of 395s instead of 410s.
#
# `coldbuild.py`'s ~200s narrow figure is not wrong, it is a SINGLE-SCAN measurement, and nothing
# in this file had ever exercised seven at once. Read every build-time number here as conditional
# on how many banks the layout implies.
#
# HEK293 stays on "union" — but see the note in `conjunction_cut_seeds`: with per-hotkey band
# windows it now implies 7 banks there too, and that is UNMEASURED.
# 2026-09-21, operator request: every cell to "narrow" -- a 300-seed joined cut window, replacing
# "wide" (the full 900) on the erythroid cells and "union" on HEK293.
#
# **"narrow" means the 300-seed window containing THIS hotkey's own band, not a single global 300.**
# Taking it as "the plan's predicted 300 for everybody" would break `band ⊆ cut` for group A, whose
# band space is the full 900 -- most of their bands sit outside the predicted classes entirely. The
# per-tier reading gives every hotkey a 300-wide cut AND keeps the band inside it:
#
#     group C (h1-h3)      the predicted 300 -- their band space already IS that, at width 100
#     group A (h0, h4-h9)  their own width-300 band sub-window, i.e. `band_candidates`
#
# Group A's case is byte-identical to what "union" already returned, so this is a rename for them
# and a real change only for the erythroid cells, which come off the 900-seed cut.
#
# What that costs, from the measurements already on file: `widecut_price.py` prices wide-vs-narrow
# as a WASH on score (+1.1%, 4W/2L/2T, sign p = 0.688), and `coldbuild.py` measures narrow at
# roughly HALF the cold build (K562 199s against 410s, CD34+ 182s against 386s) because the whole
# difference is the Cas12a bank scan. Against a prefetch-lead p10 of 395s that is the entry's own
# recommendation -- it calls dropping the wide cut on the erythroid cells the right trade and
# records the status quo as operator preference, so this moves TOWARD what was measured.
#
# Two consequences to expect. `conjunction.CELL_CONFIG`'s k/group rows were tuned at the 900-seed
# cut, and the wall moves with the bank: `conj_wall.py` measured the erythroid wall at 12 on BOTH
# spans but HEK293's at 9 narrow against 8 wide, so HEK293 gains a step of slack at k=8 rather than
# losing one. And `CONJUNCTION_MIN_BUDGET_S` of 700 on K562/CD34+ was raised FOR the wide cut; at
# ~200s cold it is now very conservative, which costs availability only at the margin.
# **2026-09-21, operator request: narrow on CD34+/HUDEP-2, wide on HEK293/K562.**
#
# "narrow" is now the shared 300-seed window itself (`shared_window()`), so on those two cells the
# cut and the band space COINCIDE -- and because all ten hotkeys share that window, narrow costs
# ONE bank per (contract, cell), not ten. That is the constraint that sank the 2026-09-20 per-hotkey
# narrow cut (seven banks, seven concurrent 566-580s scans, a lost submission); it does not apply
# here.
#
# **What the narrow cut buys, measured 2026-09-21** (`loop_narrow.py` against `loop_axis.py` +
# `loop_clean.py`, 40 tasks per arm, 36 shared, same oracle windows, cut the only variable):
#
#   | cell       | clean density  wide -> narrow | clean-hit rate wide -> narrow |
#   | HUDEP-2    | 15.9% of 900 -> 35.6% of 300  |  0/17  ->  9/30  (30.0%)      |
#   | K562       | 15.2%        -> 34.0%         |  2/28  ->  4/14  (28.6%)      |
#   | CD34+_HSPC | 15.9%        -> 35.1%         |  2/20  ->  7/20  (35.0%)      |
#   | HEK293     |  1.9%        ->  7.8%         |  0/19  ->  1/8   (12.5%)      |
#   | TOTAL      |                               |  4/84 (4.8%) -> 21/72 (29.2%) |
#
# Fisher exact p = 6.8e-05; on the 36 shared tasks alone 18/65 against 4/74, so it is not the task
# lists drifting. Two multiplied effects: the clean set roughly doubles in density (every narrow
# bank returns at the 300,000 `bank_keep` cap, where HEK293's wide banks collapse to 140k-182k),
# AND the wide clean set is DEPLETED inside the band window while the narrow one is uniform over it
# -- a min-union over 900 has no reason to spend its freedom on the 300 seeds that matter.
#
# **The measurement favours narrow on all four cells, so wide on HEK293/K562 is the operator's
# call and this is the number it costs.** On HEK293 the case for wide is weak either way: 7.8%
# density is still too thin for the floor to pay (mean cons went 0.4186 wide -> 0.3988 narrow on a
# small sample), and `conj_wall.py` puts the narrow wall at 9 against wide's 8, so wide is also the
# span where `band_k` 9 is unreachable. On K562 the narrow arm measured a 4x better clean-hit rate
# (7.1% -> 28.6%) and `coldbuild.py` puts narrow at roughly half the cold build (199s vs 410s).
# Both are recorded here as evidence, not as a change.
# **2026-09-22: HEK293 -> `narrow`, on the 50-task six-arm sweep** (`loop_sweep.py`, k 6/8/9 x
# {wide 900, narrow 300}, 11 loops, the real seed classes as the window). The cut span is not a
# wash on this cell -- it decides availability:
#
#   | arm       | loops formed of 11 | 11/11 tasks | pool margin | clean seeds (abs) |
#   | narrow/k8 |               11.0 |       50/50 |       2.75x |              23.2 |
#   | wide/k8   |               11.0 |       50/50 |       1.61x |              17.2 |
#   | narrow/k9 |               11.0 |       50/50 |       1.29x |              14.4 |
#   | wide/k9   |                0.3 |        0/50 |           - |      never built  |
#
# **`wide/k9` formed ZERO loops on 45 of 50 tasks** -- the band-formation wall drops from 9 to 8
# once the cut is 900, reproducing `conj_wall.py`'s n=4 result at n=50. At the shipped k=8 the wide
# cut still costs a full step of pool margin (1.61x against 2.75x) and ~25% of the absolute clean
# set, and buys only build time this cell's 600s gate does not need.
#
# What it does NOT buy is band reach: narrow and wide found 47 and 46 seeds of 150 respectively, so
# the cut decides availability and floor, not how often a band lands. K562 keeps `wide` -- its wall
# is 12 on both spans, so the same argument does not transfer, and its own arm is still being swept.
# 2026-09-23: NARROW on all four. K562 was the last cell still on `wide`; with the band tiling a
# 100-seed class, a 900-seed cut is the arm `conj_tile.py` measured WORST -- band hits fall 7 -> 3
# as the cut widens 200 -> 400, and CLAUDE.md's own wide-cut pricing was a wash (+1.1%, p = 0.688)
# while costing roughly double the cold build.
CUT_MODE: dict[str, str] = {
    "HEK293": "narrow",
    "K562": "narrow",
    "CD34+_HSPC": "narrow",
    "HUDEP-2": "narrow",
}

# Kept as a derived view so callers that predate `CUT_MODE` keep working. Note a cell on "union" is
# NOT in here -- a research script testing `cell in WIDE_CUT_CELLS` will correctly see HEK293 as no
# longer taking the full 900.
WIDE_CUT_CELLS: frozenset = frozenset(c for c, m in CUT_MODE.items() if m == "wide")


def conjunction_cut_seeds(instance, cell, predicted=None, band_candidates=None):
    """This hotkey's CUT window for one cell, or None to keep the band space.

    `None` means "build the cut over whatever space the band is drawn from" -- the behaviour every
    hotkey had before any of this existed, and what any hotkey outside both band groups, or any cell
    outside `CUT_MODE`, still gets.

    `cell` is REQUIRED rather than defaulted on purpose. A caller that forgets it gets a TypeError
    instead of silently building the wrong cut span on a cell measured to be hurt by it, which is
    the same reasoning `miner.sh`'s `assert_disjoint_windows` uses when it exits rather than
    warning. `predicted` and `band_candidates` are required for the same reason, but only on a
    "union" cell -- a caller that omits them there raises rather than quietly falling back to a span
    that was measured as the worse arm. Research scripts written against the old two-argument
    signature therefore raise on HEK293 and keep working on the erythroid cells, which is the
    intended failure: they need updating, not a silent default.
    """
    if band_layout(instance, cell) is None:
        return None
    mode = CUT_MODE.get(cell)
    if mode == "wide":
        return list(FULL_SPACE)
    if mode == "union":
        if band_candidates is None:
            raise TypeError(f"cut mode 'union' on {cell} needs `band_candidates`")
        # The narrowest span satisfying band ⊆ cut, which is what "union" has always MEANT. Until
        # 2026-09-20 the band came out of the plan's predicted classes, so that span was
        # `predicted ∪ band` (300 seeds for a predicted-group hotkey, 450 for a complement one).
        # The band space is now the full 900 and the plan no longer steers it, so the predicted term
        # is vestigial and dropping it leaves the band's own 300-seed window. That is not a new arm:
        # it is exactly `hek_cutspan.py`'s A/narrow, the best-measured HEK293 option on every
        # score-adjacent term (wall 9 so k=9 stays reachable, clean 24.0/900, margin 2.71x).
        return sorted(set(int(s) for s in band_candidates))
    if mode == "narrow":
        # The shared 300-seed window: cut == band space, so `band ⊆ cut` holds trivially and every
        # hotkey folds the SAME seed list into `bank_key` -- one bank per (contract, cell).
        #
        # `band_candidates` is still required rather than defaulted, and is cross-checked against
        # `shared_window()`. A caller that drifts from the shared window would silently fragment
        # the bank ten ways, which is the exact failure the 2026-09-20 layout hit; raising is how
        # that gets caught at the call site instead of in a GPU queue.
        if band_candidates is None:
            raise TypeError(f"cut mode 'narrow' on {cell} needs `band_candidates` for {instance}")
        cand = sorted(set(int(s) for s in band_candidates))
        cut = cut_window()
        # The cut is the SAME 200 seeds on every hotkey, so `bank_key` folds in one seed list and
        # the fleet pays one bank scan per (contract, cell) rather than seven. What is per-hotkey
        # is the BAND, and it has to sit inside the cut or `band <= cut` fails and a band seed
        # stops being a clean seed -- the property the whole construction rests on.
        if not set(cand) <= set(cut):
            raise ValueError(
                f"cut mode 'narrow': {instance}'s band candidates are not inside the cut; "
                f"{len(set(cand) - set(cut))} of {len(cand)} seeds fall outside the "
                f"{len(cut)}-seed cut window. band <= cut would not hold.")
        return cut
    return None

# The width a FULL_HK hotkey takes out of the joined space. `None` (or anything >= the joined span)
# means the WHOLE space -- 300 seeds at three width-100 classes -- with no rotation offset, which is
# the "full width" tier CLAUDE.md prices at band 11 on K562 and 7 on HEK293.
#
# It is a separate knob from `sub_width` on purpose. The rotating slices must take the cell's
# VALIDATED width, because `group_size` is tuned at that width; a full-width hotkey is by definition
# not taking the cell width, so deriving its width from CELL_CONFIG only ever narrowed it back to a
# half-stride slice. Both `fallback_for` here and `window_plan._joined_full` read this one value so
# the plan and the fallback cannot drift.
#
# What full width actually costs, measured on the live fleet rather than extrapolated: the band is
# **11 at 300 against 12 at 225 on K562, and flat at 7 on HEK293** -- one seed and none. That is far
# cheaper than the contiguous width sweep (13/12/11/9 at 100/150/200/300) predicts, and CLAUDE.md's
# `fidelity_window.py` row measures `weighted x fidelity` FLAT at 221.6-223.4 across contiguous 100
# / contiguous 225 / joined 225 / joined 300, so the width is free on the product too. What it buys
# is reach: the whole joined space instead of 75% of it, which is the only thing that matters while
# one hotkey is carrying the fleet.
#
# `main_max_fail` follows the span automatically -- `all_hdr._scaled_max_fail` holds the z-score of
# the tuned value, giving 135 on K562/HUDEP-2, 156 on CD34+_HSPC and 165 on HEK293 at span 300 --
# and `all_hdr.WIDE_WINDOW_VARIANTS` caps `variants` at 44000. Set this to a number to go back to a
# half-stride slice of that width.
FULL_SUB_WIDTH: int | None = None

# Used only if `all_hdr.CELL_CONFIG` cannot be read, so a cron run never dies on an import.
DEFAULT_WIDTH = 225

# The fallback's joined space, as width-100 class indices (0 -> 100-199 ... 8 -> 900-999).
#
# The plan's three classes come from `JOINED_SOURCE` (`repeat_last`, shipped). The fallback has no
# prediction to read, so it takes a FIXED triple — and that costs nothing, because band position is
# free under a uniform generator and the generator is measured uniform (seeds uniform: min-gap
# median 94 observed against 93 simulated). Any three classes are worth as many as any other three.
#
# What is NOT free is the choice of joined-vs-spread, and it goes the other way here than in the
# plan: eleven bands of ~12 cannot sit distinctly in 300 seeds, so the union is ~104 of 900 (30.8%
# of rounds) against ~126 (36.4%) for the same eleven spread over the full 900 (`band_overlap.py`).
# The plan pays that ~5.6 points to concentrate on a predicted space. The fallback has nothing to
# concentrate ON, so set `FALLBACK_SPREAD = True` to take the wider union instead.
FALLBACK_CLASSES = [0, 2, 7]          # 100-199, 300-399, 800-899
FALLBACK_SPREAD = False               # True -> rotate over the whole 100-999 instead


def sub_width(cell, default=DEFAULT_WIDTH):
    """The cell's validated band width, from `all_hdr.CELL_CONFIG`.

    100 for K562 and HUDEP-2, 150 for CD34+_HSPC, 75 for HEK293 — validated over five contracts
    each. Read from CELL_CONFIG rather than restated so the two cannot drift: `group_size` is tuned
    at the width, so a layout that hands a cell the wrong width also hands it the wrong group.
    """
    try:
        from niome_subnet.genomics import all_hdr as _AH
        lo, hi = _AH.CELL_CONFIG[cell]["hdr_range"]
        return hi - lo + 1
    except Exception:
        return default


def expand(classes):
    """Width-100 class indices -> the sorted joined seed list."""
    return sorted(s for w in classes for s in range(w * 100 + 100, w * 100 + 200))


def spans_of(seeds):
    """A sorted seed list -> ascending, non-overlapping [lo, hi] ranges."""
    seeds = sorted(seeds)
    out, lo, prev = [], seeds[0], seeds[0]
    for x in seeds[1:]:
        if x != prev + 1:
            out.append([lo, prev])
            lo = x
        prev = x
    out.append([lo, prev])
    return out


def slice_at(seeds, offset, width):
    """The circular width-`width` slice starting at index `offset`, as [lo, hi] spans."""
    n = len(seeds)
    width = min(width or DEFAULT_WIDTH, n)
    return spans_of({seeds[(offset + k) % n] for k in range(width)})


def rotated(seeds, width, count=None):
    """One slice per rotating hotkey: hotkey h takes indices (h*STRIDE + k) mod n for k < width."""
    n = max(1, len(seeds))
    return [slice_at(seeds, (h * STRIDE) % n, width)
            for h in range(count if count is not None else len(ROTATE_HK))]


def half_stride(seeds, width):
    """h0's slice: `width` seeds at a half-stride offset, or the whole space if it is not wider."""
    if not width or width >= len(seeds):
        return spans_of(seeds)
    return slice_at(seeds, (STRIDE // 2) % len(seeds), width)


# Three width-34 windows at stride 34 over the class. 3 * 34 = 102 against 100 seeds, so this is
# NOT an exact partition: the third window WRAPS and overlaps the first by 2 seeds (h9 takes
# 400-401 + 468-499). Union is still 100/100 and 98 of the 100 sit in exactly one window.
#
# Handing all three the whole class would not have worked -- `conjunction.py` names it directly:
# hotkeys given the same window and the same bank build the SAME band, the correlation that sank
# all-cut at fleet level. Slicing costs almost nothing here because all-HDR's clean band barely
# moves with window width (measured 14-15 at width 20, 14 at 50, 13 at 100), so three width-34
# windows hold roughly three independent 13-14 seed bands instead of one.
HDR_ONLY_SUB_WIDTH = 34
HDR_ONLY_STRIDE = 34
# 2026-09-23: h7's band space is 400-499, NOT the conjunction's 200-299. That is the point of it.
# The conjunction seven band inside 200-299 and reach 400-499 only through their CUT, so no band
# in the fleet covers the second class at all. h7 fills exactly that gap, which also makes its
# band disjoint from all seven conjunction bands by construction rather than by luck -- the
# overlap the 200-299 placement would have had is gone.
# Moved 400-499 -> 500-599 with the rest of the layout; still a class no conjunction band reaches,
# so the three spike hotkeys stay disjoint from the seven by construction.
HDR_ONLY_CLASSES: tuple[int, ...] = (1,)

assert not (set(HDR_ONLY_HK) & set(LOOP_HK)), \
    "a hotkey cannot both play the conjunction layout and be all-HDR-only"
assert all(0 <= c <= 8 for c in HDR_ONLY_CLASSES), "an all-HDR-only class is outside 100-999"
assert len(set(HDR_ONLY_HK)) == len(HDR_ONLY_HK), "a hotkey is listed twice in HDR_ONLY_HK"
# Coverage, not an exact partition: `slice_at` wraps, so an overshoot is fine and only an
# UNDERSHOOT would leave seeds no hotkey can band. A full revolution would put two hotkeys on the
# identical slice, which the second assert rules out.
assert len(HDR_ONLY_HK) * HDR_ONLY_STRIDE >= len(expand(HDR_ONLY_CLASSES)), \
    (f"{len(HDR_ONLY_HK)} hotkeys x stride {HDR_ONLY_STRIDE} does not reach "
     f"{len(expand(HDR_ONLY_CLASSES))} seeds: some seed sits in no hotkey's window")
assert (len(HDR_ONLY_HK) - 1) * HDR_ONLY_STRIDE < len(expand(HDR_ONLY_CLASSES)), \
    "the last hotkey's offset wraps a full revolution onto another hotkey's slice"
assert HDR_ONLY_SUB_WIDTH <= len(expand(HDR_ONLY_CLASSES)), "the sub window is wider than the class"
assert set(HDR_ONLY_CLASSES) <= set(CUT_CLASSES), \
    ("the all-HDR-only band space must sit inside the cut classes, or its band seeds fall outside "
     "every cut the fleet builds")

HDR_HK: list[str] = LOOP_HK
HDR_CLASSES: tuple[int, ...] = BAND_CLASSES
HDR_SUB_WIDTH = BAND_SUB_WIDTH
HDR_STRIDE = BAND_STRIDE

assert len(HDR_HK) * HDR_STRIDE >= len(expand(HDR_CLASSES)), \
    (f"{len(HDR_HK)} hotkeys x stride {HDR_STRIDE} does not reach "
     f"{len(expand(HDR_CLASSES))} seeds: some seed sits in no hotkey's window")
assert len(set(HDR_HK)) == len(HDR_HK), "a hotkey is listed twice in HDR_HK"
assert HDR_SUB_WIDTH <= len(expand(HDR_CLASSES)), "the sub window is wider than the space"


def _sub_window(seeds, width, frac):
    """`conjunction.sub_window`, imported lazily so this module stays import-cycle free."""
    from niome_subnet.genomics.conjunction import sub_window
    return sub_window(seeds, width, frac)


def hdr_space() -> list[int]:
    """The 200 seeds every all-HDR hotkey draws its window from."""
    return expand(HDR_CLASSES)


def hdr_window_for(instance, cell=None):
    """This hotkey's all-HDR band space as [[lo, hi], ...] spans, or None if it plays none.

    None is the signal for "not in this layout", which keeps the caller's own plan/env fallback
    rather than handing back a guess -- the same contract `fallback_for` has.
    """
    lay = layout_for(cell)
    group = hdr_hk(cell)
    if instance in group:
        # A slice of `HDR_ONLY_CLASSES`, a class no conjunction band reaches. The slices partition
        # the class, so these hotkeys are decorrelated from each other by WINDOW and from the
        # conjunction group by both window and construction. At one hotkey the slice is the whole
        # class.
        only = expand(HDR_ONLY_CLASSES)
        w = min(lay["hdr_width"], len(only))
        return slice_at(only, (group.index(instance) * lay["hdr_stride"]) % len(only), w)
    if instance in conj_hk(cell):
        # all-HDR is the FALLBACK under the conjunction, and it must land on the SAME window or a
        # decline silently re-correlates the siblings onto a different tiling.
        return spans_of(_sub_window(band_window(),
                                    band_sub_width(instance, cell) or BAND_SUB_WIDTH,
                                    band_offset_frac(instance, cell) or 0.0))
    return None


def fallback_for(cell, instance):
    """This hotkey's window when no usable plan exists -> [[lo, hi], ...], or None if unknown.

    Deterministic and offline: same construction as the plan, fixed classes, width from the cell.
    Returns None for a hotkey that is in neither list, so the caller keeps its own last resort
    rather than building on a guess.
    """
    seeds = (expand(range(9)) if FALLBACK_SPREAD else expand(FALLBACK_CLASSES))
    width = sub_width(cell)
    if instance in FULL_HK:
        return half_stride(seeds, FULL_SUB_WIDTH)
    if instance in ROTATE_HK:
        return slice_at(seeds, (ROTATE_HK.index(instance) * STRIDE) % len(seeds), width)
    return None
