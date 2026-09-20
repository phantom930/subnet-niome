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
BAND_HK: list[str] = ["niome_hotkey", "niome_hotkey4", "niome_hotkey5", "niome_hotkey6",
                      "niome_hotkey7", "niome_hotkey8", "niome_hotkey9"]
BAND_STRIDE = 100
# The band space is now the FULL 100-999, not the plan's three predicted classes. `band_offset_frac`
# is expressed against this, so the fractions run 0, 1/9, ... 6/9 and land the windows exactly where
# the comment above says.
BAND_SPAN = 900

# Group B -- empty since 2026-09-19. The complement tier is meaningless now that group A's own band
# space IS the whole 900; there is no complement left to cover.
REST_HK: list[str] = []
REST_STRIDE = 100
REST_SPAN = 0

# The whole seed space a round can draw from. The band split above and the wide cut window below are
# both expressed against this one list so they cannot drift apart.
FULL_SPACE = list(range(100, 1000))

# Per-hotkey offset overrides. EMPTY and should stay so: the sequence below already places seven
# windows across 900 with no wrap, and an override would pull one off that.
BAND_OFFSET_OVERRIDES: dict[str, int] = {}

assert (len(BAND_HK) - 1) * BAND_STRIDE + BAND_SUB_WIDTH == BAND_SPAN, \
    "group A's windows do not span 100-999 exactly (overlapping, non-wrapping tiling)"
assert len(REST_HK) * REST_STRIDE == REST_SPAN, "group B does not tile its span exactly"
assert BAND_SPAN <= len(FULL_SPACE), "the band span claims more than 100-999"
assert not (set(BAND_HK) & set(REST_HK)), "a hotkey cannot be in both band groups"


def band_layout(instance):
    """``(hotkeys, stride, span, on_predicted)`` for a band hotkey, else None.

    ``on_predicted`` is what separates the two groups: True means this hotkey's sub-window is cut
    out of the plan's three predicted classes, False means out of the 600 seeds the plan left over.
    """
    if instance in BAND_HK:
        return BAND_HK, BAND_STRIDE, BAND_SPAN, True
    if instance in REST_HK:
        return REST_HK, REST_STRIDE, REST_SPAN, False
    return None


def band_offset(instance):
    """This hotkey's band sub-window offset, in seeds of ITS OWN space, or None."""
    layout = band_layout(instance)
    if layout is None:
        return None
    hotkeys, stride, span, _ = layout
    off = BAND_OFFSET_OVERRIDES.get(instance, hotkeys.index(instance) * stride)
    # No modulo. Since 2026-09-20 the windows are PLACED across 100-999 rather than rotated inside
    # a smaller space, so an offset past `span - width` would wrap `sub_window` back to the front
    # and collide with the first hotkey instead of erroring. The import-time assertion above is what
    # guarantees it cannot happen; this raises rather than wrapping if it ever does.
    if off + BAND_SUB_WIDTH > span:
        raise ValueError(f"band offset {off} + width {BAND_SUB_WIDTH} runs past the {span}-seed "
                         f"band space for {instance}; the layout no longer tiles")
    return off


def band_offset_frac(instance):
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
    layout = band_layout(instance)
    if layout is None:
        return None
    return band_offset(instance) / float(layout[2])


def band_space(instance, predicted):
    """The seeds this hotkey cuts its band sub-window out of, or None if it plays no band.

    ``predicted`` is the plan's joined space — the three predicted width-100 classes, 300 seeds, as
    a sorted seed list. A group A hotkey draws from exactly that; a group B hotkey draws from its
    complement in 100-999, the 600 seeds no group A hotkey can reach.

    Returning the complement here rather than at the call site is what keeps the two groups from
    drifting: the 600 is DEFINED as "whatever the plan did not predict", so a plan that moves its
    classes moves both groups together and the fleet still covers 100-999 exactly once per group.
    """
    layout = band_layout(instance)
    if layout is None:
        return None
    # 2026-09-20: the band space is the FULL 100-999 for every band hotkey. `predicted` is accepted
    # and ignored -- the plan's three classes no longer steer the band at all, only the all-HDR rung
    # below the conjunction still reads them. Kept in the signature so callers do not have to know
    # which layout is live.
    if layout[3]:
        return list(FULL_SPACE)
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
CUT_MODE: dict[str, str] = {
    "HEK293": "union",
    "CD34+_HSPC": "union",
    "K562": "union",
    "HUDEP-2": "union",
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
    if band_layout(instance) is None:
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
