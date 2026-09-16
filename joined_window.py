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
FULL_HK: list[str] = ["niome_hotkey", "niome_hotkey1", "niome_hotkey2", "niome_hotkey3"]

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
BAND_SUB_WIDTH = 150
BAND_STRIDE = 75
BAND_HK: list[str] = ["niome_hotkey", "niome_hotkey1", "niome_hotkey2", "niome_hotkey3"]
# The joined space these offsets are expressed against: 3 classes x 100 seeds. The offset is carried
# as a FRACTION so it stays proportional if a plan ever yields a space of another size.
BAND_SPAN = 300


def band_offset_frac(instance):
    """This hotkey's band sub-window offset, as a fraction of the joined space, or None.

    None means "not a band hotkey here" and the caller keeps `ConjunctionConfig`'s own default,
    which is the single-hotkey value the 12-contract replication was measured at.
    """
    if instance not in BAND_HK:
        return None
    return ((BAND_HK.index(instance) * BAND_STRIDE) % BAND_SPAN) / float(BAND_SPAN)

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
