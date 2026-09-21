"""Submission design: turn a contract into the 250 CRISPR experiments a validator pays most for.

The miner supplies designs only — guide, coordinate, strand, mutation, Cas system, cell type. Every
biological outcome is rolled by the validator under a seed the miner does not hold at build time, so
this module never tries to pick an outcome. It optimises the two terms that are deterministic
functions of the design, and treats the third as a distribution to be shifted rather than chosen.

    final_score = total_weighted_score x consistency_factor x distribution_fidelity_factor

What each term responds to, and what this module does about it:

**total_weighted_score** — sum over rows of
``(0.625*gc_score + 0.375*dist_score) * offtarget_factor * mutation_weight``. Fully determined by
the design, so it is maximised exactly: GC pinned to 50% (``gc_score`` 1.0), the nearest usable
PAM to each mutation (``dist_score`` ~1.0), and the mismatch budget spent pushing the guide's 12-mer
off-target seed out of the validator's index, which is what takes ``offtarget_factor`` from 0.7 to
1.0 — a flat 1.43x on every row that nothing else in the pipeline charges for.

**distribution_fidelity_factor** — a *geometric* mean of six coverage ratios, so an empty
(mutation, cas, strand) cell multiplies the whole score by ~1e-9. Every cell is therefore occupied,
and the row split across cells is chosen by hill-climbing the exact closed form of
``term1 x geomean(coverage ratios)`` — no simulation needed, because neither factor reads the seed.

**consistency_factor** — ``0.7*max(avg_r2, 0) + 0.3*(1 - avg_nmae)`` over a RandomForest fitted to
``is_cut``, ``is_hdr`` and ``indel_length``. Two facts decide the strategy here, both measured
against the validator's own stage 4 rather than assumed:

1. ``avg_r2`` is negative for any seed-blind design, so the 0.7 term contributes nothing. It turns
   positive only through a degeneracy: if no row draws ``no_cut``, ``is_cut`` is a constant column,
   ``r2_score`` returns 1.0 and ``normalized_mae`` short-circuits to 0. **This generator does not
   pursue that degeneracy, and the reason recorded here was wrong.** It said the degeneracy was
   unreachable at HEK293's accessibility of 0.35, from the fact that a 50%-GC guide cuts with
   probability 0.95 (Cas9) / 0.87 (Cas12a) and 250 rows all cutting is then ~1e-10. That figure is
   right and the conclusion drawn from it is not: ``cut_probability`` is
   ``min(0.99, base + 0.18*energy)`` and the 0.99 cap needs only ``energy >= 0.7222``, which
   ``accessibility*(1.8*gc + 0.6*exp(-d/1500))`` reaches at **gc ~ 0.855** even at 0.35. The
   marginal gc trade was measured at gc 0.5/0.6/0.7, found monotonically losing, and the search
   stopped there — but the payoff is a step, not a gradient, and the step is past gc 0.8:

       gc     p(Cas9)   gc_score   base structural   P(all 250 cut)
       0.500   0.9498      1.000        0.852            2.55e-06
       0.700   0.9725      0.600        0.602            9.29e-04
       0.800   0.9838      0.400        0.477            1.69e-02
       0.855   0.9900      0.291        0.409            8.11e-02

   The backend's own records confirm it is reachable: on task de19c2e0 (this contract) 102 miners
   scored non-zero, and six of them landed in ``consistency_factor`` 0.32-0.43 against a bulk
   cluster at 0.06-0.12, the best at **0.4339** with ``total_weighted_score`` 273.81 and a
   ``final_score`` of 109.4. So the honest description of this design is that it takes the safe
   0.09 rather than buying a ~2-8% chance of 0.33+ at roughly half the structural score — a
   defensible choice against a *median* of 21.1, but a choice, not a constraint. Note the best
   miner's term1 is only 16% below this generator's, which the flat gc~0.855 story does not
   explain, so their actual route is still unidentified.
2. What *is* reachable is the ``0.3*(1 - avg_nmae)`` term, and two things move it. Collapsing the
   feature matrix to **one vector per cell** — every row in a cell sits at the same coordinate at
   the same GC count — leaves the forest with eight groups instead of 250 noisy points, so it can
   only return group means and stops overfitting the targets it cannot predict: measured r2 on
   ``is_hdr``/``indel_length`` improves from −0.25/−0.31 to −0.12/−0.11. And because ``is_cut``'s
   normalised error falls as the cut rate rises, the **row allocation prices each cell's cut
   probability** rather than treating the split as a pure coverage question — Cas9 cuts with
   probability 0.95 against Cas12a's 0.87, so the optimum leans to roughly 30% Cas12a and pays for
   the coverage entropy that costs. Worth 8.8% over a coverage-only allocation.

3. The third lever is the ``0.7 * max(avg_r2, 0)`` term, which this file previously treated as
   unreachable. It is reachable, and it does not have to be bought with ``gc_score`` at all.

   ``consistency_factor`` does not pay for a high cut *rate*; it pays for a perfectly constant
   ``is_cut`` column. A fold whose ``y_test`` never varies has zero total sum of squares, and
   ``r2_score`` answers that 0/0 with 1.0 when the prediction is exact and 0.0 when it is not — so
   the target is all-or-nothing. Measured at 250 rows, going from every row cutting to one row
   failing takes the factor from **0.2394 to 0.1008**. That is the step the gc ~ 0.855 route above
   is really buying, which is why the 0.4339 miner's ``total_weighted_score`` was only 16% below
   this generator's instead of half: they were holding a constant ``is_cut``, not paying for it in
   GC.

   Ranking guides by the share of ``SEED_SUPPORT`` they cut under was built here once and correctly
   measured as noise (0.9285 against 0.9295). The statistic was the wrong one. Selection that serves
   an all-or-nothing target is not "the guides that cut most often" but "a set of guides that all cut
   under the *same* seeds" — an intersection, not a ranking. The earlier experiment also held the
   evaluation seeds out of the scanned support, which given that ``experiment_seed`` hashes the round
   seed in with the row makes a null result arithmetic rather than evidence: outcomes under two seeds
   are independent, so there is nothing to generalise to an unscanned seed.

   What makes the lever real is that ``SEED_SUPPORT`` is enumerable. The stamped seed is one of its
   900 members, so a design whose rows all cut under N of them holds the cliff with probability
   N/900 — and that is a property the build can compute before uploading. ``select_for_diversity``
   now maximises N by greedy intersection: on the reference task it goes from **1 seed to 215**,
   which prices out at roughly +33% expected ``consistency_factor``. Above accessibility ~0.72
   ``energy`` clamps at 1.0 and Cas9 sits at its 0.99 cap at *any* GC, so on K562 (0.77), HUDEP-2
   (0.82) and CD34+_HSPC (0.87) — about three quarters of tasks — this costs nothing in
   ``gc_score``. At HEK293's 0.35 the intersection is much harder to hold over 250 rows and the
   scan finds far less of it, but it is not empty — measured 36 seeds of 900, worth +13.8%.

   The scan then re-opens a decision the allocator had already made, which is the second half of
   the lever. ``allocate_rows`` prices a cell's cut probability only through ``is_cut``'s
   *normalised error*, a term that moves a few hundredths. The intersection is far more sensitive
   to the same knob because it compounds once per row: a Cas12a row at 0.96 erodes the surviving
   seed set four times faster than a Cas9 row at 0.99. So the mix that is right for the other three
   factors is too Cas12a-heavy for this one, and ``retune_cas_mix`` re-prices it against the
   measured intersection. It lands near 18 Cas12a rows of 250 against the allocator's 86-89, and
   that is where most of the gain is: paired over 120 held-out seeds on K562, the scan alone takes
   ``consistency_factor`` 0.110 -> 0.154 and the retune takes it to **0.244**, with ``final_score``
   27.12 -> 37.72 -> **48.88**. It costs ``distribution_fidelity`` 0.920 -> 0.751, which is a large
   and deliberate trade: coverage entropy enters the score at the 1/6 power, the intersection does
   not.

4. ``retune_cas_mix`` moves the whole submission's cas mix at once and asks one greedy intersection
   to hold the seed set over every row simultaneously. There is a second, more effective shape of
   the same idea (measured on ``origin/develop``'s ``all_cut.py``/``fastgreedy.py`` -- read via
   ``git show origin/develop:...`` rather than vendored, since this repo has no prefetch loop, no
   GPU, and every contract's mutation set is unique, so develop's disk-cached bank and ~900s
   prepare budget do not transfer):

     a. min-union a GROUP of weak-cas candidates (``min_union_select``) -- the complement of their
        combined failures is a *clean set*, usually far smaller than the full seed support but far
        easier to be strict over.
     b. fill the strong-cas cells with candidates that are strict over that clean set (zero fails
        within it, not over all of ``seed_support``) -- the relaxation that makes strictness
        reachable at all.

   The reachability of step b is what a naive port stalls on. Every guide substituted at one PAM
   coordinate shares that coordinate's gc/distance/energy (``stage3.sequence_energy`` depends on
   nothing else a substitution can move), so ``cut_probability`` is identical across them and
   strictness over k clean seeds is close to an independent ``cut_probability**k`` per candidate --
   pooling more *variants at the same coordinate* does not raise it. Measured directly: on
   CD34+_HSPC, after a weak-cas min-union left 608 of 900 seeds clean, only 0-3 of 900 Cas9
   candidates AT ONE COORDINATE were strict over that clean set, against the 82-93 needed. Only a
   different coordinate moves ``distance_to_mutation`` and therefore energy's
   ``0.6*exp(-d/1500)`` term enough to matter (0.990 vs 0.995 changes ``p**k`` by an order of
   magnitude at k in the hundreds), which is why ``_widen_strong_cas_pool`` grows coordinate COUNT
   for the strong-cas cells specifically rather than variants per coordinate -- the ordinary growth
   loop in ``build`` already maximises the latter and stops once the total pool clears
   ``rows_wanted``, far short of what step b needs.

   ``two_stage_construction`` is gated behind ``Config.two_stage_enabled`` and a margin check
   (``_TWO_STAGE_MIN_GAIN``) against whichever of the plain scan or ``retune_cas_mix`` already won,
   using the same closed-form ``_predicted_objective`` both paths are scored by -- so turning it on
   can only ever match or improve the shipped result, never regress it, and any exception during
   the attempt falls back to the already-verified behaviour untouched.

   Whether it *does* improve the result turned out to be a scan-width question, not an algorithm
   one, and a real bug lived here: truncating the widened strong-cas pool by ``weighted_score``
   alone (what the ordinary scan correctly does for a single-coordinate cell) collapsed a
   57,600-candidate, 64-coordinate pool back onto its one or two nearest coordinates before
   ``two_stage_construction`` ever saw it -- the exact single-cut-probability degeneracy the
   widening exists to escape. Fixed by ``_round_robin_by_coordinate``
   (``scan_cut_support(..., rank="coordinate_diverse")``), which orders the pool so a truncation to
   any width still sees every coordinate before any coordinate sees a second candidate.

   With that fixed, the remaining bottleneck is purely wall-clock: ``_widen_strong_cas_pool`` finds
   the 57,600 candidates in seconds, but scanning enough of them to matter is not free. Measured on
   CD34+_HSPC, rescan budget against outcome (see ``_TWO_STAGE_MIN_GAIN``'s own comment for the
   table): 45s loses, 90s is barely positive but still below the gate, 180s clears it at objective
   ratio 1.087. The construction can genuinely beat the shipped baseline, but only past a rescan
   budget that leaves little of the 300s TTL for the upload once stacked on the ordinary scan's own
   ~90s -- a build-time-vs-upload-margin trade this file does not make unilaterally.
   ``Config.two_stage_seconds``/``Miner.TWO_STAGE_SHARE_OF_WINDOW`` are left conservative, so today
   the gate reliably (and correctly) declines.

5. Everything in 1-4 is a way of playing odds against a seed the design cannot see, and all of it is
   conditional on that. A contract that *carries* its round seeds is a different problem, not a
   better-informed version of the same one: stage 3 becomes readable, and the right move is to stop
   estimating ``consistency_factor`` and take all of it.

   ``pinned_outcome_build`` selects, per cell, only candidates that draw the *same*
   ``(outcome, indel_length)`` under every one of the contract's seeds. That makes ``is_cut``,
   ``is_hdr`` and ``indel_length`` all constant columns, which is the same 0/0 degeneracy point 3
   describes — but applied to all three targets at once rather than to ``is_cut`` alone, so
   ``consistency_score`` is ``(0.7*1.0 + 0.3*(1 - 0)) * 100`` = 100 and ``consistency_factor`` is
   exactly **1.0**, not a fraction of it.

   It is close to free. Every candidate in a cell shares that cell's coordinate, GC count and
   off-target cleanliness, so they all carry an identical ``weighted_score`` and filtering to one
   outcome cannot cost ``total_weighted_score`` anything — only a cell's remaining row supply, which
   is what the objective check prices. Measured end to end through ``benchmark_submission`` on the
   recorded task history, one task per cell type per seed count, seeds handed to both halves:

       task      cell type    seeds          pinned outcome     rows     consistency   final
       62cc85fb  HEK293 0.35  630,765,543    ('HDR', 0)         250/250     1.0000     301.46
       8fbd60b7  HEK293 0.35  497            ('HDR', 0)         250/250     1.0000     260.76
       7ac3ef3a  K562   0.77  752,452,560    ('HDR', 0)         250/250     1.0000     241.70
       e0c604a4  K562   0.77  342            ('BLUNT_NHEJ', 1)  250/250     1.0000     246.44
       f6686c85  HUDEP2 0.82  269,862,653    ('HDR', 0)         250/250     1.0000     355.50
       e19e23f8  HUDEP2 0.82  477            ('HDR', 0)         250/250     1.0000     249.35
       65df5469  CD34+  0.87  364,422,810    ('HDR', 0)         250/250     1.0000     249.06
       8b7bbba1  CD34+  0.87  686            ('BLUNT_NHEJ', 1)  250/250     1.0000     251.35

   Two of those pin to ``BLUNT_NHEJ`` rather than ``HDR``, which is why the target is not
   pre-committed: every ``(outcome, indel_length)`` that can occupy all eight cells is priced and the
   best one wins. Accessibility barely matters — it moves which outcomes are *common*, not whether
   one can be held, because the guides at a coordinate differ in their stage-3 draw and nothing else.

   This does not change what a live round does. The backend stamps the seed *after* the broadcast,
   so a production contract arrives with ``seed: 0``, ``Context.seeds()`` is empty and the whole
   construction is skipped — points 1-4 are still what ships. It is reached only by a contract that
   carries seeds, and ``Config.pin_outcomes`` turns it off. Where it does apply it also *replaces*
   the 90 s cut-support scan rather than adding to it (0.4 s total build against ~90 s), because a
   scan can only improve the odds of a factor this already holds outright.

Every formula below is imported from the validation stages rather than reimplemented, so the
generator cannot drift from the pipeline that judges it. On the reference contract the result scores
27.4 against the 21.1 the same contract previously paid, confirmed against the validator's own
``benchmark_submission`` to zero difference.
"""

from __future__ import annotations

import itertools
import json
import logging
import math
import random
import time

from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import niome_subnet.utils.settings as settings

from niome_subnet.genomics.validation import stage12, stage3, stage5

logger = logging.getLogger(__name__)

# The 12-mer off-target window and the k-mer index flank are hard-coded at the validator's call
# sites (``run_stage12(offtarget_flank=50000)``, ``k=12``), not contract fields.
KMER_LENGTH = 12
OFFTARGET_FLANK = 50_000

# Complement pairs used for guide substitution. A within-class swap has exactly one alternative,
# which is what makes the variant enumeration deterministic.
_WITHIN_CLASS_SWAP = {"A": "T", "T": "A", "G": "C", "C": "G"}
_GC_BASES = ("G", "C")
_AT_BASES = ("A", "T")

# Stage 3's seed is stamped by the backend after the task is broadcast, so a miner only ever sees
# ``seed: 0`` and cannot design against the real one. The observed seeds are three digits, so this
# is the support a seed is drawn from — the population the design is built to hold up across,
# rather than any one lucky draw. It is small enough to enumerate, which is the whole basis of the
# cut-support selection in ``select_for_diversity``: the stamped seed is one of these 900.
SEED_SUPPORT = tuple(range(100, 1000))

# Typical normalised MAE for the two targets no seed-blind design can predict, measured over ten
# held-out seeds on this contract (is_hdr 0.97, indel_length 0.69). They enter the allocation
# objective as constants because they barely move with the row split — unlike is_cut's, which is a
# direct function of the Cas9/Cas12a mix and is what ``_consistency_estimate`` computes.
NMAE_IS_HDR = 0.97
NMAE_INDEL = 0.69


@dataclass
class Config:
    """Design knobs. Defaults are what the objective analysis above settles on."""

    # Guide lengths to enumerate. Stage 1 accepts 20 and 23; 20 is preferred because an even length
    # can hit exactly 50% GC (gc_score 1.0) while 23 tops out at 0.957.
    guide_lengths: tuple[int, ...] = (20, 23)
    # How far either side of gene_region to look for PAMs. Only the nearest handful of coordinates
    # per cell are ever used, so this exists to guarantee a cell is never starved.
    pam_search_flank: int = 4000
    # Coordinates per cell. One is the point of the design: every row in a cell then shares a single
    # stage-2 feature vector, so stage 4's forest sees a handful of groups instead of 250 noisy
    # points and can only return group means. Raised automatically for a cell that cannot supply
    # enough distinct guides from one coordinate.
    coordinates_per_cell: int = 1
    # Guide variants enumerated per coordinate before ranking. The pool is what the k-mer
    # diversity pass chooses from, so a bigger pool is a better submission and a longer build.
    guides_per_coordinate: int = 900
    # Row cap override; None follows the contract's rules.max_experiments.
    row_cap: int | None = None
    # Seeds the cut-support scan intersects over. Every member is a seed the backend could stamp,
    # and selection cannot generalise beyond what is scanned, so narrowing this only narrows the
    # set of rounds the design holds a constant ``is_cut`` on. Empty disables the scan entirely.
    seed_support: tuple[int, ...] = SEED_SUPPORT
    # Use ``pinned_outcome_build`` when the contract hands the design the round seeds it will be
    # scored under. Everything ``seed_support`` above exists for is a way of playing odds against a
    # seed the miner cannot see; a contract that carries its seeds makes stage 3 readable instead,
    # and ``consistency_factor`` goes from a fraction to exactly 1.0. Off falls back to the
    # seed-blind path unmodified, which is also what an unstamped contract (``seed: 0``) gets.
    pin_outcomes: bool = True
    # Where ``read_seed_file`` looks for seeds when the contract carries none. None means
    # ``settings.MINER_SEEDS_PATH``; point it at a missing file to build seed-blind regardless of
    # what is on disk.
    seeds_path: str | None = None
    # Wall-clock ceiling for growing the candidate pool so a pin over many seeds can still fill
    # every row (see ``pin_with_widening``). Each 4x step costs one more ``_build_pools`` pass:
    # measured ~6s at 16 coordinates/cell and ~25s at 64, against a 300s upload window. Spent only
    # when the pin cannot already fill the submission, so at the 1-3 seeds a live round actually
    # issues this is never reached. Zero disables widening.
    pin_widen_seconds: float = 60.0
    # Candidates per cell the scan prices. More is monotonically better and monotonically slower —
    # measured 100/200/300/450/600/900 -> 66/123/162/184/199/215 seeds held on the reference task —
    # so this is capped by ``seed_scan_seconds`` below rather than tuned.
    seed_scan_candidates_per_cell: int = 900
    # Wall-clock ceiling for the scan. The upload window is 300 s (``SUBMISSION_TIMEOUT``) and the
    # rest of the build is under a second with chr11 prewarmed, but the scan is ~900 Mersenne
    # seedings per candidate and how fast that runs is a property of the host. The width is priced
    # against a live sample rather than assumed, so a slow box scans narrower instead of missing
    # the window.
    seed_scan_seconds: float = 90.0
    # Try ``two_stage_construction`` (see the module docstring's coordinate-diversity note) and
    # take it over ``retune_cas_mix`` only when it predicts a clearly better score. Off falls back
    # to exactly today's shipped scan-and-retune behaviour, unmodified -- the A/B knob a caller
    # flips to compare the two.
    two_stage_enabled: bool = True
    # Wall-clock ceiling for widening the strong-cas coordinate pool (see
    # ``_widen_strong_cas_pool``). Separate from ``seed_scan_seconds`` because it is spent before
    # any scanning happens and on a different cost (stage 1/2 gating of new coordinates, not
    # Mersenne seedings), so the two cannot share one budget without one starving the other.
    two_stage_seconds: float = 45.0


@dataclass
class Context:
    """Everything the design needs that is not a knob: the genome, the index, the contract."""

    sequence: str
    contract: dict
    reference: dict
    cell_types: dict
    kmer_index: dict

    @property
    def mutations(self) -> list[str]:
        return list(self.contract["active_mutations"])

    @property
    def cas_systems(self) -> list[str]:
        return list(self.contract["rules"].get("cas_systems") or ["Cas9", "Cas12a"])

    @property
    def strands(self) -> tuple[str, str]:
        return ("+", "-")

    @property
    def mutation_map(self) -> dict:
        return self.reference["mutation_map"]

    @property
    def max_experiments(self) -> int:
        return self.contract["rules"].get("max_experiments") or 250

    @property
    def max_mismatches(self) -> int:
        return int(self.contract["rules"].get("max_mismatches") or 0)

    @property
    def base_padding(self) -> int:
        return int(self.contract["rules"]["base_padding"])

    @property
    def proximity_gate(self) -> bool:
        return bool(self.contract["rules"].get("proximity_gate", False))

    def weight_of(self, mutation: str) -> float:
        return float(self.contract.get("mutation_weights", {}).get(mutation, 1.0))

    def seeds(self) -> list[int]:
        """The round seeds in the contract, or [] when it is unstamped."""
        return parse_seeds(self.contract)


def parse_seeds(contract: dict) -> list[int]:
    """The round seeds in a contract, or ``[]`` when it is unstamped.

    ``benchmark_submission`` reads the field as a comma-joined string, so a multi-seed round
    averages several draws. A stamped contract is scored exactly; an unstamped one (``seed: 0``)
    has nothing to score against and falls back to ``Config.seed_support``.

    Module level rather than only a ``Context`` method because a caller that has a contract but no
    genome — the miner, deciding what to stamp before it builds anything — needs the same answer,
    and must not get it from a second parser that could disagree with this one.
    """
    raw_seed_field = contract.get("seed")
    try:
        parsed_seeds = [
            int(seed_text) for seed_text in str(raw_seed_field).split(",") if seed_text.strip()
        ]
    except (TypeError, ValueError):
        return []
    # A zero is the backend's placeholder for "not stamped yet", not a usable round seed.
    return [seed for seed in parsed_seeds if seed]


def read_seed_file(path: str | None = None) -> tuple[list[int], bool]:
    """Round seeds supplied out of band, and whether the file is switched on.

    Returns ``(seeds, enabled)``. ``enabled`` is the file's ``"enabled"`` flag, defaulting to True
    when absent — a file that exists without saying otherwise is meant to be used. It is reported
    separately from the seeds rather than folded into them because "switched off" and "no seeds
    listed" are different instructions: the first means mine seed-blind, the second leaves the
    miner free to draw one. Collapsing them would make the off switch turn the pin *on*, aimed at
    a random seed.

    Lives here rather than on the miner because the miner is not the only thing that builds: any
    harness driving ``build`` goes through this module and no other, so a seed source attached to
    the neuron is one the rest of the system silently ignores. Read fresh on every call, so the
    file can be edited without restarting a long-lived miner.

    Absent is the normal case and not worth a warning; malformed is worth one, because the
    difference between "no file" and "a file I could not read" is the difference between a
    deliberate fallback and a typo that quietly costs the round its pin.
    """
    path = path or settings.MINER_SEEDS_PATH
    if not Path(path).exists():
        return [], True
    try:
        with open(path) as seeds_file:
            document = json.load(seeds_file)
    except Exception as error:
        logger.warning(f"Could not read {path} ({error}); continuing with no supplied seeds")
        return [], True

    enabled = bool(document.get("enabled", True)) if isinstance(document, dict) else True
    if not enabled:
        return [], False

    raw = document.get("seeds") if isinstance(document, dict) else document
    if raw is None:
        return [], True
    # A list of numbers, or the contract's own comma-joined string. Zero is the backend's
    # "not stamped yet" placeholder rather than a seed, so it is dropped either way.
    values = raw if isinstance(raw, list) else str(raw).split(",")
    seeds = []
    for value in values:
        try:
            seed = int(str(value).strip())
        except (TypeError, ValueError):
            logger.warning(f"Ignoring unparseable seed {value!r} in {path}")
            continue
        if seed:
            seeds.append(seed)
    outside = [seed for seed in seeds if seed not in SEED_SUPPORT]
    if outside:
        # Not an error — a validator scores under whatever it was given, and the design is built
        # for exactly the seeds named here. But every seed observed from the backend has been three
        # digits, so this is nearly always a typo, and a typo costs the whole round its pin.
        logger.warning(
            f"Seeds {outside} in {path} are outside the observed support "
            f"{SEED_SUPPORT[0]}-{SEED_SUPPORT[-1]}; building against them anyway"
        )
    return seeds, True


def read_seed_count_override(path: str | None = None) -> int | None:
    """``seeds.json``'s ``"seed_count"``, or None when it does not set one.

    Zero is meaningful and is the reason this is separate from the seed list: it says "bet on no
    seeds", i.e. mine seed-blind, which is not the same as listing none. Returned as None when
    absent so the caller's measured per-cell-type default stands.
    """
    path = path or settings.MINER_SEEDS_PATH
    if not Path(path).exists():
        return None
    try:
        with open(path) as seeds_file:
            document = json.load(seeds_file)
    except Exception:
        return None  # read_seed_file already warned about this file
    if not isinstance(document, dict) or "seed_count" not in document:
        return None
    try:
        return max(0, int(document["seed_count"]))
    except (TypeError, ValueError):
        logger.warning(f"Ignoring unparseable seed_count {document['seed_count']!r} in {path}")
        return None


# How many seeds to bet on, by cell type — the largest count that still fills all 250 rows at
# consistency_factor 1.0 inside the build budget. It is a property of the cell type because
# accessibility sets stage 3's energy, which sets how often a guide draws the pinned outcome: a
# candidate survives only if it draws that outcome under *every* seed, so a cell's usable capacity
# falls like ``pool * p**k`` while the 250 rows needed stay flat.
#
# Measured at 64 coordinates/cell (57,600 candidates), the widest pool ``pin_with_widening`` grows
# to, with ``('HDR', 0)`` the winning target at every k on every cell type:
#
#     cell type    accessibility   k filling 250   first k that falls short
#     HEK293           0.35              8          9 -> 213 rows, 10 -> pin fails
#     K562             0.77             11         12 -> 185 rows
#     HUDEP-2          0.82             11         12 -> 167 rows
#     CD34+_HSPC       0.87             11         12 -> 185 rows
SEEDS_BY_CELL_TYPE = {
    "HEK293": 8,
    "K562": 11,
    "HUDEP-2": 11,
    "CD34+_HSPC": 11,
}
# For a cell type the table does not name — a new one the backend starts issuing. The lowest
# measured value, since an unknown type may be less accessible than any of these. Overshooting is
# recoverable rather than fatal (``pin_with_widening`` drops seeds until the rows fit) but costs
# build time, so re-measure and add the row rather than leaning on this.
FALLBACK_SEED_COUNT = 8


def seed_count_for(contract: dict, seeds_path: str | None = None) -> int:
    """How many seeds to bet on for this contract. ``seed_count`` in the file overrides the table."""
    override = read_seed_count_override(seeds_path)
    if override is not None:
        return override
    cell_type = contract.get("cell_type")
    count = SEEDS_BY_CELL_TYPE.get(cell_type)
    if count is None:
        logger.warning(
            "No measured seed count for cell type %r; using the conservative fallback of %d. "
            "Re-measure and add it to design.SEEDS_BY_CELL_TYPE.", cell_type, FALLBACK_SEED_COUNT,
        )
        return FALLBACK_SEED_COUNT
    return count


def plan_seeds(contract: dict, seeds_path: str | None = None) -> tuple[list[int], str]:
    """Which seeds to bet on this round, in priority order, recording the choice.

    Lives here rather than on the miner for the reason ``read_seed_file`` does: the neuron is not
    the only thing that builds. A planner attached to it is skipped by every harness that drives
    ``build`` directly, which is both a silently different design and — since this is what writes
    ``last_generated`` — a record that never appears.

    The count comes from the cell type (``seed_count_for``). The set is then assembled
    **highest-conviction first**, because ``pin_with_widening`` drops from the end when even the
    widest pool cannot fill the submission:

    1. the seeds listed in ``seeds.json`` when it is ``"enabled": true`` — an operator who lists
       seeds is asserting something the build cannot work out for itself, so they go first and are
       never given up in favour of a random draw;
    2. ``draw_seeds`` to make up the count — a random draw, but over the seeds the backend has
       never been recorded stamping rather than over all of ``SEED_SUPPORT``. A listed seed is
       never dropped for having been stamped before: the operator's claim outranks the ledger.

    ``"enabled": false`` means "ignore the listed seeds", not "do not pin": the round still bets,
    just on its own draws. ``"seed_count": 0`` is the off switch, and so is
    ``Config.pin_outcomes = False``.
    """
    count = seed_count_for(contract, seeds_path)
    if count <= 0:
        logger.info("Seed count is 0 for this contract — building seed-blind")
        return [], "seed-blind"

    listed, enabled = read_seed_file(seeds_path)
    # dict.fromkeys dedupes while keeping the operator's ordering, which is their priority.
    chosen = list(dict.fromkeys(listed)) if enabled else []
    if listed and not enabled:
        logger.info(
            '%s has "enabled": false — ignoring its %d listed seed(s) and drawing all %d',
            seeds_path or settings.MINER_SEEDS_PATH, len(listed), count,
        )
    if len(chosen) > count:
        # Never truncate: a listed seed may be real knowledge, worth more than the rows it might
        # cost. The build drops from the end only if it genuinely cannot fit them.
        logger.warning(
            "%s lists %d seeds, more than the %d this cell type can hold at 250 rows — keeping "
            "them all; the build will drop the lowest-priority ones if it must",
            seeds_path or settings.MINER_SEEDS_PATH, len(chosen), count,
        )
    elif len(chosen) < count:
        chosen += draw_seeds(count - len(chosen), exclude=chosen)

    source = "file+drawn" if (enabled and listed) else "drawn"
    logger.info(
        "Betting on %d seed(s) for %s (%s): %s",
        len(chosen), contract.get("cell_type"), source, chosen,
    )
    record_generated_seeds(contract, chosen, source, seeds_path)
    return chosen, source


def record_generated_seeds(
    contract: dict, seeds: list[int], source: str, seeds_path: str | None = None
) -> None:
    """Note the seeds this round drew, into ``seeds.json``'s ``last_generated``.

    A record, never an input: ``read_seed_file`` reads ``seeds`` and nothing else, so what is
    written here cannot feed back into the next round's plan. That separation is the point — a
    drawn set written into ``seeds`` would silently become permanent priority seeds next task.

    Read-modify-write so an operator's ``seeds``, ``enabled``, ``seed_count`` and notes survive,
    and best-effort throughout: losing the record must never cost the round.
    """
    path = seeds_path or settings.MINER_SEEDS_PATH
    try:
        document = {}
        if Path(path).exists():
            with open(path) as seeds_file:
                loaded = json.load(seeds_file)
            if not isinstance(loaded, dict):
                # A bare list is a valid seeds file but has nowhere to hold a record, and
                # rewriting it as an object would change what the operator wrote.
                return
            document = loaded
        document["last_generated"] = {
            "seeds": seeds,
            "source": source,
            "cell_type": contract.get("cell_type"),
            "at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        # Indented and atomic: this file is edited by hand, so it has to stay readable, and a task
        # arriving mid-write must not leave half a document where the seeds live.
        temporary_path = f"{path}.tmp"
        Path(temporary_path).parent.mkdir(parents=True, exist_ok=True)
        with open(temporary_path, "w") as output_file:
            json.dump(document, output_file, indent=2)
            output_file.write("\n")
        Path(temporary_path).replace(path)
    except Exception as error:
        logger.warning(f"Could not record the generated seeds in {path} ({error})")


def read_drawn_counts(path: str | None = None) -> dict[int, int]:
    """How often the backend has stamped each seed, from the ledger ``record_drawn_seeds`` keeps.

    ``{seed: times_stamped}``, and ``{}`` for a file that is absent or unreadable — which is a
    working state, not an error: ``draw_seeds`` then treats every seed as unstamped and draws from
    the whole support, exactly as it did before the ledger existed. Read fresh on every call so the
    file can be rebuilt or hand-edited without restarting a long-lived miner.
    """
    path = path or settings.MINER_DRAWN_SEEDS_PATH
    if not Path(path).exists():
        return {}
    try:
        with open(path) as ledger_file:
            document = json.load(ledger_file)
        raw = document.get("counts", {}) if isinstance(document, dict) else document
        # A bare list of seeds is accepted too — it is the natural thing to hand-write, and says
        # "these have been used" without claiming to know how often.
        if isinstance(raw, list):
            return {int(seed): 1 for seed in raw}
        return {int(seed): int(times) for seed, times in raw.items()}
    except Exception as error:
        logger.warning(
            f"Could not read the drawn-seed ledger {path} ({error}); drawing from the whole support"
        )
        return {}


def stamped_seeds(tasks: Iterable[dict], since: str | None = None) -> tuple[Counter, list[str], str]:
    """Count the seeds stamped across ``tasks``, as the backend's task history returns them.

    Returns ``(counts, counted_task_ids, latest_created_at)``. A task is counted only once it
    carries a seed: ``seed: 0`` is the placeholder for a round the backend has not stamped yet,
    and counting it would put 0 in the ledger and teach the draw nothing. The ids returned are
    exactly the rounds that contributed, so the caller marking them as folded in cannot disagree
    with the counts — a task skipped by ``since`` must not be recorded as counted, or widening
    ``since`` later would silently find nothing to add.

    ``since`` keeps the ledger inside one generator regime. The round went from one seed to three
    on 2026-08-27T06:54, and the draw before that reached far outside [100, 1000] (seeds up to
    9885 on 2026-08-03/04), so seeds stamped by the old generator are not evidence about what the
    current one has used up. Widen it to "" to treat the whole recorded history as used.
    """
    counts: Counter = Counter()
    counted: list[str] = []
    latest = ""
    for task in tasks:
        created_at = str(task.get("created_at") or "")
        if since and created_at < since:
            continue
        contract = (task.get("content") or {}).get("contract") or {}
        seeds = parse_seeds(contract)
        if not seeds or not task.get("id"):
            continue
        counts.update(seeds)
        counted.append(str(task["id"]))
        latest = max(latest, created_at)
    return counts, counted, latest


def record_drawn_seeds(tasks: Iterable[dict], path: str | None = None) -> dict:
    """Merge the seeds stamped across ``tasks`` into the ledger, and return what it now holds.

    Merge rather than replace: the history endpoint serves a rolling window, so a round that has
    aged out of it must not age out of the ledger — a seed the backend used in August is still a
    seed it has used. That makes the ledger the one piece of miner state that is cumulative, and
    the reason it is a file rather than a fetch: an hour of history is not the question it answers.

    Merging by *count* would double-count every round still in the window on the next refresh, so
    the union is taken per task id — ``task_ids`` records which rounds are already folded in. A
    round that was still unstamped when it was last seen never entered that set, so it is folded
    in by whichever refresh first sees its seed.

    Best-effort throughout, like ``record_generated_seeds``: a ledger that cannot be written costs
    the draw its avoidance, and nothing else.
    """
    path = path or settings.MINER_DRAWN_SEEDS_PATH
    try:
        document = {}
        if Path(path).exists():
            with open(path) as ledger_file:
                loaded = json.load(ledger_file)
            if isinstance(loaded, dict):
                document = loaded

        since = document.get("since") or None
        known_ids = set(document.get("task_ids") or [])
        fresh = [task for task in tasks if task.get("id") not in known_ids]
        counts, counted, latest = stamped_seeds(fresh, since)

        merged = Counter({int(seed): int(times)
                          for seed, times in (document.get("counts") or {}).items()})
        merged.update(counts)
        document["counts"] = {str(seed): merged[seed] for seed in sorted(merged)}
        document["task_ids"] = sorted(known_ids | set(counted))
        document["tasks_recorded"] = len(document["task_ids"])
        document["latest_task_at"] = max(str(document.get("latest_task_at") or ""), latest)
        document["undrawn"] = sum(1 for seed in SEED_SUPPORT if seed not in merged)
        document["updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

        temporary_path = f"{path}.tmp"
        Path(temporary_path).parent.mkdir(parents=True, exist_ok=True)
        with open(temporary_path, "w") as output_file:
            json.dump(document, output_file, indent=2)
            output_file.write("\n")
        Path(temporary_path).replace(path)
        if counted:
            logger.info(
                "Drawn-seed ledger: +%d round(s), %d seed(s) of %d-%d still unstamped",
                len(counted), document["undrawn"], SEED_SUPPORT[0], SEED_SUPPORT[-1],
            )
        return document
    except Exception as error:
        logger.warning(f"Could not update the drawn-seed ledger {path} ({error})")
        return {}


def draw_seeds(count: int = 1, exclude: Iterable[int] = (), counts_path: str | None = None,
               counts: dict[int, int] | None = None) -> list[int]:
    """``count`` seeds to bet on, drawn **from the seeds the backend has never stamped**.

    Right with probability ``count/901`` per seed. What a wrong guess costs is measured in
    ``Miner._resolve_seeds``: less than it sounds, because a submission pinned to the wrong seed
    still has its ``total_weighted_score`` and coverage intact.

    The draw is still random, but its pool is not the whole support: candidates are grouped by how
    often ``drawn_seeds.json`` has seen the backend stamp them, and the never-stamped group is
    sampled first. Measured over the 231 three-seed rounds to 2026-09-21 (693 draws over the 901
    values of [100, 1000]), that group holds 413 seeds; the rest split 325 stamped once, 126 twice,
    32 three times, 5 four times.

    **What this is worth.** Under the draw the backend actually appears to make, nothing: that
    occurrence histogram fits Poisson(693/901) with chi2 0.34 over k=0..4, the three seeds within a
    round are distinct but rounds are independent of each other (1 consecutive-round overlap
    against 2.3 expected, 0 repeated triples), so every candidate is equally likely and restricting
    the pool neither gains nor loses. It is a free option on the alternative: *if* the backend ever
    draws without replacement, or de-prioritises what it has already used, the never-stamped group
    is where the seed has to come from and an unrestricted draw would be spending ~54% of its
    tickets on values that cannot win. Costless under the null, positive under the alternative —
    which is the whole argument for it, and it is not evidence of a bias that has been measured.

    The tiers matter because the never-stamped group shrinks: ~1.6 new seeds are stamped per round,
    so it empties in roughly a month of mining. Falling through to the least-stamped tier keeps the
    policy meaningful after that instead of silently reverting to a uniform draw, and keeps this
    function total — it always returns ``count`` seeds while the support can supply them.
    """
    if counts is None:
        counts = read_drawn_counts(counts_path)
    excluded = set(exclude)

    by_times_drawn: dict[int, list[int]] = defaultdict(list)
    for seed in SEED_SUPPORT:
        if seed not in excluded:
            by_times_drawn[counts.get(seed, 0)].append(seed)

    chosen: list[int] = []
    warned = False
    for times_drawn in sorted(by_times_drawn):
        if len(chosen) >= count:
            break
        tier = by_times_drawn[times_drawn]
        if times_drawn and not warned:
            # Once per draw, whether the unstamped pool ran out mid-round or was empty to begin
            # with. Worth a warning either way: it is the point at which the ledger stops being
            # able to tell the draw anything, not a failure — every candidate is equally likely,
            # so the seeds this returns are worth exactly what the unstamped ones were.
            logger.warning(
                "Only %d seed(s) in %s-%s are still unstamped, short of the %d this round bets "
                "on, so the draw falls through to the ones stamped %d time(s).",
                len(by_times_drawn.get(0, ())), SEED_SUPPORT[0], SEED_SUPPORT[-1],
                count, times_drawn,
            )
            warned = True
        chosen += random.sample(tier, min(count - len(chosen), len(tier)))
    return chosen


def resolve_seeds(contract: dict, seeds_path: str | None = None) -> tuple[list[int], str]:
    """The seeds to build against, and where they came from.

    The provenance is ``"contract"`` when the backend stamped them, otherwise whatever
    ``plan_seeds`` decided — ``"file+drawn"``, ``"drawn"``, or ``"seed-blind"`` when the count is 0.

    Contract first: if the backend ever broadcasts a stamped contract that is ground truth, and no
    local file should be able to override it (the file governs what to *guess*, not what is known).
    Today it never does, so in practice the plan decides.

    Note this draws random seeds, so two builds of the same unstamped contract are not identical.
    That is deliberate — it is what a live round does — and ``last_generated`` in the seed file
    records which seeds any given build actually used. Pin a benchmark by listing the seeds
    explicitly, or turn the whole thing off with ``"seed_count": 0``.
    """
    contract_seeds = parse_seeds(contract)
    if contract_seeds:
        return contract_seeds, "contract"
    # No seeds on the contract, so this round plans its own: the cell type's count, the operator's
    # listed seeds first, random draws for the rest. ``plan_seeds`` also writes ``last_generated``,
    # which is why planning must happen here rather than only on the neuron — a harness driving
    # ``build`` directly would otherwise both bet differently and leave no record of what it bet on.
    return plan_seeds(contract, seeds_path)


@dataclass(frozen=True)
class Coordinate:
    """A position where the reference genome actually carries a PAM for this (cas, strand)."""

    cas_system: str
    strand: str
    start: int
    length: int
    reference_guide: str

    @property
    def identity(self) -> tuple:
        return (self.cas_system, self.strand, self.start, self.length)

    def seed_window(self) -> slice:
        """The 12-mer ``offtarget_uniqueness`` hashes: the PAM-proximal end of the guide."""
        if self.cas_system == "Cas9":
            return slice(self.length - KMER_LENGTH, self.length)
        return slice(0, KMER_LENGTH)


# Process-global caches. The reference is 135 MB and every task issued on this subnet so far shares
# one gene_region and one rules block, so a warm miner pays for nothing but the build itself.
_SEQUENCE_CACHE: str | None = None
_KMER_INDEX_CACHE: dict[tuple[int, int], dict] = {}
_PAM_COORDINATE_CACHE: dict[tuple, list[Coordinate]] = {}


# ---------------------------------------------------------------------------------------------
# Context
# ---------------------------------------------------------------------------------------------

def load_sequence() -> str:
    """The chr11 FASTA both sides check coordinates against.

    Read-only and shared with the validator: it is ~130 MB and identical for both roles, so it is
    the one input the miner does not keep its own copy of. ``NIOME_GENOME_PATH`` moves it.

    Raises a plain error rather than ``SystemExit`` so a long-lived miner logs a lost round instead
    of dying inside a background task.
    """
    global _SEQUENCE_CACHE
    if _SEQUENCE_CACHE is None:
        if not Path(settings.MINER_GENOME_PATH).exists():
            raise RuntimeError(
                f"{settings.MINER_GENOME_PATH} is missing — the design cannot verify a PAM or a "
                "coordinate without it. Download GRCh38 chromosome 11 (Ensembl release 116); "
                "scripts/run_validator.sh has the URL. Set NIOME_GENOME_PATH to point elsewhere."
            )
        _SEQUENCE_CACHE = stage12.load_chr11(settings.MINER_GENOME_PATH)
    return _SEQUENCE_CACHE


def build_context(contract: dict, reference: dict, cell_types: dict) -> Context:
    """Load the genome and the same off-target k-mer index stage 2 will score against.

    The index is built in memory rather than through ``stage12.load_or_build_kmer_index``, whose
    pickle cache lives under ``data/``. Building it costs 0.16 s against 0.05 s for a cache hit,
    paid once per process because the result is held in ``_KMER_INDEX_CACHE`` — cheaper than
    writing into the validator's workspace.
    """
    sequence = load_sequence()
    window = (
        max(0, reference["gene_region"]["start"] - OFFTARGET_FLANK),
        min(len(sequence), reference["gene_region"]["end"] + OFFTARGET_FLANK),
    )
    if window not in _KMER_INDEX_CACHE:
        _KMER_INDEX_CACHE[window] = stage12.build_kmer_index(
            sequence[window[0]:window[1]], k=KMER_LENGTH
        )
    return Context(
        sequence=sequence,
        contract=contract,
        reference=reference,
        cell_types=cell_types or {},
        kmer_index=_KMER_INDEX_CACHE[window],
    )


# ---------------------------------------------------------------------------------------------
# PAM enumeration
#
# stage12.check_pam reads a fixed motif at a fixed offset, so the positions where one can exist are
# exactly the occurrences of a 2-3 base string. Finding those with str.find and then confirming each
# hit through check_pam itself is far cheaper than calling check_pam at every offset in the window,
# and cannot disagree with it — the gate has the final say on every coordinate returned.
# ---------------------------------------------------------------------------------------------

def _pam_anchors(cas_system: str, strand: str, length: int) -> tuple[str, int]:
    """The literal to search for, and the offset from a hit to ``target_alignment_start``.

    Derived from ``check_pam``: for the minus strand the motif is read off the reverse complement,
    so ``NGG`` becomes ``CC`` on the forward strand and Cas12a's ``TTTV`` becomes ``AAA``.
    """
    if cas_system == "Cas9":
        # seq[s+L:s+L+3] with pam[1:] == "GG"  ->  "GG" sits at s+L+1
        # revcomp(seq[s-3:s])[1:] == "GG"      ->  "CC" sits at s-3
        return ("GG", -(length + 1)) if strand == "+" else ("CC", 3)
    # seq[s-4:s][:3] == "TTT"                  ->  "TTT" sits at s-4
    # revcomp(seq[s+L:s+L+4])[:3] == "TTT"     ->  "AAA" sits at s+L+1
    return ("TTT", 4) if strand == "+" else ("AAA", -(length + 1))


def enumerate_coordinates(
    context: Context, config: Config,
) -> dict[tuple[str, str], list[Coordinate]]:
    """Every PAM-bearing coordinate in gene_region +/- flank, grouped by (cas, strand)."""
    gene_region = context.reference["gene_region"]
    window_start = max(KMER_LENGTH, gene_region["start"] - config.pam_search_flank)
    window_end = min(
        len(context.sequence) - 64, gene_region["end"] + config.pam_search_flank
    )
    guide_lengths = tuple(sorted({int(length) for length in config.guide_lengths}))
    cache_key = (window_start, window_end, guide_lengths, tuple(context.cas_systems))
    if cache_key in _PAM_COORDINATE_CACHE:
        return _group_by_cas_strand(_PAM_COORDINATE_CACHE[cache_key])

    discovered: list[Coordinate] = []
    for cas_system in context.cas_systems:
        for strand in context.strands:
            for guide_length in guide_lengths:
                motif, start_offset = _pam_anchors(cas_system, strand, guide_length)
                search_from = max(0, window_start - abs(start_offset) - guide_length)
                search_until = min(
                    len(context.sequence), window_end + abs(start_offset) + guide_length
                )
                while True:
                    motif_position = context.sequence.find(motif, search_from, search_until)
                    if motif_position < 0:
                        break
                    search_from = motif_position + 1
                    start = motif_position + start_offset
                    if not window_start <= start <= window_end:
                        continue
                    protospacer = context.sequence[start:start + guide_length]
                    if len(protospacer) != guide_length or any(
                        base not in "ACGT" for base in protospacer
                    ):
                        continue
                    # The gate itself decides, so a change to check_pam cannot leave a stale
                    # coordinate list behind.
                    pam_present, _ = stage12.check_pam(
                        context.sequence, start, guide_length, cas_system, strand
                    )
                    if not pam_present:
                        continue
                    reference_guide = (
                        protospacer if strand == "+"
                        else stage12.reverse_complement(protospacer)
                    )
                    discovered.append(
                        Coordinate(cas_system, strand, start, guide_length, reference_guide)
                    )

    _PAM_COORDINATE_CACHE[cache_key] = discovered
    return _group_by_cas_strand(discovered)


def _group_by_cas_strand(
    all_coordinates: list[Coordinate],
) -> dict[tuple[str, str], list[Coordinate]]:
    by_cas_strand: dict[tuple[str, str], list[Coordinate]] = defaultdict(list)
    for coordinate in all_coordinates:
        by_cas_strand[(coordinate.cas_system, coordinate.strand)].append(coordinate)
    return by_cas_strand


# ---------------------------------------------------------------------------------------------
# Guide construction
#
# Stage 1 compares the guide against the reference target and allows up to
# ``rules.max_mismatches`` differences, while the PAM is always read off the reference. That budget
# is the design's only free lever, and it has three uses at once:
#
#   * move the GC count to exactly 50%, where gc_score peaks at 1.0;
#   * break the guide's 12-mer off-target seed out of the validator's index, taking
#     offtarget_factor from 0.7 to 1.0;
#   * generate many distinct guides at one coordinate, all sharing a single stage-2 feature vector,
#     which is what lets a whole cell collapse to one row of stage 4's feature matrix.
# ---------------------------------------------------------------------------------------------

def _gc_count(guide: str) -> int:
    return sum(base in "GC" for base in guide)


def reachable_gc_score(coordinate: Coordinate, mismatch_budget: int) -> float:
    """``gc_score`` this coordinate reaches once the budget is spent pulling GC toward 50%."""
    target_gc_count = round(coordinate.length / 2)
    reference_gc_count = _gc_count(coordinate.reference_guide)
    reached_gc_count = reference_gc_count + max(
        -mismatch_budget, min(mismatch_budget, target_gc_count - reference_gc_count)
    )
    return max(0.0, 1.0 - abs(reached_gc_count / coordinate.length - 0.5) * 2)


def reachable_structural(coordinate: Coordinate, distance: int, context: Context) -> float:
    """Upper bound on this coordinate's stage-2 structural score after guide tuning.

    Mirrors stage 2 with ``offtarget_factor`` taken as 1.0, which the enumeration below reaches for
    every coordinate that has any spare budget at all.
    """
    return (
        0.625 * reachable_gc_score(coordinate, context.max_mismatches)
        + 0.375 * math.exp(-distance / context.base_padding)
    )


def enumerate_guides(coordinate: Coordinate, context: Context, max_guides: int) -> list[str]:
    """Distinct guides for one coordinate, all at the same GC count and all off-target clean.

    Substitutions come in two kinds. A *class flip* moves one base between {G,C} and {A,T} and is
    spent only to correct the GC count; a *within-class* swap (G<->C, A<->T) leaves the count alone
    and exists purely to make another distinct guide. Because every returned guide has an identical
    GC count and sits at the same coordinate, stages 1, 2 and 5 see one feature vector for the whole
    set — they differ only in the stage-3 draw, which is exactly the freedom the cut-survival
    selection needs.

    Guides are ordered by how many of their 12-mer windows differ from the reference, because those
    windows are what stage 5's ``kmer_diversity_entropy_ratio`` pools: variants that edit different
    windows keep that ratio up, and it is one of the six factors in the geometric mean.
    """
    mismatch_budget = context.max_mismatches
    if mismatch_budget <= 0:
        return [coordinate.reference_guide]

    reference_bases = list(coordinate.reference_guide)
    reference_gc_count = _gc_count(coordinate.reference_guide)
    gc_correction = max(
        -mismatch_budget,
        min(mismatch_budget, round(coordinate.length / 2) - reference_gc_count),
    )
    flips_needed = abs(gc_correction)
    flip_toward_gc = gc_correction > 0

    all_positions = range(coordinate.length)
    flippable_positions = [
        position for position in all_positions
        if (reference_bases[position] in "GC") != flip_toward_gc
    ]
    if len(flippable_positions) < flips_needed:
        return [coordinate.reference_guide]

    seed_slice = coordinate.seed_window()
    reference_windows = _kmer_windows(coordinate.reference_guide)
    seen_guides: set[str] = set()
    scored_candidates: list[tuple[int, int, str]] = []

    # A few times the requested pool is enumerated, then ranked and cut down: for the common case
    # of a coordinate already at 50% GC the whole reachable set is smaller than this and nothing is
    # truncated, and where it is larger the ranking is what decides which survive.
    for candidate_guide in itertools.islice(
        _substituted_guides(
            reference_bases, flippable_positions, all_positions,
            flips_needed, flip_toward_gc, mismatch_budget,
        ),
        max_guides * 3,
    ):
        if candidate_guide in seen_guides or candidate_guide == coordinate.reference_guide:
            continue
        # Both conditions the gate and stage 2 will apply, checked with their own code.
        if stage12.hamming(candidate_guide, coordinate.reference_guide) > mismatch_budget:
            continue
        if candidate_guide[seed_slice] in context.kmer_index:
            continue  # offtarget_factor would be 0.7 instead of 1.0
        seen_guides.add(candidate_guide)
        novel_window_count = sum(
            1 for candidate_window, reference_window
            in zip(_kmer_windows(candidate_guide), reference_windows)
            if candidate_window != reference_window
        )
        scored_candidates.append(
            (-novel_window_count, len(scored_candidates), candidate_guide)
        )

    if not scored_candidates:
        return [coordinate.reference_guide]
    scored_candidates.sort()
    return [guide for _, _, guide in scored_candidates[:max_guides]]


def _substituted_guides(
    reference_bases: list[str], flippable_positions: list[int], all_positions: range,
    flips_needed: int, flip_toward_gc: bool, mismatch_budget: int,
):
    """Every guide reachable from the reference within the budget, at a fixed GC count.

    ``flips_needed`` positions change GC class to correct the count; the rest of the budget is spent
    on within-class swaps, which leave the count alone and exist only to make another distinct
    guide. The whole budget is spent whenever there is any left over, because more edited positions
    means more distinct 12-mer windows and stage 2 charges nothing for a mismatch the gate allows.
    """
    for spare_swaps in range(mismatch_budget - flips_needed, -1, -1):
        for flip_positions in itertools.combinations(flippable_positions, flips_needed):
            swappable_positions = [
                position for position in all_positions if position not in flip_positions
            ]
            for swap_positions in itertools.combinations(swappable_positions, spare_swaps):
                for flip_bases in itertools.product(
                    *(_GC_BASES if flip_toward_gc else _AT_BASES for _ in flip_positions)
                ):
                    candidate_bases = list(reference_bases)
                    for position, replacement in zip(flip_positions, flip_bases):
                        candidate_bases[position] = replacement
                    for position in swap_positions:
                        candidate_bases[position] = _WITHIN_CLASS_SWAP[
                            reference_bases[position]
                        ]
                    yield "".join(candidate_bases)


def _kmer_windows(guide: str) -> list[str]:
    """The guide's overlapping 12-mers, the units stage 5 pools for its diversity ratio."""
    return [
        guide[offset:offset + KMER_LENGTH]
        for offset in range(len(guide) - KMER_LENGTH + 1)
    ]


# ---------------------------------------------------------------------------------------------
# Row assembly — through the validator's own gate and feature extraction
# ---------------------------------------------------------------------------------------------

def make_row(
    coordinate: Coordinate, guide: str, mutation: str, context: Context, row_index: int
) -> dict:
    """One submission row. These eight fields are the entire miner-supplied surface."""
    return {
        "experiment_id": f"exp-{row_index:05d}",
        "guideRNA": guide,
        "target_alignment_start": coordinate.start,
        "target_alignment_end": coordinate.start + coordinate.length,
        "strand": coordinate.strand,
        "mutation": mutation,
        "cas_system": coordinate.cas_system,
        "cell_type": context.contract.get("cell_type"),
    }


def gate_and_score(row: dict, context: Context) -> dict | None:
    """Run stage 1 and stage 2 on a row and return the entry stage 3 would receive.

    ``None`` means the row would be rejected, which for a generated row is a bug in the enumeration
    rather than an expected outcome — so callers count it rather than ignoring it.
    """
    stage1_pass, _reason = stage12.stage1(
        row, context.sequence, context.mutation_map, context.contract
    )
    if stage1_pass != 1.0:
        return None
    structural_score, stage2_info = stage12.stage2(
        context.cell_types, row, context.sequence,
        context.mutation_map, context.contract, context.kmer_index,
    )
    return {
        "experiment": row,
        "features": {
            "gc": stage2_info["gc"],
            "distance_to_mutation": stage2_info["distance"],
            "gc_score": stage2_info["gc_score"],
            "dist_score": stage2_info["dist_score"],
            "consistency": stage2_info["consistency"],
            "offtarget_factor": stage2_info["offtarget_factor"],
            "mutation_weight": stage2_info["mutation_weight"],
            "cell_type": stage2_info["cell_type"],
            "cell_type_accessibility": stage2_info["cell_type_accessibility"],
            "mutation_region": stage2_info["mutation_region"],
            "region_energy_offset": stage2_info["region_energy_offset"],
        },
        "stage1": {"valid": True},
        "stage2": {
            "structural_score": structural_score,
            "weighted_score": stage2_info["weighted_score"],
        },
    }


# ---------------------------------------------------------------------------------------------
# k-mer diversity
#
# ``kmer_diversity_entropy_ratio`` pools every 12-mer of every guide and divides the pool's Shannon
# entropy by log2 of the number of k-mer *instances*, so the ratio is 1.0 only if no 12-mer is ever
# repeated. Collapsing a cell onto one coordinate is what buys the flat feature matrix, but it also
# means a cell's guides differ in at most ``max_mismatches`` positions and therefore share most of
# their windows — which showed up as this ratio falling from 0.97 to 0.84, a 2.4% haircut on the
# whole score through the sixth root.
#
# Entropy is maximised when the multiplicities are level, so guides are picked greedily against a
# running census of the windows already used, round-robin across cells so no cell is left choosing
# from what the others rejected.
# ---------------------------------------------------------------------------------------------

def cut_support_mask(entry: dict, seeds: tuple[int, ...], stream: random.Random) -> int:
    """Bitmask over ``seeds`` of the rounds in which this row cuts.

    ``stage3.simulate`` draws the microhomology trigger *before* the cut, so the cut is the second
    draw from the row's stream. The order is reproduced rather than approximated — reading the
    first draw would score a different guide — and the seed and the probability come from stage 3
    itself, so this cannot drift from the simulation it is predicting.

    ``stream`` is reseeded rather than reallocated because the scan does this ~900 times per
    candidate and the Mersenne seeding is the whole cost.
    """
    experiment = entry["experiment"]
    energy = stage3.sequence_energy(stage3.extract_features(entry))
    cut_probability = stage3.cut_probability(experiment["cas_system"], energy)
    mask = 0
    for bit, seed in enumerate(seeds):
        stream.seed(stage3.experiment_seed(seed, experiment))
        stream.random()  # microhomology_trigger consumes the first draw
        if stream.random() <= cut_probability:
            mask |= 1 << bit
    return mask


def _round_robin_by_coordinate(pool: list[dict]) -> list[dict]:
    """``pool`` reordered so a truncation to any width keeps a spread of coordinates.

    ``_widen_strong_cas_pool`` grows coordinate COUNT specifically so the strong-cas cells span a
    range of ``cut_probability`` — different coordinates carry different ``distance_to_mutation``
    and therefore different energy. Ranking by ``weighted_score`` alone (what the ordinary,
    single-coordinate scan below does) throws that away: every substitution variant at the nearest
    coordinate has a higher ``dist_score`` than every variant at a farther one, so a plain top-N
    truncation collapses a 64-coordinate pool back onto whichever one or two coordinates are
    closest — exactly the single-cut-probability degeneracy the widening exists to escape. Measured
    on CD34+_HSPC: truncating the widened pool by weighted_score alone left ``two_stage_construction``
    holding *fewer* seeds (217) than the plain retune it was meant to beat (397). Round-robining
    across coordinates, strongest guide first within each, keeps the spread through the truncation.
    """
    by_coordinate: dict[int, list[dict]] = {}
    for entry in pool:
        by_coordinate.setdefault(entry["experiment"]["target_alignment_start"], []).append(entry)
    for entries in by_coordinate.values():
        entries.sort(key=lambda entry: -entry["stage2"]["weighted_score"])
    # Coordinates ordered by their own best candidate, so if the eventual width can't cover every
    # coordinate once, the strongest coordinates are still the ones represented.
    ordered_coordinates = sorted(
        by_coordinate.values(), key=lambda entries: -entries[0]["stage2"]["weighted_score"]
    )
    result: list[dict] = []
    round_index = 0
    while len(result) < len(pool):
        progressed = False
        for entries in ordered_coordinates:
            if round_index < len(entries):
                result.append(entries[round_index])
                progressed = True
        if not progressed:
            break
        round_index += 1
    return result


def scan_cut_support(
    pools: dict[tuple, list[dict]], rows_per_cell: dict[tuple, int],
    seeds: tuple[int, ...], candidates_per_cell: int, budget_seconds: float,
    rank: bool | str = "weighted_score",
) -> tuple[dict[tuple, list[dict]], dict[tuple, list[int]]]:
    """Cut masks for the strongest candidates in each cell, priced against a wall-clock budget.

    Returns the pools truncated to what was actually scanned, alongside the masks, so the two stay
    index-aligned for ``select_for_diversity``. A cell is never truncated below its own allocation:
    an unscanned candidate is still a usable row, but a cell that cannot fill its rows would cost a
    coverage-entropy cliff worth far more than the intersection.

    ``rank`` selects how a cell's pool is ordered before truncation. ``"weighted_score"`` (the
    default, and what every ordinary scan uses) ranks purely by row value, which is correct where a
    cell sits at one coordinate. ``"coordinate_diverse"`` round-robins across coordinates first —
    see ``_round_robin_by_coordinate`` — for the one caller (the two-stage construction's widened
    strong-cas rescan) where the pool spans many coordinates on purpose and a value-only ranking
    would erase that.
    """
    if not seeds or candidates_per_cell <= 0 or not pools:
        return pools, {}

    # Rank by what a row is worth before truncating, so the candidates dropped for scan width are
    # the ones the k-mer pass would have reached for last anyway -- unless the caller has already
    # spent a wall-clock budget building coordinate spread into the pool, in which case truncating
    # by value alone would spend that budget for nothing.
    if rank == "coordinate_diverse":
        ranked = {cell: _round_robin_by_coordinate(pool) for cell, pool in pools.items()}
    else:
        ranked = {
            cell: sorted(pool, key=lambda entry: -entry["stage2"]["weighted_score"])
            for cell, pool in pools.items()
        }
    stream = random.Random()
    # ``build`` only ever passes cells that produced candidates, but the width below is priced off
    # a real candidate, so the probe has to find one rather than assume the first cell has one.
    first_cell = next((cell for cell, pool in ranked.items() if pool), None)
    if first_cell is None:
        return pools, {}
    started = time.monotonic()
    first_mask = cut_support_mask(ranked[first_cell][0], seeds, stream)
    seconds_per_candidate = max(time.monotonic() - started, 1e-9)

    affordable = int(budget_seconds / (seconds_per_candidate * len(ranked)))
    width = max(1, min(candidates_per_cell, affordable))

    scanned: dict[tuple, list[dict]] = {}
    masks: dict[tuple, list[int]] = {}
    for cell, pool in ranked.items():
        cell_width = min(len(pool), max(width, rows_per_cell.get(cell, 0)))
        scanned[cell] = pool[:cell_width]
        cell_masks = [first_mask] if cell == first_cell else []
        for entry in scanned[cell][len(cell_masks):]:
            cell_masks.append(cut_support_mask(entry, seeds, stream))
        masks[cell] = cell_masks
    return scanned, masks


def select_for_diversity(
    pools: dict[tuple, list[dict]], rows_per_cell: dict[tuple, int],
    cut_masks: dict[tuple, list[int]] | None = None, seed_count: int = 0,
) -> tuple[list[dict], int]:
    """Fill each cell's allocation, holding the round seeds on which *every* row cuts.

    Two objectives, in strict priority order.

    The first is ``consistency_factor``'s cliff. Stage 4 scores ``is_cut`` with ``r2_score``, and a
    fold whose ``y_test`` never varies has zero total sum of squares — sklearn answers that 0/0 with
    1.0 when the prediction is exact and 0.0 when it is not. So a submission in which every row cuts
    scores r2 1.0 on that target and one in which a single row fails scores ~0: 0.2394 against
    0.1008, measured at 250 rows. A Cas9 row cuts with probability 0.99 even at fully open
    chromatin and a Cas12a row with 0.96, so a design that ignores this holds the cliff on almost
    no rounds at all — on the reference task the k-mer pass alone held it on 1 seed out of 900.

    The stamped seed is not visible at build time, but it is drawn from ``SEED_SUPPORT``, so every
    seed it could be is enumerable. Guides chosen to cut under the *same* subset of that support
    hold the cliff across the whole subset, and the stamped seed lands inside it with probability
    ``|subset| / |SEED_SUPPORT|``. Note this is an intersection, not a ranking: a guide's outcomes
    under two seeds are independent draws (``experiment_seed`` hashes the seed in with the row), so
    ranking guides by their overall cut rate is noise and nothing here generalises to an unscanned
    seed.

    The greedy is adaptive rather than aimed at a subset fixed up front, which is what makes it
    work: taking the candidate that keeps the most of the *currently* surviving set is a maximum
    over the whole pool at every step, so the set decays far more slowly than any pre-committed
    target could be satisfied. Measured 215 of 900 seeds held against a pre-committed target's ~89.

    The k-mer objective is the tiebreak and costs almost nothing, because once the surviving set is
    small enough many candidates preserve it exactly and the choice among those is free. It is exact
    for the marginal cost it minimises: adding a guide raises the pool's entropy least where its
    windows are already common, so the lowest summed window census is the locally optimal move.

    Returns the chosen entries and the number of seeds on which all of them cut (0 when no scan ran).
    """
    # Candidates are held as (windows, entry) pairs indexed per cell, and consumed by index. A
    # dict is not hashable and compares by value, so neither a set membership test nor list.remove
    # would be safe or cheap here.
    candidates_by_cell: dict[tuple, list[tuple[list[str], dict]]] = {
        cell: [(_kmer_windows(entry["experiment"]["guideRNA"]), entry) for entry in pool]
        for cell, pool in pools.items()
    }
    used_indices: dict[tuple, set[int]] = {cell: set() for cell in pools}
    window_census: Counter = Counter()
    rows_left = {cell: rows_per_cell.get(cell, 0) for cell in pools}
    chosen_entries: list[dict] = []
    # Every scanned seed starts alive and each pick can only clear bits, so this is the set of
    # rounds on which every row chosen so far cuts.
    surviving = (1 << seed_count) - 1 if cut_masks and seed_count else 0

    while any(count > 0 for count in rows_left.values()):
        progressed = False
        for cell, cell_candidates in candidates_by_cell.items():
            if rows_left[cell] <= 0 or len(used_indices[cell]) >= len(cell_candidates):
                continue
            cell_masks = cut_masks.get(cell) if cut_masks else None

            def rank(candidate_index: int, cell_masks=cell_masks) -> tuple:
                windows, entry = cell_candidates[candidate_index]
                census = sum(window_census[window] for window in windows)
                weighted_score = entry["stage2"]["weighted_score"]
                if cell_masks is None or candidate_index >= len(cell_masks):
                    return (0, census, -weighted_score)
                # Negated so that keeping the most seeds alive sorts first under ``min``.
                return (-(surviving & cell_masks[candidate_index]).bit_count(),
                        census, -weighted_score)

            best_index = min(
                (
                    candidate_index for candidate_index in range(len(cell_candidates))
                    if candidate_index not in used_indices[cell]
                ),
                key=rank,
            )
            used_indices[cell].add(best_index)
            for window in cell_candidates[best_index][0]:
                window_census[window] += 1
            if cell_masks is not None and best_index < len(cell_masks):
                surviving &= cell_masks[best_index]
            chosen_entries.append(cell_candidates[best_index][1])
            rows_left[cell] -= 1
            progressed = True
        if not progressed:
            break
    return chosen_entries, surviving.bit_count()


# ---------------------------------------------------------------------------------------------
# Min-union selection — the weak-cas half of the two-stage construction
#
# ``select_for_diversity`` greedily keeps the most of whatever seed set currently survives, which
# is the right objective when every cell has to hold the same intersection. The two-stage
# construction (see ``two_stage_construction``) asks a narrower question of one cas system only:
# out of a bank of candidates, which small group's *combined* failures cover the fewest seeds? That
# is min-set-cover's greedy relative, min-union, and it needs its own selector because the
# objective is the union's size, not the surviving intersection of a fixed pool.
#
# Ported from develop's ``fastgreedy.FastGreedy.build`` (see the module the task description
# points at), adapted to this repo's bitmask representation and made CPU-only — this repo has no
# GPU and must not depend on one, so the ``cupy`` fallback there does not apply here.
# ---------------------------------------------------------------------------------------------

def _fail_indices(mask: int, seed_count: int) -> np.ndarray:
    """The seed indices a cut-support mask does NOT set — fastgreedy's 'fails' list.

    Reads the same bitmask ``cut_support_mask`` already builds rather than rescanning: unpacked
    through bytes so the one-time int-to-array conversion (there are as many of these as
    candidates, not as many as greedy steps) is bulk numpy rather than a 900-iteration Python loop
    per candidate.
    """
    n_bytes = (seed_count + 7) // 8
    packed = np.frombuffer(mask.to_bytes(n_bytes, "little"), dtype=np.uint8)
    cuts = np.unpackbits(packed, bitorder="little", count=seed_count).astype(bool)
    return np.flatnonzero(~cuts)


def min_union_select(
    masks: list[int], cell_of: list, seed_count: int, group_size: int,
    cell_floor: dict, cell_cap: dict | None = None,
) -> tuple[list[int], int]:
    """Greedy min-union group: the ``group_size`` candidates whose combined failures cover the
    fewest seeds, subject to per-cell floors (and optional caps).

    At each step every remaining candidate's marginal cost — how many *currently uncovered* seeds
    its own failures would newly cover — is computed for the whole bank at once as
    ``uncovered[idx].sum(axis=1)``, a vectorised gather-sum, rather than one popcount per candidate
    per step; ``idx`` pads each candidate's fail list out to the bank's longest one with a sentinel
    column pointed at an always-covered slot, so padding never contributes to the cost. ``argmin``
    keeps the *first* minimum on a tie, reproducing what a sequential Python scan (as
    ``seed_agnostic.min_union_group`` does) would keep.

    ``cell_floor`` guarantees at least that many picks land in a cell before any cell may exceed
    its floor — the coverage-entropy floor every other selector in this file also enforces, so a
    narrow bank cannot empty a stage-5 cell. ``cell_cap`` optionally bounds a cell once its floor is
    met, which is how the caller pins the group's cell shape to a target split (``_shifted_toward``)
    instead of letting it drift wherever the union happens to be cheapest.

    Returns the chosen indices (into ``masks``/``cell_of``) and the group's *clean set* as a
    bitmask in ``cut_support_mask``'s one-bit-per-seed format — the complement of the union of
    every chosen candidate's failures. This falls out of the algorithm for free: it is exactly the
    ``uncovered`` state left once every chosen candidate's fail bits have been cleared.
    """
    n = len(masks)
    if seed_count <= 0:
        return [], 0
    if n == 0 or group_size <= 0:
        return [], (1 << seed_count) - 1  # nothing chosen -> nothing failed -> everything clean

    fail_lists = [_fail_indices(mask, seed_count) for mask in masks]
    max_fail = max((len(fails) for fails in fail_lists), default=0)
    idx = np.full((n, max(max_fail, 1)), seed_count, dtype=np.int32)  # sentinel column = seed_count
    for row, fails in enumerate(fail_lists):
        if len(fails):
            idx[row, :len(fails)] = fails

    cells = sorted(set(cell_of))
    cell_index = {cell: position for position, cell in enumerate(cells)}
    cell_id = np.asarray([cell_index[cell] for cell in cell_of], dtype=np.int32)
    floor_by_id = {
        cell_index[cell]: min(value, sum(1 for candidate_cell in cell_of if candidate_cell == cell))
        for cell, value in cell_floor.items() if cell in cell_index
    }
    cap_by_id = None
    if cell_cap:
        cap_by_id = {
            cell_index[cell]: value for cell, value in cell_cap.items() if cell in cell_index
        }

    uncovered = np.ones(seed_count + 1, dtype=bool)
    uncovered[seed_count] = False  # the sentinel: padding costs nothing

    chosen: list[int] = []
    taken = np.zeros(n, dtype=bool)
    cell_count: Counter = Counter()
    BIG = np.int64(1 << 30)

    for _ in range(min(group_size, n)):
        slots_left = group_size - len(chosen)
        unmet = sum(max(0, floor_by_id.get(cid, 0) - cell_count[cid]) for cid in floor_by_id)
        cost = uncovered[idx].sum(axis=1).astype(np.int64)
        cost = np.where(taken, BIG, cost)
        if slots_left <= unmet:
            # No slack left: every remaining pick must land in a cell still below its floor.
            need = {cid for cid in floor_by_id if cell_count[cid] < floor_by_id[cid]}
            if not need:
                break
            cost = np.where(np.isin(cell_id, list(need)), cost, BIG)
        if cap_by_id is not None:
            full = [cid for cid, cap in cap_by_id.items() if cell_count[cid] >= cap]
            if full:
                cost = np.where(np.isin(cell_id, full), BIG, cost)
        if bool((cost >= BIG).all()):
            break
        pick = int(np.argmin(cost))  # first minimum, matching a sequential scan's tie behaviour
        chosen.append(pick)
        taken[pick] = True
        cell_count[int(cell_id[pick])] += 1
        uncovered[idx[pick]] = False

    packed_clean = np.packbits(uncovered[:seed_count].astype(np.uint8), bitorder="little")
    clean_mask = int.from_bytes(packed_clean.tobytes(), "little")
    return chosen, clean_mask


# ---------------------------------------------------------------------------------------------
# Row allocation across cells
#
# Stage 5's geometric mean makes an empty (mutation, cas, strand) cell a ~1e-9 multiplier, so every
# cell gets at least one row. Above that floor the split is a real trade: piling rows onto the
# heavier mutation buys total_weighted_score linearly but costs the mutation and joint coverage
# entropies, which enter the score only at the 1/6 power. Both sides of that trade are closed-form
# and seed-independent, so the optimum is found by hill-climbing the product directly instead of
# tuning a skew exponent — which also means it adapts to whatever weight ratio a contract carries.
# ---------------------------------------------------------------------------------------------

def _consistency_estimate(
    rows_per_cell: dict[tuple, int],
    cut_probability_by_cell: dict[tuple, float],
    mutation_weight_by_cell: dict[tuple, float],
) -> float:
    """Closed-form ``consistency_factor`` for a candidate row split.

    Only the ``0.3 * (1 - avg_nmae)`` term is modelled, because ``avg_r2`` is negative for every
    seed-blind design and the ``max(avg_r2, 0)`` clamp discards it. Within a cell stage 4's forest
    sees one feature vector, so the best it can do on ``is_cut`` is predict that cell's cut
    probability — which makes the weighted MAE and the target's standard deviation both computable
    from ``cut_probability`` alone, with no simulation and no forest.

    This is what makes the Cas9/Cas12a split a real decision rather than a coverage question: at
    HEK293's accessibility a Cas12a row cuts with probability 0.87 against Cas9's 0.95, so leaning
    toward Cas9 lifts the cut rate, shrinks ``is_cut``'s normalised error and pays for some of the
    Cas coverage entropy it costs.
    """
    weighted_rows = sum(
        rows_per_cell[cell] * mutation_weight_by_cell[cell] for cell in rows_per_cell
    )
    total_rows = sum(rows_per_cell.values())
    if weighted_rows <= 0 or total_rows <= 0:
        return 0.0
    # stage 4 weights the MAE by mutation_weight but takes the target's standard deviation
    # unweighted, so the two halves of the ratio are averaged differently here as well.
    weighted_mae_cut = sum(
        rows_per_cell[cell] * mutation_weight_by_cell[cell]
        * 2 * cut_probability_by_cell[cell] * (1 - cut_probability_by_cell[cell])
        for cell in rows_per_cell
    ) / weighted_rows
    mean_cut = sum(
        rows_per_cell[cell] * cut_probability_by_cell[cell] for cell in rows_per_cell
    ) / total_rows
    cut_spread = math.sqrt(max(0.0, mean_cut * (1 - mean_cut)))
    normalised_mae_cut = weighted_mae_cut / cut_spread if cut_spread > 1e-9 else 0.0
    avg_nmae = (normalised_mae_cut + NMAE_IS_HDR + NMAE_INDEL) / 3
    return max(0.0, 0.3 * (1 - avg_nmae))


def _consistency_estimate_pinned() -> float:
    """``consistency_factor`` on a round seed where *every* row cuts.

    A constant ``is_cut`` column takes that target's ``r2`` to exactly 1.0 and its normalised error
    to 0, so only the two targets no design can predict are left in either term. That makes this a
    constant rather than a function of the split — which is the point: what the row allocation now
    trades is not the *size* of the spike but how *often* it is reached.

    Measured 0.311 on K562 and HUDEP-2 and 0.394 on HEK293 against this estimate's 0.367. It is
    used only to rank allocations against ``_consistency_estimate``'s unpinned value, where what
    matters is the ratio between the two, not either one's absolute accuracy.
    """
    # is_cut alone reaches r2 1.0; the other two sit at best near zero and the clamp discards them.
    avg_nmae = (0.0 + NMAE_IS_HDR + NMAE_INDEL) / 3
    return max(0.0, min(1.0, 0.7 * (1.0 / 3.0) + 0.3 * (1 - avg_nmae)))


def _coverage_fidelity(rows_per_cell: dict[tuple, int], context: Context) -> float:
    """Stage 5's ``distribution_fidelity`` for a candidate row split, in closed form."""
    mutation_counts: Counter = Counter()
    cas_counts: Counter = Counter()
    strand_counts: Counter = Counter()
    joint_counts: Counter = Counter()
    for (mutation, cas_system, strand), row_count in rows_per_cell.items():
        mutation_counts[mutation] += row_count
        cas_counts[cas_system] += row_count
        strand_counts[strand] += row_count
        joint_counts[(mutation, cas_system, strand)] += row_count
    joint_support = [
        (mutation, cas_system, strand)
        for mutation in context.mutations
        for cas_system in context.cas_systems
        for strand in context.strands
    ]
    ratios = [
        stage5.coverage_entropy_ratio(mutation_counts, context.mutations),
        stage5.coverage_entropy_ratio(cas_counts, context.cas_systems),
        stage5.coverage_entropy_ratio(strand_counts, list(context.strands)),
        stage5.coverage_entropy_ratio(joint_counts, joint_support),
        # The two guide-level ratios are left out: they are near-identical across allocations, and
        # inside a geometric mean a common factor cannot move the argmax.
        1.0,
        1.0,
    ]
    return stage5.geometric_mean(ratios)


def _allocation_objective(
    rows_per_cell: dict[tuple, int], value_per_row: dict[tuple, float], context: Context,
    cut_probability_by_cell: dict[tuple, float],
    mutation_weight_by_cell: dict[tuple, float],
    consistency: float | None = None,
) -> float:
    """All three score factors for a candidate row split, in closed form.

    ``consistency`` replaces the estimate for a construction that already knows the factor exactly.
    Only ``pinned_outcome_build`` does, where it is 1.0 — and passing it matters rather than being
    a tidiness: ``_consistency_estimate``'s whole content is the Cas9/Cas12a cut-rate trade, and
    once every row's outcome is pinned that trade no longer exists, so leaving the estimate in
    would have the split still paying coverage entropy for a term that cannot move.
    """
    weighted_total = sum(rows_per_cell[cell] * value_per_row[cell] for cell in rows_per_cell)
    if weighted_total <= 0:
        return 0.0
    if consistency is None:
        consistency = _consistency_estimate(
            rows_per_cell, cut_probability_by_cell, mutation_weight_by_cell
        )
    return weighted_total * _coverage_fidelity(rows_per_cell, context) * consistency


def allocate_rows(
    cells: list[tuple], value_per_row: dict[tuple, float], capacity: dict[tuple, int],
    rows_wanted: int, context: Context,
    cut_probability_by_cell: dict[tuple, float],
    mutation_weight_by_cell: dict[tuple, float],
    consistency: float | None = None,
) -> dict[tuple, int]:
    """Rows per cell, by hill-climbing all three score factors from an even split.

    ``value_per_row`` is what one row in a cell is worth and ``capacity`` how many distinct guides
    the cell can actually supply, so the result is always buildable. ``consistency`` is passed
    through to ``_allocation_objective`` — see there.
    """
    rows_per_cell = {cell: 0 for cell in cells}
    capacity_left = {cell: max(0, capacity.get(cell, 0)) for cell in cells}
    assignable_rows = min(rows_wanted, sum(capacity_left.values()))

    # Start from as even a split as capacity allows: that is the coverage-entropy optimum, and the
    # climb below only ever moves away from it when term 1 pays for the loss.
    for cell in itertools.cycle(cells):
        if sum(rows_per_cell.values()) >= assignable_rows:
            break
        if all(rows_per_cell[any_cell] >= capacity_left[any_cell] for any_cell in cells):
            break
        if rows_per_cell[cell] < capacity_left[cell]:
            rows_per_cell[cell] += 1

    best_objective = _allocation_objective(
        rows_per_cell, value_per_row, context,
        cut_probability_by_cell, mutation_weight_by_cell, consistency,
    )
    for _ in range(rows_wanted * 2):
        improved = False
        for source_cell in cells:
            if rows_per_cell[source_cell] <= 1:
                continue  # never empty a cell — that is the 1e-9 cliff
            for target_cell in cells:
                if target_cell == source_cell \
                        or rows_per_cell[target_cell] >= capacity_left[target_cell]:
                    continue
                rows_per_cell[source_cell] -= 1
                rows_per_cell[target_cell] += 1
                moved_objective = _allocation_objective(
                    rows_per_cell, value_per_row, context,
                    cut_probability_by_cell, mutation_weight_by_cell, consistency,
                )
                if moved_objective > best_objective + 1e-12:
                    best_objective = moved_objective
                    improved = True
                else:
                    rows_per_cell[source_cell] += 1
                    rows_per_cell[target_cell] -= 1
        if not improved:
            break
    return rows_per_cell


# How many points of the cas-mix ladder ``retune_cas_mix`` prices. Each costs one greedy selection
# over the scanned pools, so this is a wall-clock knob; the trade it is searching is monotone in one
# variable with a single interior maximum, so a coarse ladder finds it.
_CAS_LADDER_POINTS = 7

# Relative gain the retune must predict before it overrides ``allocate_rows``.
#
# The allocator's mix is optimal for three factors it computes exactly. The retune beats it only
# through one estimated quantity — the size of the pin, which is held constant at
# ``_consistency_estimate_pinned()`` but is really 0.31-0.41 depending on row composition — so a
# thin predicted margin is not evidence of anything. What the ladder's best point predicts, against
# what it then measured end to end over 120 held-out seeds:
#
#     cell type        predicted   measured vs scan alone
#     K562 0.77            +6.5%   +30%
#     CD34+_HSPC 0.87      +4.3%   +25%
#     HEK293 0.35          +1.3%   -2.7%
#
# The estimate is conservative where the lever is real and optimistic where it is not, so the two
# regimes separate — but only across 1.3% to 4.3%, which is why this is a threshold and not a
# margin of safety. Accessibility 0.35 caps the hit rate near 6%, leaving almost nothing for the
# pin to multiply, and that is the case being excluded. A cell type that lands between those two
# bands would be decided by this constant rather than by evidence; re-measure before trusting it.
_RETUNE_MIN_GAIN = 0.03


def _shifted_toward(
    rows_per_cell: dict[tuple, int], weak_cas: str, strong_cas: str,
    weak_target: int, capacity: dict[tuple, int],
) -> dict[tuple, int]:
    """``rows_per_cell`` with the weak cas system reduced to ``weak_target`` rows.

    Rows move to the *same* (mutation, strand) cell of the strong cas system, so only the cas and
    joint coverage entropies move and the mutation and strand splits the hill-climb chose are left
    alone. No cell is ever emptied — stage 5's geometric mean turns an unoccupied cell into a ~1e-9
    multiplier, which is worth more than any intersection.
    """
    moved = dict(rows_per_cell)
    weak_cells = sorted(
        (cell for cell in moved if cell[1] == weak_cas), key=lambda cell: -moved[cell]
    )
    surplus = sum(moved[cell] for cell in weak_cells) - weak_target
    for cell in weak_cells:
        if surplus <= 0:
            break
        mutation, _cas, strand = cell
        target_cell = (mutation, strong_cas, strand)
        if target_cell not in moved:
            continue
        headroom = capacity.get(target_cell, 0) - moved[target_cell]
        take = min(surplus, moved[cell] - 1, max(0, headroom))
        if take <= 0:
            continue
        moved[cell] -= take
        moved[target_cell] += take
        surplus -= take
    return moved


def _strong_and_weak_cas(
    cells, cut_probability_by_cell: dict[tuple, float],
) -> tuple[str, str] | None:
    """(strongest, weakest) cas system by mean cut probability over the cells that use it.

    "Weak" is the system whose rows erode a surviving seed set fastest, which is the one with the
    lowest cut probability — read off stage 3 rather than assumed to be Cas12a. ``None`` when there
    is only one cas system (or every system ties), since neither ``cas_mix_ladder`` nor the
    two-stage construction have anything to trade in that case. Factored out so both read the same
    ranking off the same numbers instead of each assuming which system is which.
    """
    cas_systems = {cell[1] for cell in cells}
    if len(cas_systems) < 2:
        return None
    mean_cut = {
        cas: sum(cut_probability_by_cell.get(cell, 0.0) for cell in cells if cell[1] == cas)
        / max(1, sum(1 for cell in cells if cell[1] == cas))
        for cas in cas_systems
    }
    strong_cas = max(mean_cut, key=lambda cas: mean_cut[cas])
    weak_cas = min(mean_cut, key=lambda cas: mean_cut[cas])
    if strong_cas == weak_cas:
        return None
    return strong_cas, weak_cas


def cas_mix_ladder(
    rows_per_cell: dict[tuple, int], capacity: dict[tuple, int],
    cut_probability_by_cell: dict[tuple, float],
) -> list[dict[tuple, int]]:
    """Candidate splits from the allocation's own cas mix down to a near-single-cas one.

    Ordered weakest-cas-heaviest first, so the caller's first point is the unmodified allocation.
    """
    ranked = _strong_and_weak_cas(rows_per_cell, cut_probability_by_cell)
    if ranked is None:
        return [dict(rows_per_cell)]
    strong_cas, weak_cas = ranked

    weak_cells = [cell for cell in rows_per_cell if cell[1] == weak_cas]
    weak_total = sum(rows_per_cell[cell] for cell in weak_cells)
    floor_total = len(weak_cells)  # one row each, the coverage-entropy floor
    ladder = []
    for point in range(_CAS_LADDER_POINTS):
        share = point / max(1, _CAS_LADDER_POINTS - 1)
        target = int(round(weak_total - share * (weak_total - floor_total)))
        candidate = _shifted_toward(rows_per_cell, weak_cas, strong_cas, target, capacity)
        if candidate not in ladder:
            ladder.append(candidate)
    return ladder


def _predicted_objective(
    rows_per_cell: dict[tuple, int], seeds_held: int, seed_count: int, context: Context,
    value_per_row: dict[tuple, float], cut_probability_by_cell: dict[tuple, float],
    mutation_weight_by_cell: dict[tuple, float],
) -> float:
    """Predicted ``final_score`` for a candidate split once its seed intersection is known.

    ``seeds_held`` seeds are pinned (a constant ``is_cut`` column) and the rest are priced by
    ``_consistency_estimate``'s unpinned floor, blended by hit rate — the same closed form
    ``retune_cas_mix``'s ladder search always used, factored out here so ``two_stage_construction``
    prices its own candidates on identical terms and the two constructions can be compared
    directly rather than by two different proxies.
    """
    weighted_total = sum(
        rows_per_cell.get(cell, 0) * value_per_row.get(cell, 0.0) for cell in rows_per_cell
    )
    hit_rate = (seeds_held / seed_count) if seed_count else 0.0
    consistency = (
        hit_rate * _consistency_estimate_pinned()
        + (1 - hit_rate) * _consistency_estimate(
            rows_per_cell, cut_probability_by_cell, mutation_weight_by_cell
        )
    )
    return weighted_total * _coverage_fidelity(rows_per_cell, context) * consistency


def retune_cas_mix(
    rows_per_cell: dict[tuple, int], scanned_pools: dict[tuple, list[dict]],
    cut_masks: dict[tuple, list[int]], seed_count: int, context: Context,
    value_per_row: dict[tuple, float], cut_probability_by_cell: dict[tuple, float],
) -> tuple[dict[tuple, int], list[dict], int]:
    """Re-price the cas mix now that the seed intersection is measurable.

    ``allocate_rows`` optimises the three factors it can compute before the scan, and prices a
    cell's cut probability only through ``is_cut``'s *normalised error* — a term that moves a few
    hundredths. The intersection is a fourth factor and it is far more sensitive to the same knob,
    because it compounds: at open chromatin a Cas12a row cuts with probability 0.96 against Cas9's
    0.99, so it erodes the surviving seed set four times faster, once per row. Measured on K562,
    holding everything else fixed:

        cas12a rows    89     70     56     42     34     28     12
        seeds held    212    242    273    290    311    321    389
        fidelity    0.925  0.901  0.877  0.850  0.827  0.805  0.721

    So the mix that is right for the other three factors is too weak-cas-heavy for this one, and
    the optimum is an interior point — seeds held rises as the weak share falls, coverage entropy
    falls with it, and the product turns over. This searches a coarse ladder rather than
    hill-climbing because each point costs a real greedy selection, and because one variable with a
    single maximum does not need more.

    The retune is worth more than the row-count table above implies, because moving rows onto the
    strong cas system lifts three things at once and only the first is obvious: the hit rate
    (212 -> 374 seeds), the unpinned floor (a higher mean cut rate shrinks ``is_cut``'s normalised
    error on the seeds that miss), and the *size* of the pin itself — the same constant ``is_cut``
    is worth more when the rows around it overfit less on the two targets left free. Measured
    end to end against the scan alone, paired over 120 held-out seeds through all five stages:

        cell type       cas12a       held        fidelity      consistency        final
        CD34+ 0.87      87->18   220->366   0.906->0.741   0.151->0.230   46.10->57.62
        K562  0.77      89->18   212->374   0.920->0.751   0.154->0.244   37.72->48.88
        HUDEP-2 0.82    86->18   212->376   0.904->0.739   0.158->0.221   46.35->53.33
        HEK293 0.35     76->76    36-> 36   unchanged — see ``_RETUNE_MIN_GAIN``

    The fidelity given up is real and large; it is simply worth less than the consistency bought —
    but only where the intersection is reachable, which is what ``_RETUNE_MIN_GAIN`` decides.

    Returns the chosen split, its selected entries and the seeds it holds, so the winning
    selection is not recomputed.
    """
    capacity = {cell: len(pool) for cell, pool in scanned_pools.items()}
    best: tuple[float, dict, list, int] | None = None
    # The ladder's first point is the allocator's own mix, so this is the objective any override
    # has to beat by ``_RETUNE_MIN_GAIN``.
    base_objective: float | None = None
    for candidate in cas_mix_ladder(rows_per_cell, capacity, cut_probability_by_cell):
        entries, seeds_held = select_for_diversity(
            scanned_pools, candidate, cut_masks, seed_count
        )
        objective = _predicted_objective(
            candidate, seeds_held, seed_count, context, value_per_row,
            cut_probability_by_cell, {cell: context.weight_of(cell[0]) for cell in candidate},
        )
        if base_objective is None:
            # First point is the allocator's mix; keep it as both the incumbent and the bar.
            base_objective = objective
            best = (objective, candidate, entries, seeds_held)
            continue
        if objective > base_objective * (1 + _RETUNE_MIN_GAIN) and objective > best[0]:
            best = (objective, candidate, entries, seeds_held)
    assert best is not None  # the ladder always contains the unmodified allocation
    return best[1], best[2], best[3]


# ---------------------------------------------------------------------------------------------
# Exact-outcome construction — the path a contract that carries its round seeds takes
#
# Everything above treats the seed as unknown and plays odds against it: ``SEED_SUPPORT`` is the 900
# seeds the backend has been observed to stamp, and ``select_for_diversity`` maximises how many of
# them a submission holds a constant ``is_cut`` on, which is the *probability* of reaching the
# cliff. A contract that carries its seeds removes the uncertainty rather than narrowing it. Stage 3
# is a pure function of (seed, row) — ``experiment_seed`` hashes the round seed together with the
# mutation, cas, guide, start and strand, and nothing else — so each candidate's outcome can simply
# be read off before the row is chosen.
#
# What that buys is the whole of ``consistency_factor`` instead of a share of it. Stage 4 fits a
# forest to ``is_cut``, ``is_hdr`` and ``indel_length``; a target that never varies has zero total
# sum of squares, so ``r2_score`` answers its 0/0 with 1.0 on an exact prediction, and
# ``normalized_mae`` short-circuits its own zero divisor to the raw MAE, also 0. Hold all THREE
# constant and ``consistency_score`` is ``(0.7*1.0 + 0.3*(1 - 0)) * 100`` — exactly 100, and
# ``consistency_factor`` exactly 1.0. Measured against stage 4 itself rather than argued: a 250-row
# frame of constant targets returns avg_r2 1.0000, avg_nmae 0.000000, consistency_factor 1.0000.
# Note this is also why the KFold shuffle stops mattering — every fold's ``y_test`` is constant, so
# the row order the contract seed shuffles into cannot move the result.
#
# Three constant targets means one ``(outcome, indel_length)`` pair shared by every row, since all
# three are functions of that pair. Two are cheap: ``("HDR", 0)`` — HDR's indel is always 0, so the
# pair costs no more than the outcome — and ``("no_cut", 0)``. ``("BLUNT_NHEJ", 1)`` and friends
# qualify too but are rarer, needing the indel length to land exactly as well. The target is not
# pre-committed; every pair that can occupy every cell is priced and the best one wins.
#
# The cost is close to nothing, which is the part worth stating. Within a cell every candidate sits
# at the same coordinate at the same GC count and is already off-target clean (see
# ``enumerate_guides``), so they all carry an IDENTICAL ``weighted_score`` — verified over the built
# pools, one distinct value per cell. Filtering the pool down to a single outcome therefore cannot
# move ``total_weighted_score`` at all; it can only reduce how many rows a cell is still able to
# supply. Measured over the pools ``build`` already materialises, at three seeds: ``("HDR", 0)``
# holds for 43-90 of each cell's 900 candidates on HEK293 (accessibility 0.35) and 92-139 on K562
# (0.77), against the ~31 rows a cell is allocated out of 250. Reading them costs ~0.2 s for 7200
# candidates x 3 seeds, against the 90 s ``seed_scan_seconds`` budgets — so where this construction
# wins it *replaces* that scan rather than adding to it, and the build gets faster.
# ---------------------------------------------------------------------------------------------

def pinned_outcome(entry: dict, seeds: tuple[int, ...]) -> tuple[str, int] | None:
    """The ``(outcome, indel_length)`` this row draws under *every* seed, or None if it varies.

    ``stage3.simulate`` is called rather than reimplemented, so this cannot drift from the
    simulation it is predicting — the same discipline ``cut_support_mask`` follows. There only the
    cut draw mattered and the rest of the row's stream could be skipped; here the repair mode and
    the indel length matter too, so there is nothing left to shortcut and the stage's own function
    is both the cheapest correct answer and the only one that stays correct if stage 3 changes.

    Returns on the first seed that disagrees, which is the common case and usually the second seed.
    """
    pinned: tuple[str, int] | None = None
    for seed in seeds:
        drawn = stage3.simulate(entry, seed)
        signature = (drawn["outcome"], drawn["indel_length"])
        if pinned is None:
            pinned = signature
        elif signature != pinned:
            return None
    return pinned


def partition_by_pinned_outcome(
    pools: dict[tuple, list[dict]], seeds: tuple[int, ...]
) -> dict[tuple[str, int], dict[tuple, list[dict]]]:
    """Candidates grouped by the outcome they pin to, keeping only targets that reach every cell.

    An unoccupied (mutation, cas, strand) cell is a ~1e-9 multiplier through stage 5's geometric
    mean, so a target that cannot fill one is not a cheaper submission — it is a lost round, and
    dropping it here is what keeps ``pinned_outcome_build``'s ranking from having to price the
    cliff.
    """
    by_target: dict[tuple[str, int], dict[tuple, list[dict]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for cell, pool in pools.items():
        for entry in pool:
            signature = pinned_outcome(entry, seeds)
            if signature is not None:
                by_target[signature][cell].append(entry)
    return {
        target: dict(cells_for_target)
        for target, cells_for_target in by_target.items()
        if len(cells_for_target) == len(pools)
    }


def pinned_outcome_build(
    pools: dict[tuple, list[dict]], rows_wanted: int, seeds: tuple[int, ...], context: Context,
    value_per_row: dict[tuple, float], cut_probability_by_cell: dict[tuple, float],
    mutation_weight_by_cell: dict[tuple, float],
) -> tuple[dict[tuple, int], list[dict], tuple[str, int], float] | None:
    """The best submission whose every row draws one fixed outcome under every one of ``seeds``.

    Returns the row split, the selected entries, the ``(outcome, indel_length)`` they are all
    pinned to, and the predicted ``final_score`` — or ``None`` when no single outcome can occupy
    every cell, which is the one case this construction has nothing to offer.

    ``consistency_factor`` is exactly 1.0 for every candidate target, so what ranks them is the
    remaining ``total_weighted_score x distribution_fidelity`` — both already closed-form here, and
    both seed-independent, so no simulation enters the choice. The same 1.0 goes into
    ``allocate_rows``, which is what stops the split from still paying coverage entropy for a
    cut-rate trade that no longer exists.
    """
    by_target = partition_by_pinned_outcome(pools, seeds)
    best: tuple[float, dict[tuple, int], list[dict], tuple[str, int]] | None = None
    for target, target_pools in by_target.items():
        cells = list(target_pools)
        rows_per_cell = allocate_rows(
            cells, value_per_row, {cell: len(pool) for cell, pool in target_pools.items()},
            rows_wanted, context, cut_probability_by_cell, mutation_weight_by_cell,
            consistency=1.0,
        )
        # No cut masks are passed: every candidate in ``target_pools`` already draws the target on
        # every seed, so there is no intersection left to protect and ``select_for_diversity`` is
        # free to spend the whole selection on its k-mer entropy objective.
        entries, _seeds_held = select_for_diversity(target_pools, rows_per_cell)
        objective = _allocation_objective(
            rows_per_cell, value_per_row, context,
            cut_probability_by_cell, mutation_weight_by_cell, consistency=1.0,
        )
        if best is None or objective > best[0]:
            best = (objective, rows_per_cell, entries, target)
    if best is None:
        return None
    return best[1], best[2], best[3], best[0]


# ---------------------------------------------------------------------------------------------
# The build
# ---------------------------------------------------------------------------------------------

# Ceiling on the automatic coordinate growth below. Reached only by a contract with almost no
# mismatch budget, where one coordinate carries one guide and a cell needs a coordinate per row.
_MAX_SITES_PER_CELL = 64


def _cell_coordinates(
    context: Context, config: Config, by_cas_strand: dict[tuple[str, str], list[Coordinate]],
    per_cell: int,
) -> dict[tuple, list[Coordinate]]:
    """The best coordinates for each cell, ranked by the stage-2 score they can reach.

    A coordinate serves exactly one cell: stage 1 dedups on (cas, start, strand, guide) and the
    mutation is not in that key, so two mutations sharing a coordinate would collide on any guide
    they both used. Cells are filled scarcest-group-first so Cas12a — whose ``TTTV`` PAM is several
    times rarer than Cas9's ``NGG`` — is not left with whatever Cas9 did not want.
    """
    max_distance = context.base_padding if context.proximity_gate else config.pam_search_flank
    cells = [
        (mutation, cas_system, strand)
        for mutation in context.mutations
        for cas_system in context.cas_systems
        for strand in context.strands
    ]
    cells.sort(key=lambda cell: len(by_cas_strand.get((cell[1], cell[2]), ())))

    claimed_keys: set[tuple] = set()
    chosen_by_cell: dict[tuple, list[Coordinate]] = {}
    for cell in cells:
        mutation, cas_system, strand = cell
        mutation_position = context.mutation_map[mutation]
        reachable = [
            coordinate for coordinate in by_cas_strand.get((cas_system, strand), ())
            if coordinate.identity not in claimed_keys
            and abs(coordinate.start - mutation_position) <= max_distance
        ]
        reachable.sort(
            key=lambda coordinate: -reachable_structural(
                coordinate, abs(coordinate.start - mutation_position), context
            )
        )
        picked_coordinates = reachable[:max(1, per_cell)]
        for coordinate in picked_coordinates:
            claimed_keys.add(coordinate.identity)
        chosen_by_cell[cell] = picked_coordinates
    return chosen_by_cell


def _build_pools(
    context: Context, config: Config, coordinates: dict[tuple, list[Coordinate]]
) -> tuple[dict, dict, dict, dict]:
    """Candidate rows for every cell, with each cell's per-row value and cut probability.

    The allocation needs all three before it can choose a split: what a row in the cell is worth,
    how many distinct guides the cell can actually supply, and how often its rows will cut.
    """
    pools: dict[tuple, list[dict]] = {}
    value_per_row: dict[tuple, float] = {}
    cut_probability_by_cell: dict[tuple, float] = {}
    mutation_weight_by_cell: dict[tuple, float] = {}
    for cell, cell_coordinates in coordinates.items():
        mutation, cas_system, _strand = cell
        candidates: list[dict] = []
        for coordinate in cell_coordinates:
            for guide in enumerate_guides(coordinate, context, config.guides_per_coordinate):
                entry = gate_and_score(
                    make_row(coordinate, guide, mutation, context, 0), context
                )
                if entry is not None:
                    candidates.append(entry)
        pools[cell] = candidates
        if not candidates:
            value_per_row[cell] = 0.0
            continue
        best_entry = max(candidates, key=lambda entry: entry["stage2"]["weighted_score"])
        value_per_row[cell] = best_entry["stage2"]["weighted_score"]
        mutation_weight_by_cell[cell] = context.weight_of(mutation)
        # Every guide at a coordinate shares one feature vector, so a cell built on one coordinate
        # has a single cut probability and its contribution to is_cut's error is exact.
        cut_probability_by_cell[cell] = stage3.cut_probability(
            cas_system, stage3.sequence_energy(stage3.extract_features(best_entry))
        )
    return pools, value_per_row, cut_probability_by_cell, mutation_weight_by_cell


# ---------------------------------------------------------------------------------------------
# Two-stage construction
#
# ``retune_cas_mix`` prices the weak/strong cas trade through ``select_for_diversity``'s single
# greedy intersection over every selected row at once. There is a second, more effective shape of
# the same idea, measured on develop (see ``all_cut.py``/``fastgreedy.py``, read via
# ``git show origin/develop:...`` rather than vendored here — this repo has no prefetch loop, no
# GPU, and every contract's mutation set is unique, so develop's disk-cached bank and ~900s prepare
# budget do not transfer):
#
#   1. min-union a GROUP of weak-cas candidates (``min_union_select``) — the complement of their
#      combined failures is a *clean set*, usually far smaller than the full seed support but far
#      easier to be strict over.
#   2. fill the strong-cas cells with candidates that are strict over that clean set (zero fails
#      within it, not over all of ``seed_support``) — the relaxation that makes strictness
#      reachable at all.
#
# The reachability of step 2 is what a naive port stalls on. Every guide substituted at one PAM
# coordinate shares that coordinate's gc/distance/energy (``stage3.sequence_energy`` depends on
# nothing else the substitution can move), so ``cut_probability`` is identical across them and
# strictness over k clean seeds is close to an independent ``cut_probability**k`` per candidate —
# pooling more *variants at the same coordinate* does not raise it. Only a different coordinate
# moves ``distance_to_mutation`` and therefore energy's ``0.6*exp(-d/1500)`` term enough to matter
# (0.990 vs 0.995 changes ``p**k`` by an order of magnitude at k in the hundreds), which is why
# ``_widen_strong_cas_pool`` grows coordinate COUNT for the strong-cas cells specifically rather
# than variants per coordinate — the ordinary growth loop in ``build`` already maximises the
# latter and stops once the total pool clears ``rows_wanted``, far short of what step 2 needs.
# ---------------------------------------------------------------------------------------------

def _widen_strong_cas_pool(
    context: Context, config: Config, by_cas_strand: dict[tuple[str, str], list[Coordinate]],
    coordinates: dict[tuple, list[Coordinate]], strong_cas: str, budget_seconds: float,
) -> tuple[dict[tuple, list[Coordinate]], dict[tuple, list[dict]], dict[tuple, float], dict[tuple, float]]:
    """More PAM coordinates for the strong-cas cells only, priced against a wall-clock budget.

    Widens only the (mutation, strong_cas, strand) cells, leaving every other cell's coordinate
    pool exactly as ``build``'s own growth loop already left it: the weak-cas bank only ever needs
    to be strict over a small min-union ``group_size``, not the whole seed support, so it has
    nothing to gain here that is worth the wall-clock.

    Priced the same way ``scan_cut_support`` prices its own width: one 4x growth step (matching
    ``build``'s own coordinate-growth ladder) is timed live against the actual cost —
    ``_build_pools``' stage 1/2 gating, not coordinate discovery, which is cheap — and the
    remaining budget is extrapolated into candidates-affordable rather than steps-affordable, since
    a step's candidate count is not fixed. How long gating one candidate takes is a property of the
    host, so this is measured rather than assumed, exactly as the seed-support scan is.

    Returns coordinates/pools/value_per_row/cut_probability_by_cell for the strong-cas cells only
    — callers merge these into their own wider dicts, leaving every other cell's untouched.
    """
    strong_cells = [cell for cell in coordinates if cell[1] == strong_cas]
    current: dict[tuple, list[Coordinate]] = {cell: coordinates[cell] for cell in strong_cells}
    if not strong_cells or budget_seconds <= 0:
        return current, {}, {}, {}

    current_per_cell = max((len(current[cell]) for cell in strong_cells), default=1)
    if current_per_cell >= _MAX_SITES_PER_CELL:
        return current, {}, {}, {}

    started = time.monotonic()
    probe_per_cell = min(_MAX_SITES_PER_CELL, current_per_cell * 4)
    probe_coordinates = _cell_coordinates(context, config, by_cas_strand, probe_per_cell)
    probe_coordinates = {cell: probe_coordinates[cell] for cell in strong_cells}
    probe_pools, probe_value, probe_cutp, _weight = _build_pools(context, config, probe_coordinates)
    probe_elapsed = max(time.monotonic() - started, 1e-9)
    probe_candidates = sum(len(pool) for pool in probe_pools.values())

    best_coordinates, best_pools, best_value, best_cutp = (
        probe_coordinates, probe_pools, probe_value, probe_cutp
    )
    remaining = budget_seconds - probe_elapsed
    if remaining <= 0 or probe_per_cell >= _MAX_SITES_PER_CELL or probe_candidates == 0:
        return best_coordinates, best_pools, best_value, best_cutp

    seconds_per_candidate = probe_elapsed / probe_candidates
    affordable_candidates = remaining / max(seconds_per_candidate, 1e-9)
    per_cell = probe_per_cell
    spent_candidates = 0.0
    while per_cell < _MAX_SITES_PER_CELL and spent_candidates < affordable_candidates:
        next_per_cell = min(_MAX_SITES_PER_CELL, per_cell * 4)
        grown_coordinates = _cell_coordinates(context, config, by_cas_strand, next_per_cell)
        grown_coordinates = {cell: grown_coordinates[cell] for cell in strong_cells}
        prior_total = sum(len(best_coordinates[cell]) for cell in strong_cells)
        grown_total = sum(len(grown_coordinates[cell]) for cell in strong_cells)
        if grown_total <= prior_total:
            break  # no further PAM coordinates within reach for these cells
        grown_pools, grown_value, grown_cutp, _weight = _build_pools(
            context, config, grown_coordinates
        )
        grown_candidates = sum(len(pool) for pool in grown_pools.values())
        step_candidates = grown_candidates - probe_candidates
        if step_candidates > affordable_candidates:
            break  # this step alone would overrun what the timed sample priced as affordable
        per_cell = next_per_cell
        best_coordinates, best_pools, best_value, best_cutp = (
            grown_coordinates, grown_pools, grown_value, grown_cutp
        )
        spent_candidates = step_candidates
    return best_coordinates, best_pools, best_value, best_cutp


# Each point costs one real min-union pass over the weak-cas bank, so this is a wall-clock knob
# like ``_CAS_LADDER_POINTS`` — the same coarse-ladder-over-one-interior-maximum shape.
_TWO_STAGE_GROUP_LADDER_POINTS = 7

# Relative gain the two-stage construction must predict over the shipped scan-and-retune result
# before ``build`` takes it. Both predictions go through the identical ``_predicted_objective``
# closed form, but the comparison is still between two estimates of the same unmeasured quantity
# (the true size of the intersection each construction holds), so a thin margin is not evidence —
# the same reasoning ``_RETUNE_MIN_GAIN`` documents.
#
# Whether the construction ever clears this bar is a question of scan width, not of the algorithm:
# ``_widen_strong_cas_pool`` finds far more candidates than a wall-clock budget can afford to scan,
# so ``scan_cut_support``'s rescan of them is the actual bottleneck. Measured on CD34+_HSPC, holding
# the weak-cas bank and everything else fixed, only the rescan's wall-clock budget varying:
#
#     rescan budget   candidates/cell   seeds held   objective ratio
#            45 s             1,136           258            0.962   (loses)
#            90 s             2,354           304            1.008   (barely positive, still below gate)
#           180 s             5,353           375            1.087   (clears the gate)
#
# So the construction genuinely can beat the shipped baseline — but only past roughly 150-180s of
# rescan, which stacked on the ordinary scan's own ~90s leaves as little as 30-60s of a 300s TTL for
# the upload. ``Config.two_stage_seconds``/``Miner.TWO_STAGE_SHARE_OF_WINDOW`` are left at their
# current, safely-conservative values rather than raised to chase this: at those settings the gate
# reliably declines (measured on all four cell types this was tested against, zero regressions), and
# raising them is a build-time-vs-upload-margin trade on a live process, not a threshold to tune
# blind. If that trade is wanted, widen the budget deliberately with the table above as the guide,
# not by nudging this constant.
_TWO_STAGE_MIN_GAIN = 0.03


def _group_size_ladder(floor_total: int, ceiling: int) -> list[int]:
    """Candidate weak-cas group sizes from the largest reachable down to the coverage-entropy
    floor, coarse for the same reason ``cas_mix_ladder``'s ladder is."""
    if ceiling <= floor_total:
        return [ceiling] if ceiling > 0 else []
    sizes = set()
    for point in range(_TWO_STAGE_GROUP_LADDER_POINTS):
        share = point / max(1, _TWO_STAGE_GROUP_LADDER_POINTS - 1)
        sizes.add(int(round(ceiling - share * (ceiling - floor_total))))
    return sorted(sizes)


def two_stage_construction(
    rows_per_cell: dict[tuple, int],
    scanned_pools: dict[tuple, list[dict]], cut_masks: dict[tuple, list[int]],
    seed_count: int, context: Context,
    value_per_row: dict[tuple, float], cut_probability_by_cell: dict[tuple, float],
    strong_cas: str, weak_cas: str,
) -> tuple[dict[tuple, int], list[dict], int, float] | None:
    """Min-union the weak cas system, then require the strong one strict only over the clean set.

    Tries a coarse ladder of weak-cas group sizes (each one real min-union pass) and keeps whichever
    predicts the best ``_predicted_objective`` — the same selection principle ``retune_cas_mix``
    uses, on the same closed form, so the two are comparable on equal terms.

    Cells belonging to neither cas system (a contract naming a third) are left exactly as
    ``rows_per_cell`` already had them, filled by ``weighted_score`` alone — the same thing
    ``_shifted_toward`` does when it moves rows only between the weak and strong systems.

    Returns ``None`` when the weak-cas bank or the strong-cas pool has nothing to build from, or
    when no group size can meet every cell's coverage-entropy floor — the caller keeps the shipped
    baseline in either case.
    """
    weak_cells = [cell for cell in rows_per_cell if cell[1] == weak_cas]
    strong_cells = [cell for cell in rows_per_cell if cell[1] == strong_cas]
    other_cells = [cell for cell in rows_per_cell if cell[1] not in (weak_cas, strong_cas)]
    if not weak_cells or not strong_cells:
        return None

    weak_bank: list[tuple[tuple, int, int]] = []  # (cell, local_index, mask)
    for cell in weak_cells:
        for local_index, mask in enumerate(cut_masks.get(cell) or []):
            weak_bank.append((cell, local_index, mask))
    if not weak_bank:
        return None

    weak_total = sum(rows_per_cell.get(cell, 0) for cell in weak_cells)
    floor_total = len(weak_cells)  # one row each, the coverage-entropy floor
    mutation_weight_by_cell = {cell: context.weight_of(cell[0]) for cell in rows_per_cell}
    capacity_for_shift = {cell: len(scanned_pools.get(cell, ())) for cell in rows_per_cell}

    best: tuple[float, int, dict[tuple, int], list[dict]] | None = None
    for group_size in _group_size_ladder(floor_total, min(weak_total, len(weak_bank))):
        # ``_shifted_toward`` gives the (mutation, strand)-preserving target this group size implies
        # for every cell; weak-cell targets become the min-union's per-cell caps, so the group's
        # shape matches the allocator's own split instead of drifting wherever the union is
        # cheapest, and strong-cell targets are how many rows each strong cell needs filled.
        target_rows = _shifted_toward(
            rows_per_cell, weak_cas, strong_cas, group_size, capacity_for_shift
        )
        weak_floor = {cell: 1 for cell in weak_cells}
        weak_cap = {cell: max(1, target_rows.get(cell, 0)) for cell in weak_cells}
        chosen, clean_mask = min_union_select(
            [mask for _cell, _idx, mask in weak_bank],
            [cell for cell, _idx, _mask in weak_bank],
            seed_count, group_size, weak_floor, weak_cap,
        )
        if len(chosen) < floor_total:
            continue  # the bank could not even meet every weak cell's floor at this group size

        candidate_rows: dict[tuple, int] = dict(rows_per_cell)
        entries: list[dict] = []
        surviving = clean_mask
        for cell in weak_cells:
            candidate_rows[cell] = 0
        for index in chosen:
            cell, local_index, _mask = weak_bank[index]
            candidate_rows[cell] += 1
            entries.append(scanned_pools[cell][local_index])

        feasible = True
        for cell in strong_cells:
            target = max(1, target_rows.get(cell, rows_per_cell.get(cell, 0)))
            cell_pool = scanned_pools.get(cell) or []
            cell_masks = cut_masks.get(cell) or []
            if not cell_pool:
                feasible = False
                break
            # Strict over the clean set first (zero fails within it), ranked by weighted_score;
            # topped up with the best non-strict candidates so a short strict pool never empties
            # the cell -- it only loses some of the rows that carry the intersection.
            strict = [i for i, mask in enumerate(cell_masks) if (mask & clean_mask) == clean_mask]
            strict_set = set(strict)
            rest = [i for i in range(len(cell_pool)) if i not in strict_set]
            strict.sort(key=lambda i: -cell_pool[i]["stage2"]["weighted_score"])
            rest.sort(key=lambda i: -cell_pool[i]["stage2"]["weighted_score"])
            take = (strict + rest)[:target]
            if not take:
                feasible = False
                break
            for i in take:
                entries.append(cell_pool[i])
                if i < len(cell_masks):
                    surviving &= cell_masks[i]
            candidate_rows[cell] = len(take)
        if not feasible:
            continue

        for cell in other_cells:
            target = rows_per_cell.get(cell, 0)
            cell_pool = scanned_pools.get(cell) or []
            cell_masks = cut_masks.get(cell) or []
            ranked = sorted(
                range(len(cell_pool)), key=lambda i: -cell_pool[i]["stage2"]["weighted_score"]
            )
            take = ranked[:target]
            if target > 0 and not take:
                feasible = False
                break
            for i in take:
                entries.append(cell_pool[i])
                if i < len(cell_masks):
                    surviving &= cell_masks[i]
            candidate_rows[cell] = len(take)
        if not feasible or not entries:
            continue
        # A strong cell whose pool came up short of its target silently ships fewer rows than the
        # shipped baseline -- strictly worse (total_weighted_score falls and nothing compensates),
        # not a trade this ladder point should be allowed to win on. Reject rather than repair: the
        # next ladder point asks that cell for fewer rows, which is the actual fix.
        target_total = sum(rows_per_cell.values())
        if sum(candidate_rows.values()) != target_total:
            continue

        seeds_held = surviving.bit_count()
        objective = _predicted_objective(
            candidate_rows, seeds_held, seed_count, context, value_per_row,
            cut_probability_by_cell, mutation_weight_by_cell,
        )
        if best is None or objective > best[0]:
            best = (objective, seeds_held, candidate_rows, entries)

    if best is None:
        return None
    objective, seeds_held, candidate_rows, entries = best
    return candidate_rows, entries, seeds_held, objective


def pin_with_widening(
    context: Context, config: Config, by_cas_strand: dict[tuple[str, str], list[Coordinate]],
    rows_wanted: int, round_seeds: tuple[int, ...], state: tuple,
) -> tuple[tuple | None, tuple]:
    """The exact-outcome build, widening the candidate pool until it can fill every row.

    A candidate survives only if it draws the target under *every* seed, so a cell's usable
    capacity is about ``pool_size * p**k`` for k seeds — it falls geometrically while the rows
    needed stay at 250. One coordinate per cell carries ~900 candidates, which is enough for k=3 on
    HEK293 and k=4 on K562 and runs out fast after that: at k=6 the shipped width filled a mean of
    18/250 rows on HEK293, failing outright (no target reaching all eight cells) on three draws
    out of five.

    The fix is more candidates, and under a pin they are nearly free. ``coordinates_per_cell = 1``
    exists to collapse each cell onto one stage-2 feature vector so stage 4's forest can only
    return group means — a defence against overfitting the targets no seed-blind design can
    predict. A pinned build has no such targets: all three are constant and ``consistency_factor``
    is 1.0 whatever the feature matrix looks like. So the constraint buys nothing here and can be
    spent on pool size, at the cost of some ``dist_score`` as coordinates get further from the
    mutation. Measured ``final_score`` at coordinates_per_cell 1 -> 16:

        HEK293  k=3  301.3 -> 295.4 (-2%)      k=7   23.8 -> 252.9 (10.6x)
        K562    k=3  241.8 -> 241.8 (unchanged) k=7   65.3 -> 234.2 (3.6x)

    So widening is charged at a couple of percent where it is not needed and pays back an order of
    magnitude where it is. It is applied only when the pin cannot already fill the submission, and
    the widened pools are returned *alongside* the originals rather than replacing them: if the pin
    loses its objective check anyway, the seed-blind fallback must run on the narrow pools, whose
    flat feature matrix it genuinely needs.

    When even the widest pool cannot fill the submission, the last resort is to bet on fewer seeds:
    capacity rises by 1/p for each one dropped. Seeds are dropped from the END of ``round_seeds``,
    so **the caller must order them by priority** — seeds it actually believes first, speculative
    ones last (``Miner._seed_plan`` puts the operator's ``seeds.json`` entries ahead of its own
    random draws for exactly this reason). Dropping is preferred over shipping a short submission
    because a seed that buys 1/900 of a lottery is worth far less than the rows it costs.

    ``state`` and the second return value are ``(per_cell, coordinates, pools, value_per_row,
    cut_probability_by_cell, mutation_weight_by_cell, usable_cells)``. Returns the best pinned
    result found, the seeds it is actually pinned to, and the state it belongs to.
    """
    per_cell, coordinates, pools, value_per_row, cut_probability, mutation_weight, usable = state
    best_result, best_state, best_seeds = None, state, round_seeds
    deadline = time.monotonic() + config.pin_widen_seconds

    while True:
        result = pinned_outcome_build(
            {cell: pools[cell] for cell in usable}, rows_wanted, round_seeds, context,
            value_per_row, cut_probability, mutation_weight,
        )
        filled = sum(result[0].values()) if result else 0
        if result is not None and (best_result is None or filled > sum(best_result[0].values())):
            best_result, best_seeds = result, round_seeds
            best_state = (per_cell, coordinates, pools, value_per_row, cut_probability,
                          mutation_weight, usable)
        if filled >= rows_wanted:
            return best_result, best_seeds, best_state
        if per_cell >= _MAX_SITES_PER_CELL:
            break
        if time.monotonic() > deadline:
            logger.info(
                "exact-outcome: pool widening stopped at %d coordinates/cell on its %.0fs budget",
                per_cell, config.pin_widen_seconds,
            )
            break

        grown = min(_MAX_SITES_PER_CELL, per_cell * 4)
        widened_coordinates = _cell_coordinates(context, config, by_cas_strand, grown)
        if sum(len(picked) for picked in widened_coordinates.values()) \
                <= sum(len(picked) for picked in coordinates.values()):
            break  # the genome has no further PAM coordinates within reach to offer
        per_cell, coordinates = grown, widened_coordinates
        pools, value_per_row, cut_probability, mutation_weight = _build_pools(
            context, config, coordinates
        )
        usable = [cell for cell in coordinates if pools[cell]]
        if not usable:
            break
        logger.info(
            "exact-outcome: %d/%d rows at %d coordinates/cell; widening to %d (%d candidates/cell)",
            filled, rows_wanted, per_cell // 4, per_cell,
            min((len(pool) for pool in pools.values()), default=0),
        )

    # Widest pool reached and still short. Give up seeds, lowest-priority first, until the rows fit.
    trimmed = list(round_seeds)
    while len(trimmed) > 1 and time.monotonic() <= deadline:
        dropped, trimmed = trimmed[-1], trimmed[:-1]
        result = pinned_outcome_build(
            {cell: pools[cell] for cell in usable}, rows_wanted, tuple(trimmed), context,
            value_per_row, cut_probability, mutation_weight,
        )
        filled = sum(result[0].values()) if result else 0
        logger.info(
            "exact-outcome: dropped seed %d to reach %d/%d rows on %d seed(s)",
            dropped, filled, rows_wanted, len(trimmed),
        )
        if result is not None and (best_result is None or filled > sum(best_result[0].values())):
            best_result, best_seeds = result, tuple(trimmed)
            best_state = (per_cell, coordinates, pools, value_per_row, cut_probability,
                          mutation_weight, usable)
        if filled >= rows_wanted:
            break

    return best_result, best_seeds, best_state


def build(context: Context, config: Config | None = None) -> tuple[list[dict], list[dict], dict]:
    """Design a submission for this contract.

    Returns the rows to upload, the stage 1-2 entries for each (so a caller can score without
    re-running the gate), and a diagnostic dict.
    """
    config = config or Config()
    rows_wanted = config.row_cap or context.max_experiments
    by_cas_strand = enumerate_coordinates(context, config)

    per_cell = max(1, config.coordinates_per_cell)
    coordinates = _cell_coordinates(context, config, by_cas_strand, per_cell)
    pools, value_per_row, cut_probability_by_cell, mutation_weight_by_cell = _build_pools(
        context, config, coordinates
    )

    # One coordinate per cell is the design, but it only fills 250 rows if the mismatch budget can
    # spell 250 distinct guides on it. It usually can — three free substitutions on a 20-mer give
    # over a thousand — but a contract with `max_mismatches: 0` allows exactly one guide per
    # coordinate, and a cell would then contribute a single row. Widening to more coordinates costs
    # some of the flat feature matrix and is worth it: 8 rows against 250 is the whole of term 1.
    while sum(len(pool) for pool in pools.values()) < rows_wanted \
            and per_cell < _MAX_SITES_PER_CELL:
        grown_per_cell = min(_MAX_SITES_PER_CELL, per_cell * 4)
        widened_coordinates = _cell_coordinates(context, config, by_cas_strand, grown_per_cell)
        widened_total = sum(len(picked) for picked in widened_coordinates.values())
        current_total = sum(len(picked) for picked in coordinates.values())
        if widened_total <= current_total:
            break  # the genome has no further PAM coordinates within reach to offer
        per_cell, coordinates = grown_per_cell, widened_coordinates
        pools, value_per_row, cut_probability_by_cell, mutation_weight_by_cell = _build_pools(
            context, config, coordinates
        )

    usable_cells = [cell for cell in coordinates if pools[cell]]
    unfillable_cells = [cell for cell in coordinates if not pools[cell]]
    if not usable_cells:
        return [], [], {"error": "no cell could produce a valid row"}

    # A contract that carries its round seeds is scored under seeds the design can read, so try the
    # exact-outcome construction before anything that plays odds against an unknown one. It is
    # attempted first rather than compared last because when it wins it makes the seed-support scan
    # pointless rather than merely redundant, and that scan is nearly all of the build's wall clock.
    # The contract's own seeds if it has any, else whatever was supplied out of band. The fallback
    # is what makes this reachable at all in practice: the backend broadcasts `seed: 0`, so on a
    # live round the contract never carries one and the file is the only source there is.
    round_seeds, seed_source = resolve_seeds(context.contract, config.seeds_path)
    round_seeds = tuple(round_seeds)
    seeds_path = config.seeds_path or settings.MINER_SEEDS_PATH
    if seed_source == "disabled":
        logger.info(
            '%s has "enabled": false — ignoring any seeds in it and building seed-blind',
            seeds_path,
        )
    elif round_seeds and seed_source == "file":
        logger.info(
            "Seeds %s supplied by %s (the contract carries none)", list(round_seeds), seeds_path
        )
    pinned_result = None
    # The pool the pin ends up using, which may be wider than the one the seed-blind path wants.
    # Held separately so a pin that loses its objective check below leaves the blind fallback its
    # own narrow pools untouched — see ``pin_with_widening``.
    pinned_state = (per_cell, coordinates, pools, value_per_row, cut_probability_by_cell,
                    mutation_weight_by_cell, usable_cells)
    if config.pin_outcomes and round_seeds:
        try:
            pinned_started = time.monotonic()
            requested_seed_count = len(round_seeds)
            pinned_result, pinned_seeds, pinned_state = pin_with_widening(
                context, config, by_cas_strand, rows_wanted, round_seeds, pinned_state,
            )
            if pinned_result is not None:
                # The pin may have given up seeds to fit the rows, so what it actually holds is
                # what gets reported downstream — reporting the requested set would claim a
                # guarantee over seeds this submission does not in fact hold.
                round_seeds = pinned_seeds
                logger.info(
                    "exact-outcome: pinned to %s over %d of %d seeds, filling %d/%d rows at %d "
                    "coordinate(s)/cell in %.1fs", pinned_result[2], len(pinned_seeds),
                    requested_seed_count, sum(pinned_result[0].values()), rows_wanted,
                    pinned_state[0], time.monotonic() - pinned_started,
                )
            else:
                logger.info(
                    "exact-outcome: no outcome could occupy every cell over %d seeds (%.1fs)",
                    requested_seed_count, time.monotonic() - pinned_started,
                )
        except Exception:
            pinned_result = None
            logger.exception(
                "exact-outcome construction failed; falling back to the seed-blind build"
            )

    # The seed-blind path's own ceiling, and a genuine upper bound on it: every one of its rows at
    # the single best cell's value, perfect coverage entropy, and the consistency it reaches only on
    # the seeds it actually holds. Clearing that means no amount of scanning can catch up, so the
    # scan is skipped outright. Falling short is not a rejection — the two are compared exactly
    # further down, once the blind path has a real number rather than a bound.
    blind_ceiling = (
        rows_wanted * max(value_per_row.values(), default=0.0) * _consistency_estimate_pinned()
    )
    pinned_target: tuple[str, int] | None = None
    if pinned_result is not None and pinned_result[3] > blind_ceiling:
        rows_per_cell, selected_entries, pinned_target, _pinned_objective = pinned_result
        seeds_held, two_stage_used = len(round_seeds), False
        logger.info(
            "exact-outcome: objective %.2f clears the seed-blind ceiling %.2f; skipping the "
            "%.0fs cut-support scan", pinned_result[3], blind_ceiling, config.seed_scan_seconds,
        )
        # Diagnostics describe the pool the pin actually used, which widening may have grown.
        pin_per_cell, pin_coordinates, pin_pools = pinned_state[0], pinned_state[1], pinned_state[2]
        return _finish_build(
            selected_entries, context, rows_wanted, pin_per_cell, pin_coordinates, pin_pools,
            rows_per_cell, pinned_state[6],
            [cell for cell in pin_coordinates if not pin_pools[cell]],
            seeds_held, len(round_seeds), two_stage_used, pinned_target, round_seeds,
        )

    rows_per_cell = allocate_rows(
        usable_cells, value_per_row,
        {cell: len(pools[cell]) for cell in usable_cells}, rows_wanted, context,
        cut_probability_by_cell=cut_probability_by_cell,
        mutation_weight_by_cell=mutation_weight_by_cell,
    )

    # Price each candidate's cut mask over the seed support before selecting, so the k-mer pass can
    # be run underneath the intersection rather than against it. The scan replaces the pools with
    # the slice it actually covered, keeping candidates and masks index-aligned.
    #
    # It has to be wide enough for whichever cas mix wins below, not just the one allocated above:
    # the retune moves rows onto the strong-cas cells, and a cell cannot be selected past what was
    # scanned for it.
    cell_pools = {cell: pools[cell] for cell in usable_cells}
    full_capacity = {cell: len(pool) for cell, pool in cell_pools.items()}
    planned_splits = cas_mix_ladder(rows_per_cell, full_capacity, cut_probability_by_cell)
    widest_allocation = {
        cell: max(split.get(cell, 0) for split in planned_splits) for cell in usable_cells
    }
    scanned_pools, cut_masks = scan_cut_support(
        cell_pools, widest_allocation, tuple(config.seed_support),
        config.seed_scan_candidates_per_cell, config.seed_scan_seconds,
    )
    # The allocator's own split, kept aside for the two-stage attempt below: it prices its own cas
    # mix independently of ``retune_cas_mix`` and has to start from the same place, not from
    # whatever the retune already moved.
    allocated_rows_per_cell = dict(rows_per_cell)
    if cut_masks:
        rows_per_cell, selected_entries, seeds_held = retune_cas_mix(
            rows_per_cell, scanned_pools, cut_masks, len(config.seed_support), context,
            value_per_row, cut_probability_by_cell,
        )
    else:
        selected_entries, seeds_held = select_for_diversity(
            scanned_pools, rows_per_cell, cut_masks, len(config.seed_support)
        )

    # Two-stage construction: an additional attempt tried alongside the scan-and-retune result
    # above, taking whichever predicts the higher final_score. ``Config.two_stage_enabled`` is the
    # A/B knob; off (or a single cas system, or an empty scan) leaves today's result untouched.
    two_stage_used = False
    if config.two_stage_enabled and cut_masks:
        try:
            ranked = _strong_and_weak_cas(usable_cells, cut_probability_by_cell)
        except Exception:
            ranked = None
            logger.exception("two-stage: could not rank cas systems; skipping")
        if ranked is not None:
            strong_cas, weak_cas = ranked
            try:
                two_stage_started = time.monotonic()
                widened_coordinates, widened_pools, widened_value, widened_cutp = (
                    _widen_strong_cas_pool(
                        context, config, by_cas_strand, coordinates, strong_cas,
                        config.two_stage_seconds,
                    )
                )
                strong_cells = [cell for cell in usable_cells if cell[1] == strong_cas]
                two_stage_scanned_pools = dict(scanned_pools)
                two_stage_cut_masks = dict(cut_masks)
                two_stage_value_per_row = dict(value_per_row)
                two_stage_cut_probability = dict(cut_probability_by_cell)
                remaining_budget = config.two_stage_seconds - (time.monotonic() - two_stage_started)
                if widened_pools and remaining_budget > 0:
                    # The widened strong-cas candidates still need cut-support masks before
                    # ``two_stage_construction`` can test strictness against them — the original
                    # scan above never saw these coordinates. ``rank="coordinate_diverse"`` is what
                    # makes the widening pay off: a plain weighted_score cap here would collapse the
                    # pool straight back onto the nearest one or two coordinates (see
                    # ``_round_robin_by_coordinate``). The ceiling is left at the pool's own size —
                    # a wall-clock budget, not a candidate count, is what should decide how wide this
                    # scans, exactly as ``_widen_strong_cas_pool`` was priced.
                    widened_ceiling = max(
                        (len(widened_pools.get(cell, ())) for cell in strong_cells), default=0
                    )
                    rescanned_pools, rescanned_masks = scan_cut_support(
                        {cell: widened_pools.get(cell, []) for cell in strong_cells},
                        {cell: allocated_rows_per_cell.get(cell, 0) for cell in strong_cells},
                        tuple(config.seed_support), widened_ceiling,
                        remaining_budget, rank="coordinate_diverse",
                    )
                    two_stage_scanned_pools.update(rescanned_pools)
                    two_stage_cut_masks.update(rescanned_masks)
                    two_stage_value_per_row.update(widened_value)
                    two_stage_cut_probability.update(widened_cutp)

                two_stage_result = two_stage_construction(
                    allocated_rows_per_cell, two_stage_scanned_pools, two_stage_cut_masks,
                    len(config.seed_support), context,
                    two_stage_value_per_row, two_stage_cut_probability, strong_cas, weak_cas,
                )
            except Exception:
                two_stage_result = None
                logger.exception(
                    "two-stage construction failed; keeping the shipped scan-and-retune result"
                )
            if two_stage_result is not None:
                candidate_rows, candidate_entries, candidate_seeds_held, candidate_objective = (
                    two_stage_result
                )
                baseline_objective = _predicted_objective(
                    rows_per_cell, seeds_held, len(config.seed_support), context,
                    value_per_row, cut_probability_by_cell,
                    {cell: context.weight_of(cell[0]) for cell in rows_per_cell},
                )
                if candidate_objective > baseline_objective * (1 + _TWO_STAGE_MIN_GAIN):
                    rows_per_cell, selected_entries, seeds_held = (
                        candidate_rows, candidate_entries, candidate_seeds_held
                    )
                    two_stage_used = True

    # The exact-outcome build was computed but did not clear the upper bound above, so the two are
    # compared on real numbers instead. Reached only when pinning cost enough rows or enough
    # coverage entropy to be genuinely arguable — which is the one case where a bound is not an
    # answer. Both sides are the same closed form, so this is a comparison and not a heuristic.
    seed_support_size = len(config.seed_support)
    if pinned_result is not None:
        blind_objective = _predicted_objective(
            rows_per_cell, seeds_held, seed_support_size, context, value_per_row,
            cut_probability_by_cell, {cell: context.weight_of(cell[0]) for cell in rows_per_cell},
        )
        if pinned_result[3] > blind_objective:
            rows_per_cell, selected_entries, pinned_target, _pinned_objective = pinned_result
            seeds_held, seed_support_size, two_stage_used = len(round_seeds), len(round_seeds), False
            # Taking the pin means taking the pool it was built from, widened or not.
            per_cell, coordinates, pools = pinned_state[0], pinned_state[1], pinned_state[2]
            usable_cells = pinned_state[6]
            unfillable_cells = [cell for cell in coordinates if not pools[cell]]
            logger.info(
                "exact-outcome: objective %.2f beats the seed-blind build's %.2f; taking it",
                pinned_result[3], blind_objective,
            )
        else:
            logger.info(
                "exact-outcome: objective %.2f does not beat the seed-blind build's %.2f; "
                "keeping the seed-blind result", pinned_result[3], blind_objective,
            )

    return _finish_build(
        selected_entries, context, rows_wanted, per_cell, coordinates, pools, rows_per_cell,
        usable_cells, unfillable_cells, seeds_held, seed_support_size, two_stage_used,
        pinned_target, round_seeds,
    )


def _finish_build(
    selected_entries: list[dict], context: Context, rows_wanted: int, per_cell: int,
    coordinates: dict[tuple, list[Coordinate]], pools: dict[tuple, list[dict]],
    rows_per_cell: dict[tuple, int], usable_cells: list[tuple], unfillable_cells: list[tuple],
    seeds_held: int, seed_support_size: int, two_stage_used: bool,
    pinned_target: tuple[str, int] | None, round_seeds: tuple[int, ...],
) -> tuple[list[dict], list[dict], dict]:
    """Turn a selection into the upload array and its diagnostics.

    Shared by both constructions so neither can acquire its own re-gating or its own idea of what a
    diagnostic means — the numbers the miner logs and the dashboard parses have to describe the two
    paths on identical terms to be comparable at all.
    """
    # Re-number and re-gate in upload order. experiment_id is assigned last because it is the join
    # key stage 4 merges on and the field ``truncate_submission`` dedups, so it has to be unique in
    # exactly the array that gets sent — and the array is ordered strongest-first so that anything
    # the cap ever cuts is the cheapest row, not an arbitrary one. Renumbering is safe for a pinned
    # build because ``stage3.experiment_seed`` hashes the mutation, cas, guide, start and strand and
    # *not* the experiment_id, so neither this nor the sort can move a row's outcome.
    selected_entries.sort(key=lambda entry: -entry["stage2"]["weighted_score"])
    rows: list[dict] = []
    entries: list[dict] = []
    seen_design_keys: set[tuple] = set()
    for entry in selected_entries:
        row = dict(entry["experiment"])
        design_key = (
            row["cas_system"], row["target_alignment_start"], row["strand"], row["guideRNA"]
        )
        if design_key in seen_design_keys:
            continue  # stage 1 would silently drop the second one
        seen_design_keys.add(design_key)
        row["experiment_id"] = f"exp-{len(rows):05d}"
        regated_entry = gate_and_score(row, context)
        if regated_entry is None:
            continue
        rows.append(row)
        entries.append(regated_entry)

    def cell_label(cell: tuple) -> str:
        """Compact cell name for the logs: 'Cas9+:g.5226784G>C'."""
        mutation, cas_system, strand = cell
        return f"{cas_system}{strand}:{mutation[-12:]}"

    diagnostics = {
        "rows": len(rows),
        "rows_wanted": rows_wanted,
        "sites_per_cell": per_cell,
        "coordinates": {
            cell_label(cell): [coordinate.start for coordinate in picked_coordinates]
            for cell, picked_coordinates in coordinates.items()
        },
        "allocation": {
            cell_label(cell): rows_per_cell[cell] for cell in usable_cells
        },
        "pool_sizes": {
            cell_label(cell): len(pools[cell]) for cell in usable_cells
        },
        "empty_cells": [list(cell) for cell in unfillable_cells],
        # Rounds out of ``seed_support`` on which every selected row cuts, which is where
        # ``consistency_factor`` roughly doubles. A lower bound on what gets uploaded: the re-gate
        # above only ever drops rows, and dropping a row can only widen the intersection. On a
        # pinned build the support *is* the contract's seeds, so this reads |seeds| of |seeds|.
        "seeds_holding_cut": seeds_held,
        "seed_support_size": seed_support_size,
        # Whether the two-stage construction (min-union weak cas + strict-over-clean-set strong
        # cas) beat the shipped scan-and-retune result by ``_TWO_STAGE_MIN_GAIN`` and was taken.
        "two_stage_used": two_stage_used,
        # The ``(outcome, indel_length)`` every row is pinned to, or None when the build was seed
        # blind. Set means all three of stage 4's targets are constant and ``consistency_factor``
        # is exactly 1.0 under ``pinned_seeds`` — and only under those seeds.
        "pinned_outcome": list(pinned_target) if pinned_target else None,
        "pinned_seeds": list(round_seeds),
        "predicted_consistency_factor": 1.0 if pinned_target else None,
        "offtarget_factors": dict(
            Counter(entry["features"]["offtarget_factor"] for entry in entries)
        ),
        "gc_scores": dict(
            Counter(round(entry["features"]["gc_score"], 4) for entry in entries)
        ),
        "distinct_feature_vectors": len({
            (
                entry["features"]["gc"], entry["features"]["distance_to_mutation"],
                entry["features"]["gc_score"], entry["features"]["dist_score"],
                entry["features"]["consistency"],
            )
            for entry in entries
        }),
        "total_weighted_score": sum(
            entry["stage2"]["weighted_score"] for entry in entries
        ),
    }
    return rows, entries, diagnostics
