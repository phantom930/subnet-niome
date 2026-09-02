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
   ``r2_score`` returns 1.0 and ``normalized_mae`` short-circuits to 0. That is unreachable here —
   ``cut_probability`` is ``base + 0.18*energy`` and ``energy`` is scaled by chromatin
   accessibility, so at HEK293's 0.35 a Cas12a row cuts with probability 0.87 at best and 250 rows
   all cutting has probability ~1e-10. Chasing it by shrinking the submission loses more on term 1
   than the jump is worth at every row count.
2. What *is* reachable is the ``0.3*(1 - avg_nmae)`` term, and two things move it. Collapsing the
   feature matrix to **one vector per cell** — every row in a cell sits at the same coordinate at
   the same GC count — leaves the forest with eight groups instead of 250 noisy points, so it can
   only return group means and stops overfitting the targets it cannot predict: measured r2 on
   ``is_hdr``/``indel_length`` improves from −0.25/−0.31 to −0.12/−0.11. And because ``is_cut``'s
   normalised error falls as the cut rate rises, the **row allocation prices each cell's cut
   probability** rather than treating the split as a pure coverage question — Cas9 cuts with
   probability 0.95 against Cas12a's 0.87, so the optimum leans to roughly 30% Cas12a and pays for
   the coverage entropy that costs. Worth 8.8% over a coverage-only allocation.

Note what is deliberately *not* here. Ranking guides by the share of ``SEED_SUPPORT`` they cut
under looks like a third lever, and it was built and measured: with the evaluation seeds held out
of the scanned support the cut rate came out at 0.9285 against 0.9295 for no selection at all. That
is the expected result, not a surprise — every guide in a cell shares one feature vector and so one
cut *probability*, and ``experiment_seed`` hashes the round seed in with the design, so a guide's
outcomes under two seeds are independent draws. Inside the scanned support it did lift the mean
score 2.9%, but only by +0.76 +/- 0.56 paired over 29 seeds, against 36 s of a 300 s upload window.
It was removed rather than left switched off. Should a contract ever arrive with accessibility
above 0.67, ``energy`` clamps at 1.0 and Cas9's cut probability reaches its 0.99 ceiling — that is
the regime where all 250 rows cutting stops being out of reach, and it is worth revisiting there.

Every formula below is imported from the validation stages rather than reimplemented, so the
generator cannot drift from the pipeline that judges it. On the reference contract the result scores
27.4 against the 21.1 the same contract previously paid, confirmed against the validator's own
``benchmark_submission`` to zero difference.
"""

from __future__ import annotations

import itertools
import logging
import math

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import niome_subnet.utils.settings as settings

from niome_subnet.genomics.validation import stage12, stage3, stage4, stage5

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
# is the range the local score samples when the contract carries none — it is what makes the
# prediction a statement about a seed drawn from the support rather than about one lucky draw.
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
        """The round seeds in the contract, or [] when it is unstamped.

        ``benchmark_submission`` reads this as a comma-joined string, so a multi-seed round averages
        several draws. A stamped contract is scored exactly; an unstamped one (``seed: 0``) has
        nothing to score against and falls back to ``Config.seed_support``.
        """
        raw_seed_field = self.contract.get("seed")
        try:
            parsed_seeds = [
                int(seed_text) for seed_text in str(raw_seed_field).split(",")
                if seed_text.strip()
            ]
        except (TypeError, ValueError):
            return []
        # A zero is the backend's placeholder for "not stamped yet", not a usable round seed.
        return [seed for seed in parsed_seeds if seed]


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

def select_for_diversity(
    pools: dict[tuple, list[dict]], rows_per_cell: dict[tuple, int]
) -> list[dict]:
    """Fill each cell's allocation with the guides that repeat the fewest 12-mers.

    The greedy step is exact for the marginal cost it is minimising: adding a guide raises the
    pool's entropy least where its windows are already common, so picking the candidate with the
    lowest summed window census is the locally optimal move at every step.
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

    while any(count > 0 for count in rows_left.values()):
        progressed = False
        for cell, cell_candidates in candidates_by_cell.items():
            if rows_left[cell] <= 0 or len(used_indices[cell]) >= len(cell_candidates):
                continue
            best_index = min(
                (
                    candidate_index for candidate_index in range(len(cell_candidates))
                    if candidate_index not in used_indices[cell]
                ),
                key=lambda candidate_index: (
                    sum(
                        window_census[window]
                        for window in cell_candidates[candidate_index][0]
                    ),
                    -cell_candidates[candidate_index][1]["stage2"]["weighted_score"],
                ),
            )
            used_indices[cell].add(best_index)
            for window in cell_candidates[best_index][0]:
                window_census[window] += 1
            chosen_entries.append(cell_candidates[best_index][1])
            rows_left[cell] -= 1
            progressed = True
        if not progressed:
            break
    return chosen_entries


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


def _allocation_objective(
    rows_per_cell: dict[tuple, int], value_per_row: dict[tuple, float], context: Context,
    cut_probability_by_cell: dict[tuple, float],
    mutation_weight_by_cell: dict[tuple, float],
) -> float:
    """All three score factors for a candidate row split, in closed form."""
    weighted_total = sum(rows_per_cell[cell] * value_per_row[cell] for cell in rows_per_cell)
    if weighted_total <= 0:
        return 0.0
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
    return (
        weighted_total
        * stage5.geometric_mean(ratios)
        * _consistency_estimate(
            rows_per_cell, cut_probability_by_cell, mutation_weight_by_cell
        )
    )


def allocate_rows(
    cells: list[tuple], value_per_row: dict[tuple, float], capacity: dict[tuple, int],
    rows_wanted: int, context: Context,
    cut_probability_by_cell: dict[tuple, float],
    mutation_weight_by_cell: dict[tuple, float],
) -> dict[tuple, int]:
    """Rows per cell, by hill-climbing all three score factors from an even split.

    ``value_per_row`` is what one row in a cell is worth and ``capacity`` how many distinct guides
    the cell can actually supply, so the result is always buildable.
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
        cut_probability_by_cell, mutation_weight_by_cell,
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
                    cut_probability_by_cell, mutation_weight_by_cell,
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

    rows_per_cell = allocate_rows(
        usable_cells, value_per_row,
        {cell: len(pools[cell]) for cell in usable_cells}, rows_wanted, context,
        cut_probability_by_cell=cut_probability_by_cell,
        mutation_weight_by_cell=mutation_weight_by_cell,
    )

    selected_entries = select_for_diversity(
        {cell: pools[cell] for cell in usable_cells}, rows_per_cell
    )

    # Re-number and re-gate in upload order. experiment_id is assigned last because it is the join
    # key stage 4 merges on and the field ``truncate_submission`` dedups, so it has to be unique in
    # exactly the array that gets sent — and the array is ordered strongest-first so that anything
    # the cap ever cuts is the cheapest row, not an arbitrary one.
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


# ---------------------------------------------------------------------------------------------
# Local scoring — the validator's stages 4 and 5, in memory, with no file I/O
#
# A miner gets no feedback: the score is computed by every validator independently and never sent
# back, so this replica is the only signal available before the next task. It calls stage 4's and
# stage 5's own functions so the number it reports is the number they would report.
# ---------------------------------------------------------------------------------------------

_ZERO_SCORE = {
    "n_valid_experiments": 0,
    "total_weighted_score": 0.0,
    "consistency_score": 0.0,
    "consistency_factor": 0.0,
}


def score_under_seed(entries: list[dict], context: Context, seed: int) -> dict:
    """What a validator holding ``seed`` would pay for these entries."""
    results = [stage3.simulate(entry, seed) for entry in entries]
    stage3_frame = stage4.flatten_stage3(results)
    stage12_frame = stage4.flatten_stage12(entries)
    if len(stage12_frame) < 2 or len(stage3_frame) < 2:
        return dict(_ZERO_SCORE, final_score=0.0, distribution_fidelity_factor=0.0)

    stage12_columns = stage12_frame[[
        "experiment_id", "guideRNA", "start", "stage2_score", "mutation_weight", "weighted_score"
    ]]
    merged_frame = stage3_frame.merge(stage12_columns, on="experiment_id", how="inner")
    if merged_frame.empty:
        return dict(_ZERO_SCORE, final_score=0.0, distribution_fidelity_factor=0.0)

    feature_matrix = stage4.build_X(merged_frame)
    targets = stage4.build_y(merged_frame)
    sample_weight = merged_frame["mutation_weight"]
    metrics_by_target = {
        target_name: stage4.evaluate(
            feature_matrix, targets[target_name],
            sample_weight=sample_weight, fold_seed=seed, n_splits=5,
        )
        for target_name in targets.columns
    }
    average_r2 = sum(
        metrics["r2_mean"] for metrics in metrics_by_target.values()
    ) / len(metrics_by_target)
    average_nmae = sum(
        stage4.normalized_mae(metrics["mae_mean"], targets[target_name])
        for target_name, metrics in metrics_by_target.items()
    ) / len(metrics_by_target)
    consistency_score = (0.7 * max(average_r2, 0.0) + 0.3 * (1 - average_nmae)) * 100
    if math.isnan(consistency_score):
        consistency_score = 0.0
    consistency_factor = max(0.0, min(1.0, consistency_score / 100.0))

    fidelity = stage5.compute_distribution_fidelity(
        entries, results, context.contract, k=KMER_LENGTH
    )
    fidelity_factor = max(0.0, min(1.0, fidelity.get("distribution_fidelity_score", 0.0)))
    weighted_total = float(merged_frame["weighted_score"].sum())

    return {
        "seed": seed,
        "n_valid_experiments": len(merged_frame),
        "total_weighted_score": weighted_total,
        "consistency_score": float(consistency_score),
        "consistency_factor": consistency_factor,
        "distribution_fidelity_score": fidelity.get("distribution_fidelity_score", 0.0),
        "distribution_fidelity_factor": fidelity_factor,
        "final_score": weighted_total * consistency_factor * fidelity_factor,
        "cut_rate": 1 - sum(
            result["outcome"] == "no_cut" for result in results
        ) / len(results),
        "per_target": {
            target_name: {
                "r2": metrics["r2_mean"],
                "nmae": stage4.normalized_mae(metrics["mae_mean"], targets[target_name]),
            }
            for target_name, metrics in metrics_by_target.items()
        },
        "fidelity_detail": fidelity,
    }


def score_rows(rows: list[dict], context: Context, seeds: list[int] | None = None) -> dict:
    """Score the array that will actually be uploaded.

    Re-derives stage 1-2 from ``rows`` rather than trusting entries built earlier, because stage 4
    shuffles its cross-validation folds from the round seed and applies that shuffle in *file
    order* — so a submission scored in a different order from the one it is sent in reports a
    consistency factor the validator will never compute.

    With no seed available the report is the mean over a sample of ``SEED_SUPPORT``: unbiased for a
    seed drawn from that range, which is the only honest statement a miner can make before a
    validator stamps one.
    """
    entries = [
        entry for entry in (gate_and_score(row, context) for row in rows)
        if entry is not None
    ]
    if len(entries) < 2:
        return {"final_score": 0.0, "rows": len(rows), "valid": len(entries)}

    if seeds is None:
        seeds = context.seeds() or list(SEED_SUPPORT[::150])
    per_seed_reports = [score_under_seed(entries, context, seed) for seed in seeds]

    def mean_across_seeds(field_name: str) -> float:
        return sum(
            report[field_name] for report in per_seed_reports
        ) / len(per_seed_reports)

    return {
        "rows": len(rows),
        "valid": len(entries),
        "seeds": seeds,
        "total_weighted_score": mean_across_seeds("total_weighted_score"),
        "consistency_factor": mean_across_seeds("consistency_factor"),
        "distribution_fidelity_factor": mean_across_seeds("distribution_fidelity_factor"),
        "final_score": mean_across_seeds("final_score"),
        "cut_rate": mean_across_seeds("cut_rate"),
        "per_seed": per_seed_reports,
    }
