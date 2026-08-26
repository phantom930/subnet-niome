"""Seed-agnostic clustered experiment generation for HEK293 contracts.

HEK293 tasks are broadcast with a provisional seed and may be scored with a different seed.  The
normal generator deliberately searches guide variants against design-keyed Stage-3 outcomes, so
its learned construction does not survive that change. This module never consults such an outcome
while choosing rows:

* materialise and rank Cas9 and Cas12a candidate pools with independent deterministic RNG streams;
* allocate exact global 70/30 Cas and 50/50 strand totals, retaining the existing adaptive
  seed-free candidate as a guaranteed fallback;
* independently build a typed, mutation/Cas-correlated candidate (Cas9 20-mer at 10/20 GC,
  Cas12a 23-mer at 12/23 GC) from exhaustive same-locus Hamming variants;
* discard non-canonical Cas12a TTTT PAMs and sites that cannot reach 40--60% guide GC;
* rank Cas9 by validator distance and Cas12a by its nominal staggered-cut midpoint;
* preserve same-locus feature clusters while selecting globally unique variants by the combined
  submission's marginal 12-mer entropy.

Repeated feature vectors let Stage 4 estimate conditional outcome rates across folds instead of
overfitting one stochastic outcome per feature vector. A paired anonymous Monte Carlo gate chooses
between the two complete candidates from the published Stage-3 probability laws and exact Stage-4
feature groups. Its streams depend only on a fixed version, replicate and canonical row ordinal --
never the contract seed or a guide/design hash. ``stage3.simulate`` is called only after the gate has
irreversibly selected the submitted rows.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import random
from collections import Counter
from dataclasses import dataclass, replace

import genExp as G
from niome_subnet.genomics.validation import stage3, stage5, stage12


@dataclass(frozen=True)
class Hek293ClusterConfig:
    """Fixed knobs measured over the historical HEK293 task set."""

    cas9_share: float = 0.70
    target_gc: float = 0.50
    min_guide_gc: float = 0.40
    max_guide_gc: float = 0.60
    reject_cas12a_tttt: bool = True
    # Baseline mutation allocation. The optimizer also tests 68% and the configured cap, retaining
    # this power-derived allocation when it is better for the task's actual sites and guides.
    mutation_balance_power: float = 1.25
    max_mutation_share: float = 0.72
    adaptive_mutation_share: float = 0.68
    # Generate a wider same-locus candidate set, then select guides by marginal 12-mer entropy.
    # Locus counts stay unchanged, preserving the feature clusters Stage 4 relies on.
    diversity_candidate_multiplier: int = 8
    maximum_variant_request: int = 512
    # Deliberately constant and unrelated to contract.seed.  The share component preserves the
    # exact deterministic configuration used by the 70/30 historical benchmark.
    guide_rng_seed: int = 0xC1A57E + 700
    minimum_variant_request: int = 64
    variant_oversample: int = 2
    # Provisional typed-candidate gate. These fields deliberately live in one immutable config so
    # held-out validation can disable the candidate or adjust the material-gain threshold without
    # changing construction code.
    typed_candidate_enabled: bool = True
    anonymous_gate_replicates: int = 256
    anonymous_gate_folds: int = 3
    anonymous_gate_standard_errors: float = 1.0
    anonymous_gate_minimum_lcb_gain: float = 0.40

    def __post_init__(self) -> None:
        if not 0.0 < self.cas9_share < 1.0:
            raise ValueError("cas9_share must be strictly between 0 and 1")
        if not 0.0 <= self.target_gc <= 1.0:
            raise ValueError("target_gc must be between 0 and 1")
        if not 0.0 <= self.min_guide_gc <= self.max_guide_gc <= 1.0:
            raise ValueError("guide GC bounds must satisfy 0 <= min <= max <= 1")
        if not self.min_guide_gc <= self.target_gc <= self.max_guide_gc:
            raise ValueError("target_gc must fall inside the accepted guide GC range")
        if self.mutation_balance_power < 0.0:
            raise ValueError("mutation_balance_power must be non-negative")
        if not 0.5 <= self.max_mutation_share < 1.0:
            raise ValueError("max_mutation_share must satisfy 0.5 <= share < 1")
        if not 0.5 <= self.adaptive_mutation_share <= self.max_mutation_share:
            raise ValueError(
                "adaptive_mutation_share must be between 0.5 and max_mutation_share"
            )
        if self.diversity_candidate_multiplier < 1 or self.maximum_variant_request < 1:
            raise ValueError("diversity candidate sizes must be positive")
        if self.minimum_variant_request < 1 or self.variant_oversample < 1:
            raise ValueError("variant request sizes must be positive")
        if self.anonymous_gate_replicates < 2:
            raise ValueError("anonymous_gate_replicates must be at least 2")
        if self.anonymous_gate_folds < 2:
            raise ValueError("anonymous_gate_folds must be at least 2")
        if self.anonymous_gate_standard_errors < 0.0:
            raise ValueError("anonymous_gate_standard_errors must be non-negative")
        if self.anonymous_gate_minimum_lcb_gain < 0.0:
            raise ValueError("anonymous_gate_minimum_lcb_gain must be non-negative")


DEFAULT_CONFIG = Hek293ClusterConfig()


@dataclass(frozen=True)
class Hek293GateDiagnostics:
    """Seed-free evidence used to choose the final HEK293 candidate."""

    estimator: str
    replicates: int
    folds: int
    standard_error_multiplier: float
    minimum_lcb_gain: float
    fallback_mean_consistency: float | None
    typed_mean_consistency: float | None
    fallback_mean_final: float | None
    typed_mean_final: float | None
    paired_final_delta_mean: float | None
    paired_final_delta_standard_error: float | None
    paired_final_delta_lcb: float | None
    selected_candidate: str
    typed_candidate_failure: str | None = None


@dataclass
class Hek293ClusteredBuild:
    """A complete build, aligned in the exact order submitted to the validator."""

    rows: list[dict]
    valid: list[dict]
    results: list[dict]
    weight_skew: float
    quotas: dict[tuple[str, str, str], int]
    site_clusters: int
    feature_clusters: int
    mutation_totals: dict[str, int]
    deterministic_tws: float
    deterministic_fidelity: float
    deterministic_objective: float
    quota_candidates_evaluated: int
    quota_candidate_failures: tuple[str, ...]
    # Defaults keep the public return type backward-compatible while an optional typed candidate
    # is disabled or infeasible.
    selected_candidate: str = "adaptive_fallback"
    gate_diagnostics: Hek293GateDiagnostics | None = None

    @property
    def cas_counts(self) -> dict[str, int]:
        return dict(Counter(row["cas_system"] for row in self.rows))

    @property
    def strand_counts(self) -> dict[str, int]:
        return dict(Counter(row["strand"] for row in self.rows))


def _validate_cas_roster(ctx: G.Context) -> None:
    allowed = set(ctx.cas_systems)
    expected = {"Cas9", "Cas12a"}
    if allowed != expected:
        raise ValueError(
            "HEK293 clustered generation requires exactly Cas9 and Cas12a; "
            f"contract allows {sorted(allowed)}"
        )


def _oriented_pam(
    ctx: G.Context,
    cas: str,
    strand: str,
    start: int,
    length: int,
) -> str:
    """Return the reference PAM in guide orientation, mirroring Stage 1 exactly."""

    if cas == "Cas9":
        raw = (
            ctx.seq[start + length:start + length + 3]
            if strand == "+"
            else ctx.seq[start - 3:start]
        )
    elif cas == "Cas12a":
        raw = (
            ctx.seq[start - 4:start]
            if strand == "+"
            else ctx.seq[start + length:start + length + 4]
        )
    else:
        return ""
    return raw if strand == "+" else stage12.reverse_complement(raw)


def _pam_allowed(cas: str, pam: str, config: Hek293ClusterConfig) -> bool:
    """Apply the stricter PAM policy on top of the validator's binary PAM gate.

    The validator treats every accepted PAM equally: Cas9 NGG and Cas12a TTTN. Cas9 therefore has
    no additional modeled weak subgroup. For Cas12a, TTTT is the permissive borderline case; the
    strict canonical set is TTTV, where V is A, C or G.
    """

    if cas == "Cas9":
        return len(pam) == 3 and pam[1:] == "GG"
    if cas == "Cas12a":
        if len(pam) != 4 or pam[:3] != "TTT":
            return False
        return not config.reject_cas12a_tttt or pam[3] in "ACG"
    return False


def _gc_count_bounds(length: int, config: Hek293ClusterConfig) -> tuple[int, int]:
    """Inclusive integer GC-count bounds, avoiding floating-point boundary surprises."""

    minimum = math.ceil(config.min_guide_gc * length - 1e-12)
    maximum = math.floor(config.max_guide_gc * length + 1e-12)
    return minimum, maximum


def _site_can_reach_gc_band(
    site: G.Site,
    ctx: G.Context,
    config: Hek293ClusterConfig,
) -> bool:
    """Whether the mismatch budget can move this site's guide into the accepted GC band."""

    minimum, maximum = _gc_count_bounds(site.length, config)
    reference_gc = sum(base in "GC" for base in site.ref_guide)
    reachable_minimum = max(0, reference_gc - ctx.max_mismatches)
    reachable_maximum = min(site.length, reference_gc + ctx.max_mismatches)
    return reachable_maximum >= minimum and reachable_minimum <= maximum


def _guide_gc_allowed(guide: str, config: Hek293ClusterConfig) -> bool:
    minimum, maximum = _gc_count_bounds(len(guide), config)
    gc_count = sum(base in "GC" for base in guide)
    return minimum <= gc_count <= maximum


def _cas12a_cut_distance_twice(site: G.Site, mutation_position: int) -> int:
    """Twice the distance to the nominal midpoint of Cas12a's staggered 18/23-nt cuts.

    Stage 2 scores alignment-start distance, not biochemical cut distance. This is therefore a
    Cas12a-only selection policy: the scored distance remains available as the next tie-breaker.
    Integer doubled coordinates avoid rounding the 20.5-nt cut midpoint.
    """

    if site.strand == "+":
        cut_low, cut_high = site.start + 18, site.start + 23
    else:
        cut_low, cut_high = site.start + site.length - 23, site.start + site.length - 18
    return abs(2 * mutation_position - (cut_low + cut_high))


def _ranked_pool(
    ctx: G.Context,
    sites: list[G.Site],
    generation_config: G.GenConfig,
    config: Hek293ClusterConfig,
    mutation: str,
    cas: str,
    strand: str,
) -> list[G.Site]:
    position = ctx.mutation_map[mutation]
    max_distance = generation_config.max_distance
    if ctx.contract["rules"].get("proximity_gate", False):
        max_distance = min(max_distance, ctx.base_padding)

    pool = [
        site
        for site in sites
        if site.cas == cas
        and site.strand == strand
        and site.length in generation_config.lengths
        and abs(site.start - position) <= max_distance
        and _pam_allowed(
            site.cas,
            _oriented_pam(ctx, site.cas, site.strand, site.start, site.length),
            config,
        )
        and _site_can_reach_gc_band(site, ctx, config)
    ]
    # Cas9 follows the validator's alignment-start distance. Cas12a has a 5' PAM and staggered
    # distal cuts, so it follows the nominal cut midpoint first and retains validator distance as
    # its next key. The remaining fields only choose between equally close eligible sites.
    def rank(site: G.Site) -> tuple:
        validator_distance = abs(site.start - position)
        quality = (
            G.gc_miss(site, ctx.max_mismatches),
            -stage12.offtarget_uniqueness(site.ref_guide, site.cas, ctx.kmer_index),
            -G.achievable_structural(site, validator_distance, ctx),
            site.start,
            site.length,
            site.ref_guide,
        )
        if site.cas == "Cas12a":
            return (_cas12a_cut_distance_twice(site, position), validator_distance, *quality)
        return (validator_distance, *quality)

    pool.sort(
        key=rank
    )
    return pool


def _largest_remainder(total: int, labels: list[str], weights: list[float]) -> dict[str, int]:
    """Deterministic largest-remainder apportionment over named labels."""

    weight_total = sum(max(weight, 0.0) for weight in weights)
    if weight_total <= 0.0:
        raise ValueError("allocation weights must contain a positive value")
    exact = [total * max(weight, 0.0) / weight_total for weight in weights]
    counts = [int(value) for value in exact]
    remainder = total - sum(counts)
    order = sorted(range(len(labels)), key=lambda index: (-(exact[index] - counts[index]), index))
    for index in order[:remainder]:
        counts[index] += 1
    return dict(zip(labels, counts))


def _power_mutation_totals(
    ctx: G.Context,
    config: Hek293ClusterConfig,
) -> dict[str, int]:
    """The existing power-derived mutation allocation, before its concentration cap."""

    rows_wanted = ctx.max_experiments
    mutation_weights = ctx.contract.get("mutation_weights", {})
    return _largest_remainder(
        rows_wanted,
        ctx.mutations,
        [
            max(float(mutation_weights.get(mutation, 1.0)), 1e-9)
            ** config.mutation_balance_power
            for mutation in ctx.mutations
        ],
    )


def _mutation_total_candidates(
    ctx: G.Context,
    config: Hek293ClusterConfig,
) -> list[dict[str, int]]:
    """Small deterministic mutation frontier measured to cover the useful trade-off.

    HEK293 contracts in the task history have exactly two active mutations. Validating that
    invariant keeps the concentration cap and the 68/72 frontier unambiguous instead of silently
    inventing a redistribution policy for an unseen multi-mutation task.
    """

    if len(ctx.mutations) != 2:
        raise ValueError(
            "HEK293 adaptive mutation allocation requires exactly two active mutations; "
            f"got {len(ctx.mutations)}"
        )
    rows_wanted = ctx.max_experiments
    maximum = math.floor(rows_wanted * config.max_mutation_share + 1e-12)
    minimum_majority = (rows_wanted + 1) // 2
    if maximum < minimum_majority:
        raise ValueError(
            f"max_mutation_share={config.max_mutation_share} cannot allocate {rows_wanted} rows"
        )

    mutation_weights = ctx.contract.get("mutation_weights", {})
    heavy = min(
        ctx.mutations,
        key=lambda mutation: (
            -float(mutation_weights.get(mutation, 1.0)),
            ctx.mutations.index(mutation),
        ),
    )
    light = next(mutation for mutation in ctx.mutations if mutation != heavy)
    power_totals = _power_mutation_totals(ctx, config)
    power_heavy = max(minimum_majority, min(maximum, power_totals[heavy]))
    adaptive_heavy = max(
        minimum_majority,
        min(maximum, int(round(rows_wanted * config.adaptive_mutation_share))),
    )

    candidates: list[dict[str, int]] = []
    seen: set[tuple[int, int]] = set()
    for heavy_total in (power_heavy, adaptive_heavy, maximum):
        totals = {heavy: heavy_total, light: rows_wanted - heavy_total}
        key = tuple(totals[mutation] for mutation in ctx.mutations)
        if key not in seen:
            seen.add(key)
            candidates.append(totals)
    return candidates


def _quotas_from_mutation_totals(
    ctx: G.Context,
    config: Hek293ClusterConfig,
    mutation_totals: dict[str, int],
) -> dict[tuple[str, str, str], int]:
    """Allocate mutation totals over Cas and strand with exact global margins.

    Cas rounding happens once globally rather than independently inside every mutation.
    Mutation/Cas cells are as close to independent as integer margins allow, maximising joint
    entropy for the proposed mutation totals.
    """

    if len(ctx.strands) != 2:
        raise ValueError(f"HEK293 balanced allocation requires two strands, got {ctx.strands}")
    rows_wanted = ctx.max_experiments
    if set(mutation_totals) != set(ctx.mutations) or sum(mutation_totals.values()) != rows_wanted:
        raise ValueError(
            f"mutation totals must cover {ctx.mutations} and sum to {rows_wanted}: "
            f"{mutation_totals}"
        )
    maximum = math.floor(rows_wanted * config.max_mutation_share + 1e-12)
    if max(mutation_totals.values()) > maximum:
        raise ValueError(
            f"mutation allocation exceeds {config.max_mutation_share:.0%} cap: {mutation_totals}"
        )

    cas_totals = _largest_remainder(
        rows_wanted,
        ctx.cas_systems,
        [config.cas9_share if cas == "Cas9" else 1.0 - config.cas9_share
         for cas in ctx.cas_systems],
    )

    # Integer transportation table for mutation x Cas with the exact row/column margins above.
    exact = {
        (mutation, cas): mutation_totals[mutation] * cas_totals[cas] / rows_wanted
        for mutation in ctx.mutations
        for cas in ctx.cas_systems
    }
    mutation_cas = {cell: int(value) for cell, value in exact.items()}
    mutation_missing = {
        mutation: mutation_totals[mutation]
        - sum(mutation_cas[(mutation, cas)] for cas in ctx.cas_systems)
        for mutation in ctx.mutations
    }
    cas_missing = {
        cas: cas_totals[cas]
        - sum(mutation_cas[(mutation, cas)] for mutation in ctx.mutations)
        for cas in ctx.cas_systems
    }
    candidates = sorted(
        exact,
        key=lambda cell: (-(exact[cell] - mutation_cas[cell]),
                          ctx.mutations.index(cell[0]), ctx.cas_systems.index(cell[1])),
    )
    while any(value > 0 for value in mutation_missing.values()):
        chosen = next(
            (cell for cell in candidates
             if mutation_missing[cell[0]] > 0 and cas_missing[cell[1]] > 0),
            None,
        )
        if chosen is None:
            raise RuntimeError("could not reconcile HEK293 mutation/Cas allocation margins")
        mutation_cas[chosen] += 1
        mutation_missing[chosen[0]] -= 1
        cas_missing[chosen[1]] -= 1

    # Split joint cells across strands while enforcing an exact global strand balance.
    first_strand, second_strand = ctx.strands
    first_target = (rows_wanted + 1) // 2
    quotas: dict[tuple[str, str, str], int] = {}
    odd_cells: list[tuple[str, str]] = []
    first_count = 0
    for mutation in ctx.mutations:
        for cas in ctx.cas_systems:
            joint_total = mutation_cas[(mutation, cas)]
            half = joint_total // 2
            quotas[(mutation, cas, first_strand)] = half
            quotas[(mutation, cas, second_strand)] = half
            first_count += half
            if joint_total % 2:
                odd_cells.append((mutation, cas))
    extras_for_first = first_target - first_count
    for index, (mutation, cas) in enumerate(odd_cells):
        strand = first_strand if index < extras_for_first else second_strand
        quotas[(mutation, cas, strand)] += 1

    if sum(quotas.values()) != rows_wanted or any(quota <= 0 for quota in quotas.values()):
        raise RuntimeError(f"invalid HEK293 balanced quotas: {quotas}")
    return quotas


def _balanced_quotas(
    ctx: G.Context,
    config: Hek293ClusterConfig,
) -> dict[tuple[str, str, str], int]:
    """Compatibility helper returning the capped power-derived candidate's quotas."""

    mutation_totals = _mutation_total_candidates(ctx, config)[0]
    return _quotas_from_mutation_totals(ctx, config, mutation_totals)


@dataclass
class _KmerEntropyState:
    """Incremental form of Stage 5's exact rolling-12-mer entropy calculation."""

    counts: Counter
    total: int = 0
    count_log_mass: float = 0.0

    @staticmethod
    def _mass(count: int) -> float:
        return count * math.log2(count) if count > 0 else 0.0

    def ratio_after(self, guide: str) -> float:
        additions = Counter(stage5.extract_kmers(guide, 12))
        new_total = self.total + sum(additions.values())
        if new_total <= 1:
            return 0.0
        new_mass = self.count_log_mass
        for kmer, addition in additions.items():
            old = self.counts[kmer]
            new_mass += self._mass(old + addition) - self._mass(old)
        entropy = math.log2(new_total) - new_mass / new_total
        return entropy / math.log2(new_total)

    def add(self, guide: str) -> None:
        additions = Counter(stage5.extract_kmers(guide, 12))
        for kmer, addition in additions.items():
            old = self.counts[kmer]
            self.count_log_mass += self._mass(old + addition) - self._mass(old)
            self.counts[kmer] += addition
            self.total += addition


def _variant_seed(
    config: Hek293ClusterConfig,
    cas: str,
    mutation: str,
    strand: str,
    site: G.Site,
) -> int:
    """Stable per-site RNG seed; no Cas or cell can perturb another pool's candidate stream."""

    cas_seed = config.guide_rng_seed if cas == "Cas9" else config.guide_rng_seed ^ 0x12A
    key = "|".join(str(value) for value in (
        cas,
        mutation,
        strand,
        site.start,
        site.length,
        site.ref_guide,
    ))
    digest = int(hashlib.sha256(key.encode()).hexdigest()[:16], 16)
    return cas_seed ^ digest


def _variant_candidates(
    ctx: G.Context,
    site: G.Site,
    mutation: str,
    config: Hek293ClusterConfig,
    request: int,
    *,
    combinatorial: bool = False,
) -> list[tuple[dict, dict]]:
    """Materialise a seed-independent same-locus guide pool in Stage-1-valid form.

    The fallback keeps its independently seeded historical tuner. The typed candidate enumerates
    every Hamming-0..3 substitution, retains the one achievable GC count shared by that site's
    variants, ranks by validator off-target uniqueness and then caps the materialised pool. No
    outcome RNG or guide-derived outcome seed is consulted by either path.
    """

    minimum_gc, maximum_gc = _gc_count_bounds(site.length, config)
    target_gc_count = max(
        minimum_gc,
        min(maximum_gc, int(round(config.target_gc * site.length))),
    )
    if combinatorial:
        reference = site.ref_guide
        reference_gc = sum(base in "GC" for base in reference)
        gc_delta = max(
            -ctx.max_mismatches,
            min(ctx.max_mismatches, target_gc_count - reference_gc),
        )
        achieved_gc = reference_gc + gc_delta
        guide_set: set[str] = set()
        if achieved_gc == reference_gc:
            guide_set.add(reference)
        for mismatch_count in range(1, ctx.max_mismatches + 1):
            for positions in itertools.combinations(range(site.length), mismatch_count):
                replacements = [
                    tuple(base for base in "ACGT" if base != reference[position])
                    for position in positions
                ]
                for bases in itertools.product(*replacements):
                    candidate = list(reference)
                    for position, base in zip(positions, bases):
                        candidate[position] = base
                    guide = "".join(candidate)
                    if sum(base in "GC" for base in guide) == achieved_gc:
                        guide_set.add(guide)
        guides = sorted(
            guide_set,
            key=lambda guide: (
                -stage12.offtarget_uniqueness(guide, site.cas, ctx.kmer_index),
                guide,
            ),
        )[:request]
    else:
        rng = random.Random(_variant_seed(config, site.cas, mutation, site.strand, site))
        guides = G.tune_variants(site, target_gc_count, ctx, request, rng)
    candidates: list[tuple[dict, dict]] = []
    for candidate_index, guide in enumerate(guides):
        if not _guide_gc_allowed(guide, config):
            continue
        experiment = G.make_experiment(
            site,
            guide,
            mutation,
            ctx,
            f"candidate-{site.cas}-{site.strand}-{site.start}-{candidate_index}",
        )
        entry = G.build_valid_entry(experiment, ctx)
        if entry is not None:
            candidates.append((experiment, entry))
    return candidates


def _select_diverse_guides(
    candidates: list[tuple[dict, dict]],
    take: int,
    config: Hek293ClusterConfig,
    entropy: _KmerEntropyState,
    used_designs: set[tuple],
    used_guides: set[str],
    experiment_index: int,
) -> tuple[list[tuple[dict, dict]], int]:
    """Choose same-locus variants by quality, then exact marginal 12-mer entropy."""

    remaining = list(candidates)
    selected: list[tuple[dict, dict]] = []
    while len(selected) < take:
        eligible: list[tuple[int, tuple]] = []
        for index, (experiment, entry) in enumerate(remaining):
            guide = experiment["guideRNA"]
            design = G.design_key(experiment)
            if design in used_designs or guide in used_guides:
                continue
            gc_miss = abs(entry["features"]["gc"] - config.target_gc)
            quality = (
                entry["features"]["offtarget_factor"],
                -gc_miss,
                entropy.ratio_after(guide),
                -index,
            )
            eligible.append((index, quality))
        if not eligible:
            break
        best_index, _quality = max(eligible, key=lambda item: item[1])
        experiment, entry = remaining.pop(best_index)
        experiment_id = f"exp-{experiment_index:05d}"
        experiment["experiment_id"] = experiment_id
        # build_valid_entry retains the same experiment dictionary, but assign explicitly to keep
        # that alignment true if its representation ever changes.
        entry["experiment"]["experiment_id"] = experiment_id
        guide = experiment["guideRNA"]
        used_designs.add(G.design_key(experiment))
        used_guides.add(guide)
        entropy.add(guide)
        selected.append((experiment, entry))
        experiment_index += 1
    return selected, experiment_index


class _CandidateCapacityError(RuntimeError):
    """One quota proposal cannot be filled from its eligible guide pools."""


@dataclass
class _SeedFreeCandidate:
    rows: list[dict]
    valid: list[dict]
    quotas: dict[tuple[str, str, str], int]
    mutation_totals: dict[str, int]
    site_clusters: int
    deterministic_tws: float
    deterministic_fidelity: float
    deterministic_objective: float


def _construct_seed_free_candidate(
    ctx: G.Context,
    cells: list[tuple[str, str, str]],
    site_pools: dict[tuple[str, str, str], list[G.Site]],
    quotas: dict[tuple[str, str, str], int],
    mutation_totals: dict[str, int],
    config: Hek293ClusterConfig,
    *,
    cell_configs: dict[tuple[str, str, str], Hek293ClusterConfig] | None = None,
    combinatorial: bool = False,
    exhaustive_cluster: bool = False,
) -> _SeedFreeCandidate:
    """Construct and score one quota proposal without sampling any Stage-3 outcome."""

    built_by_cas: dict[str, list[tuple[dict, dict]]] = {
        cas: [] for cas in ctx.cas_systems
    }
    used_designs_by_cas: dict[str, set[tuple]] = {
        cas: set() for cas in ctx.cas_systems
    }
    # Site ranking and variant RNG remain independent per Cas. Full-guide uniqueness and rolling
    # 12-mer entropy are global because Stage 5 computes both over the combined submission.
    used_guides: set[str] = set()
    entropy = _KmerEntropyState(Counter())
    site_clusters: Counter = Counter()
    experiment_index = 0

    for cas in ctx.cas_systems:
        cas_cells = [cell for cell in cells if cell[1] == cas]
        cas_cells.sort(key=lambda cell: len(site_pools[cell]))
        for mutation, _cas, strand in cas_cells:
            cell = (mutation, cas, strand)
            cell_config = cell_configs.get(cell, config) if cell_configs else config
            wanted = quotas[(mutation, cas, strand)]
            made = 0
            for site in site_pools[(mutation, cas, strand)]:
                if made >= wanted:
                    break
                remaining = wanted - made
                requested = max(
                    cell_config.minimum_variant_request,
                    remaining * cell_config.variant_oversample,
                )
                baseline = _variant_candidates(
                    ctx,
                    site,
                    mutation,
                    cell_config,
                    requested,
                    combinatorial=combinatorial,
                )
                baseline_available = [
                    candidate for candidate in baseline
                    if G.design_key(candidate[0]) not in used_designs_by_cas[cas]
                    and candidate[0]["guideRNA"] not in used_guides
                ]
                cluster_take = min(remaining, len(baseline_available))
                if cluster_take <= 0 and not exhaustive_cluster:
                    continue

                expanded_request = min(
                    cell_config.maximum_variant_request,
                    requested * cell_config.diversity_candidate_multiplier,
                )
                expanded = _variant_candidates(
                    ctx,
                    site,
                    mutation,
                    cell_config,
                    expanded_request,
                    combinatorial=combinatorial,
                )
                candidate_by_design: dict[tuple, tuple[dict, dict]] = {}
                for candidate in expanded + baseline_available:
                    candidate_by_design.setdefault(G.design_key(candidate[0]), candidate)
                if exhaustive_cluster:
                    globally_available = [
                        candidate for candidate in candidate_by_design.values()
                        if G.design_key(candidate[0]) not in used_designs_by_cas[cas]
                        and candidate[0]["guideRNA"] not in used_guides
                    ]
                    cluster_take = min(remaining, len(globally_available))
                selected, experiment_index = _select_diverse_guides(
                    list(candidate_by_design.values()),
                    cluster_take,
                    cell_config,
                    entropy,
                    used_designs_by_cas[cas],
                    used_guides,
                    experiment_index,
                )
                built_by_cas[cas].extend(selected)
                selected_count = len(selected)
                if selected_count:
                    site_clusters[(mutation, cas, strand, site.start, site.length)] += selected_count
                made += selected_count

            if made != wanted:
                raise _CandidateCapacityError(
                    f"cell {(mutation, cas, strand)} built {made}/{wanted} rows"
                )

    merged = [item for cas in ctx.cas_systems for item in built_by_cas[cas]]
    rows = [experiment for experiment, _entry in merged]
    valid = [entry for _experiment, entry in merged]
    rows_wanted = ctx.max_experiments
    used_designs = {G.design_key(row) for row in rows}

    if len(rows) != rows_wanted or len(valid) != rows_wanted:
        raise RuntimeError(
            f"HEK293 clustered build produced {len(rows)} rows; expected {rows_wanted}"
        )
    if len({row["experiment_id"] for row in rows}) != len(rows):
        raise RuntimeError("HEK293 clustered build produced duplicate experiment_id values")
    if len(used_designs) != len(rows):
        raise RuntimeError("HEK293 clustered build produced duplicate validator design keys")
    if len({row["guideRNA"] for row in rows}) != len(rows):
        raise RuntimeError("HEK293 clustered build produced duplicate full guide strings")
    observed_quotas = Counter(
        (row["mutation"], row["cas_system"], row["strand"])
        for row in rows
    )
    if dict(observed_quotas) != quotas:
        raise RuntimeError(
            "HEK293 final joint distribution differs from allocation: "
            f"observed={dict(observed_quotas)}, expected={quotas}"
        )

    ordered_rows = G.order_rows(rows, valid)
    valid_by_id = {entry["experiment"]["experiment_id"]: entry for entry in valid}
    ordered_valid = [valid_by_id[row["experiment_id"]] for row in ordered_rows]

    bad_pam_ids: list[str] = []
    bad_gc_ids: list[str] = []
    for row in ordered_rows:
        pam = _oriented_pam(
            ctx,
            row["cas_system"],
            row["strand"],
            row["target_alignment_start"],
            len(row["guideRNA"]),
        )
        if not _pam_allowed(row["cas_system"], pam, config):
            bad_pam_ids.append(row["experiment_id"])
        if not _guide_gc_allowed(row["guideRNA"], config):
            bad_gc_ids.append(row["experiment_id"])
    if bad_pam_ids or bad_gc_ids:
        raise RuntimeError(
            "HEK293 final quality audit failed: "
            f"bad PAM rows={bad_pam_ids}, bad GC rows={bad_gc_ids}"
        )

    deterministic_tws = sum(
        entry["stage2"]["weighted_score"] for entry in ordered_valid
    )
    # Stage 5 needs only support labels and rules. Strip the seed before this pre-choice call so
    # neither candidate construction nor its deterministic objective can even observe it.
    seed_free_contract = {
        key: value for key, value in ctx.contract.items() if key != "seed"
    }
    fidelity_detail = stage5.compute_distribution_fidelity(
        ordered_valid,
        [],
        seed_free_contract,
        k=12,
    )
    deterministic_fidelity = max(
        0.0,
        min(1.0, fidelity_detail.get("distribution_fidelity_score", 0.0)),
    )
    return _SeedFreeCandidate(
        rows=ordered_rows,
        valid=ordered_valid,
        quotas=quotas,
        mutation_totals=dict(mutation_totals),
        site_clusters=len(site_clusters),
        deterministic_tws=deterministic_tws,
        deterministic_fidelity=deterministic_fidelity,
        deterministic_objective=deterministic_tws * deterministic_fidelity,
    )


def _heavy_and_light_mutations(ctx: G.Context) -> tuple[str, str]:
    mutation_weights = ctx.contract.get("mutation_weights", {})
    heavy = min(
        ctx.mutations,
        key=lambda mutation: (
            -float(mutation_weights.get(mutation, 1.0)),
            ctx.mutations.index(mutation),
        ),
    )
    light = next(mutation for mutation in ctx.mutations if mutation != heavy)
    return heavy, light


def _typed_x170_quotas(ctx: G.Context) -> tuple[
    dict[tuple[str, str, str], int], dict[str, int]
]:
    """The measured 250-row typed allocation with exact Cas and strand margins."""

    if ctx.max_experiments != 250 or len(ctx.mutations) != 2:
        raise _CandidateCapacityError(
            "typed x170 requires exactly two mutations and max_experiments=250"
        )
    if set(ctx.strands) != {"+", "-"}:
        raise _CandidateCapacityError(
            f"typed x170 requires '+'/'-' strands, got {ctx.strands}"
        )
    heavy, light = _heavy_and_light_mutations(ctx)
    quotas = {
        (heavy, "Cas9", "+"): 85,
        (heavy, "Cas9", "-"): 85,
        (heavy, "Cas12a", "+"): 5,
        (heavy, "Cas12a", "-"): 5,
        (light, "Cas9", "+"): 3,
        (light, "Cas9", "-"): 2,
        (light, "Cas12a", "+"): 32,
        (light, "Cas12a", "-"): 33,
    }
    mutation_totals = {heavy: 180, light: 70}
    cas_counts = Counter()
    strand_counts = Counter()
    mutation_cas_counts = Counter()
    for (mutation, cas, strand), count in quotas.items():
        cas_counts[cas] += count
        strand_counts[strand] += count
        mutation_cas_counts[(mutation, cas)] += count
    expected_pairs = {
        (heavy, "Cas9"): 170,
        (heavy, "Cas12a"): 10,
        (light, "Cas9"): 5,
        (light, "Cas12a"): 65,
    }
    if (
        cas_counts != Counter({"Cas9": 175, "Cas12a": 75})
        or strand_counts != Counter({"+": 125, "-": 125})
        or dict(mutation_cas_counts) != expected_pairs
    ):
        raise RuntimeError("internal typed x170 quota margins are inconsistent")
    return quotas, mutation_totals


def _construct_typed_x170_candidate(
    ctx: G.Context,
    sites: list[G.Site],
    generation_config: G.GenConfig,
    config: Hek293ClusterConfig,
) -> _SeedFreeCandidate:
    """Build the Cas-typed large-cluster candidate without any outcome sampling."""

    # Exhaustive same-GC variant materialisation is intentionally calibrated for the task's
    # three-mismatch budget. A larger future budget would make that search combinatorial; retain
    # the bounded adaptive candidate instead of risking an unexpectedly expensive build.
    if ctx.max_mismatches != 3:
        raise _CandidateCapacityError(
            "typed x170 requires the calibrated max_mismatches=3 contract; "
            f"received {ctx.max_mismatches}"
        )

    quotas, mutation_totals = _typed_x170_quotas(ctx)
    cells = list(quotas)
    cas_configs = {
        "Cas9": replace(config, target_gc=0.50),
        "Cas12a": replace(config, target_gc=12 / 23),
    }
    cas_lengths = {"Cas9": (20,), "Cas12a": (23,)}
    cell_configs = {
        cell: cas_configs[cell[1]]
        for cell in cells
    }
    site_pools: dict[tuple[str, str, str], list[G.Site]] = {}
    for cell in cells:
        mutation, cas, strand = cell
        typed_generation_config = replace(
            generation_config,
            lengths=cas_lengths[cas],
        )
        site_pools[cell] = _ranked_pool(
            ctx,
            sites,
            typed_generation_config,
            cell_configs[cell],
            mutation,
            cas,
            strand,
        )
        target_gc_count = int(round(cell_configs[cell].target_gc * cas_lengths[cas][0]))
        # The general quality band accepts 40--60% GC, but the typed marker is stricter: every
        # selected locus must be able to reach exactly 10/20 or 12/23 within the mismatch budget.
        site_pools[cell] = [
            site for site in site_pools[cell]
            if abs(
                sum(base in "GC" for base in site.ref_guide) - target_gc_count
            ) <= ctx.max_mismatches
        ]

    candidate = _construct_seed_free_candidate(
        ctx,
        cells,
        site_pools,
        quotas,
        mutation_totals,
        config,
        cell_configs=cell_configs,
        combinatorial=True,
        exhaustive_cluster=True,
    )
    bad_markers: list[tuple[str, str, int, int, int, int]] = []
    for row in candidate.rows:
        expected_length = cas_lengths[row["cas_system"]][0]
        expected_gc = 10 if row["cas_system"] == "Cas9" else 12
        observed_gc = sum(base in "GC" for base in row["guideRNA"])
        if len(row["guideRNA"]) != expected_length or observed_gc != expected_gc:
            bad_markers.append((
                row["experiment_id"], row["cas_system"], len(row["guideRNA"]),
                observed_gc, expected_length, expected_gc,
            ))
    if bad_markers:
        raise _CandidateCapacityError(
            "typed x170 could not preserve exact Cas length/GC markers for rows "
            f"{bad_markers[:8]}"
        )
    return candidate


@dataclass(frozen=True)
class _AnonymousRow:
    x: tuple[float, float, float, float, float, float, int]
    targets: tuple[float, float, float]
    weight: float


_ANONYMOUS_GATE_VERSION = 1
_ANONYMOUS_OUTCOME_BASE = 0xA110CE
_ANONYMOUS_FOLD_BASE = 0xF01D


def _canonical_seed_free_entries(valid: list[dict]) -> list[dict]:
    """Canonicalise rows using biological labels/features only, never guide sequence or hash."""

    decorated = []
    for original_index, entry in enumerate(valid):
        experiment = entry["experiment"]
        features = entry["features"]
        key = (
            -float(features.get("mutation_weight", 1.0)),
            experiment["mutation"],
            experiment["cas_system"],
            experiment.get("strand", ""),
            float(features["gc"]),
            float(features["distance_to_mutation"]),
            float(features["gc_score"]),
            float(features["dist_score"]),
            float(features["consistency"]),
            original_index,
        )
        decorated.append((key, entry))
    decorated.sort(key=lambda item: item[0])
    return [entry for _key, entry in decorated]


def _anonymous_row_rng(replicate: int, ordinal: int) -> random.Random:
    # Arithmetic mixing is deliberate: only fixed version/replicate/ordinal enter this stream.
    seed = (
        _ANONYMOUS_OUTCOME_BASE
        + _ANONYMOUS_GATE_VERSION * 1_000_003
        + replicate * 104_729
        + ordinal * 15_485_863
    ) & 0xFFFFFFFF
    return random.Random(seed)


def _anonymous_rows(valid: list[dict], replicate: int) -> list[_AnonymousRow]:
    """Sample the public Stage-3 laws without ``stage3.simulate`` or a design-derived seed."""

    rows: list[_AnonymousRow] = []
    for ordinal, entry in enumerate(_canonical_seed_free_entries(valid)):
        features = stage3.extract_features(entry)
        cas = entry["experiment"]["cas_system"]
        rng = _anonymous_row_rng(replicate, ordinal)
        energy = stage3.sequence_energy(features)
        mh = stage3.microhomology_trigger(features, rng)
        if rng.random() > stage3.cut_probability(cas, energy):
            outcome = "no_cut"
            indel = 0
        else:
            outcome = stage3.repair_mode(cas, energy, mh, rng)
            indel = stage3.sample_indel_length(outcome, rng)
        rows.append(_AnonymousRow(
            x=(
                float(features["gc"]),
                float(features["distance"]),
                float(features["gc_score"]),
                float(features["dist_score"]),
                float(features["consistency"]),
                float(energy),
                int(mh),
            ),
            targets=(
                0.0 if outcome == "no_cut" else 1.0,
                1.0 if outcome == "HDR" else 0.0,
                float(indel),
            ),
            weight=float(features["mutation_weight"]),
        ))
    return rows


def _anonymous_fold_ids(row_count: int, folds: int, replicate: int) -> list[int]:
    indices = list(range(row_count))
    seed = (
        _ANONYMOUS_FOLD_BASE
        + _ANONYMOUS_GATE_VERSION * 1_000_033
        + replicate * 130_363
    ) & 0xFFFFFFFF
    random.Random(seed).shuffle(indices)
    fold_ids = [0] * row_count
    for rank, index in enumerate(indices):
        fold_ids[index] = rank % folds
    return fold_ids


def _weighted_r2(actual: list[float], predicted: list[float], weights: list[float]) -> float:
    total_weight = sum(weights)
    if total_weight <= 0.0:
        return 0.0
    mean = sum(weight * value for value, weight in zip(actual, weights)) / total_weight
    residual = sum(
        weight * (value - estimate) ** 2
        for value, estimate, weight in zip(actual, predicted, weights)
    )
    total = sum(weight * (value - mean) ** 2 for value, weight in zip(actual, weights))
    if total <= 1e-15:
        return 1.0 if residual <= 1e-15 else 0.0
    return 1.0 - residual / total


def _group_mean_consistency(rows: list[_AnonymousRow], folds: int, replicate: int) -> float:
    """Approximate Stage 4 by mutation-weighted exact-X train-group means."""

    fold_ids = _anonymous_fold_ids(len(rows), folds, replicate)
    target_stats: list[tuple[float, float]] = []
    for target_index in range(3):
        all_targets = [row.targets[target_index] for row in rows]
        all_mean = sum(all_targets) / len(all_targets)
        scale = math.sqrt(
            sum((value - all_mean) ** 2 for value in all_targets) / len(all_targets)
        )
        fold_r2: list[float] = []
        fold_mae: list[float] = []
        for fold in range(folds):
            group_weight: Counter = Counter()
            group_weighted_target: Counter = Counter()
            global_weight = 0.0
            global_weighted_target = 0.0
            for index, row in enumerate(rows):
                if fold_ids[index] == fold:
                    continue
                value = row.targets[target_index]
                group_weight[row.x] += row.weight
                group_weighted_target[row.x] += row.weight * value
                global_weight += row.weight
                global_weighted_target += row.weight * value
            fallback = global_weighted_target / global_weight if global_weight > 0.0 else 0.0
            actual: list[float] = []
            predicted: list[float] = []
            weights: list[float] = []
            for index, row in enumerate(rows):
                if fold_ids[index] != fold:
                    continue
                actual.append(row.targets[target_index])
                if group_weight[row.x] > 0.0:
                    predicted.append(group_weighted_target[row.x] / group_weight[row.x])
                else:
                    predicted.append(fallback)
                weights.append(row.weight)
            fold_r2.append(_weighted_r2(actual, predicted, weights))
            weight_sum = sum(weights)
            fold_mae.append(
                sum(
                    weight * abs(value - estimate)
                    for value, estimate, weight in zip(actual, predicted, weights)
                ) / weight_sum
                if weight_sum > 0.0 else 0.0
            )
        r2_mean = sum(fold_r2) / len(fold_r2)
        mae_mean = sum(fold_mae) / len(fold_mae)
        normalized_mae = mae_mean if scale < 1e-9 else mae_mean / scale
        target_stats.append((r2_mean, normalized_mae))
    average_r2 = sum(value[0] for value in target_stats) / len(target_stats)
    average_nmae = sum(value[1] for value in target_stats) / len(target_stats)
    return max(
        0.0,
        min(1.0, 0.7 * max(average_r2, 0.0) + 0.3 * (1.0 - average_nmae)),
    )


def _gate_candidates(
    fallback: _SeedFreeCandidate,
    typed: _SeedFreeCandidate,
    config: Hek293ClusterConfig,
) -> tuple[_SeedFreeCandidate, Hek293GateDiagnostics]:
    fallback_consistencies: list[float] = []
    typed_consistencies: list[float] = []
    paired_final_deltas: list[float] = []
    for replicate in range(config.anonymous_gate_replicates):
        fallback_consistency = _group_mean_consistency(
            _anonymous_rows(fallback.valid, replicate),
            config.anonymous_gate_folds,
            replicate,
        )
        typed_consistency = _group_mean_consistency(
            _anonymous_rows(typed.valid, replicate),
            config.anonymous_gate_folds,
            replicate,
        )
        fallback_consistencies.append(fallback_consistency)
        typed_consistencies.append(typed_consistency)
        paired_final_deltas.append(
            typed.deterministic_objective * typed_consistency
            - fallback.deterministic_objective * fallback_consistency
        )
    replicate_count = len(paired_final_deltas)
    fallback_mean_consistency = sum(fallback_consistencies) / replicate_count
    typed_mean_consistency = sum(typed_consistencies) / replicate_count
    delta_mean = sum(paired_final_deltas) / replicate_count
    delta_variance = sum(
        (value - delta_mean) ** 2 for value in paired_final_deltas
    ) / (replicate_count - 1)
    delta_standard_error = math.sqrt(delta_variance / replicate_count)
    delta_lcb = (
        delta_mean
        - config.anonymous_gate_standard_errors * delta_standard_error
    )
    choose_typed = delta_lcb >= config.anonymous_gate_minimum_lcb_gain
    selected_name = "typed_x170" if choose_typed else "adaptive_fallback"
    diagnostics = Hek293GateDiagnostics(
        estimator=f"anonymous_exact_x_group_mean_v{_ANONYMOUS_GATE_VERSION}",
        replicates=replicate_count,
        folds=config.anonymous_gate_folds,
        standard_error_multiplier=config.anonymous_gate_standard_errors,
        minimum_lcb_gain=config.anonymous_gate_minimum_lcb_gain,
        fallback_mean_consistency=fallback_mean_consistency,
        typed_mean_consistency=typed_mean_consistency,
        fallback_mean_final=fallback.deterministic_objective * fallback_mean_consistency,
        typed_mean_final=typed.deterministic_objective * typed_mean_consistency,
        paired_final_delta_mean=delta_mean,
        paired_final_delta_standard_error=delta_standard_error,
        paired_final_delta_lcb=delta_lcb,
        selected_candidate=selected_name,
    )
    return (typed if choose_typed else fallback), diagnostics


def generate_seed_agnostic_clustered(
    ctx: G.Context,
    sites: list[G.Site],
    generation_config: G.GenConfig,
    config: Hek293ClusterConfig = DEFAULT_CONFIG,
) -> Hek293ClusteredBuild:
    """Build HEK293 rows using only seed-free evidence, then simulate the chosen rows once.

    All Cas-specific site pools are materialised before quota proposals are constructed. The
    existing adaptive candidate remains the fallback. If the typed x170 candidate is feasible, a
    fixed anonymous exact-X group-mean Monte Carlo gate selects between the two using the public
    Stage-3 probability laws. Neither the contract seed nor a guide/design-derived outcome stream
    is read, and ``stage3.simulate`` is called only for the already-selected rows.
    """

    if ctx.contract.get("cell_type") != "HEK293":
        raise ValueError(
            "HEK293 clustered generation received cell_type "
            f"{ctx.contract.get('cell_type')!r}"
        )

    _validate_cas_roster(ctx)
    cells = [
        (mutation, cas, strand)
        for mutation in ctx.mutations
        for cas in ctx.cas_systems
        for strand in ctx.strands
    ]
    site_pools = {
        cell: _ranked_pool(
            ctx,
            sites,
            generation_config,
            config,
            cell[0],
            cell[1],
            cell[2],
        )
        for cell in cells
    }

    successful: list[_SeedFreeCandidate] = []
    failures: list[str] = []
    for mutation_totals in _mutation_total_candidates(ctx, config):
        quotas = _quotas_from_mutation_totals(ctx, config, mutation_totals)
        try:
            candidate = _construct_seed_free_candidate(
                ctx,
                cells,
                site_pools,
                quotas,
                mutation_totals,
                config,
            )
        except _CandidateCapacityError as exc:
            failures.append(f"{mutation_totals}: {exc}")
            continue
        successful.append(candidate)

    if not successful:
        detail = "; ".join(failures) if failures else "no quota proposals were generated"
        raise RuntimeError(f"HEK293 adaptive quota allocation found no feasible candidate: {detail}")

    fallback = max(
        successful,
        key=lambda candidate: (
            candidate.deterministic_objective,
            candidate.deterministic_fidelity,
            -max(candidate.mutation_totals.values()),
            tuple(candidate.mutation_totals[mutation] for mutation in ctx.mutations),
        ),
    )

    typed_failure: str | None = None
    typed: _SeedFreeCandidate | None = None
    if config.typed_candidate_enabled:
        try:
            typed = _construct_typed_x170_candidate(
                ctx,
                sites,
                generation_config,
                config,
            )
        except _CandidateCapacityError as exc:
            typed_failure = str(exc)
    else:
        typed_failure = "typed candidate disabled by Hek293ClusterConfig"

    if typed is not None:
        chosen, gate_diagnostics = _gate_candidates(fallback, typed, config)
    else:
        chosen = fallback
        gate_diagnostics = Hek293GateDiagnostics(
            estimator=f"anonymous_exact_x_group_mean_v{_ANONYMOUS_GATE_VERSION}",
            replicates=0,
            folds=config.anonymous_gate_folds,
            standard_error_multiplier=config.anonymous_gate_standard_errors,
            minimum_lcb_gain=config.anonymous_gate_minimum_lcb_gain,
            fallback_mean_consistency=None,
            typed_mean_consistency=None,
            fallback_mean_final=None,
            typed_mean_final=None,
            paired_final_delta_mean=None,
            paired_final_delta_standard_error=None,
            paired_final_delta_lcb=None,
            selected_candidate="adaptive_fallback",
            typed_candidate_failure=typed_failure,
        )

    # This is intentionally the first design-keyed outcome simulation. Candidate construction and
    # the anonymous gate are complete, and chosen rows cannot change after here.
    feature_clusters: Counter = Counter()
    for entry in chosen.valid:
        features = entry["features"]
        energy = stage3.sequence_energy(stage3.extract_features(entry))
        feature_clusters[(
            features["gc"],
            features["distance_to_mutation"],
            features["gc_score"],
            features["dist_score"],
            features["consistency"],
            energy,
        )] += 1
    diagnostic_results = [G.simulate(entry, ctx) for entry in chosen.valid]
    return Hek293ClusteredBuild(
        rows=chosen.rows,
        valid=chosen.valid,
        results=diagnostic_results,
        weight_skew=config.mutation_balance_power,
        quotas=chosen.quotas,
        site_clusters=chosen.site_clusters,
        feature_clusters=len(feature_clusters),
        mutation_totals=chosen.mutation_totals,
        deterministic_tws=chosen.deterministic_tws,
        deterministic_fidelity=chosen.deterministic_fidelity,
        deterministic_objective=chosen.deterministic_objective,
        quota_candidates_evaluated=len(successful),
        quota_candidate_failures=tuple(failures),
        selected_candidate=gate_diagnostics.selected_candidate,
        gate_diagnostics=gate_diagnostics,
    )
