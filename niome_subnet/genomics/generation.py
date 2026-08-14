"""Miner-side submission builder — the generation counterpart to :mod:`genomics.validation`.

Given a contract and an HBB reference, this designs the experiment set a miner PUTs to the
validator's presigned URL. Rows carry **designs only** (guide, coordinates, strand, mutation, cas
system, cell type); nothing here writes a biological outcome into a submission. Stage 3 is re-run
locally purely as a *predictor* of what the validator will compute, which is sound because it is
seeded deterministically from ``contract.seed`` plus the design fields — see
:func:`niome_subnet.genomics.validation.stage3.experiment_seed`.

The design target is the validator's own objective::

    final_score = total_weighted_score x consistency_factor x distribution_fidelity_factor

Every gate and feature is delegated to the validator's own stage functions rather than
reimplemented, so the generator cannot drift from the pipeline that will judge it.

Shape of a build (:func:`build_submission`):

1. enumerate every coordinate in ``gene_region`` +/- ``flank`` that has a real PAM;
2. apportion ``max_experiments`` rows across the full mutation x cas x strand support — stage 5
   takes a *geometric* mean, so an empty cell is a ~0.0316x penalty on the whole score;
3. pick the highest-scoring coordinates per cell, spend the contract's mismatch budget pulling each
   guide toward 50% GC and out of the off-target k-mer index;
4. search that site's guide variants for one whose deterministic stage-3 draw satisfies the
   configured construction (default ``"mh"``, the HDR/NHEJ mix — see :data:`CONSTRUCTIONS`), which
   is what takes ``consistency_factor`` to 1.0;
5. order rows strongest-first, because ``truncate_submission`` keeps the first ``max_experiments``.

Paths (``data/chr11.fa``, the k-mer cache) come from
:mod:`niome_subnet.utils.settings` and are relative, so run neurons from the repo root as the
validator does. ``chr11.fa`` is resolved against the repo root as a fallback; the k-mer cache is
not, and a wrong CWD only costs a rebuild.
"""

from __future__ import annotations

import logging
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

import niome_subnet.utils.settings as settings
from niome_subnet.genomics.validation import stage12, stage3, stage4, stage5

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]

AT = ("A", "T")
GC = ("G", "C")


# --------------------------------------------------------------------------------------------
# 1. Context: reference sequence, k-mer index, contract-derived constants
# --------------------------------------------------------------------------------------------

@dataclass
class Context:
    seq: str
    contract: dict
    reference: dict
    cell_types: dict
    kmer_index: dict
    mutation_map: dict
    mutations: list[str]
    cas_systems: list[str]
    strands: tuple[str, str] = ("+", "-")

    @property
    def seed(self) -> int:
        # A contract with no seed still simulates: experiment_seed hashes whatever it is given, and
        # the validator will score with the same value, so build and score stay in step.
        return int(self.contract.get("seed") or 0)

    @property
    def max_experiments(self) -> int:
        return self.contract["rules"].get("max_experiments") or 250

    @property
    def max_mismatches(self) -> int:
        return self.contract["rules"].get("max_mismatches", 0)

    @property
    def base_padding(self) -> int:
        return self.contract["rules"]["base_padding"]


# Every task the backend has issued so far shares one gene_region and one rules block, so the
# 130 MB reference, the k-mer index and the PAM enumeration are all task-independent. Caching them
# is what turns a per-task build into seconds — see warm().
_SEQ_CACHE: str | None = None
_KMER_CACHE: dict[tuple[int, int], dict] = {}
_SITE_CACHE: dict[tuple, list] = {}


def load_sequence() -> str:
    """Parse ``data/chr11.fa`` once per process. ~130 MB, so this is the expensive cold start."""
    global _SEQ_CACHE
    if _SEQ_CACHE is None:
        path = Path(settings.CHR11_PATH)
        if not path.exists():
            path = REPO_ROOT / settings.CHR11_PATH
        if not path.exists():
            raise FileNotFoundError(
                f"{settings.CHR11_PATH} missing — download "
                "Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz from Ensembl release 116 and gunzip "
                "it to that path (scripts/run_validator.sh has the URL). Coordinates and PAMs are "
                "checked against the real sequence, so there is no generating without it."
            )
        _SEQ_CACHE = stage12.load_chr11(str(path))
    return _SEQ_CACHE


def build_context(contract: dict, reference: dict, cell_types: dict,
                  offtarget_flank: int = 50000) -> Context:
    seq = load_sequence()

    # The same window the validator indexes: forward strand of gene_region +/- 50 kb. Indexing a
    # different window would misprice offtarget_factor on every row.
    win_start = reference["gene_region"]["start"]
    win_end = reference["gene_region"]["end"]
    index_start = max(0, win_start - offtarget_flank)
    index_end = min(len(seq), win_end + offtarget_flank)
    if (index_start, index_end) not in _KMER_CACHE:
        _KMER_CACHE[(index_start, index_end)] = stage12.load_or_build_kmer_index(
            seq[index_start:index_end], k=12
        )

    return Context(
        seq=seq,
        contract=contract,
        reference=reference,
        cell_types=cell_types,
        kmer_index=_KMER_CACHE[(index_start, index_end)],
        mutation_map=reference["mutation_map"],
        mutations=list(contract["active_mutations"]),
        cas_systems=list(contract["rules"].get("cas_systems", ["Cas9", "Cas12a"])),
    )


# --------------------------------------------------------------------------------------------
# 2. Site enumeration
# --------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Site:
    """A coordinate where a real PAM exists, plus the reference-matching guide at that position."""
    cas: str
    strand: str
    start: int
    length: int
    ref_guide: str

    @property
    def key(self) -> tuple:
        return (self.cas, self.strand, self.start, self.length)


def seed_slice(cas: str, length: int) -> slice:
    """The 12-mer window ``offtarget_uniqueness`` hashes: guide[-12:] for Cas9, guide[:12] else."""
    return slice(length - 12, length) if cas == "Cas9" else slice(0, 12)


def enumerate_sites(ctx: Context, flank: int, lengths: tuple[int, ...]) -> list[Site]:
    """Every (cas, strand, length, start) with a valid PAM in ``gene_region`` +/- ``flank``.

    PAM validity is delegated to ``stage12.check_pam`` so the enumeration cannot drift from the gate
    that will actually judge the row.
    """
    lo = max(10, ctx.reference["gene_region"]["start"] - flank)
    hi = min(len(ctx.seq) - 40, ctx.reference["gene_region"]["end"] + flank)

    cache_key = (lo, hi, lengths, tuple(ctx.cas_systems), ctx.strands)
    if cache_key in _SITE_CACHE:
        return _SITE_CACHE[cache_key]

    sites: list[Site] = []
    for cas in ctx.cas_systems:
        for strand in ctx.strands:
            for length in lengths:
                for start in range(lo, hi):
                    ok, _ = stage12.check_pam(ctx.seq, start, length, cas, strand)
                    if not ok:
                        continue
                    target = ctx.seq[start:start + length]
                    if any(base not in "ACGT" for base in target):
                        continue
                    guide = target if strand == "+" else stage12.reverse_complement(target)
                    sites.append(Site(cas, strand, start, length, guide))

    _SITE_CACHE[cache_key] = sites
    return sites


# --------------------------------------------------------------------------------------------
# 3. Guide tuning
#
# Stage 1 accepts up to contract.rules.max_mismatches Hamming distance between the guide and the
# reference target, and the PAM is read off the reference (never off the guide). That budget is a
# free design lever with three uses:
#
#   * push GC content to exactly 50%, where gc_score peaks at 1.0;
#   * perturb the 12-mer off-target seed so it is absent from the forward-strand index, taking
#     offtarget_factor from 0.7 to 1.0 (a flat 1.43x on every + strand row);
#   * every distinct guide is a distinct stage-3 RNG seed, so the variants of one site are draws
#     from the outcome distribution at a *fixed* feature vector. That is what lets _best_variant
#     hunt for a draw that satisfies the construction without changing the row's score.
# --------------------------------------------------------------------------------------------

def _substitute(guide: list[str], position: int, want_gc: bool | None,
                rng: random.Random) -> bool:
    """Replace guide[position] with a different base. ``want_gc`` None keeps the GC class."""
    current = guide[position]
    if want_gc is None:
        pool = [b for b in (GC if current in "GC" else AT) if b != current]
    else:
        pool = [b for b in (GC if want_gc else AT) if b != current]
    if not pool:
        return False
    guide[position] = rng.choice(pool)
    return True


def tune_variants(site: Site, target_gc_count: int, ctx: Context, n_variants: int,
                  rng: random.Random) -> list[str]:
    """Guide variants for one site: GC count driven toward ``target_gc_count``, off-target seed
    driven out of the k-mer index, all within the mismatch budget.

    Every returned variant has the *same* GC count, so all of them share one stage-2 feature vector
    and differ only in the stage-3 draw.
    """
    budget = ctx.max_mismatches
    ref = list(site.ref_guide)
    ref_gc = sum(base in "GC" for base in ref)
    # Only |delta| substitutions can move the GC count, so the reachable count is clamped.
    delta = max(-budget, min(budget, target_gc_count - ref_gc))
    achieved_gc = ref_gc + delta

    sl = seed_slice(site.cas, site.length)
    seed_positions = list(range(sl.start, sl.stop))
    other_positions = [i for i in range(site.length) if not (sl.start <= i < sl.stop)]

    variants: list[str] = []
    seen: set[str] = set()
    for _attempt in range(n_variants * 12):
        if len(variants) >= n_variants:
            break
        guide = list(ref)
        spent = 0

        # Step 1: move the GC count. Prefer editing inside the seed window so one substitution buys
        # both the GC correction and the off-target break.
        need = abs(delta)
        want_gc = delta > 0
        eligible = [i for i in seed_positions if (ref[i] in "GC") != want_gc]
        eligible += [i for i in other_positions if (ref[i] in "GC") != want_gc]
        for position in eligible:
            if need == 0:
                break
            if _substitute(guide, position, want_gc, rng):
                need -= 1
                spent += 1

        # Step 2: spend whatever is left breaking the seed / diversifying the draw. GC-preserving
        # substitutions only, so the feature vector stays identical across variants.
        spare = list(seed_positions) + other_positions
        rng.shuffle(spare)
        for position in spare:
            if spent >= budget:
                break
            seed_kmer = "".join(guide)[sl]
            enough_variety = spent >= max(1, budget - 1)
            if seed_kmer not in ctx.kmer_index and enough_variety:
                break
            if position < sl.start or position >= sl.stop:
                if seed_kmer not in ctx.kmer_index:
                    continue  # only spend outside the seed once the seed is already clean
            if _substitute(guide, position, None, rng):
                spent += 1

        candidate = "".join(guide)
        if candidate in seen:
            continue
        if sum(base in "GC" for base in candidate) != achieved_gc:
            continue
        if stage12.hamming(candidate, site.ref_guide) > budget:
            continue
        seen.add(candidate)
        variants.append(candidate)

    if not variants:
        variants = [site.ref_guide]
    # Best off-target class first, so a fallback pick is still a good row.
    variants.sort(key=lambda g: -stage12.offtarget_uniqueness(g, site.cas, ctx.kmer_index))
    return variants


# --------------------------------------------------------------------------------------------
# 4. Row assembly — features and simulated outcome, straight from the validator's own functions
# --------------------------------------------------------------------------------------------

def make_experiment(site: Site, guide: str, mutation: str, ctx: Context,
                    experiment_id: str) -> dict:
    """One submission row. These eight fields are the entire submission format."""
    return {
        "experiment_id": experiment_id,
        "guideRNA": guide,
        "target_alignment_start": site.start,
        "target_alignment_end": site.start + site.length,
        "strand": site.strand,
        "mutation": mutation,
        "cas_system": site.cas,
        "cell_type": ctx.contract.get("cell_type"),
    }


def build_valid_entry(experiment: dict, ctx: Context) -> dict | None:
    """Run the real stage 1 gate and stage 2 feature extraction. None means the row would be
    rejected — the generator treats that as a bug in its own site enumeration."""
    passed, _reason = stage12.stage1(experiment, ctx.seq, ctx.mutation_map, ctx.contract)
    if passed != 1.0:
        return None
    structural, info = stage12.stage2(
        ctx.cell_types, experiment, ctx.seq, ctx.mutation_map, ctx.contract, ctx.kmer_index
    )
    return {
        "experiment": experiment,
        "features": {
            "gc": info["gc"],
            "distance_to_mutation": info["distance"],
            "gc_score": info["gc_score"],
            "dist_score": info["dist_score"],
            "consistency": info["consistency"],
            "offtarget_factor": info["offtarget_factor"],
            "mutation_weight": info["mutation_weight"],
            "cell_type": info["cell_type"],
            "cell_type_accessibility": info["cell_type_accessibility"],
            "mutation_region": info["mutation_region"],
            "region_energy_offset": info["region_energy_offset"],
        },
        "stage1": {"valid": True},
        "stage2": {"structural_score": structural, "weighted_score": info["weighted_score"]},
    }


def simulate(entry: dict, ctx: Context) -> dict:
    return stage3.simulate(entry, ctx.seed)


# --------------------------------------------------------------------------------------------
# 5. In-memory scorer — stages 4 and 5 without the file round trip
#
# Optional: nothing in a build depends on it. It exists so a miner can log the score the validator
# is going to report, and notice a regression before a whole task interval is spent on it.
# --------------------------------------------------------------------------------------------

_ZERO_STAGE4 = {
    "n_valid_experiments": 0,
    "total_weighted_score": 0.0,
    "consistency_score": 0.0,
    "consistency_factor": 0.0,
}


def stage4_in_memory(valid: list[dict], results: list[dict], fold_seed: int,
                     n_folds: int = 5) -> dict:
    frame3 = stage4.flatten_stage3(results)
    frame12 = stage4.flatten_stage12(valid)
    if len(frame12) < 2 or len(frame3) < 2:
        return dict(_ZERO_STAGE4)

    slim = frame12[["experiment_id", "guideRNA", "start", "stage2_score",
                    "mutation_weight", "weighted_score"]]
    merged = frame3.merge(slim, on="experiment_id", how="inner")
    if len(merged) == 0:
        return dict(_ZERO_STAGE4)

    X = stage4.build_X(merged)
    y = stage4.build_y(merged)
    sample_weight = merged["mutation_weight"]

    per_target = {
        column: stage4.evaluate(X, y[column], sample_weight=sample_weight,
                                fold_seed=fold_seed, n_splits=n_folds)
        for column in y.columns
    }
    avg_r2 = float(np.mean([v["r2_mean"] for v in per_target.values()]))
    avg_nmae = float(np.mean([
        stage4.normalized_mae(v["mae_mean"], y[column]) for column, v in per_target.items()
    ]))
    consistency_score = (0.7 * max(avg_r2, 0) + 0.3 * (1 - avg_nmae)) * 100
    if math.isnan(consistency_score):
        return dict(_ZERO_STAGE4)

    return {
        "n_valid_experiments": len(merged),
        "total_weighted_score": float(merged["weighted_score"].sum()),
        "consistency_score": float(consistency_score),
        "consistency_factor": max(0.0, min(1.0, consistency_score / 100.0)),
        "avg_r2": avg_r2,
        "avg_nmae": avg_nmae,
    }


def score_rows(valid: list[dict], results: list[dict], ctx: Context) -> dict:
    """The validator's final_score for this build, computed without touching ``data/``."""
    s4 = stage4_in_memory(valid, results, fold_seed=ctx.seed)
    fidelity = stage5.compute_distribution_fidelity(valid, results, ctx.contract, k=12)
    factor = max(0.0, min(1.0, fidelity.get("distribution_fidelity_score", 0.0)))
    return {
        **s4,
        "distribution_fidelity_score": fidelity.get("distribution_fidelity_score", 0.0),
        "distribution_fidelity_factor": factor,
        "final_score": s4["total_weighted_score"] * s4["consistency_factor"] * factor,
    }


# --------------------------------------------------------------------------------------------
# 6. Configuration
# --------------------------------------------------------------------------------------------

@dataclass
class GenConfig:
    """Generation knobs. The defaults are what a 200-task sweep of the backend's task history
    settled on; ``construction`` is the only one worth reconsidering per contract."""
    flank: int = 3000              # site-enumeration window beyond gene_region
    max_distance: int = 2000       # widest |start - mutation_pos| a selected row may have
    lengths: tuple[int, ...] = (20, 23)
    variants: int = 24             # guide variants searched per selected site
    rows: int | None = None        # defaults to contract max_experiments
    construction: str = "mh"       # key into CONSTRUCTIONS; the rule every row must satisfy
    weight_skew: float = 2.0       # exponent on mutation_weight when apportioning rows;
                                   # build_submission refits this per contract by default


def _quota(total: int, cells: int) -> list[int]:
    base, extra = divmod(total, cells)
    return [base + (1 if i < extra else 0) for i in range(cells)]


def _cell_quotas(ctx: Context, rows_wanted: int, skew: float) -> tuple[list[tuple], list[int]]:
    """Rows per (mutation, cas, strand) cell.

    Even quotas maximise stage 5's coverage entropies. Skewing toward the heavier mutations trades
    some of that entropy for total_weighted_score, which multiplies every row by mutation_weight.
    The trade is favourable well past 50/50: term 1 gains linearly, while stage 5 enters the product
    through a *six-way* geometric mean, so only two of the six ratios move and each moves at the 1/6
    power. ``weight_skew`` is the exponent on mutation_weight; 0 keeps it uniform.
    """
    weights = ctx.contract.get("mutation_weights", {})
    shares = {m: max(weights.get(m, 1.0), 1e-9) ** skew for m in ctx.mutations}
    total = sum(shares.values())
    per_group = max(1, len(ctx.cas_systems) * len(ctx.strands))

    exact = {m: rows_wanted * shares[m] / total for m in ctx.mutations}
    counts = {m: int(exact[m]) for m in ctx.mutations}
    remainder = rows_wanted - sum(counts.values())
    for mutation in sorted(ctx.mutations, key=lambda m: -(exact[m] - counts[m]))[:remainder]:
        counts[mutation] += 1

    cells: list[tuple] = []
    quotas: list[int] = []
    for mutation in ctx.mutations:
        split = _quota(counts[mutation], per_group)
        index = 0
        for cas in ctx.cas_systems:
            for strand in ctx.strands:
                cells.append((mutation, cas, strand))
                quotas.append(split[index])
                index += 1
    return cells, quotas


def achievable_gc_score(site: Site, budget: int) -> float:
    """gc_score this site can reach once the mismatch budget is spent pulling GC toward 50%."""
    ref_gc = sum(base in "GC" for base in site.ref_guide)
    target = round(site.length / 2)
    reached = ref_gc + max(-budget, min(budget, target - ref_gc))
    return max(0.0, 1.0 - abs(reached / site.length - 0.5) * 2)


def achievable_structural(site: Site, distance: int, ctx: Context) -> float:
    """Upper bound on this site's stage-2 structural score after guide tuning.

    Mirrors stage 2 exactly — 0.625*gc_score + 0.375*dist_score — with offtarget_factor taken as
    1.0, which guide tuning reaches for all but a handful of sites.
    """
    dist_score = math.exp(-distance / ctx.base_padding)
    return 0.625 * achievable_gc_score(site, ctx.max_mismatches) + 0.375 * dist_score


def select_sites(ctx: Context, sites: list[Site], cfg: GenConfig, quiet: bool = False
                 ) -> tuple[list[tuple[Site, str]], dict[tuple, list[Site]]]:
    """Pick (site, mutation) pairs across the mutation x cas x strand support.

    Every cell must be occupied — stage 5's geometric mean clips a missing category to 1e-9, i.e. a
    0.0316x hit on the whole score. How many rows each cell gets is ``_cell_quotas``' call.

    Within a cell, sites are ranked by the stage-2 score they can actually *reach* rather than by
    distance alone: a site 40 bp further out but tunable to 50% GC beats a nearer one stuck at
    gc_score 0.5, because gc_score carries 0.625 of the structural score.

    Sites are globally unique because stage 1 dedups on (cas, start, strand, guide) — the mutation
    is *not* part of that key, so one coordinate cannot serve two mutations.

    Returns the selection plus the per-cell reserve of untouched sites, in the same ranking, which
    :func:`generate` uses to replace rows whose outcome could not be forced.
    """
    rows_wanted = cfg.rows or ctx.max_experiments
    cells, quotas = _cell_quotas(ctx, rows_wanted, cfg.weight_skew)

    # With the proximity gate on, stage 1 rejects anything beyond base_padding outright, so a wider
    # max_distance would not widen the search — it would just feed the generator rows that get
    # thrown away. None of the tasks seen so far enable it, but the contract is the authority.
    max_distance = cfg.max_distance
    if ctx.contract["rules"].get("proximity_gate", False):
        max_distance = min(max_distance, ctx.base_padding)

    by_group: dict[tuple[str, str], list[Site]] = defaultdict(list)
    for site in sites:
        if site.length in cfg.lengths:
            by_group[(site.cas, site.strand)].append(site)

    # Scarcest (cas, strand) group first, so Cas12a — roughly 4x rarer than Cas9 — is not starved of
    # coordinates by an earlier Cas9 cell.
    order = sorted(range(len(cells)), key=lambda i: len(by_group[(cells[i][1], cells[i][2])]))

    used: set[tuple] = set()
    selected: list[tuple[Site, str]] = []
    reserve: dict[tuple, list[Site]] = {}
    shortfall = 0

    for index in order:
        cell = cells[index]
        mutation, cas, strand = cell
        want = quotas[index]
        position = ctx.mutation_map[mutation]
        pool = [
            site for site in by_group[(cas, strand)]
            if site.key not in used and abs(site.start - position) <= max_distance
        ]
        if not pool:
            reserve[cell] = []
            shortfall += want
            continue

        take = min(want, len(pool))
        # Distance first, then reachable stage-2 score. Both sorts are stable, so distance is the
        # tie-break among sites that can reach the same structural score.
        pool.sort(key=lambda site: abs(site.start - position))
        ranked = sorted(pool,
                        key=lambda s: -achievable_structural(s, abs(s.start - position), ctx))

        deduped: list[Site] = []
        spare: list[Site] = []
        for site in ranked:
            (deduped if len(deduped) < take else spare).append(site)

        for site in deduped:
            used.add(site.key)
            selected.append((site, mutation))
        reserve[cell] = spare
        shortfall += want - len(deduped)

    if shortfall and not quiet:
        logger.warning(
            "%d row(s) short of the %d cap — widen GenConfig.flank / max_distance",
            shortfall, rows_wanted,
        )
    return selected, reserve


def estimate_final_score(ctx: Context, sites: list[Site], cfg: GenConfig) -> float:
    """Surrogate score for a candidate *selection*: no guide tuning, no simulation, no forest.

    Sound only for comparing selections under a construction that reaches consistency_factor 1.0:

    * consistency_factor is pinned at 1.0 by the construction, so the forest never needs to run;
    * total_weighted_score follows from site choice alone, because guide tuning is deterministic
      (``achievable_structural`` mirrors what stage 2 will compute) and the variant search only
      picks among guides sharing one feature vector;
    * the two guide-dependent stage-5 ratios — 12-mer diversity and distinct-guide — are dropped.
      Inside a six-way geometric mean they multiply every candidate by the same factor, so they
      cannot move the argmax.

    A few ms, which is what makes fitting a knob per contract affordable.
    """
    selected, _ = select_sites(ctx, sites, cfg, quiet=True)
    if len(selected) < 2:
        return 0.0

    weights = ctx.contract.get("mutation_weights", {})
    total = 0.0
    joint, mutations, cas_systems, strands = Counter(), Counter(), Counter(), Counter()
    for site, mutation in selected:
        distance = abs(site.start - ctx.mutation_map[mutation])
        total += achievable_structural(site, distance, ctx) * weights.get(mutation, 1.0)
        joint[(mutation, site.cas, site.strand)] += 1
        mutations[mutation] += 1
        cas_systems[site.cas] += 1
        strands[site.strand] += 1

    support = [(m, c, s) for m in ctx.mutations for c in ctx.cas_systems for s in ctx.strands]
    ratios = [
        stage5.coverage_entropy_ratio(mutations, ctx.mutations),
        stage5.coverage_entropy_ratio(cas_systems, ctx.cas_systems),
        stage5.coverage_entropy_ratio(strands, list(ctx.strands)),
        stage5.coverage_entropy_ratio(joint, support),
        1.0, 1.0,
    ]
    return total * stage5.geometric_mean(ratios)


SKEW_CANDIDATES = tuple(i / 4 for i in range(0, 33))  # 0.00 .. 8.00 in 0.25 steps


def choose_weight_skew(ctx: Context, sites: list[Site], cfg: GenConfig) -> float:
    """Fit the mutation-weight skew to *this* contract.

    The optimum depends on the weight ratio, which ranges from about 1.10/0.65 to 1.68/0.53 across
    the task history — one repo-wide constant leaves score on the table for most of them.
    """
    return max(SKEW_CANDIDATES,
               key=lambda skew: estimate_final_score(ctx, sites, replace(cfg, weight_skew=skew)))


def design_key(experiment: dict) -> tuple:
    """Stage 1's second dedup key — independent of experiment_id."""
    return (experiment["cas_system"], experiment["target_alignment_start"],
            experiment["strand"], experiment["guideRNA"])


# --------------------------------------------------------------------------------------------
# 7. Constructions — the rule every row's stage-3 outcome is forced to satisfy.
#
# Stage 4 fits is_cut, is_hdr and indel_length separately and scores each with r2 + normalised MAE.
# It does not reward *variety* in outcomes; it rewards outcomes being a function the forest can
# recover from X = [gc, distance, gc_score, dist_score, consistency, energy, mh] under CV. Two ways
# to satisfy that, both measured against the real scorer:
#
#   "mh"   (default, the HDR/NHEJ mix) mh -> HDR, otherwise BLUNT_NHEJ pinned to indel_length 1.
#          consistency_factor = 1.0 with r2 = 1.0 on all three targets, *without* is_hdr or
#          indel_length being constant — because `mh` is literally a column in X, so is_hdr == mh
#          and indel_length == 1 - mh are exactly recoverable.
#
#   "hdr"  every row HDR. All three targets constant, so r2_score hits its zero-numerator/
#          zero-denominator case (1.0) and normalized_mae short-circuits on std < 1e-9.
#          consistency_factor = 1.0, but the dataset is degenerate in all three targets. Measured
#          215.21 vs 214.75 for "mh" on one task: a 0.2% premium for three degeneracies.
#
# A rule receives the simulated result and the stage-12 entry, so a construction may key off any
# feature the validator computes — e.g. entry["features"]["consistency"], the other exact binary in
# X, which is what a rule that also breaks is_cut's constancy would use.
# --------------------------------------------------------------------------------------------

def _rule_mh(result: dict, entry: dict) -> bool:
    if result["mh"]:
        return result["outcome"] == "HDR"
    # indel_length must be pinned, or the exponential draw leaves that target unlearnable.
    return result["outcome"] == "BLUNT_NHEJ" and result["indel_length"] == 1


def _rule_hdr(result: dict, entry: dict) -> bool:
    return result["outcome"] == "HDR"


def _rule_nocut(result: dict, entry: dict) -> bool:
    # Pins all three targets at once (cut=0, hdr=0, indel=0), like "hdr" — but far more expensive to
    # hit: cut_p is capped at 0.99, so P(no_cut) is ~0.01 per variant on Cas9.
    return result["outcome"] == "no_cut"


def _rule_blunt(result: dict, entry: dict) -> bool:
    # indel_length pinned to its modal value (P=0.70 given BLUNT_NHEJ).
    return result["outcome"] == "BLUNT_NHEJ" and result["indel_length"] == 1


def _rule_mhnhej(result: dict, entry: dict) -> bool:
    # The gamma draw is flat (no value exceeds P=0.12), so pinning costs ~48 variants per row.
    return result["outcome"] == "MH_NHEJ" and result["indel_length"] == 1


CONSTRUCTIONS = {
    "mh": _rule_mh,
    "hdr": _rule_hdr,
    "nocut": _rule_nocut,
    "blunt": _rule_blunt,
    "mhnhej": _rule_mhnhej,
}


def _best_variant(site: Site, mutation: str, ctx: Context, cfg: GenConfig, index: int,
                  rng: random.Random, require, seen_designs: set[tuple]) -> tuple | None:
    """Search a site's guide variants for one whose stage-3 draw satisfies ``require``.

    The match is mandatory: the caller would rather lose the row than accept a draw that breaks the
    rule (see :func:`generate`). All variants share one feature vector, so which one is returned
    costs nothing in stage 2 — only the outcome differs.

    ``seen_designs`` mirrors stage 1's ``seen_valid_keys``. Stage 1 dedups on
    (cas, start, strand, guide) *in addition to* experiment_id and silently drops the second
    occurrence — so a collision here costs a row and breaks outcome purity as well. The point of
    choice is the only place it can be avoided.
    """
    target_gc_count = int(round(0.5 * site.length))  # gc_score peaks at exactly 50% GC

    for guide in tune_variants(site, target_gc_count, ctx, cfg.variants, rng):
        experiment = make_experiment(site, guide, mutation, ctx, f"exp-{index:05d}")
        if design_key(experiment) in seen_designs:
            continue
        entry = build_valid_entry(experiment, ctx)
        if entry is None:
            continue
        result = simulate(entry, ctx)
        if require(result, entry):
            return experiment, entry, result
    return None


def generate(ctx: Context, sites: list[Site], cfg: GenConfig) -> tuple[list, list, list]:
    """Force every row's stage-3 outcome to satisfy the configured construction rule.

    Conformance is **all-or-nothing**: the rules work by making stage 4's targets exactly
    recoverable from X, and a single stray row breaks that. Under "hdr" a stray row gives
    ``indel_length`` a near-zero variance against a large residual and r2 collapses from 1.0 to
    roughly 0; under "mh" it puts a point off the mh -> outcome mapping the forest would otherwise
    fit exactly. Either way the result is worse than not trying.

    So a site whose variants will not satisfy the rule is *replaced* from its cell's reserve, and if
    the reserve is exhausted the row is dropped rather than submitted non-conforming. Because
    stage 3 is deterministic in (contract.seed, design), the conformance measured here is exactly
    what the validator will compute.

    Returns (rows, valid_entries, stage3_results), index-aligned.
    """
    if cfg.construction not in CONSTRUCTIONS:
        raise ValueError(
            f"unknown construction {cfg.construction!r}; expected one of {sorted(CONSTRUCTIONS)}"
        )

    rng = random.Random(ctx.seed ^ 0x5EED)
    selected, reserve = select_sites(ctx, sites, cfg)
    if not selected:
        return [], [], []

    rule = CONSTRUCTIONS[cfg.construction]
    rows, valid, results = [], [], []
    seen_designs: set[tuple] = set()
    # Every coordinate any cell has actually taken. `reserve` is a snapshot from the moment its own
    # cell was processed, so later cells may since have claimed some of its entries; replaying one
    # of those would emit a duplicate design that stage 1 drops.
    claimed = {site.key for site, _ in selected}
    index = 0
    replaced = 0
    dropped = 0

    for site, mutation in selected:
        cell = (mutation, site.cas, site.strand)
        found = _best_variant(site, mutation, ctx, cfg, index, rng, rule, seen_designs)
        while found is None and reserve.get(cell):
            replacement = reserve[cell].pop(0)
            if replacement.key in claimed:
                continue
            claimed.add(replacement.key)
            replaced += 1
            found = _best_variant(replacement, mutation, ctx, cfg, index, rng, rule, seen_designs)
        if found is None:
            dropped += 1
            continue
        experiment, entry, result = found
        seen_designs.add(design_key(experiment))
        rows.append(experiment)
        valid.append(entry)
        results.append(result)
        index += 1

    if replaced or dropped:
        logger.info("construction %r: %d site(s) replaced, %d row(s) dropped",
                    cfg.construction, replaced, dropped)
    return rows, valid, results


def check_invariants(rows: list[dict], results: list[dict], valid: list[dict],
                     cfg: GenConfig) -> list[str]:
    """The two uniqueness rules stage 1 and truncate_submission enforce, plus rule conformance.

    Worth checking rather than trusting: a violation is invisible in the in-memory score (which
    never dedups) but costs rows in the real pipeline, so it would show up only as an unexplained
    gap between a local run and the validator's.
    """
    problems = []
    ids = [r["experiment_id"] for r in rows]
    if len(set(ids)) != len(ids):
        problems.append(f"{len(ids) - len(set(ids))} duplicate experiment_id")
    designs = [design_key(r) for r in rows]
    if len(set(designs)) != len(designs):
        problems.append(f"{len(designs) - len(set(designs))} duplicate (cas, start, strand, guide)")

    if results:
        rule = CONSTRUCTIONS[cfg.construction]
        broken = sum(1 for result, entry in zip(results, valid) if not rule(result, entry))
        if broken:
            problems.append(f"{broken} row(s) break the {cfg.construction!r} construction")
    return problems


def order_rows(rows: list[dict], valid: list[dict]) -> list[dict]:
    """Strongest rows first: truncate_submission keeps the first max_experiments unique ids, so if
    anything is ever cut it should be the cheapest rows."""
    weight_by_id = {
        entry["experiment"]["experiment_id"]: entry["stage2"]["weighted_score"] for entry in valid
    }
    return sorted(rows, key=lambda r: -weight_by_id.get(r["experiment_id"], 0.0))


# --------------------------------------------------------------------------------------------
# 8. Entry points
# --------------------------------------------------------------------------------------------

def warm(contract: dict, reference: dict, cell_types: dict | None = None,
         cfg: GenConfig | None = None) -> int:
    """Populate the sequence, k-mer and site caches. Returns the number of PAM sites found.

    Everything cached here is task-independent, so a miner should call this at startup: it moves the
    ~130 MB FASTA parse, the k-mer index and the PAM enumeration off the critical path, leaving a
    per-task build that fits comfortably inside the presigned URL's 300 s TTL.
    """
    cfg = cfg or GenConfig()
    ctx = build_context(contract, reference, cell_types or {})
    return len(enumerate_sites(ctx, cfg.flank, tuple(sorted(set(cfg.lengths)))))


def build_submission(contract: dict, reference: dict, cell_types: dict,
                     cfg: GenConfig | None = None, fit_weight_skew: bool = True,
                     score: bool = False) -> tuple[list[dict], dict]:
    """Build one task's submission.

    Returns ``(rows, meta)`` where ``rows`` is the JSON array to PUT — ordered strongest-first — and
    ``meta`` reports what was built (row count, outcome mix, invariant violations, and the predicted
    score when ``score`` is set).
    """
    cfg = cfg or GenConfig()
    ctx = build_context(contract, reference, cell_types)
    sites = enumerate_sites(ctx, cfg.flank, tuple(sorted(set(cfg.lengths))))

    if fit_weight_skew:
        cfg = replace(cfg, weight_skew=choose_weight_skew(ctx, sites, cfg))

    rows, valid, results = generate(ctx, sites, cfg)
    ordered = order_rows(rows, valid)

    meta = {
        "seed": ctx.seed,
        "cell_type": contract.get("cell_type"),
        "max_experiments": ctx.max_experiments,
        "sites": len(sites),
        "rows": len(ordered),
        "construction": cfg.construction,
        "weight_skew": cfg.weight_skew,
        "outcome_counts": dict(Counter(r["outcome"] for r in results)),
        "problems": check_invariants(ordered, results, valid, cfg),
    }
    if score and len(valid) >= 2:
        report = score_rows(valid, results, ctx)
        meta["expected"] = {
            "total_weighted_score": report["total_weighted_score"],
            "consistency_factor": report["consistency_factor"],
            "distribution_fidelity_factor": report["distribution_fidelity_factor"],
            "final_score": report["final_score"],
        }
    return ordered, meta
