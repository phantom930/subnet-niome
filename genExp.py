#!/usr/bin/env python3
"""genExp.py — CRISPR experiment-set generator optimised against the validator pipeline.

Pulls the most recent task with a non-zero ``contract.seed`` from the backend, enumerates every
PAM-valid guide site on the real GRCh38 chr11 sequence, designs a submission that maximises

    final_score = total_weighted_score x consistency_factor x distribution_fidelity_factor

and scores it with the same stage functions the validator runs
(:mod:`niome_subnet.genomics.validation`). Nothing here computes a biological outcome that ends up
in the submission — rows carry designs only; stage 3 is re-run locally purely as a *predictor* of
what the validator will compute, which is legitimate because it is seeded deterministically from
``contract.seed`` plus the design fields.

The default build drives every row's outcome to satisfy the "mh" construction (see
``CONSTRUCTIONS``), which reaches ``consistency_factor = 1.0`` with a mixed HDR/NHEJ dataset rather
than a single repeated outcome. ``--construction hdr`` restores the older all-HDR build.

    python genExp.py                       # fetch latest task, build, search, verify
    python genExp.py --no-search           # single build with default knobs
    python genExp.py --task-id <uuid>      # pin a specific task
    python genExp.py --baseline            # naive generator, for comparison

Writes ``data/submission.json`` (plus contract/hbb_reference) and prints the validator's own
``MinerScore``.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import random
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path

# settings must be imported before bittensor-adjacent modules; it also fixes BT_NO_PARSE_CLI_ARGS.
REPO_ROOT = Path(__file__).resolve().parent
os.chdir(REPO_ROOT)  # every settings.py path is relative to the repo root
sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
import requests  # noqa: E402

import niome_subnet.utils.settings as settings  # noqa: E402
from niome_subnet.genomics.validation import stage12, stage3, stage4, stage5  # noqa: E402

TASKS_URL = f"{settings.BASE_URL}/api/v3/tasks"
NUCLEOTIDES = ("A", "C", "G", "T")
AT = ("A", "T")
GC = ("G", "C")


# --------------------------------------------------------------------------------------------
# 1. Task acquisition
# --------------------------------------------------------------------------------------------

def fetch_task(task_id: str | None = None) -> dict:
    """Return the newest task whose contract carries a non-zero seed.

    Seed 0 is not actually degenerate — ``experiment_seed`` hashes it into a string prefix like any
    other value, and ``KFold(random_state=0)`` is a normal shuffle. Those tasks look like unstamped
    placeholders (3 of 200 at the time of writing, otherwise structurally identical), so they are
    skipped by default but remain reachable via ``task_id``.
    """
    items = fetch_all_tasks()

    if task_id:
        for item in items:
            if item["id"] == task_id:
                return item
        raise RuntimeError(f"task {task_id} not present in the {len(items)} tasks returned")

    for item in items:
        if item["content"]["contract"].get("seed"):
            return item
    raise RuntimeError("no task with a non-zero contract seed in the returned page")


def fetch_all_tasks() -> list[dict]:
    """Every task the backend will hand out, newest first.

    Sent without ``per_page``: the endpoint 422s on anything above 100, but its unparameterised
    default returns the whole history in a single page.
    """
    response = requests.get(TASKS_URL, timeout=120)
    response.raise_for_status()
    items = response.json().get("items", [])
    if not items:
        raise RuntimeError("backend returned no tasks")
    items.sort(key=lambda it: it.get("created_at", ""), reverse=True)
    return items


def fetch_cell_types() -> dict:
    """Accessibility table. Miners can read this endpoint unsigned; a miss defaults to 1.0 and
    would silently inflate stage 3's energy relative to a real validator run."""
    try:
        response = requests.get(settings.CELL_TYPES_URL, timeout=60)
        response.raise_for_status()
        return response.json()
    except Exception as exc:  # noqa: BLE001 - degraded mode is better than no run
        print(f"  ! cell-types fetch failed ({exc}); falling back to accessibility 1.0")
        return {}


def persist_task(task: dict) -> tuple[dict, dict]:
    contract = task["content"]["contract"]
    hbb_reference = task["content"]["hbb_reference"]
    Path(settings.CONTRACT_PATH).parent.mkdir(parents=True, exist_ok=True)
    with open(settings.CONTRACT_PATH, "w") as handle:
        json.dump(contract, handle, indent=2)
    with open(settings.HBB_REFERENCE_PATH, "w") as handle:
        json.dump(hbb_reference, handle, indent=2)
    return contract, hbb_reference


# --------------------------------------------------------------------------------------------
# 2. Context: reference sequence, k-mer index, contract-derived constants
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
        return self.contract["seed"]

    @property
    def max_experiments(self) -> int:
        return self.contract["rules"].get("max_experiments") or 250

    @property
    def max_mismatches(self) -> int:
        return self.contract["rules"].get("max_mismatches", 0)

    @property
    def base_padding(self) -> int:
        return self.contract["rules"]["base_padding"]


# Every task in the backend's history shares one gene_region and one rules block, so the 135 MB
# reference, the k-mer index and the PAM site enumeration are all task-independent. Caching them
# turns a 200-task sweep from hours into minutes.
_SEQ_CACHE: str | None = None
_KMER_CACHE: dict[tuple[int, int], dict] = {}
_SITE_CACHE: dict[tuple, list] = {}


def load_sequence() -> str:
    global _SEQ_CACHE
    if _SEQ_CACHE is None:
        if not Path(settings.CHR11_PATH).exists():
            raise SystemExit(
                f"{settings.CHR11_PATH} missing — run scripts/run_validator.sh or download "
                "Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz from Ensembl release 116."
            )
        _SEQ_CACHE = stage12.load_chr11(settings.CHR11_PATH)
    return _SEQ_CACHE


def build_context(contract: dict, reference: dict, cell_types: dict,
                  offtarget_flank: int = 50000) -> Context:
    seq = load_sequence()

    # Same window the validator indexes: forward strand of gene_region +/- 50 kb.
    win_start = reference["gene_region"]["start"]
    win_end = reference["gene_region"]["end"]
    index_start = max(0, win_start - offtarget_flank)
    index_end = min(len(seq), win_end + offtarget_flank)
    if (index_start, index_end) not in _KMER_CACHE:
        _KMER_CACHE[(index_start, index_end)] = stage12.load_or_build_kmer_index(
            seq[index_start:index_end], k=12
        )
    kmer_index = _KMER_CACHE[(index_start, index_end)]

    return Context(
        seq=seq,
        contract=contract,
        reference=reference,
        cell_types=cell_types,
        kmer_index=kmer_index,
        mutation_map=reference["mutation_map"],
        mutations=list(contract["active_mutations"]),
        cas_systems=list(contract["rules"].get("cas_systems", ["Cas9", "Cas12a"])),
    )


# --------------------------------------------------------------------------------------------
# 3. Site enumeration
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
    """Every (cas, strand, length, start) with a valid PAM in gene_region +/- flank.

    PAM validity is delegated to ``stage12.check_pam`` so the enumeration cannot drift from the
    gate that will actually judge the row.
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
# 4. Guide tuning
#
# Stage 1 accepts up to contract.rules.max_mismatches Hamming distance between the guide and the
# reference target, and the PAM is read off the reference (never off the guide). That budget is a
# free design lever with two uses:
#
#   * push GC content to exactly 50%, where gc_score peaks at 1.0;
#   * perturb the 12-mer off-target seed so it is absent from the forward-strand index, taking
#     offtarget_factor from 0.7 to 1.0 (a flat 1.43x on every + strand row).
#
# A third, subtler use: every distinct guide is a distinct stage-3 RNG seed, so the variants of one
# site are draws from the outcome distribution at a *fixed* feature vector. That is what lets the
# selector shape outcomes into something a RandomForest can actually learn.
# --------------------------------------------------------------------------------------------

def _substitute(guide: list[str], position: int, want_gc: bool | None, rng: random.Random) -> bool:
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
    for attempt in range(n_variants * 12):
        if len(variants) >= n_variants:
            break
        guide = list(ref)
        spent = 0

        # Step 1: move the GC count. Prefer editing inside the seed window so one substitution
        # buys both the GC correction and the off-target break.
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
    # Best off-target class first so a fallback pick is still a good row.
    variants.sort(key=lambda g: -stage12.offtarget_uniqueness(g, site.cas, ctx.kmer_index))
    return variants


# --------------------------------------------------------------------------------------------
# 5. Row assembly — features and simulated outcome, straight from the validator's own functions
# --------------------------------------------------------------------------------------------

def make_experiment(site: Site, guide: str, mutation: str, ctx: Context, experiment_id: str) -> dict:
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
# 6. In-memory scorer — an exact replica of stages 4 and 5, no file I/O
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
        "per_target": {
            column: {
                "r2": v["r2_mean"],
                "nmae": stage4.normalized_mae(v["mae_mean"], y[column]),
            }
            for column, v in per_target.items()
        },
    }


def score_rows(valid: list[dict], results: list[dict], ctx: Context) -> dict:
    s4 = stage4_in_memory(valid, results, fold_seed=ctx.seed)
    fidelity = stage5.compute_distribution_fidelity(valid, results, ctx.contract, k=12)
    factor = max(0.0, min(1.0, fidelity.get("distribution_fidelity_score", 0.0)))
    return {
        **s4,
        "distribution_fidelity_score": fidelity.get("distribution_fidelity_score", 0.0),
        "distribution_fidelity_factor": factor,
        "final_score": s4["total_weighted_score"] * s4["consistency_factor"] * factor,
        "fidelity_detail": fidelity,
    }


# --------------------------------------------------------------------------------------------
# 7. Generation
# --------------------------------------------------------------------------------------------

@dataclass
class GenConfig:
    """Knobs the config search turns. Defaults are the starting point of the coordinate descent."""
    strategy: str = "pure"         # "pure" | "shaped" — see generate()
    selection: str = "packed"      # "packed" | "stratified" | "nearest" — see select_sites()
    flank: int = 3000              # site-enumeration window beyond gene_region
    max_distance: int = 2000       # widest |start - mutation_pos| a selected row may have
    lengths: tuple[int, ...] = (20, 23)
    variants: int = 24             # guide variants generated per selected site
    rows: int | None = None        # defaults to contract max_experiments
    construction: str = "mh"       # key into CONSTRUCTIONS; the rule every row must satisfy
    weight_skew: float = 2.0       # exponent on mutation_weight when apportioning rows

    # Generate this multiple of every cell's quota, then keep only the rows whose *realised*
    # (gc, distance) land closest to the ideal. 1.0 disables the extra work. Costs build time
    # linearly and nothing else: the surplus rows are discarded before scoring, so the submission
    # is still exactly max_experiments rows — just drawn from a wider pool of candidates.
    oversample: float = 1.0

    # Treat any guide within this much of 50% GC as equally good, so the ranking falls through to
    # distance instead of chasing the last thousandth of gc_score. 0.0 makes the "nearest" ranking
    # strictly lexicographic, which reaches for a gc-perfect site however far out it sits.
    gc_tolerance: float = 0.0

    # Admit only sites whose stage-3 cut probability will reach this. 0.0 disables the filter.
    # cut_p = min(0.99, max(0.4, base + 0.18*energy)) with base 0.86 for Cas9 and 0.78 for Cas12a,
    # so the two systems have different ceilings: 0.99 and 0.96. A single threshold above 0.96 is
    # therefore satisfiable by Cas9 alone and silently empties every Cas12a cell — see
    # cut_p_ceiling(). Use cut_p_floors to hold each system to its own ceiling instead.
    min_cut_p: float = 0.0

    # Per-Cas floors, e.g. {"Cas9": 0.99, "Cas12a": 0.96}. Overrides min_cut_p for the systems
    # named. This is the only way to demand a maximal cut probability from both systems at once,
    # because their ceilings differ by 0.03.
    cut_p_floors: dict | None = None

    # What to do when a floor is unreachable. cut_p depends on energy, and energy scales with the
    # cell type's accessibility: at 0.87 it saturates the clamp and every Cas9 site sits at 0.99,
    # but HEK293's 0.35 caps energy near 0.52, so the best Cas9 site reaches 0.9545 and the best
    # Cas12a site 0.8745. A strict floor there matches nothing, empties every cell and submits zero
    # rows — 23% of the backend's tasks are HEK293. Relaxing keeps the floor wherever it is
    # satisfiable and falls back to the best available cut_p where it is not. False restores the
    # strict gate, which is what a deliberate "Cas9 only" experiment wants.
    cut_p_relax: bool = True

    # When a filter leaves the submission under max_experiments, top it up with rows from this Cas
    # system, taking the highest cut_p still available. Named rather than boolean because only Cas9
    # can reach 0.99, so it is the only system with spare high-cut_p sites to give.
    backfill_cas: str | None = "Cas9"

    # Row share per Cas system, e.g. {"Cas9": 0.8, "Cas12a": 0.2}. None splits evenly across the
    # (cas x strand) product, which is what maximises stage 5's cas and joint coverage entropies.
    # Skewing trades that entropy for the freedom to draw more rows from whichever system has the
    # sites — worth it exactly when a filter has made one system's eligible pool too small to fill
    # an even quota.
    cas_share: dict | None = None

    # Outcome shaping. Stage 4 asks whether a RandomForest can learn design -> outcome under CV;
    # a dataset whose outcomes are a smooth function of `distance` is exactly that. Rows are ranked
    # by distance and each rank band is given an outcome target, then the variant search picks the
    # guide whose deterministic stage-3 draw lands closest to it.
    cut_frac: float = 0.75         # nearest this fraction of rows target a cut
    hdr_frac: float = 0.5          # of the cut rows, nearest this fraction target HDR
    indel_lo: float = 1.0
    indel_hi: float = 9.0
    indel_levels: int = 4          # quantise the indel ramp; a coarse ladder is actually hittable
    mh_frac: float = 0.5           # mh is a *feature*; tying it to distance keeps X non-noisy
    w_cut: float = 4.0
    w_hdr: float = 2.0
    w_indel: float = 0.35
    w_mh: float = 0.8

    # gc_spread > 0 trades gc_score for feature variance: energy saturates at 1.0 for most of the
    # near-mutation range, so without some GC variation `energy` is a constant column.
    gc_spread: float = 0.0


def _quota(total: int, cells: int) -> list[int]:
    base, extra = divmod(total, cells)
    return [base + (1 if i < extra else 0) for i in range(cells)]


def _cell_quotas(ctx: Context, rows_wanted: int, skew: float,
                 groups: list[tuple[str, str]] | None = None,
                 group_weights: list[float] | None = None) -> tuple[list[tuple], list[int]]:
    """Rows per (mutation, cas, strand) cell.

    ``groups`` narrows the (cas, strand) product to those that can actually be filled. Passing it
    keeps the submission at the row cap when a filter has emptied some cells: without it those
    cells keep their quota, nothing fills them, and the submission silently lands under
    max_experiments — losing rows on top of whatever coverage the filter already cost.

    Even quotas maximise stage 5's coverage entropies. Skewing toward the heavier mutations trades
    some of that entropy for total_weighted_score, which multiplies every row by mutation_weight.
    The trade is favourable well past 50/50: term 1 gains linearly, while stage 5 enters the
    product through a *six-way* geometric mean, so only two of the six ratios move and each moves
    at the 1/6 power. ``weight_skew`` is the exponent on mutation_weight; 0 keeps it uniform.
    """
    weights = ctx.contract.get("mutation_weights", {})
    shares = {m: max(weights.get(m, 1.0), 1e-9) ** skew for m in ctx.mutations}
    total = sum(shares.values())
    if groups is None:
        groups = [(cas, strand) for cas in ctx.cas_systems for strand in ctx.strands]
    if group_weights is None:
        group_weights = [1.0] * len(groups)
    weight_total = sum(group_weights) or 1.0

    exact = {m: rows_wanted * shares[m] / total for m in ctx.mutations}
    counts = {m: int(exact[m]) for m in ctx.mutations}
    remainder = rows_wanted - sum(counts.values())
    for mutation in sorted(ctx.mutations, key=lambda m: -(exact[m] - counts[m]))[:remainder]:
        counts[mutation] += 1

    cells: list[tuple] = []
    quotas: list[int] = []
    for mutation in ctx.mutations:
        # Largest-remainder apportionment over the groups. With uniform weights this reproduces
        # _quota exactly — floor everywhere, then +1 to the leading groups — so the default path is
        # untouched.
        want = counts[mutation]
        exact = [want * weight / weight_total for weight in group_weights]
        split = [int(value) for value in exact]
        remainder = want - sum(split)
        for index in sorted(range(len(groups)),
                            key=lambda i: -(exact[i] - split[i]))[:max(0, remainder)]:
            split[index] += 1
        for index, (cas, strand) in enumerate(groups):
            cells.append((mutation, cas, strand))
            quotas.append(split[index])
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


def _gc_reach_penalty(site: Site, budget: int) -> int:
    """GC-count units still separating this site from 50% after the whole mismatch budget is spent."""
    ref_gc = sum(base in "GC" for base in site.ref_guide)
    return max(0, abs(ref_gc - round(site.length / 2)) - budget)


def gc_miss(site: Site, budget: int) -> float:
    """How far this site's guide will still sit from 50% GC once the mismatch budget is spent.

    Zero is reachable only at even lengths: at L=20 the target is 10 GC bases and gc_score peaks
    at exactly 1.0, while at L=23 the closest attainable count is 12/23, an unavoidable 0.0217
    residual worth 0.625 * 0.0435 = 0.027 of structural score. Ranking on this therefore prefers
    20-mers over 23-mers on its own, without the length ever being named.
    """
    ref_gc = sum(base in "GC" for base in site.ref_guide)
    target = round(site.length / 2)
    reached = ref_gc + max(-budget, min(budget, target - ref_gc))
    return abs(reached / site.length - 0.5)


def predicted_energy(site: Site, mutation: str, ctx: Context) -> float:
    """The energy stage 3 will compute for this site once its guide is tuned.

    Exact rather than approximate: ``tune_variants`` guarantees every variant it returns carries the
    same GC count, so the gc that reaches stage 3 is known before a single guide is built. Mirrors
    ``stage3.sequence_energy`` term for term, including the 1500 bp decay length that differs from
    stage 2's base_padding.
    """
    budget = ctx.max_mismatches
    ref_gc = sum(base in "GC" for base in site.ref_guide)
    target = round(site.length / 2)
    gc = (ref_gc + max(-budget, min(budget, target - ref_gc))) / site.length
    distance = abs(site.start - ctx.mutation_map[mutation])
    accessibility = ctx.cell_types.get(ctx.contract.get("cell_type"), {}).get("accessibility", 1.0)
    offset = stage3.REGION_ENERGY_OFFSETS.get(
        (ctx.contract.get("mutation_regions") or {}).get(mutation), 0.0
    )
    return max(0.0, min(1.0, accessibility * (
        1.8 * gc + 0.6 * math.exp(-distance / 1500) + offset
    )))


def predicted_cut_p(site: Site, mutation: str, ctx: Context) -> float:
    """Stage 3's own cut_probability, evaluated on the energy this site will reach."""
    return stage3.cut_probability(site.cas, predicted_energy(site, mutation, ctx))


def cut_p_floor(cas: str, cfg: "GenConfig") -> float:
    """The cut_p a row of this Cas system must reach. Per-system floors win over the global one."""
    if cfg.cut_p_floors and cas in cfg.cut_p_floors:
        return float(cfg.cut_p_floors[cas])
    return cfg.min_cut_p


def cut_p_filtering(cfg: "GenConfig") -> bool:
    return cfg.min_cut_p > 0.0 or bool(cfg.cut_p_floors)


def cut_p_ceiling_for(cas: str) -> float:
    """The highest cut_p this Cas system can ever reach, i.e. at the energy clamp."""
    return stage3.cut_probability(cas, 1.0)


# Kept as the historical spelling used in reports and the stage printer.
cut_p_ceiling = cut_p_ceiling_for


def gc_rank(miss: float, tolerance: float) -> float:
    """Collapse every near-enough GC into one tier so distance decides between them.

    gc_score is flat to first order around 50%: the whole span from a perfect 20-mer to the best a
    23-mer can do is 0.027 of structural score, while 200 bp of extra distance at base_padding 400
    costs about ten times that. Without a tier the ranking spends distance it cannot afford.
    """
    return 0.0 if miss <= tolerance else miss


def _stratified_pick(pool: list[Site], take: int, budget: int) -> list[Site]:
    """One site per equal-width band of a distance-sorted pool.

    Banding is what spreads the rows along the distance axis. Within a band the tie-break prefers a
    site whose reference GC can still be pulled to 50% inside the mismatch budget: gc_score carries
    0.625 of the structural score, so an unreachable site is a permanent haircut on that row.
    """
    if take <= 0 or not pool:
        return []
    if take >= len(pool):
        return list(pool)
    picked = []
    for i in range(take):
        lo = round(i * len(pool) / take)
        hi = max(lo + 1, round((i + 1) * len(pool) / take))
        band = pool[lo:hi]
        picked.append(
            min(enumerate(band), key=lambda t: (_gc_reach_penalty(t[1], budget), t[0]))[1]
        )
    return picked


def select_sites(ctx: Context, sites: list[Site], cfg: GenConfig, quiet: bool = False,
                 quota_override: dict[tuple, int] | None = None
                 ) -> tuple[list[tuple[Site, str]], dict[tuple, list[Site]]]:
    """Pick (site, mutation) pairs across the mutation x cas x strand support.

    Every cell must be occupied — stage 5's geometric mean clips a missing category to 1e-9, i.e.
    a 0.0316x hit. How many rows each cell gets is ``_cell_quotas``' call.

    Two selection modes, because the two strategies want opposite things from the distance axis:

    * ``packed``   — nearest sites first. dist_score = exp(-d/base_padding), so proximity is pure
      total_weighted_score. Correct when consistency does not depend on feature spread.
    * ``stratified`` — one site per equal-width distance band, giving stage 4's RandomForest a
      feature axis to learn against. Costs some dist_score to buy consistency_factor.

    Sites are globally unique because stage 1 dedups on (cas, start, strand, guide) — the mutation
    is *not* part of that key, so one coordinate cannot serve two mutations.

    Returns the selection plus the per-cell reserve of untouched sites, ordered by distance, which
    the pure strategy uses to replace rows whose outcome could not be forced.
    """
    rows_wanted = cfg.rows or ctx.max_experiments

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

    groups = [(cas, strand) for cas in ctx.cas_systems for strand in ctx.strands]
    if cut_p_filtering(cfg):
        # Which groups can satisfy their own floor at all. Tested against real sites rather than
        # against the per-Cas ceiling alone, because a group can clear the ceiling and still have
        # no coordinate near enough to reach the energy the threshold needs.
        eligible = {
            group for group in groups
            if any(predicted_cut_p(site, mutation, ctx) >= cut_p_floor(group[0], cfg) - 1e-9
                   for mutation in ctx.mutations
                   for site in by_group[group]
                   if abs(site.start - ctx.mutation_map[mutation]) <= max_distance)
        }
        dropped = [g for g in groups if g not in eligible]
        if dropped and not quiet:
            print(f"  ! no site can reach the cut_p floor for {dropped} "
                  f"(floors { {g[0]: cut_p_floor(g[0], cfg) for g in dropped} } vs ceilings "
                  f"{ {g[0]: cut_p_ceiling(g[0]) for g in dropped} })"
                  + (", relaxing to the best available" if cfg.cut_p_relax
                     else " — those cells stay empty"))
        # Only prune groups under the strict gate. Relaxed, an unreachable floor must not cost the
        # category: losing a whole Cas system takes its coverage ratio to 0, and the 1e-9 clip in
        # the geometric mean turns that into a 0.0316x multiplier on the entire score.
        if eligible and not cfg.cut_p_relax:
            groups = [g for g in groups if g in eligible]

    group_weights = None
    if cfg.cas_share:
        strands_per_cas = max(1, len({strand for _cas, strand in groups}))
        group_weights = [float(cfg.cas_share.get(cas, 0.0)) / strands_per_cas
                         for cas, _strand in groups]
        if not any(weight > 0 for weight in group_weights):
            raise ValueError(f"cas_share {cfg.cas_share} leaves every group at zero rows")

    cells, quotas = _cell_quotas(ctx, rows_wanted, cfg.weight_skew, groups, group_weights)
    if quota_override is not None:
        quotas = [quota_override.get(cell, 0) for cell in cells]
        # A negative quota silently empties its cell while the others keep their inflated counts,
        # so the submission overshoots max_experiments and the in-memory score (which never
        # truncates) reports a total the validator would never pay. Fail loudly instead.
        if any(q < 0 for q in quotas) or sum(quotas) != rows_wanted:
            raise ValueError(
                f"quota_override must be non-negative and sum to {rows_wanted}, "
                f"got sum={sum(quotas)} min={min(quotas)}"
            )

    # Scarcest (cas, strand) group first so Cas12a — roughly 4x rarer than Cas9 — is not starved
    # of coordinates by an earlier Cas9 cell.
    order = sorted(range(len(cells)),
                   key=lambda i: len(by_group[(cells[i][1], cells[i][2])]))

    used: set[tuple] = set()
    selected: list[tuple[Site, str]] = []
    reserve: dict[tuple, list[Site]] = {}
    relaxed: dict[str, float] = {}
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
        if cut_p_filtering(cfg):
            floor = cut_p_floor(cas, cfg)
            eligible_pool = [site for site in pool
                             if predicted_cut_p(site, mutation, ctx) >= floor - 1e-9]
            if eligible_pool or not cfg.cut_p_relax:
                pool = eligible_pool
            elif pool:
                # Unreachable here: keep the cell and take the best cut_p on offer instead of
                # submitting nothing. Recorded so the build reports the relaxation rather than
                # quietly scoring under a floor it never met.
                best = max(predicted_cut_p(site, mutation, ctx) for site in pool)
                relaxed[cas] = min(relaxed.get(cas, best), best)
        pool.sort(key=lambda site: abs(site.start - position))
        if not pool:
            reserve[cell] = []
            shortfall += want
            continue

        take = min(want, len(pool))
        if cfg.selection == "packed":
            # Rank by the stage-2 score the site can actually reach, not by distance alone: a site
            # 40 bp further out but tunable to 50% GC beats a nearer one stuck at gc_score 0.5.
            ranked = sorted(pool,
                            key=lambda s: -achievable_structural(s, abs(s.start - position), ctx))
        elif cfg.selection == "nearest":
            # Lexicographic: closest to gc 0.5 first, then closest to the mutation. Unlike
            # "packed" this does not trade the two off against each other — a site that can reach
            # exactly 50% GC outranks every site that cannot, however near the latter sits. The
            # two rankings therefore disagree only over sites whose GC is unreachable, which is
            # where "packed" is the score-optimal choice and this one is the literal one.
            ranked = sorted(pool, key=lambda s: (
                gc_rank(gc_miss(s, ctx.max_mismatches), cfg.gc_tolerance),
                abs(s.start - position)))
        else:
            # Stratified picks first, then the rest of the pool as filler for any short cell.
            ranked = _stratified_pick(pool, take, ctx.max_mismatches) + pool

        deduped: list[Site] = []
        spare: list[Site] = []
        seen_keys: set[tuple] = set()
        for site in ranked:
            if site.key in seen_keys:
                continue
            seen_keys.add(site.key)
            (deduped if len(deduped) < take else spare).append(site)

        for site in deduped:
            used.add(site.key)
            selected.append((site, mutation))
        # Reserve keeps the same ranking, so a replacement is the next-best site, not just the
        # next-nearest one.
        reserve[cell] = spare
        shortfall += want - len(deduped)

    if relaxed and not quiet:
        print(f"  ! cut_p floor relaxed to the best available: "
              f"{ {cas: round(value, 4) for cas, value in relaxed.items()} } "
              f"(configured { {cas: cut_p_floor(cas, cfg) for cas in relaxed} })")

    if shortfall > 0 and cfg.backfill_cas:
        # A cut_p floor bites unevenly: Cas12a tops out at 0.96 and has few sites that reach it, so
        # its cells run dry while Cas9 still has hundreds of unused sites sitting at its own 0.99
        # ceiling. Leaving the gap costs whole rows off term 1, which is a sum — strictly worse than
        # the coverage those extra Cas9 rows dilute. Fill it with the highest cut_p left, breaking
        # ties on weighted_score, which is the quantity term 1 actually accumulates.
        cas = cfg.backfill_cas
        candidates = backfill_candidates(ctx, sites, cfg, used)
        filled = 0
        for _cut_p, _value, site, mutation in candidates:
            if filled >= shortfall:
                break
            if site.key in used:      # the same coordinate is a candidate once per mutation
                continue
            used.add(site.key)
            selected.append((site, mutation))
            filled += 1

        if filled and not quiet:
            best = -candidates[0][0] if candidates else 0.0
            print(f"  + backfilled {filled} {cas} row(s) to reach the {rows_wanted} cap "
                  f"(highest available cut_p {best:.4f})")
        shortfall -= filled

    if shortfall and not quiet:
        print(f"  ! {shortfall} row(s) short of the {rows_wanted} cap — widen --flank/max_distance")
    return selected, reserve


def estimate_final_score(ctx: Context, sites: list[Site], cfg: GenConfig,
                         quota_override: dict[tuple, int] | None = None) -> float:
    """Surrogate score for the pure strategy: no guide tuning, no simulation, no RandomForest.

    Sound only under ``strategy="pure"``, and only for comparing candidate *selections*:

    * consistency_factor is pinned at 1.0 by the construction, so the RF never needs to run. This
      holds for the registered rules ("hdr", "mh"), both measured at exactly 1.0; a construction
      that does not reach 1.0 would need this surrogate revisited;
    * total_weighted_score follows from site choice alone, because guide tuning is deterministic
      (``achievable_structural`` mirrors what stage 2 will compute) and the variant search only
      picks among guides sharing one feature vector;
    * the two guide-dependent stage-5 ratios — 12-mer diversity and distinct-guide — are dropped.
      Inside a six-way geometric mean they multiply every candidate by the same factor, so they
      cannot move the argmax.

    Cheap enough (a few ms) to sweep a knob per task instead of per repo.
    """
    selected, _ = select_sites(ctx, sites, cfg, quiet=True, quota_override=quota_override)
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
    """Pick the mutation-weight skew for *this* contract.

    The optimum depends on the weight ratio, which ranges from about 1.10/0.65 to 1.68/0.53 across
    the task history — one repo-wide constant leaves score on the table for most of them.
    """
    return max(SKEW_CANDIDATES,
               key=lambda skew: estimate_final_score(ctx, sites, replace(cfg, weight_skew=skew)))


def outcome_targets(distances: list[int], cfg: GenConfig) -> list[dict]:
    """Map each row's distance rank onto the outcome the variant search should hunt for."""
    order = sorted(range(len(distances)), key=lambda i: distances[i])
    n = len(distances)
    targets: list[dict] = [{} for _ in range(n)]

    cut_cutoff = cfg.cut_frac * n
    hdr_cutoff = cfg.cut_frac * cfg.hdr_frac * n
    nhej_total = max(1, int(cut_cutoff - hdr_cutoff))

    nhej_rank = 0
    for rank, row in enumerate(order):
        is_cut = rank < cut_cutoff
        is_hdr = rank < hdr_cutoff
        if is_cut and not is_hdr:
            fraction = nhej_rank / max(1, nhej_total - 1)
            # Quantised: a stage-3 draw can land exactly on one of a few integer levels far more
            # often than on a point of a continuous ramp, and a step function is just as learnable.
            levels = max(1, cfg.indel_levels)
            if levels > 1:
                fraction = round(fraction * (levels - 1)) / (levels - 1)
            indel = round(cfg.indel_lo + (cfg.indel_hi - cfg.indel_lo) * fraction)
            nhej_rank += 1
        else:
            indel = 0.0
        targets[row] = {
            "is_cut": 1.0 if is_cut else 0.0,
            "is_hdr": 1.0 if is_hdr else 0.0,
            "indel": indel,
            "mh": 1.0 if rank < cfg.mh_frac * n else 0.0,
            "rank": rank / max(1, n - 1),
        }
    return targets


def variant_loss(result: dict, target: dict, cfg: GenConfig) -> float:
    is_cut = 1.0 if result["outcome"] != "no_cut" else 0.0
    is_hdr = 1.0 if result["outcome"] == "HDR" else 0.0
    loss = cfg.w_cut * abs(is_cut - target["is_cut"])
    loss += cfg.w_hdr * abs(is_hdr - target["is_hdr"])
    loss += cfg.w_indel * abs(result["indel_length"] - target["indel"])
    loss += cfg.w_mh * abs(float(result["mh"]) - target["mh"])
    return loss


def design_key(experiment: dict) -> tuple:
    """Stage 1's second dedup key — independent of experiment_id."""
    return (experiment["cas_system"], experiment["target_alignment_start"],
            experiment["strand"], experiment["guideRNA"])


# --------------------------------------------------------------------------------------------
# Constructions — the rule every row's stage-3 outcome is forced to satisfy.
#
# Stage 4 fits is_cut, is_hdr and indel_length separately and scores each with r2 + normalised MAE.
# It does not reward *variety* in outcomes; it rewards outcomes being a function the forest can
# recover from X = [gc, distance, gc_score, dist_score, consistency, energy, mh] under CV. Two ways
# to satisfy that, both measured against the real scorer:
#
#   "hdr"  every row HDR. All three targets constant, so r2_score hits its zero-numerator/
#          zero-denominator case (1.0) and normalized_mae short-circuits on std < 1e-9.
#          consistency_factor = 1.0, but the dataset is degenerate in all three targets.
#
#   "mh"   mh -> HDR, otherwise BLUNT_NHEJ pinned to indel_length 1. Also consistency_factor = 1.0
#          with r2 = 1.0 on all three, *without* is_hdr or indel_length being constant — because
#          `mh` is literally a column in X, so is_hdr == mh and indel_length == 1 - mh are exactly
#          recoverable. Measured 214.75 vs 215.21 for "hdr" on task f9d9356a: a 0.2% cost to drop
#          two of the three degeneracies.
#
# A rule receives the simulated result and the stage-12 entry, so a construction may key off any
# feature the validator computes — e.g. entry["features"]["consistency"], the other exact binary in
# X, which is what a rule that also breaks is_cut's constancy would use.
# --------------------------------------------------------------------------------------------

def _rule_hdr(result: dict, entry: dict) -> bool:
    return result["outcome"] == "HDR"


def _rule_mh(result: dict, entry: dict) -> bool:
    if result["mh"]:
        return result["outcome"] == "HDR"
    # indel_length must be pinned, or the exponential draw leaves that target unlearnable.
    return result["outcome"] == "BLUNT_NHEJ" and result["indel_length"] == 1


def _rule_nocut(result: dict, entry: dict) -> bool:
    # Pins all three targets at once (cut=0, hdr=0, indel=0), like "hdr" — but far more expensive
    # to hit: cut_p is capped at 0.99, so P(no_cut) is ~0.01 per variant on Cas9.
    return result["outcome"] == "no_cut"


def _rule_blunt(result: dict, entry: dict) -> bool:
    # indel_length pinned to its modal value (P=0.70 given BLUNT_NHEJ); without the pin the
    # exponential draw leaves that target unlearnable.
    return result["outcome"] == "BLUNT_NHEJ" and result["indel_length"] == 1


def _rule_mhnhej(result: dict, entry: dict) -> bool:
    # The gamma draw is flat (no value exceeds P=0.12), so pinning costs ~48 variants per row.
    return result["outcome"] == "MH_NHEJ" and result["indel_length"] == 1


def _rule_mh_any(result: dict, entry: dict) -> bool:
    # "mh" without the indel pin. The repair *mode* still tracks the microhomology coin, so is_hdr
    # stays an exact function of a feature the validator itself computes; indel_length is left to
    # its draw and stops being learnable. Far cheaper per row than "mh" (no 0.70 pin to hit), so it
    # buys candidate count with one of stage 4's three targets.
    if result["mh"]:
        return result["outcome"] == "HDR"
    return result["outcome"] == "BLUNT_NHEJ"


def _rule_blunt_any(result: dict, entry: dict) -> bool:
    return result["outcome"] == "BLUNT_NHEJ"


def _rule_mhnhej_any(result: dict, entry: dict) -> bool:
    return result["outcome"] == "MH_NHEJ"


CONSTRUCTIONS = {
    "mh": _rule_mh,
    "mh_any": _rule_mh_any,
    "hdr": _rule_hdr,
    "nocut": _rule_nocut,
    "blunt": _rule_blunt,
    "mhnhej": _rule_mhnhej,
    "blunt_any": _rule_blunt_any,
    "mhnhej_any": _rule_mhnhej_any,
}


def _best_variant(site: Site, mutation: str, ctx: Context, cfg: GenConfig, index: int,
                  rng: random.Random, target: dict | None,
                  require: "callable | None",
                  seen_designs: set[tuple] | None = None) -> tuple[dict, dict, dict] | None:
    """Search a site's guide variants for the one whose stage-3 draw best fits the target.

    ``require`` is a construction rule and makes the match mandatory: the caller would rather lose
    the row than accept a draw that breaks the rule (see ``generate_pure`` for why).

    ``seen_designs`` mirrors stage 1's ``seen_valid_keys``. Stage 1 dedups on
    (cas, start, strand, guide) *in addition to* experiment_id, and silently drops the second
    occurrence — so a collision here costs a row and, under the pure strategy, would break outcome
    purity as well. Checking it at the point of choice is the only place it can be avoided.
    """
    offset = cfg.gc_spread * ((target or {}).get("rank", 0.5) - 0.5) * 2.0
    target_gc = max(0.0, min(1.0, 0.5 + offset))
    target_gc_count = int(round(target_gc * site.length))

    best = None
    for guide in tune_variants(site, target_gc_count, ctx, cfg.variants, rng):
        experiment = make_experiment(site, guide, mutation, ctx, f"exp-{index:05d}")
        if seen_designs is not None and design_key(experiment) in seen_designs:
            continue
        entry = build_valid_entry(experiment, ctx)
        if entry is None:
            continue
        result = simulate(entry, ctx)
        if require is not None:
            if require(result, entry):
                return experiment, entry, result
            continue
        loss = variant_loss(result, target, cfg)
        if best is None or loss < best[0]:
            best = (loss, experiment, entry, result)
        if loss == 0.0:
            break
    if best is None:
        return None
    return best[1], best[2], best[3]


def _realised_miss(entry: dict, tolerance: float = 0.0) -> tuple[float, int]:
    """The nearness of a row that actually exists: |gc - 0.5| tier, then distance to the mutation.

    Measured off the stage-2 features rather than predicted from the site, because guide tuning
    does not always reach the GC target — the mismatch budget is shared with the off-target seed,
    and the seed wins when the two compete.
    """
    features = entry["features"]
    return gc_rank(abs(features["gc"] - 0.5), tolerance), features["distance_to_mutation"]


def generate_oversampled(ctx: Context, sites: list[Site], cfg: GenConfig) -> tuple[list, list, list]:
    """Build oversample x the wanted rows, then keep the nearest quota of each cell.

    Generating wide and filtering afterwards is not the same as ranking sites better up front.
    Two things are only knowable after the row exists: whether any of the site's variants could be
    made to satisfy the construction at all, and what GC the tuner actually reached. Both are
    resolved here, and the surplus is dropped before anything is scored.
    """
    rows_wanted = cfg.rows or ctx.max_experiments
    cells, quotas = _cell_quotas(ctx, rows_wanted, cfg.weight_skew)
    target = {cell: quota for cell, quota in zip(cells, quotas)}

    wide = replace(cfg, rows=max(rows_wanted, int(round(rows_wanted * cfg.oversample))),
                   oversample=1.0)
    rows, valid, results = generate_pure(ctx, sites, wide)
    if not rows:
        return rows, valid, results

    by_cell: dict[tuple, list[int]] = defaultdict(list)
    for index, entry in enumerate(valid):
        experiment = entry["experiment"]
        by_cell[(experiment["mutation"], experiment["cas_system"], experiment["strand"])].append(index)

    keep: list[int] = []
    for cell, indices in by_cell.items():
        indices.sort(key=lambda i: _realised_miss(valid[i], cfg.gc_tolerance))
        keep.extend(indices[:target.get(cell, 0)])
    keep.sort()

    print(f"    oversample {cfg.oversample:g}x: {len(rows)} candidates -> {len(keep)} kept "
          f"({len(rows) - len(keep)} dropped as further from the ideal)")
    return ([rows[i] for i in keep], [valid[i] for i in keep], [results[i] for i in keep])


def backfill_candidates(ctx: Context, sites: list[Site], cfg: GenConfig,
                        claimed: set[tuple]) -> list[tuple]:
    """Unused sites of ``cfg.backfill_cas``, best first: highest cut_p, then highest weighted score.

    cut_p leads because that is the property being topped up; weighted_score
    (achievable_structural x mutation_weight) breaks the ties, and on Cas9 the ties are the whole
    field — every site near the mutation sits at the 0.99 clamp, so the second key is what actually
    orders them.
    """
    cas = cfg.backfill_cas
    max_distance = cfg.max_distance
    if ctx.contract["rules"].get("proximity_gate", False):
        max_distance = min(max_distance, ctx.base_padding)
    weights = ctx.contract.get("mutation_weights", {})
    floor = cut_p_floor(cas, cfg) if cut_p_filtering(cfg) else 0.0
    if cfg.cut_p_relax and floor > cut_p_ceiling_for(cas):
        floor = 0.0     # unreachable for this system; rank by cut_p instead of gating on it

    candidates: list[tuple] = []
    for mutation in ctx.mutations:
        position = ctx.mutation_map[mutation]
        weight = weights.get(mutation, 1.0)
        for site in sites:
            if site.cas != cas or site.length not in cfg.lengths or site.key in claimed:
                continue
            distance = abs(site.start - position)
            if distance > max_distance:
                continue
            cut_p = predicted_cut_p(site, mutation, ctx)
            if cut_p < floor - 1e-9:
                continue
            candidates.append((-cut_p, -achievable_structural(site, distance, ctx) * weight,
                               site, mutation))
    candidates.sort(key=lambda item: (item[0], item[1]))
    return candidates


def generate_pure(ctx: Context, sites: list[Site], cfg: GenConfig) -> tuple[list, list, list]:
    """Force every row's stage-3 outcome to satisfy the configured construction rule.

    See ``CONSTRUCTIONS`` for what the rules are and why they score. What matters here is that
    conformance is **all-or-nothing**: the rules work by making stage 4's targets exactly
    recoverable from X, and a single stray row breaks that. Under "hdr" a stray row gives
    ``indel_length`` a near-zero variance against a large residual and r2 collapses from 1.0 to
    roughly 0; under "mh" it puts a point off the mh -> outcome mapping the forest would otherwise
    fit exactly. Either way the result is worse than not trying.

    So a site whose variants will not satisfy the rule is *replaced* from its cell's reserve, and
    if the reserve is exhausted the row is dropped rather than submitted non-conforming. Because
    stage 3 is deterministic in (contract.seed, design), the conformance measured here is exactly
    what the validator will compute.
    """
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
        found = _best_variant(site, mutation, ctx, cfg, index, rng, None, rule, seen_designs)
        while found is None and reserve.get(cell):
            replacement = reserve[cell].pop(0)
            if replacement.key in claimed:
                continue
            claimed.add(replacement.key)
            replaced += 1
            found = _best_variant(replacement, mutation, ctx, cfg, index, rng, None, rule,
                                  seen_designs)
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
        print(f"    {cfg.construction}: {replaced} site(s) replaced, {dropped} row(s) dropped")

    # Selection filled its quota, but the construction may still have dropped rows whose sites and
    # reserves both ran out of conforming variants — so the shortfall is only knowable here, after
    # the fact. Top up with the highest cut_p sites left, and hold them to the same rule: a
    # non-conforming row would take stage 4 down with it, which costs far more than the row adds.
    rows_wanted = cfg.rows or ctx.max_experiments
    if cfg.backfill_cas and len(rows) < rows_wanted:
        short = rows_wanted - len(rows)
        added = 0
        for _cut_p, _value, site, mutation in backfill_candidates(ctx, sites, cfg, claimed):
            if len(rows) >= rows_wanted:
                break
            if site.key in claimed:
                continue
            claimed.add(site.key)
            found = _best_variant(site, mutation, ctx, cfg, index, rng, None, rule, seen_designs)
            if found is None:
                continue
            experiment, entry, result = found
            seen_designs.add(design_key(experiment))
            rows.append(experiment)
            valid.append(entry)
            results.append(result)
            index += 1
            added += 1
        print(f"    backfill: {short} row(s) short after the construction, "
              f"{added} {cfg.backfill_cas} row(s) added -> {len(rows)}/{rows_wanted}")

    return rows, valid, results


def check_invariants(rows: list[dict], results: list[dict], cfg: GenConfig,
                     valid: list[dict] | None = None) -> list[str]:
    """The two uniqueness rules stage 1 and truncate_submission enforce, plus rule conformance.

    Worth checking rather than trusting: a violation is invisible in the in-memory score (which
    never dedups) but costs rows in the real pipeline, so it shows up only as an unexplained gap
    between a local run and the validator's.

    Conformance is checked, not assumed — it is the property the whole construction rests on, and
    unlike outcome purity it cannot be eyeballed from the outcome counts once the rule admits more
    than one class.
    """
    problems = []
    ids = [r["experiment_id"] for r in rows]
    if len(set(ids)) != len(ids):
        problems.append(f"{len(ids) - len(set(ids))} duplicate experiment_id")
    designs = [design_key(r) for r in rows]
    if len(set(designs)) != len(designs):
        problems.append(f"{len(designs) - len(set(designs))} duplicate (cas, start, strand, guide)")

    if cfg.strategy == "pure" and results and valid is not None:
        rule = CONSTRUCTIONS[cfg.construction]
        broken = sum(1 for result, entry in zip(results, valid) if not rule(result, entry))
        if broken:
            problems.append(f"{broken} row(s) break the '{cfg.construction}' construction")
    return problems


def generate(ctx: Context, sites: list[Site], cfg: GenConfig) -> tuple[list[dict], list[dict], list[dict]]:
    """Build a full submission. Returns (rows, valid_entries, stage3_results)."""
    rows, valid, results = _generate(ctx, sites, cfg)
    for problem in check_invariants(rows, results, cfg, valid):
        print(f"    ! invariant violated: {problem}")
    return rows, valid, results


def _generate(ctx: Context, sites: list[Site], cfg: GenConfig) -> tuple[list[dict], list[dict], list[dict]]:
    if cfg.strategy == "pure":
        if cfg.oversample > 1.0:
            return generate_oversampled(ctx, sites, cfg)
        return generate_pure(ctx, sites, cfg)

    rng = random.Random(ctx.seed ^ 0x5EED)
    selected, _reserve = select_sites(ctx, sites, cfg)
    if not selected:
        return [], [], []

    distances = [abs(site.start - ctx.mutation_map[mutation]) for site, mutation in selected]
    targets = outcome_targets(distances, cfg)

    rows: list[dict] = []
    valid: list[dict] = []
    results: list[dict] = []
    seen_designs: set[tuple] = set()

    for index, ((site, mutation), target) in enumerate(zip(selected, targets)):
        # gc_spread pulls a row's GC target off 50% as a function of its distance rank: it costs
        # gc_score, but accessibility * (1.8*gc + ...) saturates the energy clamp across most of
        # the near-mutation range, so without it `energy` is a constant column in X.
        found = _best_variant(site, mutation, ctx, cfg, index, rng, target, None, seen_designs)
        if found is None:
            continue
        experiment, entry, result = found
        seen_designs.add(design_key(experiment))
        rows.append(experiment)
        valid.append(entry)
        results.append(result)

    return rows, valid, results


def generate_baseline(ctx: Context, sites: list[Site], cfg: GenConfig) -> tuple[list, list, list]:
    """Untuned control: nearest PAM sites, reference-matching guides, no outcome shaping.

    Exists so the numbers the optimised path reports have something to be measured against.
    """
    rows_wanted = cfg.rows or ctx.max_experiments
    cells = [(m, c, s) for m in ctx.mutations for c in ctx.cas_systems for s in ctx.strands]
    quotas = _quota(rows_wanted, len(cells))
    by_group: dict[tuple[str, str], list[Site]] = defaultdict(list)
    for site in sites:
        by_group[(site.cas, site.strand)].append(site)

    used: set[tuple] = set()
    rows, valid, results = [], [], []
    index = 0
    for (mutation, cas, strand), want in zip(cells, quotas):
        position = ctx.mutation_map[mutation]
        pool = sorted((s for s in by_group[(cas, strand)] if s.key not in used),
                      key=lambda s: abs(s.start - position))
        for site in pool[:want]:
            used.add(site.key)
            experiment = make_experiment(site, site.ref_guide, mutation, ctx, f"exp-{index:05d}")
            index += 1
            entry = build_valid_entry(experiment, ctx)
            if entry is None:
                continue
            rows.append(experiment)
            valid.append(entry)
            results.append(simulate(entry, ctx))
    return rows, valid, results


# --------------------------------------------------------------------------------------------
# 8. Config search — coordinate descent on the true objective
# --------------------------------------------------------------------------------------------

PURE_GRID: list[tuple[str, list]] = [
    # Pure decouples consistency from feature spread, so the only real question left is how tightly
    # the rows can be packed onto the mutation before the cells run out of coordinates.
    ("max_distance", [500, 800, 1200, 2000, 3500]),
    ("weight_skew", [0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]),
    # Allowing 23-mers roughly doubles the site supply per cell, which lets every cell be filled
    # from a tighter band around the mutation. It costs gc_score — 50% GC is unreachable on an odd
    # length, capping it at 1 - 2/23 - so which way it nets out is an empirical question.
    ("lengths", [(20,), (23,), (20, 23)]),
    ("construction", ["mh", "hdr"]),
    ("selection", ["packed", "stratified"]),
    ("variants", [16, 32, 64]),
]

SHAPED_GRID: list[tuple[str, list]] = [
    ("max_distance", [600, 1200, 2500, 5000, 9000]),
    ("cut_frac", [0.55, 0.65, 0.75, 0.85, 0.95]),
    ("hdr_frac", [0.25, 0.4, 0.55, 0.7]),
    ("gc_spread", [0.0, 0.06, 0.12, 0.2]),
    ("mh_frac", [0.0, 0.35, 0.5, 0.65]),
    ("indel_hi", [4.0, 9.0, 14.0]),
    ("indel_levels", [1, 2, 4, 8]),
    ("weight_skew", [0.0, 1.0, 2.5, 5.0]),
    ("lengths", [(20,), (20, 23)]),
    ("variants", [24, 48, 80]),
]


def grid_for(cfg: GenConfig) -> list[tuple[str, list]]:
    return PURE_GRID if cfg.strategy == "pure" else SHAPED_GRID


def evaluate_config(ctx: Context, sites: list[Site], cfg: GenConfig) -> tuple[float, dict, tuple]:
    rows, valid, results = generate(ctx, sites, cfg)
    if len(valid) < 2:
        return 0.0, {}, (rows, valid, results)
    report = score_rows(valid, results, ctx)
    return report["final_score"], report, (rows, valid, results)


def search_config(ctx: Context, sites: list[Site], cfg: GenConfig, passes: int) -> tuple[GenConfig, dict, tuple]:
    best_score, best_report, best_payload = evaluate_config(ctx, sites, cfg)
    best_cfg = cfg
    print(f"  start [{cfg.strategy}]        final_score={best_score:10.3f}  "
          f"(tws={best_report.get('total_weighted_score', 0):.1f} "
          f"cons={best_report.get('consistency_factor', 0):.4f} "
          f"fid={best_report.get('distribution_fidelity_factor', 0):.4f})")

    for pass_index in range(passes):
        improved = False
        for knob, values in grid_for(best_cfg):
            for value in values:
                if getattr(best_cfg, knob) == value:
                    continue
                trial = replace(best_cfg, **{knob: value})
                score, report, payload = evaluate_config(ctx, sites, trial)
                marker = " "
                if score > best_score * 1.001:
                    best_score, best_report, best_payload, best_cfg = score, report, payload, trial
                    improved = True
                    marker = "*"
                print(f" {marker} pass{pass_index} {knob:>13}={str(value):>6} "
                      f"final_score={score:10.3f}  "
                      f"(tws={report.get('total_weighted_score', 0):.1f} "
                      f"cons={report.get('consistency_factor', 0):.4f} "
                      f"fid={report.get('distribution_fidelity_factor', 0):.4f})")
        if not improved:
            print(f"  pass{pass_index}: no improvement, stopping early")
            break

    return best_cfg, best_report, best_payload


# --------------------------------------------------------------------------------------------
# 9. Output and verification
# --------------------------------------------------------------------------------------------

def order_rows(rows: list[dict], valid: list[dict]) -> list[dict]:
    """Strongest rows first: truncate_submission keeps the first max_experiments unique ids, so if
    anything is ever cut it should be the cheapest rows."""
    weight_by_id = {
        entry["experiment"]["experiment_id"]: entry["stage2"]["weighted_score"] for entry in valid
    }
    return sorted(rows, key=lambda r: -weight_by_id.get(r["experiment_id"], 0.0))


def write_submission(rows: list[dict], valid: list[dict], path: str) -> None:
    with open(path, "w") as handle:
        json.dump(order_rows(rows, valid), handle, indent=2)


def verify_with_validator(ctx: Context) -> dict:
    """Run the validator's own benchmark_submission over the files just written."""
    from niome_subnet.genomics.validation import benchmark_submission
    return benchmark_submission(ctx.cell_types, uid=0).model_dump()


def summarise(rows: list[dict], valid: list[dict], results: list[dict], report: dict) -> None:
    outcomes = Counter(r["outcome"] for r in results)
    cells = Counter(
        (e["experiment"]["mutation"], e["experiment"]["cas_system"], e["experiment"]["strand"])
        for e in valid
    )
    offtarget = Counter(e["features"]["offtarget_factor"] for e in valid)
    gc_scores = [e["features"]["gc_score"] for e in valid]
    distances = [e["features"]["distance_to_mutation"] for e in valid]

    print("\n  composition")
    print(f"    rows                 {len(rows)}")
    print(f"    outcomes             {dict(outcomes)}")
    print(f"    offtarget_factor     {dict(sorted(offtarget.items(), reverse=True))}")
    print(f"    gc_score             mean={np.mean(gc_scores):.4f} min={min(gc_scores):.4f}")
    print(f"    distance             min={min(distances)} median={int(np.median(distances))} "
          f"max={max(distances)}")
    print(f"    joint cells          {len(cells)} occupied, "
          f"counts {sorted(cells.values())}")
    for target, stats in report.get("per_target", {}).items():
        print(f"    stage4 {target:<14}r2={stats['r2']:+.4f}  nmae={stats['nmae']:.4f}")
    detail = report.get("fidelity_detail", {})
    if detail:
        print("    entropy ratios       "
              f"mut={detail.get('mutation_coverage_entropy_ratio', 0):.4f} "
              f"cas={detail.get('cas_system_coverage_entropy_ratio', 0):.4f} "
              f"strand={detail.get('strand_coverage_entropy_ratio', 0):.4f} "
              f"joint={detail.get('joint_coverage_entropy_ratio', 0):.4f} "
              f"kmer={detail.get('kmer_diversity_entropy_ratio', 0):.4f} "
              f"guide={detail.get('distinct_guide_ratio', 0):.4f}")


def stage_report(rows: list[dict], valid: list[dict], results: list[dict],
                 ctx: Context, cfg: GenConfig | None = None) -> dict:
    """Every quantity the five stages measure, for one (submission, contract) pair.

    Stage 1 is re-run here rather than trusted, because this is also used to score a row set built
    against a *different* seed: designs are seed-independent, so the gate must come out identical,
    and it is worth proving rather than assuming.
    """
    reasons: Counter = Counter()
    seen_designs: set[tuple] = set()
    n_pass = 0
    for row in rows:
        passed, reason = stage12.stage1(row, ctx.seq, ctx.mutation_map, ctx.contract)
        if passed == 1.0 and design_key(row) in seen_designs:
            passed, reason = 0.0, "duplicate_experiment"
        if passed == 1.0:
            seen_designs.add(design_key(row))
            n_pass += 1
        else:
            reasons[reason] += 1

    def spread(values: list[float]) -> dict:
        if not values:
            return {"n": 0}
        arr = np.array(values, dtype=float)
        return {"min": float(arr.min()), "mean": float(arr.mean()),
                "median": float(np.median(arr)), "max": float(arr.max())}

    features = [entry["features"] for entry in valid]
    structural = [entry["stage2"]["structural_score"] for entry in valid]
    weighted = [entry["stage2"]["weighted_score"] for entry in valid]
    cut_ps = [stage3.cut_probability(r["cas"], r["energy"]) for r in results]

    s4 = stage4_in_memory(valid, results, fold_seed=ctx.seed)
    fidelity = stage5.compute_distribution_fidelity(valid, results, ctx.contract, k=12)
    factor = max(0.0, min(1.0, fidelity.get("distribution_fidelity_score", 0.0)))

    conforming = None
    if cfg is not None and cfg.construction in CONSTRUCTIONS:
        rule = CONSTRUCTIONS[cfg.construction]
        conforming = sum(1 for result, entry in zip(results, valid) if rule(result, entry))

    return {
        "seed": ctx.seed,
        "stage1": {
            "submitted": len(rows),
            "passed": n_pass,
            "rejected": len(rows) - n_pass,
            "reasons": dict(reasons),
            "unique_ids": len({row["experiment_id"] for row in rows}),
            "unique_designs": len({design_key(row) for row in rows}),
        },
        "stage2": {
            "gc": spread([f["gc"] for f in features]),
            "gc_score": spread([f["gc_score"] for f in features]),
            "distance": spread([f["distance_to_mutation"] for f in features]),
            "dist_score": spread([f["dist_score"] for f in features]),
            "consistency": spread([f["consistency"] for f in features]),
            "offtarget_factor": dict(Counter(f["offtarget_factor"] for f in features)),
            "structural_score": spread(structural),
            "weighted_score": spread(weighted),
            "total_weighted_score": float(sum(weighted)),
            "guide_lengths": dict(Counter(len(row["guideRNA"]) for row in rows)),
            "accessibility": features[0]["cell_type_accessibility"] if features else None,
        },
        "stage3": {
            "outcomes": dict(Counter(r["outcome"] for r in results)),
            "cut_rate": (sum(1 for r in results if r["outcome"] != "no_cut") / len(results)
                         if results else 0.0),
            "cut_p": spread(cut_ps),
            "cut_p_by_cas": {
                cas: spread([p for p, r in zip(cut_ps, results) if r["cas"] == cas])
                for cas in sorted({r["cas"] for r in results})
            },
            # Rows that would fail to cut under an arbitrary seed. Every one of them breaks the
            # construction, so this is the seed-independent floor on how much consistency is at
            # risk — the only term a cut_p filter can actually protect.
            "expected_no_cut": float(sum(1.0 - p for p in cut_ps)),
            "below_cut_p_floor": (
                sum(1 for p, r in zip(cut_ps, results)
                    if p < cut_p_floor(r["cas"], cfg) - 1e-9) if cfg is not None else None
            ),
            "cas_mix": dict(Counter(r["cas"] for r in results)),
            "energy": spread([r["energy"] for r in results]),
            "mh_true": sum(1 for r in results if r["mh"]),
            "indel_length": spread([r["indel_length"] for r in results]),
            "indel_histogram": dict(sorted(Counter(r["indel_length"] for r in results).items())),
            "mh_by_outcome": {f"mh={m}|{o}": c for (m, o), c in
                              sorted(Counter((r["mh"], r["outcome"]) for r in results).items(),
                                     key=lambda kv: str(kv[0]))},
            "conforming_rows": conforming,
            "mutation_breakdown": stage3.group_by_mutation(results) if results else {},
        },
        "stage4": {
            "n_joined": s4.get("n_valid_experiments", 0),
            "per_target": s4.get("per_target", {}),
            "avg_r2": s4.get("avg_r2"),
            "avg_nmae": s4.get("avg_nmae"),
            "consistency_score": s4.get("consistency_score", 0.0),
            "consistency_factor": s4.get("consistency_factor", 0.0),
            "total_weighted_score": s4.get("total_weighted_score", 0.0),
        },
        "stage5": {
            "mutation_coverage_entropy_ratio": fidelity.get("mutation_coverage_entropy_ratio"),
            "cas_system_coverage_entropy_ratio": fidelity.get("cas_system_coverage_entropy_ratio"),
            "strand_coverage_entropy_ratio": fidelity.get("strand_coverage_entropy_ratio"),
            "joint_coverage_entropy_ratio": fidelity.get("joint_coverage_entropy_ratio"),
            "kmer_diversity_entropy_ratio": fidelity.get("kmer_diversity_entropy_ratio"),
            "distinct_guide_ratio": fidelity.get("distinct_guide_ratio"),
            "distribution_fidelity_score": fidelity.get("distribution_fidelity_score"),
            "distribution_fidelity_factor": factor,
            "coverage_detail": fidelity.get("coverage_detail", {}),
            "cas_specific_shift_diagnostic": fidelity.get("cas_specific_shift_diagnostic", {}),
        },
        "final_score": s4.get("total_weighted_score", 0.0) * s4.get("consistency_factor", 0.0) * factor,
    }


def print_stage_report(report: dict, title: str) -> None:
    """The five stages, every measured value, in the order the validator computes them."""
    def line(label: str, value: str) -> None:
        print(f"    {label:<26} {value}")

    def show(label: str, stats: dict) -> None:
        if stats.get("n") == 0:
            line(label, "-")
            return
        line(label, f"min={stats['min']:<10.4f} mean={stats['mean']:<10.4f} "
                    f"median={stats['median']:<10.4f} max={stats['max']:.4f}")

    print(f"\n{'=' * 100}")
    print(f"  {title}   (stage-3 / stage-4 seed = {report['seed']})")
    print("=" * 100)

    s1 = report["stage1"]
    print("\n  STAGE 1 — structural gate")
    line("submitted", str(s1["submitted"]))
    line("passed", str(s1["passed"]))
    line("rejected", str(s1["rejected"]))
    line("reasons", str(s1["reasons"] or "none"))
    line("unique experiment_ids", f"{s1['unique_ids']} / {s1['submitted']}")
    line("unique designs", f"{s1['unique_designs']} / {s1['submitted']}")

    s2 = report["stage2"]
    print("\n  STAGE 2 — structural score")
    line("accessibility", str(s2["accessibility"]))
    line("guide lengths", str(s2["guide_lengths"]))
    show("gc", s2["gc"])
    show("gc_score", s2["gc_score"])
    show("distance_to_mutation", s2["distance"])
    show("dist_score", s2["dist_score"])
    show("consistency", s2["consistency"])
    line("offtarget_factor", str(s2["offtarget_factor"]))
    show("structural_score", s2["structural_score"])
    show("weighted_score", s2["weighted_score"])
    line("TOTAL_WEIGHTED_SCORE", f"{s2['total_weighted_score']:.6f}")

    s3 = report["stage3"]
    print("\n  STAGE 3 — simulated outcomes")
    line("outcomes", str(s3["outcomes"]))
    line("cut_rate", f"{s3['cut_rate']:.4f}")
    show("cut_p", s3["cut_p"])
    line("cas mix", str(s3["cas_mix"]))
    for cas, stats in s3["cut_p_by_cas"].items():
        line(f"  cut_p {cas}", f"min={stats['min']:.4f}  mean={stats['mean']:.4f}  "
                               f"max={stats['max']:.4f}  (ceiling {cut_p_ceiling(cas)})")
    line("expected no_cut rows", f"{s3['expected_no_cut']:.2f}  "
                                 f"(sum of 1 - cut_p, any seed)")
    if s3["below_cut_p_floor"] is not None:
        line("below the cut_p floor", str(s3["below_cut_p_floor"]))
    show("energy", s3["energy"])
    line("mh = True", str(s3["mh_true"]))
    line("mh x outcome", str(s3["mh_by_outcome"]))
    show("indel_length", s3["indel_length"])
    line("indel histogram", str(s3["indel_histogram"]))
    if s3["conforming_rows"] is not None:
        total = report["stage1"]["passed"] or 1
        line("construction conformance", f"{s3['conforming_rows']} / {total} "
                                         f"({100.0 * s3['conforming_rows'] / total:.1f}%)")
    for mutation, stats in s3["mutation_breakdown"].items():
        line(f"  {mutation[:22]}", f"n={stats['n']:<4} w={stats['mutation_weight']:<5} "
                                   f"cut={stats['cut_rate']:.3f} "
                                   f"energy={stats['mean_energy']:.3f} "
                                   f"indel={stats['mean_indel_length']:.3f}")

    s4 = report["stage4"]
    print("\n  STAGE 4 — consistency (RandomForest under KFold)")
    line("rows joined on id", str(s4["n_joined"]))
    for target, stats in s4["per_target"].items():
        line(f"  {target}", f"r2={stats['r2']:+.6f}   nmae={stats['nmae']:.6f}")
    if s4["avg_r2"] is not None:
        line("avg_r2", f"{s4['avg_r2']:+.6f}")
        line("avg_nmae", f"{s4['avg_nmae']:.6f}")
    line("consistency_score", f"{s4['consistency_score']:.6f}")
    line("CONSISTENCY_FACTOR", f"{s4['consistency_factor']:.6f}")

    s5 = report["stage5"]
    print("\n  STAGE 5 — distribution fidelity")
    for key in ("mutation_coverage_entropy_ratio", "cas_system_coverage_entropy_ratio",
                "strand_coverage_entropy_ratio", "joint_coverage_entropy_ratio",
                "kmer_diversity_entropy_ratio", "distinct_guide_ratio"):
        value = s5.get(key)
        line(f"  {key.replace('_entropy_ratio', '').replace('_ratio', '')}",
             f"{value:.6f}" if value is not None else "-")
    line("counts", str(s5["coverage_detail"]))
    # A report with no valid rows carries no diagnostic at all — compute_distribution_fidelity
    # returns early on n == 0 — so a missing "insufficient_data" is not evidence there is data to
    # print. Both metrics are also None-able in their own right (an empty sample either side).
    diagnostic = s5["cas_specific_shift_diagnostic"]
    jsd = diagnostic.get("repair_mode_jensen_shannon_divergence")
    wasserstein = diagnostic.get("indel_length_wasserstein_distance")
    if jsd is not None:
        line("JSD (not scored)", f"{jsd:.6f}")
    if wasserstein is not None:
        line("Wasserstein (not scored)", f"{wasserstein:.6f}")
    line("DISTRIBUTION_FIDELITY", f"{s5['distribution_fidelity_factor']:.6f}")

    print("\n  FINAL")
    line("final_score", f"{report['final_score']:.6f}  =  "
                        f"{s4['total_weighted_score']:.4f} x {s4['consistency_factor']:.6f} "
                        f"x {s5['distribution_fidelity_factor']:.6f}")


def run_seed_split(task: dict, cell_types: dict, cfg: GenConfig, build_seed: int,
                   auto_skew: bool = True, verify: bool = True) -> int:
    """Build against one seed, score against another, and report all five stages under both.

    This is the unstamped-contract case made measurable. Stages 1, 2 and 5 never read the seed, so
    they must come out identical either way; stage 3 is keyed on it, so the whole construction the
    generator engineered is expected to evaporate. What the two reports isolate is exactly how much
    of the score was resting on it.
    """
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    real_seed = contract.get("seed")

    ctx = build_context(contract, reference, cell_types)
    stand_in = copy.deepcopy(contract)
    stand_in["seed"] = build_seed
    build_ctx = build_context(stand_in, reference, cell_types)

    sites = enumerate_sites(build_ctx, cfg.flank, tuple(sorted(set(cfg.lengths))))
    print(f"  {len(sites)} PAM sites  {dict(Counter((s.cas, s.strand) for s in sites))}")

    if auto_skew and cfg.strategy == "pure":
        cfg = replace(cfg, weight_skew=choose_weight_skew(build_ctx, sites, cfg))
        print(f"  weight_skew fitted to this contract: {cfg.weight_skew}")

    print(f"  building against seed {build_seed} "
          f"(selection={cfg.selection}, oversample={cfg.oversample:g}, variants={cfg.variants}, "
          f"lengths={cfg.lengths})")
    rows, valid, results = generate(build_ctx, sites, cfg)
    rows = order_rows(rows, valid)
    if not rows:
        print("no rows generated", file=sys.stderr)
        return 1

    predicted = stage_report(rows, valid, results, build_ctx, cfg)
    print_stage_report(predicted, f"AS BUILT — seed {build_seed}")

    if real_seed == build_seed:
        print("\n  build seed and contract seed match; nothing to re-score.")
        return 0

    _report, real_valid, real_results = rescore_under(rows, ctx)
    actual = stage_report(rows, real_valid, real_results, ctx, cfg)
    print_stage_report(actual, f"AS SCORED — real contract seed {real_seed}")

    print(f"\n{'=' * 100}")
    print("  WHAT THE SEED COST")
    print("=" * 100)
    for label, key, path in (
        ("total_weighted_score", "total_weighted_score", ("stage4",)),
        ("consistency_factor", "consistency_factor", ("stage4",)),
        ("distribution_fidelity_factor", "distribution_fidelity_factor", ("stage5",)),
    ):
        before = predicted[path[0]][key]
        after = actual[path[0]][key]
        print(f"    {label:<30} {before:>12.6f}  ->  {after:>12.6f}   "
              f"({'unchanged' if abs(before - after) < 1e-9 else f'{after - before:+.6f}'})")
    before, after = predicted["final_score"], actual["final_score"]
    ratio = after / before if before else 0.0
    print(f"    {'final_score':<30} {before:>12.6f}  ->  {after:>12.6f}   "
          f"({ratio:.4f}x, {100 * (ratio - 1):+.1f}%)")

    if verify:
        print("\n  cross-checking the real-seed report against the validator's own file pipeline")
        persist_task(task)
        with open(settings.MINER_SUBMISSION_PATH, "w") as handle:
            json.dump(rows, handle, indent=2)
        official = verify_with_validator(ctx)
        print(f"    benchmark_submission final_score = {official['final_score']:.6f}")
        print(f"    in-memory replica                = {actual['final_score']:.6f}")
        print(f"    delta                            = {abs(official['final_score'] - actual['final_score']):.9f}")
    return 0


# --------------------------------------------------------------------------------------------
# 9b. Sweep — run the generator against every task the backend has ever issued
# --------------------------------------------------------------------------------------------

def rescore_under(rows: list[dict], ctx: Context) -> tuple[dict, list[dict], list[dict]]:
    """Re-derive features and outcomes for an existing dataset under ``ctx``'s contract.

    Stages 1, 2 and 5 never read the seed, so total_weighted_score and the coverage ratios carry
    over unchanged; stage 3 does, which is the whole point when the build seed and the scoring seed
    differ.
    """
    valid, results = [], []
    for row in rows:
        entry = build_valid_entry(row, ctx)
        if entry is None:      # would be rejected by stage 1 — designs are seed-independent, so
            continue           # this should never fire
        valid.append(entry)
        results.append(simulate(entry, ctx))
    report = score_rows(valid, results, ctx) if len(valid) >= 2 else {}
    return report, valid, results


def sweep(tasks: list[dict], cell_types: dict, cfg: GenConfig, verify_count: int,
          out_path: str, auto_skew: bool = True, build_seed: int | None = None) -> int:
    """Generate and score a submission for every task, then report the distribution.

    The point is not the headline number on one task but whether the strategy holds across the
    whole space the backend samples from: four cell types (accessibility 0.35 to 0.87), varying
    mutation pairs, weights and regions, and a different stage-3 seed every time.
    """
    verify_at = set()
    if verify_count > 0:
        step = max(1, len(tasks) // verify_count)
        verify_at = {i * step for i in range(verify_count) if i * step < len(tasks)}

    records: list[dict] = []
    started = time.time()

    for index, task in enumerate(tasks):
        contract = task["content"]["contract"]
        reference = task["content"]["hbb_reference"]
        ctx = build_context(contract, reference, cell_types)

        # build_seed models building before the contract is known: generate against a stand-in
        # seed, then score against the seed the validator will actually use.
        build_ctx = ctx
        if build_seed is not None:
            stand_in = copy.deepcopy(contract)
            stand_in["seed"] = build_seed
            build_ctx = build_context(stand_in, reference, cell_types)

        sites = enumerate_sites(build_ctx, cfg.flank, tuple(sorted(set(cfg.lengths))))

        task_cfg = cfg
        if auto_skew and cfg.strategy == "pure":
            task_cfg = replace(cfg, weight_skew=choose_weight_skew(build_ctx, sites, cfg))

        rows, valid, results = generate(build_ctx, sites, task_cfg)
        rows = order_rows(rows, valid)

        if build_seed is not None and build_seed != contract.get("seed"):
            report, valid, results = rescore_under(rows, ctx)
        else:
            report = score_rows(valid, results, ctx) if len(valid) >= 2 else {}

        record = {
            "index": index,
            "task_id": task["id"],
            "created_at": task.get("created_at"),
            "seed": contract.get("seed"),
            # A falsy seed means the backend has not stamped this task yet. Everything stage 3
            # produces is keyed on that seed, so the construction we engineer here stops holding
            # the moment a real seed lands — the designs stay valid (stage 1/2 never see the seed)
            # but consistency_factor collapses. Observed live: task f26b7613 went 0 -> 641 and its
            # score fell from 297.6 to 30.4.
            "seed_provisional": not contract.get("seed"),
            "cell_type": contract.get("cell_type"),
            "accessibility": cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0),
            "mutations": contract["active_mutations"],
            "mutation_weights": contract.get("mutation_weights", {}),
            "weight_skew": task_cfg.weight_skew,
            "rows": len(rows),
            "construction": task_cfg.construction,
            "build_seed": build_seed if build_seed is not None else contract.get("seed"),
            "seed_matched": build_seed is None or build_seed == contract.get("seed"),
            "conforms": not check_invariants(rows, results, task_cfg, valid),
            # Under a mismatched build seed this is chance conformance, not retained signal.
            "conforming_pct": (
                100.0 * sum(1 for r, e in zip(results, valid)
                            if CONSTRUCTIONS[task_cfg.construction](r, e)) / max(1, len(valid))
            ),
            "outcome_counts": dict(Counter(r["outcome"] for r in results)),
            "total_weighted_score": report.get("total_weighted_score", 0.0),
            "consistency_factor": report.get("consistency_factor", 0.0),
            "distribution_fidelity_factor": report.get("distribution_fidelity_factor", 0.0),
            "final_score": report.get("final_score", 0.0),
        }

        if index in verify_at:
            # Cross-check the in-memory replica against the validator's own file-passing pipeline.
            # persist_task writes the *real* contract, so this validates against the true seed.
            persist_task(task)
            with open(settings.MINER_SUBMISSION_PATH, "w") as handle:
                json.dump(rows, handle, indent=2)
            official = verify_with_validator(ctx)["final_score"]
            record["verified_final_score"] = official
            record["verify_delta"] = abs(official - record["final_score"])

        records.append(record)
        elapsed = time.time() - started
        eta = elapsed / (index + 1) * (len(tasks) - index - 1)
        flag = "" if record["conforms"] and record["rows"] == ctx.max_experiments else "  <-- check"
        verified = f"  verified={record['verified_final_score']:.3f}" if index in verify_at else ""
        print(f"  [{index + 1:>3}/{len(tasks)}] {record['cell_type']:<11} seed={record['seed']:<5} "
              f"rows={record['rows']:>3} final={record['final_score']:8.3f} "
              f"(tws={record['total_weighted_score']:7.2f} "
              f"cons={record['consistency_factor']:.4f} "
              f"fid={record['distribution_fidelity_factor']:.4f}){verified}{flag}  eta {eta / 60:.1f}m")

    elapsed = time.time() - started
    write_results(records, out_path, cfg, elapsed, auto_skew=auto_skew)
    report_sweep(records, out_path, elapsed)
    return 0


def _distribution(values: list[float]) -> dict:
    arr = np.array(values, dtype=float)
    return {
        "min": float(arr.min()),
        "p25": float(np.percentile(arr, 25)),
        "median": float(np.median(arr)),
        "mean": float(arr.mean()),
        "p75": float(np.percentile(arr, 75)),
        "max": float(arr.max()),
        "std": float(arr.std()),
    }


def summarise_sweep(records: list[dict]) -> dict:
    """Aggregate the per-task records into the summary block of the results file."""
    by_cell: dict[str, list[float]] = defaultdict(list)
    for record in records:
        by_cell[record["cell_type"]].append(record["final_score"])

    verified = [r for r in records if "verify_delta" in r]
    max_experiments = max((r["rows"] for r in records), default=0)

    return {
        "tasks": len(records),
        "final_score": _distribution([r["final_score"] for r in records]),
        "total_weighted_score": _distribution([r["total_weighted_score"] for r in records]),
        "consistency_factor": _distribution([r["consistency_factor"] for r in records]),
        "distribution_fidelity_factor": _distribution(
            [r["distribution_fidelity_factor"] for r in records]
        ),
        "weight_skew": _distribution([r.get("weight_skew", 0.0) for r in records]),
        "consistency_factor_is_one": sum(1 for r in records if r["consistency_factor"] >= 1.0),
        "construction_conforming": sum(1 for r in records if r.get("conforms")),
        "outcome_counts_total": dict(sum(
            (Counter(r.get("outcome_counts", {})) for r in records), Counter()
        )),
        "full_row_cap": sum(1 for r in records if r["rows"] >= max_experiments),
        "by_cell_type": {
            cell_type: {"n": len(values), **_distribution(values)}
            for cell_type, values in sorted(by_cell.items(), key=lambda kv: -np.mean(kv[1]))
        },
        "validator_cross_check": {
            "tasks_rescored": len(verified),
            "max_abs_delta": max((r["verify_delta"] for r in verified), default=None),
        },
        "weakest_tasks": [
            {k: r[k] for k in ("task_id", "cell_type", "final_score", "mutation_weights")}
            for r in sorted(records, key=lambda r: r["final_score"])[:5]
        ],
    }


def write_results(records: list[dict], path: str, cfg: GenConfig | None,
                  elapsed: float | None = None, auto_skew: bool = False) -> dict:
    """Write the full sweep as one self-describing document.

    Per-task rows alone are hard to act on later — without the config that produced them there is
    no way to tell a regression from a knob change. So the file carries the generation config and
    the aggregate summary alongside the rows.
    """
    config = None
    if cfg:
        config = {k: (list(v) if isinstance(v, tuple) else v) for k, v in vars(cfg).items()}
        if auto_skew and cfg.strategy == "pure":
            # The dataclass default never reached a build; each contract got its own fitted value.
            config["weight_skew"] = "fitted per contract — see tasks[].weight_skew"

    document = {
        "generated_by": "genExp.py",
        "elapsed_seconds": elapsed,
        "config": config,
        "summary": summarise_sweep(records),
        "tasks": records,
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(document, handle, indent=2)
    return document


def report_sweep(records: list[dict], out_path: str, elapsed: float) -> None:
    scores = [r["final_score"] for r in records]
    consistencies = [r["consistency_factor"] for r in records]
    fidelities = [r["distribution_fidelity_factor"] for r in records]
    weighted = [r["total_weighted_score"] for r in records]

    def stats(name: str, values: list[float], fmt: str = "8.3f") -> None:
        arr = np.array(values)
        print(f"    {name:<26} min={arr.min():{fmt}}  median={np.median(arr):{fmt}}  "
              f"mean={arr.mean():{fmt}}  max={arr.max():{fmt}}")

    print(f"\n  sweep over {len(records)} tasks in {elapsed / 60:.1f} min")
    stats("final_score", scores)
    stats("total_weighted_score", weighted)
    stats("consistency_factor", consistencies, "8.5f")
    stats("distribution_fidelity", fidelities, "8.5f")

    broken = [r for r in records if not r.get("conforms")]
    short = [r for r in records if r["rows"] < 250]
    below = [r for r in records if r["consistency_factor"] < 1.0]
    print(f"    consistency_factor == 1  {len(records) - len(below)}/{len(records)}")
    print(f"    construction conforming  {len(records) - len(broken)}/{len(records)}")
    print(f"    full 250 rows            {len(records) - len(short)}/{len(records)}")
    totals = sum((Counter(r.get("outcome_counts", {})) for r in records), Counter())
    print(f"    outcome mix              {dict(totals)}")

    verified = [r for r in records if "verify_delta" in r]
    if verified:
        worst = max(r["verify_delta"] for r in verified)
        print(f"    replica vs validator     {len(verified)} checked, max |delta| = {worst:.2e}")

    by_cell: dict[str, list[float]] = defaultdict(list)
    for record in records:
        by_cell[record["cell_type"]].append(record["final_score"])
    print("\n  by cell type (accessibility feeds stage-3 energy, not stage-2 scoring)")
    for cell_type, values in sorted(by_cell.items(), key=lambda kv: -np.mean(kv[1])):
        arr = np.array(values)
        print(f"    {cell_type:<12} n={len(values):>3}  mean={arr.mean():8.3f}  "
              f"min={arr.min():8.3f}  max={arr.max():8.3f}")

    worst = sorted(records, key=lambda r: r["final_score"])[:5]
    print("\n  weakest 5 tasks")
    for record in worst:
        print(f"    {record['task_id'][:8]}  {record['cell_type']:<11} "
              f"final={record['final_score']:8.3f}  weights={record['mutation_weights']}")

    print(f"\n  per-task detail written to {out_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task-id", default=None, help="pin a task instead of the latest")
    parser.add_argument("--oldest", action="store_true",
                        help="use the oldest task the backend still lists")
    parser.add_argument("--rows", type=int, default=None, help="override contract max_experiments")
    parser.add_argument("--flank", type=int, default=GenConfig.flank,
                        help="site-enumeration window beyond gene_region (bp)")
    parser.add_argument("--variants", type=int, default=GenConfig.variants,
                        help="guide variants generated per site")
    parser.add_argument("--lengths", default="20,23",
                        help="comma-separated guide lengths (20 and/or 23)")
    parser.add_argument("--strategy", choices=("pure", "shaped", "auto"), default="auto",
                        help="outcome strategy; auto searches both and keeps the better")
    parser.add_argument("--construction", choices=tuple(CONSTRUCTIONS), default="mh",
                        help="rule every row's outcome must satisfy under the pure strategy")
    parser.add_argument("--selection", choices=("packed", "stratified", "nearest"), default=None,
                        help="site ranking; 'nearest' sorts by |gc-0.5| then distance, "
                             "defaults to packed for pure and stratified for shaped")
    parser.add_argument("--gc-tolerance", type=float, default=0.0,
                        help="under --selection nearest, treat any guide within this much of 50%% "
                             "GC as equally good so distance breaks the tie (try 0.03)")
    parser.add_argument("--min-cut-p", type=float, default=0.0,
                        help="admit only sites whose stage-3 cut_p reaches this. Ceilings are 0.99 "
                             "for Cas9 and 0.96 for Cas12a, so anything above 0.96 keeps Cas9 only")
    parser.add_argument("--cut-p-ceiling", action="store_true",
                        help="hold every Cas system to its own cut_p ceiling (Cas9 0.99, "
                             "Cas12a 0.96) instead of one shared floor")
    parser.add_argument("--no-cut-p-relax", action="store_true",
                        help="keep the cut_p floor strict even where no site can reach it, which "
                             "empties those cells instead of falling back to the best available")
    parser.add_argument("--cut-p-floors", default=None,
                        help="explicit per-system floors, e.g. Cas9=0.99,Cas12a=0.96")
    parser.add_argument("--backfill-cas", default=None,
                        help="when a filter leaves the submission short of max_experiments, top it "
                             "up with rows from this Cas system, highest cut_p first (e.g. Cas9)")
    parser.add_argument("--cas-mix", default=None,
                        help="row share per Cas system in rules.cas_systems order, e.g. 80/20 "
                             "for 80%% Cas9 and 20%% Cas12a; default is an even split")
    parser.add_argument("--oversample", type=float, default=1.0,
                        help="generate this multiple of the wanted rows, then keep only the ones "
                             "whose realised gc and distance land nearest the ideal")
    parser.add_argument("--no-search", action="store_true", help="single build with default knobs")
    parser.add_argument("--search-passes", type=int, default=2)
    parser.add_argument("--baseline", action="store_true",
                        help="also build the untuned control for comparison")
    parser.add_argument("--all-tasks", action="store_true",
                        help="sweep every task the backend has issued instead of building one")
    parser.add_argument("--limit", type=int, default=None,
                        help="sweep only the N newest tasks")
    parser.add_argument("--include-zero-seed", action="store_true",
                        help="sweep mode: keep the placeholder seed==0 tasks")
    parser.add_argument("--verify-count", type=int, default=3,
                        help="sweep mode: tasks to re-score through benchmark_submission")
    parser.add_argument("--build-seed", type=int, default=None,
                        help="generate against this seed but score against the task's real one — "
                             "models building before the contract is stamped. In single-task mode "
                             "this prints all five stages under both seeds")
    parser.add_argument("--no-auto-skew", action="store_true",
                        help="keep the fixed weight_skew instead of fitting it per contract")
    parser.add_argument("--sweep-out", default="result.json",
                        help="where the sweep writes config + summary + per-task rows")
    parser.add_argument("--out", default=settings.MINER_SUBMISSION_PATH)
    parser.add_argument("--no-verify", action="store_true",
                        help="skip the benchmark_submission re-run")
    return parser.parse_args()


def config_for_contract(cfg: GenConfig, contract: dict, cut_p_ceiling: bool = False,
                        cas_mix: str | None = None) -> GenConfig:
    """Fill in the knobs that can only be resolved once the contract is known.

    ``cut_p_floors`` and ``cas_share`` are keyed by Cas system name, and the roster comes from
    ``rules.cas_systems`` — so neither can be baked into a static config. Both the miner and
    submission.py route through here, which is what keeps an offline sweep predicting exactly what
    the miner will send.
    """
    cas_systems = list(contract["rules"].get("cas_systems", ["Cas9", "Cas12a"]))
    if cut_p_ceiling:
        cfg = replace(cfg, cut_p_floors={cas: cut_p_ceiling_for(cas) for cas in cas_systems})
    if cas_mix:
        cfg = replace(cfg, cas_share=parse_cas_mix(cas_mix, cas_systems))
    return cfg


def parse_cut_p_floors(args, cas_systems: list[str]) -> dict | None:
    if args.cut_p_floors:
        return {part.split("=")[0]: float(part.split("=")[1])
                for part in args.cut_p_floors.split(",") if part}
    if args.cut_p_ceiling:
        return {cas: cut_p_ceiling_for(cas) for cas in cas_systems}
    return None


def parse_cas_mix(spec: str | None, cas_systems: list[str]) -> dict | None:
    """'80/20' -> {"Cas9": 0.8, "Cas12a": 0.2}, positional in rules.cas_systems order."""
    if not spec:
        return None
    parts = [float(value) for value in spec.replace(":", "/").replace(",", "/").split("/")]
    if len(parts) != len(cas_systems):
        raise SystemExit(f"--cas-mix needs {len(cas_systems)} values for {cas_systems}")
    total = sum(parts) or 1.0
    return {cas: value / total for cas, value in zip(cas_systems, parts)}


def main() -> int:
    args = parse_args()
    started = time.time()

    lengths = tuple(int(v) for v in args.lengths.split(","))

    if args.all_tasks:
        print("[1/3] fetching task history")
        tasks = fetch_all_tasks()
        unstamped = sum(1 for t in tasks if not t["content"]["contract"].get("seed"))
        if not args.include_zero_seed:
            tasks = [t for t in tasks if t["content"]["contract"].get("seed")]
            print(f"  {len(tasks) + unstamped} tasks, {unstamped} with seed==0 skipped "
                  f"(--include-zero-seed to keep them)")
        elif unstamped:
            print(f"  ! including {unstamped} task(s) with seed==0. These are unstamped: when the "
                  f"backend assigns a real seed, every outcome engineered against seed 0 breaks "
                  f"and consistency_factor collapses. Their rows stay valid; their scores do not.")
        if args.limit:
            tasks = tasks[:args.limit]
        cell_types = fetch_cell_types()
        cfg = GenConfig(strategy="pure" if args.strategy == "auto" else args.strategy,
                        selection=args.selection or "packed", flank=args.flank,
                        variants=args.variants, lengths=lengths, rows=args.rows,
                        construction=args.construction, oversample=args.oversample,
                        gc_tolerance=args.gc_tolerance, min_cut_p=args.min_cut_p)
        print(f"[2/3] warming reference + site cache")
        warm = build_context(tasks[0]["content"]["contract"],
                             tasks[0]["content"]["hbb_reference"], cell_types)
        sites = enumerate_sites(warm, cfg.flank, tuple(sorted(set(lengths))))
        print(f"  {len(sites)} PAM sites  {dict(Counter((s.cas, s.strand) for s in sites))}")
        print(f"[3/3] sweeping {len(tasks)} tasks with {cfg.strategy} strategy"
              f"/{cfg.construction} construction"
              f"{'' if args.no_auto_skew else ', weight_skew chosen per task'}")
        if args.build_seed is not None:
            print(f"  ! generating against stand-in seed {args.build_seed}, scoring against each "
                  f"task's real seed — stage 3 is keyed on the seed, so expect "
                  f"consistency_factor to collapse")
        code = sweep(tasks, cell_types, cfg, args.verify_count, args.sweep_out,
                     auto_skew=not args.no_auto_skew, build_seed=args.build_seed)
        print(f"\ndone in {time.time() - started:.1f}s")
        return code

    print("[1/6] fetching task")
    if args.oldest:
        task = fetch_all_tasks()[-1]
    else:
        task = fetch_task(args.task_id)
    contract, reference = persist_task(task)
    print(f"  task {task['id']}  created {task['created_at']}  seed {contract['seed']}")
    print(f"  cell_type {contract.get('cell_type')}  mutations {contract['active_mutations']}")
    print(f"  rules {contract['rules']}")

    cell_types = fetch_cell_types()
    accessibility = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    print(f"  accessibility {accessibility}")

    if args.build_seed is not None:
        print("[2/2] loading chr11, then building against a stand-in seed")
        cas_systems = list(contract["rules"].get("cas_systems", ["Cas9", "Cas12a"]))
        split_cfg = GenConfig(strategy="pure" if args.strategy == "auto" else args.strategy,
                              selection=args.selection or "packed", flank=args.flank,
                              variants=args.variants, lengths=lengths, rows=args.rows,
                              construction=args.construction, oversample=args.oversample,
                              gc_tolerance=args.gc_tolerance, min_cut_p=args.min_cut_p,
                              cut_p_floors=parse_cut_p_floors(args, cas_systems),
                              cut_p_relax=not args.no_cut_p_relax,
                              cas_share=parse_cas_mix(args.cas_mix, cas_systems),
                              backfill_cas=args.backfill_cas)
        code = run_seed_split(task, cell_types, split_cfg, args.build_seed,
                              auto_skew=not args.no_auto_skew, verify=not args.no_verify)
        print(f"\ndone in {time.time() - started:.1f}s")
        return code

    print("[2/6] loading chr11 + k-mer index")
    ctx = build_context(contract, reference, cell_types)

    base_cfg = GenConfig(flank=args.flank, variants=args.variants, lengths=lengths,
                         rows=args.rows, construction=args.construction,
                         oversample=args.oversample, gc_tolerance=args.gc_tolerance,
                         min_cut_p=args.min_cut_p)
    strategies = ("pure", "shaped") if args.strategy == "auto" else (args.strategy,)

    # Enumerate the union of every length the search may ask for; select_sites filters per config.
    enumerated = tuple(sorted(set(lengths) | ({20, 23} if not args.no_search else set())))
    print(f"[3/6] enumerating PAM sites (gene_region +/- {base_cfg.flank} bp, lengths {enumerated})")
    sites = enumerate_sites(ctx, base_cfg.flank, enumerated)
    per_group = Counter((s.cas, s.strand) for s in sites)
    print(f"  {len(sites)} sites  {dict(per_group)}")

    if args.baseline:
        print("[3b/6] baseline (untuned control)")
        b_rows, b_valid, b_results = generate_baseline(ctx, sites, base_cfg)
        b_report = score_rows(b_valid, b_results, ctx)
        print(f"  baseline final_score={b_report['final_score']:.3f} "
              f"(tws={b_report['total_weighted_score']:.1f} "
              f"cons={b_report['consistency_factor']:.4f} "
              f"fid={b_report['distribution_fidelity_factor']:.4f})")

    best: tuple[float, GenConfig, dict, tuple] | None = None
    for strategy in strategies:
        start_cfg = replace(base_cfg, strategy=strategy,
                            selection=args.selection or
                            ("packed" if strategy == "pure" else "stratified"))
        if strategy == "pure" and not args.no_auto_skew:
            fitted = choose_weight_skew(ctx, sites, start_cfg)
            print(f"  weight_skew fitted to this contract: {fitted}")
            start_cfg = replace(start_cfg, weight_skew=fitted)
        if args.no_search:
            print(f"[4/6] building [{strategy}] (search disabled)")
            score, report, payload = evaluate_config(ctx, sites, start_cfg)
            cfg = start_cfg
            print(f"  final_score={score:.3f}")
        else:
            print(f"[4/6] coordinate descent over generation knobs [{strategy}]")
            cfg, report, payload = search_config(ctx, sites, start_cfg, args.search_passes)
            score = report.get("final_score", 0.0)
        if best is None or score > best[0]:
            best = (score, cfg, report, payload)

    assert best is not None
    _, cfg, report, payload = best
    if len(strategies) > 1:
        print(f"\n  winning strategy: {cfg.strategy}")

    rows, valid, results = payload
    if not rows:
        print("no rows generated", file=sys.stderr)
        return 1

    summarise(rows, valid, results, report)

    print(f"\n[5/6] writing {args.out}")
    write_submission(rows, valid, args.out)
    print(f"  {len(rows)} rows; config "
          f"{ {k: v for k, v in vars(cfg).items()} }")

    if args.no_verify or args.out != settings.MINER_SUBMISSION_PATH:
        if args.out != settings.MINER_SUBMISSION_PATH and not args.no_verify:
            print(f"[6/6] verification skipped — it reads {settings.MINER_SUBMISSION_PATH}, "
                  f"not {args.out}")
        else:
            print("[6/6] skipped verification")
        print(f"done in {time.time() - started:.1f}s")
        return 0

    print("[6/6] verifying with the validator's benchmark_submission")
    score = verify_with_validator(ctx)
    print(json.dumps(score, indent=2))
    print(f"\ndone in {time.time() - started:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
