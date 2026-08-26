"""seed_agnostic.py — build a submission whose cut survives an unknown restamped seed.

The miner sees ``contract.seed == 0`` at broadcast and the validator scores under a real seed
assigned later, so the outcome the miner engineers against seed 0 does not hold at scoring time.
Only one of stage 4's three targets can be made seed-independent: ``is_cut``. A guide cuts iff the
second draw of its per-row RNG stream is <= ``cut_p``; on a high-accessibility cell type Cas9's
``cut_p`` is 0.99, so a guide can cut under *every* seed in a window, and a group of Cas12a guides
(``cut_p`` 0.96, no single strict guide possible) can be chosen so their combined no_cut seeds are
few — leaving many seeds on which the whole submission cuts and ``is_cut`` is constant.

This module builds that submission: a strict-all-window Cas9 bank plus a min-union Cas12a group,
assembled at the score-optimal cas mix. It is used only when ``cut_p`` is high enough for strict Cas9
to exist — i.e. not on HEK293 (accessibility 0.35, where even Cas9 tops out at 0.95 and no strict
guide exists); the caller routes HEK293 to the ordinary construction.

The seed-test — sha256 -> 32-bit seed -> MT19937 draws — is the whole cost and is isolated in
``cut_fail_seeds``/``_scan_site`` so a GPU kernel can replace the CPU multiprocessing path without
touching the assembly logic above it.
"""

from __future__ import annotations

import hashlib
import itertools
import logging
import math
import multiprocessing as mp
import random
import time
from collections import Counter
from dataclasses import dataclass, field
from multiprocessing import Pool

import numpy as np

import genExp as G
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.validation import stage3

logger = logging.getLogger(__name__)

NUCLEOTIDES = ("A", "C", "G", "T")


@dataclass
class SeedAgnosticConfig:
    start_seed: int = 100
    end_seed: int = 999               # the seed window the hedge covers
    gc_min: float = 0.40
    gc_max: float = 0.60
    max_distance: int = 200           # near-mutation sites: cut_p at the clamp, high stage-2 score
    variants: int = 20000             # per-target enumeration cap (needs to be large for strict Cas9)
    # Time budget. The scan stops when this is exhausted even if the pool cap is unmet, so a build
    # always returns inside the presigned URL's TTL rather than overrunning it and losing the round.
    time_budget_s: float = 210.0
    cas_mix: str = "85/15"            # score-optimal on the robustness sweep
    cas12a_max_fail: int = 22         # accept Cas12a guides that no_cut under <= this in the window
    per_cell_min: int = 8             # floor per (mutation, cas, strand) for stage-5 coverage
    group_restarts: int = 12
    jobs: int = 12
    # Stop scanning once this many candidates are banked per Cas system. The scan exists to feed a
    # 250-row submission, so enumerating the whole reachable space is pure waste: Cas12a's deep
    # max_fail scan gets no early-out and cost 450s to bank 7,969 candidates when the min-union
    # group needs ~38. Oversampling ~20x the rows keeps the group's choice wide.
    # Cas9 strict guides contribute *zero* failed seeds by definition, so a bigger Cas9 pool buys
    # only stage-2 quality, never clean seeds. Cas12a is the opposite: every clean seed the
    # submission has comes from the min-union group, and the union shrinks as the pool it is drawn
    # from grows (measured 146 -> 447 union, 920 -> 403, 2490 -> 362, 7969 -> 329). So the budget
    # belongs to Cas12a, and Cas9 is capped at just above what the assembly consumes.
    pool_cap_cas9: int = 300      # ~1.4x the 212 Cas9 rows the assembly takes
    pool_cap_cas12a: int = 9000   # effectively the whole reachable pool
    use_gpu: bool = True          # falls back to CPU-only automatically if no device is usable
    # Cas9 must finish. Its strict guides contribute zero failed seeds, but falling short of the
    # ~212 rows the assembly wants forces extra Cas12a rows in from *outside* the min-union group,
    # and those unhedged rows wreck the union (measured: 134 strict -> union 876 -> 2.7% clean,
    # against 350 strict -> union 330 -> 63.3%). So Cas9 gets a floor in seconds, not a share, and
    # Cas12a spends whatever is left.
    cas9_budget_s: float = 130.0
    # The deadline can only be honoured between targets, and one 20,000-variant target takes tens of
    # seconds — so a naive budget overruns by roughly one target's cost (measured 60s -> 93s). Stop
    # claiming this much early so the real finish lands inside the budget rather than past it.
    claim_margin_s: float = 35.0
    max_experiments: int | None = None   # defaults to contract rules

    @property
    def seeds(self) -> range:
        return range(self.start_seed, self.end_seed + 1)


# --------------------------------------------------------------------------------------------
# The seed-test seam. Everything below build_* calls only these two; a GPU implementation
# swaps them out wholesale. The CPU path fans site scans across a process pool.
# --------------------------------------------------------------------------------------------

def _energy_at_gc(gc: float, distance: int, accessibility: float, region_offset: float) -> float:
    return max(0.0, min(1.0, accessibility * (
        1.8 * gc + 0.6 * math.exp(-distance / 1500) + region_offset)))


def cut_fail_seeds(mutation: str, cas: str, guide: str, start: int, strand: str,
                   cut_p: float, seeds: range, max_fail: int) -> frozenset | None:
    """Reference (scalar) seed-test for one guide. Kept as the oracle ``mt19937.screen_guides`` is
    checked against; the hot path batches along the guide axis instead.

    Seeds in the window this guide no_cuts under, or None once it exceeds ``max_fail``.

    The cut coin is the second draw of ``random.Random(experiment_seed)``; ``experiment_seed`` is the
    low 32 bits of sha256(round_seed|design). Returns None as an early-out so a hopeless guide is
    abandoned after ``max_fail + 1`` failures rather than scanning the whole window.
    """
    base = f"|{mutation}|{cas}|{guide}|{start}|{strand}"
    fails: list[int] = []
    for sd in seeds:
        digest = hashlib.sha256((str(sd) + base).encode()).digest()
        rng = random.Random(int.from_bytes(digest[-4:], "big"))
        rng.random()                       # microhomology coin
        if rng.random() > cut_p:           # cut coin
            fails.append(sd)
            if len(fails) > max_fail:
                return None
    return frozenset(fails)


STATE: dict = {}


def _worker_state(payload: dict) -> None:
    """Build the read-only per-worker scan state. Safe to call in the host process."""
    task = payload["task"]
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    ctx = G.build_context(contract, reference, payload["cell_types"])
    STATE.update(
        ctx=ctx,
        sites=G.enumerate_sites(ctx, 3000, (20, 23)),
        cfg=payload["cfg"],
        accessibility=payload["cell_types"].get(contract.get("cell_type"), {})
        .get("accessibility", 1.0),
        region_offset={m: stage3.REGION_ENERGY_OFFSETS.get((contract.get("mutation_regions") or {})
                                                            .get(m), 0.0) for m in ctx.mutations},
        max_fail=payload["max_fail"],
        seed_array=np.fromiter(payload["cfg"].seeds, dtype=np.int64),
    )


def _worker_init(payload: dict) -> None:
    """Forked-child initializer: quiet the inherited logger, then build the scan state.

    The setLevel is confined to children on purpose. _claim_loop also runs in a *thread* of the
    host process for the GPU lane, and raising the root level there silences the host — a miner
    whose heartbeat and upload logs vanish mid-round looks indistinguishable from a dead one.
    """
    logging.getLogger().setLevel(logging.ERROR)
    _worker_state(payload)


def _scan_site(job: tuple[int, str]) -> list[dict]:
    """Enumerate one (site, mutation) target and keep guides that cut under >= (window - max_fail).

    This is the unit a GPU kernel would replace: it is pure seed-testing over an enumerated guide
    list, with no shared state beyond the read-only context.
    """
    site_index, mutation = job
    ctx, cfg = STATE["ctx"], STATE["cfg"]
    site = STATE["sites"][site_index]
    distance = abs(site.start - ctx.mutation_map[mutation])
    if distance > cfg.max_distance:
        return []

    accessibility = STATE["accessibility"]
    region_offset = STATE["region_offset"][mutation]
    guides = enumerate_variants(site, ctx, cfg.gc_min, cfg.gc_max, ctx.max_mismatches,
                                True, cfg.variants)
    if not guides:
        return []
    max_fail = STATE["max_fail"]

    # cut_p varies with GC across a target's variants, so the screen is handed a per-guide lookup.
    cut_p_by_gc: dict[int, float] = {}

    def cut_p_of(guide: str) -> float:
        gc_count = sum(b in "GC" for b in guide)
        value = cut_p_by_gc.get(gc_count)
        if value is None:
            energy = _energy_at_gc(gc_count / site.length, distance, accessibility, region_offset)
            value = cut_p_by_gc[gc_count] = stage3.cut_probability(site.cas, energy)
        return value

    # One batched screen for the whole target instead of a scalar pass per guide: the kernel is only
    # efficient on wide batches (11 us/pair at 900, 1.3 us/pair at 450k), and dropping busted guides
    # after each seed slice keeps the strict search's early-out.
    screen = MT.screen_guides_gpu if STATE.get("use_gpu") else MT.screen_guides
    survivors = screen(guides, STATE["seed_array"], mutation, site.cas,
                       site.start, site.strand, cut_p_of, max_fail)

    kept = []
    for guide, fails in survivors.items():
        entry = G.build_valid_entry(G.make_experiment(site, guide, mutation, ctx, "cand"), ctx)
        if entry is None:
            continue
        kept.append({
            "guide": guide, "mutation": mutation, "cas_system": site.cas, "strand": site.strand,
            "target_alignment_start": site.start, "length": site.length,
            "cut_p": cut_p_of(guide), "weighted_score": entry["stage2"]["weighted_score"],
            "n_fail": int(fails.size), "fails": [int(x) for x in fails],
        })
    return kept


def enumerate_variants(site, ctx, gc_lo, gc_hi, budget, require_clean, cap):
    """Every guide within Hamming <= budget of the target whose GC lands in the band (validated)."""
    ref = list(site.ref_guide)
    L = site.length
    lo = math.ceil(gc_lo * L - 1e-9)
    hi = math.floor(gc_hi * L + 1e-9)
    sl = G.seed_slice(site.cas, L)
    ref_gc = sum(b in "GC" for b in ref)
    positions = list(range(L))
    seen: set[str] = set()
    out: list[str] = []
    for k in range(0, budget + 1):
        for combo in itertools.combinations(positions, k):
            choices = [[b for b in NUCLEOTIDES if b != ref[i]] for i in combo]
            for repl in itertools.product(*choices):
                gc = ref_gc
                for i, b in zip(combo, repl):
                    gc += (1 if b in "GC" else 0) - (1 if ref[i] in "GC" else 0)
                if not (lo <= gc <= hi):
                    continue
                guide = ref[:]
                for i, b in zip(combo, repl):
                    guide[i] = b
                gs = "".join(guide)
                if gs in seen:
                    continue
                if require_clean and gs[sl] in ctx.kmer_index:
                    continue
                seen.add(gs)
                out.append(gs)
                if len(out) >= cap:
                    return out
    return out


# --------------------------------------------------------------------------------------------
# Bank scan (parallel), min-union group, assembly.
# --------------------------------------------------------------------------------------------

def scan_bank(task, cell_types, cfg: SeedAgnosticConfig, ctx, cas, max_fail, budget_s=None):
    """All guides (across near targets) that cut under >= (window - max_fail) for the given Cas."""
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == cas for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    # Interleave the (mutation, strand) cells so an early stop still leaves every cell represented;
    # a plain site-ordered walk would fill one cell and starve the rest.
    by_cell: dict[tuple, list] = {}
    for job in jobs:
        site = sites[job[0]]
        by_cell.setdefault((job[1], site.strand), []).append(job)
    interleaved = [job for group in itertools.zip_longest(*by_cell.values())
                   for job in group if job is not None]

    payload = {"task": task, "cell_types": cell_types, "cfg": cfg, "max_fail": max_fail}
    cap = cfg.pool_cap_cas9 if cas == "Cas9" else cfg.pool_cap_cas12a
    budget_s = cfg.time_budget_s if budget_s is None else budget_s

    banked: list[dict] = []
    scanned = 0
    deadline = time.monotonic() + budget_s
    stopped = "exhausted"
    with Pool(cfg.jobs, initializer=_worker_init, initargs=(payload,)) as pool:
        results = pool.imap_unordered(_scan_site, interleaved, chunksize=1)
        for kept in results:
            banked.extend(kept)
            scanned += 1
            if len(banked) >= cap:
                stopped = "cap"
                pool.terminate()      # cap met; the remaining targets add nothing the group needs
                break
            if time.monotonic() > deadline:
                # Overrunning the TTL loses the whole round, so a partial bank beats a late one.
                stopped = "time"
                pool.terminate()
                break
    logger.info("%s scan: %d/%d targets, %d candidates (%s; cap %d, budget %.0fs)",
                cas, scanned, len(interleaved), len(banked), stopped, cap, budget_s)
    return banked


def _claim_loop(payload: dict, jobs: list, counter, lock, out_q, use_gpu: bool,
                deadline: float, stop) -> None:
    """One worker: claim the next unclaimed target, scan it, repeat until the list or time runs out.

    Work is claimed from a shared counter rather than pre-split, because per-target cost varies by
    an order of magnitude (a target yields 0 to 20,000 variants) and the GPU and a CPU core differ
    ~6x in rate. A static split would leave one side idle; claiming keeps both busy to the end.
    """
    try:
        # Only a forked child may touch the root logger; the GPU lane shares the host's.
        (_worker_init if mp.parent_process() is not None else _worker_state)(payload)
        STATE["use_gpu"] = use_gpu
        batch: list[dict] = []
        # Reserve the margin: a target claimed just before the deadline still has to finish.
        claim_until = deadline - payload["claim_margin_s"]
        while not stop.value and time.monotonic() < claim_until:
            with lock:
                index = counter.value
                if index >= len(jobs):
                    break
                counter.value = index + 1
            batch.extend(_scan_site(jobs[index]))
            if len(batch) >= 32:
                out_q.put(batch)
                batch = []
        if batch:
            out_q.put(batch)
    except Exception as exc:                       # a dead worker must not hang the join
        logger.warning("scan worker failed (%s): %s", "gpu" if use_gpu else "cpu", exc)
    finally:
        out_q.put(None)                            # this worker is done


def scan_bank_hybrid(task, cell_types, cfg: SeedAgnosticConfig, ctx, cas, max_fail,
                     budget_s=None) -> list[dict]:
    """Scan targets on the GPU and every spare CPU core at once.

    Ordering matters: the CPU children are forked *before* the parent touches cupy, because a CUDA
    context does not survive fork — a child inheriting an initialised context cannot use the device,
    and in the worst case wedges it. Forking first keeps the children pure-numpy and lets the parent
    own the GPU.
    """
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == cas for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    by_cell: dict[tuple, list] = {}
    for job in jobs:
        by_cell.setdefault((job[1], sites[job[0]].strand), []).append(job)
    jobs = [j for group in itertools.zip_longest(*by_cell.values())
            for j in group if j is not None]

    cap = cfg.pool_cap_cas9 if cas == "Cas9" else cfg.pool_cap_cas12a
    budget_s = cfg.time_budget_s if budget_s is None else budget_s
    deadline = time.monotonic() + budget_s
    payload = {"task": task, "cell_types": cell_types, "cfg": cfg, "max_fail": max_fail,
               "claim_margin_s": cfg.claim_margin_s}

    context = mp.get_context("fork")
    counter = context.Value("i", 0)
    stop = context.Value("i", 0)
    lock = context.Lock()
    out_q = context.Queue()

    workers = [context.Process(target=_claim_loop,
                               args=(payload, jobs, counter, lock, out_q, False, deadline, stop),
                               daemon=True)
               for _ in range(max(1, cfg.jobs))]
    for worker in workers:
        worker.start()                             # forked before any CUDA call in this process

    gpu_thread = None
    if cfg.use_gpu:
        import threading
        gpu_thread = threading.Thread(
            target=_claim_loop,
            args=(payload, jobs, counter, lock, out_q, True, deadline, stop), daemon=True)
        gpu_thread.start()

    expected = len(workers) + (1 if gpu_thread else 0)
    banked: list[dict] = []
    finished = 0
    while finished < expected:
        item = out_q.get()
        if item is None:
            finished += 1
            continue
        banked.extend(item)
        if len(banked) >= cap:
            stop.value = 1                         # cap met; let the workers wind down

    for worker in workers:
        worker.join(timeout=5)
        if worker.is_alive():
            worker.terminate()
    if gpu_thread:
        gpu_thread.join(timeout=5)

    why = ("cap" if len(banked) >= cap
           else "time" if time.monotonic() >= deadline else "exhausted")
    logger.info("%s hybrid scan: %d/%d targets claimed, %d candidates (%s)",
                cas, counter.value, len(jobs), len(banked), why)
    return banked


def min_union_group(candidates: list[dict], group_size: int, cfg: SeedAgnosticConfig,
                    seeds_n: int) -> list[dict]:
    """Forward-greedy + swap min-union over each candidate's no_cut seeds, with per-cell floors."""
    fails = [frozenset(c["fails"]) for c in candidates]
    cell = [(c["mutation"], c["cas_system"], c["strand"]) for c in candidates]
    available = Counter(cell)
    # The floors must be jointly satisfiable: at a small group size (a 90/10 mix leaves only ~25
    # Cas12a rows over 4 cells) a per-cell minimum of 8 would need 32, and the selector could then
    # never clear its unmet-floor bookkeeping. Clamp to an even share so every cell still gets rows —
    # an empty cell zeroes its stage-5 coverage ratio, which is the failure worth preventing.
    per_cell = min(cfg.per_cell_min, max(1, group_size // max(1, len(available))))
    floor = {k: min(per_cell, available[k]) for k in available}
    rng = random.Random(0xC0FFEE)

    def build(jitter):
        chosen, cover, cell_count = set(), Counter(), Counter()
        while len(chosen) < min(group_size, len(candidates)):
            slots = group_size - len(chosen)
            unmet = sum(max(0, floor.get(c, 0) - cell_count[c]) for c in floor)
            free = slots > unmet
            best, best_cost = None, 10 ** 9
            for j in range(len(candidates)):
                if j in chosen:
                    continue
                if not (cell_count[cell[j]] < floor.get(cell[j], 0) or free):
                    continue
                cost = sum(1 for sd in fails[j] if cover[sd] == 0)
                if cost < best_cost or (jitter and cost == best_cost and rng.random() < 0.5):
                    best_cost, best = cost, j
            if best is None:
                break
            chosen.add(best)
            cell_count[cell[best]] += 1
            for sd in fails[best]:
                cover[sd] += 1
        return chosen

    def union(chosen):
        u = set()
        for i in chosen:
            u |= fails[i]
        return len(u)

    best, best_u = None, 10 ** 9
    for r in range(cfg.group_restarts):
        chosen = build(jitter=(r > 0))
        u = union(chosen)
        if u < best_u:
            best, best_u = set(chosen), u
    logger.info("Cas12a min-union group: %d guides, union %d -> %d clean seeds of %d",
                len(best), best_u, seeds_n - best_u, seeds_n)
    return [candidates[i] for i in sorted(best)]


def _to_row(rec: dict, cell_type: str, exp_id: str) -> dict:
    start = rec["target_alignment_start"]
    return {
        "experiment_id": exp_id, "guideRNA": rec["guide"],
        "target_alignment_start": start, "target_alignment_end": start + rec["length"],
        "strand": rec["strand"], "mutation": rec["mutation"],
        "cas_system": rec["cas_system"], "cell_type": cell_type,
    }


def _repair_coverage(chosen: list[dict], ctx, contract: dict, cell_types: dict,
                     cfg: SeedAgnosticConfig, n: int) -> list[dict]:
    """Guarantee every cell is occupied and the row cap is filled, seed-agnostic or not.

    The seed-agnostic scan is best-effort under a time budget, so it can come back short or miss a
    cell entirely. Both are worse failures than a row that no_cuts under some seeds: a short
    submission loses term 1 linearly, and an empty cell costs a ~0.03x multiplier via stage 5's
    geometric mean. Filling from plain near-mutation guides trades a little seed robustness for
    keeping the score's structure intact.
    """
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    have = Counter((r["mutation"], r["cas_system"], r["strand"]) for r in chosen)
    used = {(r["cas_system"], r["target_alignment_start"], r["strand"], r["guide"])
            for r in chosen}
    wanted_cells = [(m, cas, st) for m in ctx.mutations
                    for cas in ctx.cas_systems for st in ("+", "-")]

    def fillers(cell, limit):
        """Best untaken near-mutation guides for one cell, strongest stage-2 score first."""
        mutation, cas, strand = cell
        position = ctx.mutation_map[mutation]
        out = []
        for site in sorted((s for s in sites if s.cas == cas and s.strand == strand),
                           key=lambda s: abs(s.start - position)):
            if abs(site.start - position) > cfg.max_distance * 4:
                break
            for guide in enumerate_variants(site, ctx, cfg.gc_min, cfg.gc_max,
                                            ctx.max_mismatches, True, 40):
                key = (cas, site.start, strand, guide)
                if key in used:
                    continue
                entry = G.build_valid_entry(
                    G.make_experiment(site, guide, mutation, ctx, "fill"), ctx)
                if entry is None:
                    continue
                used.add(key)
                out.append({"guide": guide, "mutation": mutation, "cas_system": cas,
                            "strand": strand, "target_alignment_start": site.start,
                            "length": site.length, "cut_p": 0.0,
                            "weighted_score": entry["stage2"]["weighted_score"],
                            "n_fail": None, "fails": [], "seed_agnostic": False})
                if len(out) >= limit:
                    return out
        return out

    # 1. every cell must reach the coverage floor
    for cell in wanted_cells:
        short = cfg.per_cell_min - have[cell]
        if short > 0:
            added = fillers(cell, short)
            chosen.extend(added)
            have[cell] += len(added)
            if added:
                logger.info("coverage repair: %s +%d ordinary guide(s)", cell, len(added))

    # 2. top the submission up to the row cap, spreading across the thinnest cells first
    while len(chosen) < n:
        cell = min(wanted_cells, key=lambda c: have[c])
        added = fillers(cell, 1)
        if not added:
            break
        chosen.extend(added)
        have[cell] += 1
    return chosen


def build_submission(contract: dict, reference: dict, cell_types: dict, task: dict,
                     cfg: SeedAgnosticConfig | None = None,
                     banks: tuple[list, list] | None = None) -> tuple[list[dict], dict]:
    """Assemble a seed-agnostic submission: strict Cas9 + min-union Cas12a at the optimal cas mix.

    ``banks`` reuses an already-scanned (cas9, cas12a_pool) pair. The scan is the whole cost and it
    depends only on the contract, not on the cas mix or the row cap — so one scan can serve repeated
    broadcasts of the same task, or a sweep over mixes. ``meta["banks"]`` returns them for reuse.
    """
    try:
        return _assemble(contract, reference, cell_types, task, cfg, banks)
    finally:
        # Freed at the build boundary, not between the two scans: within one build the Cas12a scan
        # reuses the Cas9 scan's blocks, and every exit path has to release — the fallbacks (short
        # rows, too much backfill, an exception) leave the device just as occupied as a success.
        freed = MT.free_gpu_memory()
        if freed:
            logger.info("released %.2f GB of pooled GPU memory", freed / 1e9)


def _assemble(contract: dict, reference: dict, cell_types: dict, task: dict,
              cfg: SeedAgnosticConfig | None = None,
              banks: tuple[list, list] | None = None) -> tuple[list[dict], dict]:
    """Body of ``build_submission``; see there. Split out so the GPU pool is freed on every path."""
    cfg = cfg or SeedAgnosticConfig()
    ctx = G.build_context(contract, reference, cell_types)
    n = cfg.max_experiments or ctx.max_experiments
    cas9_frac = int(cfg.cas_mix.split("/")[0]) / 100.0
    n_cas9 = round(n * cas9_frac)
    n_cas12a = n - n_cas9
    seeds_n = len(cfg.seeds)

    if banks is not None:
        cas9, cas12a_pool = banks
    else:
        # Cas9: strict over the window (max_fail 0). is_cut stays constant only if every Cas9 row
        # cuts under every window seed, so no fail budget here.
        started = time.monotonic()
        # Cas9 needs ~93s to finish, but never more than a majority of the budget: on a short
        # budget a fixed floor would consume everything and leave Cas12a — which earns every clean
        # seed — with nothing.
        cas9 = scan_bank_hybrid(task, cell_types, cfg, ctx, "Cas9", max_fail=0,
                                budget_s=min(cfg.cas9_budget_s, cfg.time_budget_s * 0.6))
        # Cas12a: no strict guide exists; keep the <= max_fail pool, then min-union n_cas12a of it.
        cas12a_pool = scan_bank_hybrid(
            task, cell_types, cfg, ctx, "Cas12a", max_fail=cfg.cas12a_max_fail,
            budget_s=max(20.0, cfg.time_budget_s - (time.monotonic() - started)))
    # Size the group to the rows that will really be Cas12a: if the Cas9 bank came up short, the
    # shortfall becomes Cas12a rows, and every one of them must come from the min-union selection
    # rather than be backfilled outside it.
    n_cas12a = max(n_cas12a, n - min(len(cas9), n_cas9))
    group = min_union_group(cas12a_pool, n_cas12a, cfg, seeds_n) if cas12a_pool else []

    # Cas9 rows: mutation-weighted like the miner, even across strands, strongest first per cell.
    by_cell: dict[tuple, list[dict]] = {}
    for r in cas9:
        by_cell.setdefault((r["mutation"], r["strand"]), []).append(r)
    for rows in by_cell.values():
        rows.sort(key=lambda r: -r["weighted_score"])
    weights = contract.get("mutation_weights", {})
    shares = {m: max(weights.get(m, 1.0), 1e-9) ** 1.25 for m in ctx.mutations}
    tot = sum(shares.values())
    chosen = list(group)
    for m in ctx.mutations:
        per_mut = round(n_cas9 * shares[m] / tot)
        for strand, want in (("+", per_mut // 2), ("-", per_mut - per_mut // 2)):
            chosen.extend(by_cell.get((m, strand), [])[:want])

    # Top up any shortfall from the strongest unused Cas9 rows.
    if len(chosen) < n:
        used = {(r["cas_system"], r["target_alignment_start"], r["strand"], r["guide"])
                for r in chosen}
        spare = sorted((r for r in cas9 if (r["cas_system"], r["target_alignment_start"],
                        r["strand"], r["guide"]) not in used), key=lambda r: -r["weighted_score"])
        chosen.extend(spare[:n - len(chosen)])

    # A cell with no rows zeroes its stage-5 coverage ratio, and the geometric mean's 1e-9 clip
    # turns that into a ~0.03x multiplier on the whole score — worse than any seed fragility. So
    # every (mutation, cas, strand) cell is backfilled from ordinary near-mutation guides when the
    # seed-agnostic pool could not supply it, and the row count is filled to the cap the same way.
    chosen = _repair_coverage(chosen, ctx, contract, cell_types, cfg, n)

    rows, seen = [], set()
    for i, rec in enumerate(chosen[:n]):
        key = (rec["cas_system"], rec["target_alignment_start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(_to_row(rec, contract.get("cell_type"), f"exp-{i:05d}"))

    # Rows outside the hedge and the window seeds the hedge fails on: the caller needs both to
    # judge whether this build is worth shipping over the ordinary construction. A submission with
    # many unhedged rows measured *worse* than the ordinary build (66 backfilled -> 0% clean ->
    # final 25.65 against ~29), because the clean-seed property is all-or-nothing across all rows.
    backfilled = sum(1 for rec in chosen[:n] if not rec.get("seed_agnostic", True))
    union: set[int] = set()
    for rec in group:
        union.update(rec.get("fails", ()))
    window = len(cfg.seeds)

    meta = {
        "backfilled": backfilled,
        "hedged_rows": len(rows) - backfilled,
        "group_failed_seed_union": len(union),
        "clean_seed_estimate": window - len(union),
        "clean_fraction_estimate": (window - len(union)) / window if window else 0.0,
        "strict_cas9_available": len(cas9),
        "cas12a_pool": len(cas12a_pool),
        "cas12a_group": len(group),
        "rows": len(rows),
        "cas_mix": dict(Counter(r["cas_system"] for r in rows)),
        "cell_counts": dict(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
        "window": [cfg.start_seed, cfg.end_seed],
        "banks": (cas9, cas12a_pool),
    }
    return rows, meta
