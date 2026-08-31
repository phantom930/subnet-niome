"""mt19937.py — batched, bit-exact replica of ``random.Random(int).random()``'s first draws.

Stage 3 seeds one RNG per row with ``random.Random(experiment_seed(round_seed, design))`` and then
takes two draws: the microhomology coin, then the cut coin. Deciding whether a guide cuts under a
seed therefore costs one full Mersenne-Twister initialisation, and a seed-agnostic search needs
millions of them. ``random.Random`` does that one seed at a time in C; this does a whole batch at
once with numpy, so the same work runs as ~1900 array operations instead of ~1900 * batch scalar ones.

**Bit-exactness is the whole contract.** CPython seeds an integer with ``init_by_array`` (not
``init_genrand``, which is what numpy's legacy ``RandomState`` uses for a scalar — they are *not*
interchangeable), and builds a double from two tempered words as ``(a >> 5) * 2**26 + (b >> 6)``
over ``2**53``. Both are reproduced here exactly, and ``verify()`` checks a large sample against
``random.Random`` itself. A silent mismatch would not crash anything — it would quietly return the
wrong guides — so nothing here should be "optimised" without re-running that check.
"""

from __future__ import annotations

import hashlib
import random

import numpy as np

N = 624
M = 397
MATRIX_A = np.uint32(0x9908B0DF)
UPPER_MASK = np.uint32(0x80000000)
LOWER_MASK = np.uint32(0x7FFFFFFF)

_MULT1 = np.uint32(1664525)          # init_by_array's first mixing multiplier
_MULT2 = np.uint32(1566083941)       # ...and its second
_INIT_KEY = np.uint32(19650218)      # the fixed seed init_by_array starts from
_TWO26 = 67108864.0
_TWO53 = 9007199254740992.0


def _init_genrand_constant() -> np.ndarray:
    """``init_genrand(19650218)`` — the state init_by_array starts from, identical for every seed."""
    mt = np.empty(N, dtype=np.uint32)
    mt[0] = _INIT_KEY
    mult = np.uint32(1812433253)
    # uint32 wraparound is the algorithm, not an accident, so the overflow notice is noise here.
    with np.errstate(over="ignore"):
        for i in range(1, N):
            prev = mt[i - 1]
            mt[i] = mult * (prev ^ (prev >> np.uint32(30))) + np.uint32(i)
    return mt


_INIT_CONST = _init_genrand_constant()


def _gpu():
    """cupy if a usable GPU is present, else None. Imported lazily so CPU-only hosts are unaffected."""
    global _GPU_MOD, _GPU_CHECKED
    if _GPU_CHECKED:
        return _GPU_MOD
    _GPU_CHECKED = True
    try:
        import cupy
        cupy.arange(1, dtype=cupy.uint32) + cupy.uint32(1)   # forces context + a JIT compile
        _GPU_MOD = cupy
    except Exception as exc:                                  # no driver, no wheel, no headers
        logger_msg = f"GPU unavailable ({type(exc).__name__}); staying on numpy"
        print(logger_msg) if False else None
        _GPU_MOD = None
    return _GPU_MOD


_GPU_MOD = None
_GPU_CHECKED = False


def free_gpu_memory() -> int:
    """Return cupy's cached device blocks to the driver; returns the bytes released.

    cupy's pool only ever grows, so a long-lived host keeps its peak allocation for the rest of the
    process — a miner that scanned once sits on ~7.4 GB of an 8 GB card forever. That is invisible
    until something else wants the device, and then it is not: the next scan's GPU lane dies on an
    OOM, the scan finishes CPU-only, and CPU-only cannot reach 250 rows inside the upload TTL, so
    the build quietly falls back to the ordinary construction. No-op if cupy was never loaded, so
    calling it on a CPU-only host does not import anything.
    """
    if _GPU_MOD is None:
        return 0
    pool = _GPU_MOD.get_default_memory_pool()
    before = pool.total_bytes()
    pool.free_all_blocks()
    _GPU_MOD.get_default_pinned_memory_pool().free_all_blocks()
    return before - pool.total_bytes()


# The state is (624, batch) uint32 = 2496 bytes per seed, so the batch is capped by device memory,
# not by anything algorithmic: 8M seeds would want ~20 GB. Chunking keeps a call of any size safe.
GPU_CHUNK = 1_500_000        # ~3.7 GB of state, comfortable on an 8 GB card


# One thread per seed, the whole init_by_array in that thread's local memory, one launch. The
# cupy-array version below does the same arithmetic but as ~1,247 device-wide kernel launches, which
# is latency-bound: it reaches ~5M pairs/s against this kernel's ~40M. Kept side by side because the
# array version is readable and is what `verify_kernel` checks this against.
_MT_SOURCE = r"""
extern "C" __global__
void mt_first_draws(const unsigned int* __restrict__ seeds, int n, int n_draws, double* out)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;

    unsigned int mt[624];                       // 2,496 B of local memory per thread
    // init_genrand(19650218): the fixed state init_by_array starts from.
    mt[0] = 19650218u;
    for (int i = 1; i < 624; ++i) {
        unsigned int prev = mt[i - 1];
        mt[i] = 1812433253u * (prev ^ (prev >> 30)) + (unsigned int)i;
    }

    unsigned int seed = seeds[tid];
    int i = 1;
    // Pass 1, 624 mixes. CPython's key array for an int < 2**32 is one word, so j is 0 whenever it
    // is read and `init_key[j] + j` collapses to the seed.
    for (int k = 0; k < 624; ++k) {
        unsigned int prev = mt[i - 1];
        mt[i] = (mt[i] ^ ((prev ^ (prev >> 30)) * 1664525u)) + seed;
        ++i;
        if (i >= 624) { mt[0] = mt[623]; i = 1; }
    }
    // Pass 2, 623 further mixes, subtracting the index.
    for (int k = 0; k < 623; ++k) {
        unsigned int prev = mt[i - 1];
        mt[i] = (mt[i] ^ ((prev ^ (prev >> 30)) * 1566083941u)) - (unsigned int)i;
        ++i;
        if (i >= 624) { mt[0] = mt[623]; i = 1; }
    }
    mt[0] = 0x80000000u;

    // The first 2*n_draws tempered words. Word kk reads mt[kk + 397]; with n_draws <= 3 that is at
    // most mt[402], still the *pre-generate* value, so no regeneration pass is needed.
    unsigned int w[6];
    int n_words = 2 * n_draws;
    for (int kk = 0; kk < n_words; ++kk) {
        unsigned int y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7FFFFFFFu);
        unsigned int v = mt[kk + 397] ^ (y >> 1) ^ ((y & 1u) ? 0x9908B0DFu : 0u);
        v ^= (v >> 11);
        v ^= (v << 7) & 0x9D2C5680u;
        v ^= (v << 15) & 0xEFC60000u;
        v ^= (v >> 18);
        w[kk] = v;
    }
    for (int k = 0; k < n_draws; ++k) {
        double a = (double)(w[2 * k] >> 5);
        double b = (double)(w[2 * k + 1] >> 6);
        out[(long long)k * n + tid] = (a * 67108864.0 + b) / 9007199254740992.0;
    }
}
"""

_MT_KERNEL = None
KERNEL_CHUNK = 4_000_000     # 4M threads x 2,496 B local = ~10 GB of backing store, streamed


def _mt_kernel(xp):
    global _MT_KERNEL
    if _MT_KERNEL is None:
        _MT_KERNEL = xp.RawKernel(_MT_SOURCE, "mt_first_draws")
    return _MT_KERNEL


def first_n_draws_kernel(seeds, xp, n_draws: int = 2, chunk: int = KERNEL_CHUNK) -> tuple:
    """``first_n_draws_gpu`` via the single-launch CUDA kernel. Same values, ~8x the throughput."""
    if not 1 <= n_draws <= 3:
        raise ValueError("n_draws must be 1..3")
    seeds = xp.asarray(seeds, dtype=xp.uint32)
    n = int(seeds.shape[0])
    if n == 0:
        return tuple(xp.empty(0, dtype=xp.float64) for _ in range(n_draws))
    if n > chunk:
        parts = [[] for _ in range(n_draws)]
        for begin in range(0, n, chunk):
            got = first_n_draws_kernel(seeds[begin:begin + chunk], xp, n_draws, chunk)
            for k, arr in enumerate(got):
                parts[k].append(arr)
        return tuple(xp.concatenate(part) for part in parts)
    out = xp.empty((n_draws, n), dtype=xp.float64)
    block = 128
    grid = (n + block - 1) // block
    _mt_kernel(xp)((grid,), (block,), (seeds, np.int32(n), np.int32(n_draws), out))
    return tuple(out[k] for k in range(n_draws))


def first_two_draws_gpu(seeds, xp, chunk: int = GPU_CHUNK) -> tuple:
    """The first two ``random()`` values on the GPU. Thin wrapper over ``first_n_draws_gpu``."""
    return first_n_draws_gpu(seeds, xp, 2, chunk)


def first_n_draws_gpu(seeds, xp, n_draws: int = 2, chunk: int = GPU_CHUNK,
                      use_kernel: bool = True) -> tuple:
    """The first ``n_draws`` ``random()`` values for each seed, on the GPU.

    Each ``random()`` consumes two tempered words, so ``n_draws`` needs ``2 * n_draws`` of them.
    Stage 3 takes three: the microhomology coin, the cut coin, then the repair-mode draw. Three is
    the practical ceiling here and it costs nothing extra — word ``k`` reads ``mt[k + M]`` with
    M = 397, so six words stay inside the 624-word state and no regeneration step is needed.
    (A fourth draw, the indel length, would be gamma/exponential variates rather than a plain
    ``random()``, so it is deliberately not covered.)

    Kept separate from the numpy version so the CPU path stays the reference implementation that
    ``verify()`` pins down: the GPU result is checked against it, not the other way round.
    """
    if not 1 <= n_draws <= 3:
        raise ValueError("n_draws must be 1..3 (see the docstring on why 3 is the ceiling)")
    if use_kernel:
        return first_n_draws_kernel(seeds, xp, n_draws)
    seeds = xp.asarray(seeds, dtype=xp.uint32)
    if seeds.shape[0] > chunk:
        # Split rather than fail: a caller should not have to know the device's memory budget.
        parts = [[] for _ in range(n_draws)]
        for begin in range(0, int(seeds.shape[0]), chunk):
            got = first_n_draws_gpu(seeds[begin:begin + chunk], xp, n_draws, chunk)
            for k, arr in enumerate(got):
                parts[k].append(arr)
        return tuple(xp.concatenate(part) for part in parts)
    n_words = 2 * n_draws
    batch = int(seeds.shape[0])
    mt = xp.repeat(xp.asarray(_INIT_CONST)[:, None], batch, axis=1)

    u30, u1 = xp.uint32(30), xp.uint32(1)
    m1, m2 = xp.uint32(1664525), xp.uint32(1566083941)
    i = 1
    for _ in range(N):
        prev = mt[i - 1]
        mt[i] = (mt[i] ^ ((prev ^ (prev >> u30)) * m1)) + seeds
        i += 1
        if i >= N:
            mt[0] = mt[N - 1]
            i = 1
    for _ in range(N - 1):
        prev = mt[i - 1]
        mt[i] = (mt[i] ^ ((prev ^ (prev >> u30)) * m2)) - xp.uint32(i)
        i += 1
        if i >= N:
            mt[0] = mt[N - 1]
            i = 1
    mt[0] = xp.uint32(0x80000000)

    out = xp.empty((n_words, batch), dtype=xp.uint32)
    for kk in range(n_words):
        y = (mt[kk] & xp.uint32(0x80000000)) | (mt[kk + 1] & xp.uint32(0x7FFFFFFF))
        out[kk] = mt[kk + M] ^ (y >> u1) ^ xp.where(
            (y & u1).astype(bool), xp.uint32(0x9908B0DF), xp.uint32(0)).astype(xp.uint32)
    for kk in range(n_words):
        y = out[kk]
        y = y ^ (y >> xp.uint32(11))
        y = y ^ ((y << xp.uint32(7)) & xp.uint32(0x9D2C5680))
        y = y ^ ((y << xp.uint32(15)) & xp.uint32(0xEFC60000))
        out[kk] = y ^ (y >> xp.uint32(18))

    return tuple(
        ((out[2 * k] >> xp.uint32(5)).astype(xp.float64) * _TWO26
         + (out[2 * k + 1] >> xp.uint32(6)).astype(xp.float64)) / _TWO53
        for k in range(n_draws)
    )


# The repair-mode weights, mirroring stage3.repair_mode. Duplicated rather than imported because
# this runs on device arrays, not scalars — but the numbers must track stage 3, and verify_rule()
# is what catches it if they drift.
_HDR_BASE = {"Cas9": 0.32, "Cas12a": 0.24}
_BLUNT_W = 0.35

# rule -> (mh_branch_target, no_mh_branch_target) over {"HDR","MH_NHEJ","BLUNT_NHEJ"}, or a set of
# acceptable outcomes when the rule does not depend on the mh coin.
RULE_SPECS = {
    "hdr": {"any": ("HDR",)},
    "not_mhnhej": {"any": ("HDR", "BLUNT_NHEJ")},
    "not_hdr": {"any": ("MH_NHEJ", "BLUNT_NHEJ")},
    # No repair condition at all — the row only has to cut. This is the production hedge's
    # criterion, expressed as a rule so the 3-draw screen can bank it with its failed-seed sets.
    "cut": {"any": ("HDR", "MH_NHEJ", "BLUNT_NHEJ")},
    "mh_any": {"mh": ("HDR",), "no_mh": ("BLUNT_NHEJ",)},
}


def outcomes_gpu(d1, d2, d3, gc, energy, cut_p, cas: str, xp):
    """Stage 3's per-row outcome for a batch, from the three draws. Returns (mh, code).

    ``code`` is 0 no_cut, 1 HDR, 2 MH_NHEJ, 3 BLUNT_NHEJ. ``gc``/``energy``/``cut_p`` broadcast
    against the draws, so they may be scalars or per-guide columns. Every comparison matches stage
    3's exactly: ``mh`` is ``draw < p_mh`` (strict), a row cuts iff ``draw <= cut_p`` (stage 3 tests
    ``> cut_p`` for the no-cut branch), and the mode is picked by ``r = draw * total`` against the
    cumulative weights in the order HDR, MH_NHEJ, BLUNT_NHEJ.
    """
    p_mh = xp.minimum(0.6, 2.2 * gc * (1.0 - gc))
    mh = d1 < p_mh
    cut = d2 <= cut_p
    hdr_w = _HDR_BASE[cas] + 0.35 * energy
    mh_nhej = xp.where(mh, 0.30, 0.12)
    r = d3 * (hdr_w + mh_nhej + _BLUNT_W)
    code = xp.where(r < hdr_w, 1, xp.where(r < hdr_w + mh_nhej, 2, 3))
    return mh, xp.where(cut, code, 0)


def rule_fails_gpu(d1, d2, d3, gc, energy, cut_p, cas: str, rule: str, xp):
    """Boolean mask of rows that BREAK ``rule``. A no_cut always breaks it (it is not an outcome)."""
    spec = RULE_SPECS[rule]
    mh, code = outcomes_gpu(d1, d2, d3, gc, energy, cut_p, cas, xp)
    ids = {"HDR": 1, "MH_NHEJ": 2, "BLUNT_NHEJ": 3}
    if "any" in spec:
        ok = xp.zeros(code.shape, dtype=bool)
        for name in spec["any"]:
            ok |= code == ids[name]
    else:
        ok_mh = xp.zeros(code.shape, dtype=bool)
        for name in spec["mh"]:
            ok_mh |= code == ids[name]
        ok_no = xp.zeros(code.shape, dtype=bool)
        for name in spec["no_mh"]:
            ok_no |= code == ids[name]
        ok = xp.where(mh, ok_mh, ok_no)
    return ~ok


def screen_guides_rule_gpu(guides: list[str], round_seeds: np.ndarray, mutation: str, cas: str,
                           start: int, strand: str, params_of, rule: str, max_fail: int,
                           target_pairs: int = 1_200_000) -> dict[str, np.ndarray]:
    """``screen_guides_gpu`` for a repair-mode rule instead of the cut gate alone.

    ``params_of(guide) -> (gc, energy, cut_p)``. Returns {guide: sorted failed seeds} for guides
    breaking the rule under at most ``max_fail`` seeds.

    Two passes rather than one, because the host work dominated otherwise. Pass 1 keeps everything
    on the device and brings back only a per-guide failure *count* per slice (a few KB), which is
    all the early-out needs. Pass 2 recovers the actual failed seeds for the survivors alone — a
    small set by definition, since anything with more than ``max_fail`` failures was dropped.
    Measured on one slice: draws 78%, the old per-guide python loop 16%, sha256 4.5%. ``params_of``
    used to be called once per guide *per slice*, rescanning every guide string; it is hoisted.
    End to end after both changes: 23.9M pairs/s at target_pairs 1.2M (17.3M at 2.4M — past the
    peak the state stops fitting), against 7.2M before and 1.87M for the 13-core CPU path.
    """
    from niome_subnet.genomics import sha256_gpu as SH

    xp = _gpu()
    if xp is None:
        raise RuntimeError("screen_guides_rule_gpu needs a GPU; no CPU fallback is implemented")
    if not guides:
        return {}

    suffix_of = {g: f"|{mutation}|{cas}|{g}|{start}|{strand}".encode() for g in guides}
    all_params = np.asarray([params_of(g) for g in guides], dtype=np.float64)   # once, not per slice
    index_of = {g: i for i, g in enumerate(guides)}
    window = int(round_seeds.size)

    # Pack the suffix table once per call: it is constant while the seed slice moves, and
    # re-packing it per slice was the dominant host cost (46% GPU utilisation on a 915-target run).
    suf_bytes, suf_lens = SH.pack([suffix_of[g] for g in guides])
    d_sufb_all = xp.asarray(suf_bytes)
    d_ulen_all = xp.asarray(suf_lens)
    d_params = xp.asarray(all_params)

    def sweep(subset_rows, collect):
        """Run the window over row indices; count failures, and collect seeds when asked."""
        alive = np.asarray(subset_rows, dtype=np.int64)
        counts = np.zeros(len(guides), dtype=np.int64)
        seeds_hit: dict[int, list[int]] = {int(i): [] for i in alive}
        begin = 0
        while begin < window and alive.size:
            span = max(32, min(window - begin, target_pairs // max(1, int(alive.size))))
            slice_seeds = round_seeds[begin:begin + span]
            digits = [str(int(sd)).encode() for sd in slice_seeds]
            rows_dev = xp.asarray(alive)
            values = SH.experiment_seeds_gpu_packed(
                digits, d_sufb_all[rows_dev], d_ulen_all[rows_dev], xp)
            d1, d2, d3 = first_n_draws_gpu(values, xp, 3)
            shape = (int(alive.size), slice_seeds.size)
            prm = d_params[rows_dev]
            bad = rule_fails_gpu(d1.reshape(shape), d2.reshape(shape), d3.reshape(shape),
                                 prm[:, 0:1], prm[:, 1:2], prm[:, 2:3], cas, rule, xp)
            per_guide = xp.asnumpy(bad.sum(axis=1))          # the only transfer in pass 1
            if collect:
                host = xp.asnumpy(bad)
                for i, row in enumerate(alive):
                    if per_guide[i]:
                        seeds_hit[int(row)].extend(int(x) for x in slice_seeds[host[i]])
            counts[alive] += per_guide
            alive = alive[counts[alive] <= max_fail]
            begin += span
        return alive, seeds_hit

    survivors, _ = sweep(np.arange(len(guides)), collect=False)
    if survivors.size == 0:
        return {}
    _, hits = sweep(survivors, collect=True)
    return {guides[int(i)]: np.array(sorted(hits[int(i)]), dtype=np.int64) for i in survivors}


def screen_guides_gpu(guides: list[str], round_seeds: np.ndarray, mutation: str, cas: str,
                      start: int, strand: str, cut_p_of, max_fail: int,
                      target_pairs: int = 600_000) -> dict[str, np.ndarray]:
    """Fully device-side screen: sha256 and the MT draws both run on the GPU.

    Same contract and same result as ``screen_guides``, and checked against it. Two differences that
    matter for throughput:

    * hashing happens on the device from a per-seed digit table and a per-guide suffix table, so the
      bus carries O(guides + seeds) bytes rather than one message per pair;
    * the seed slice grows as guides are eliminated. A fixed slice would start wide and end tiny —
      and the GPU only beats numpy above roughly 450k pairs per call — so the slice is sized to hold
      the batch near ``target_pairs`` all the way down the window.
    """
    from niome_subnet.genomics import sha256_gpu as SH

    xp = _gpu()
    if xp is None:
        return screen_guides(guides, round_seeds, mutation, cas, start, strand,
                             cut_p_of, max_fail)

    alive = list(guides)
    fails: dict[str, list[int]] = {g: [] for g in alive}
    suffix_of = {g: f"|{mutation}|{cas}|{g}|{start}|{strand}".encode() for g in alive}

    begin = 0
    window = int(round_seeds.size)
    while begin < window and alive:
        span = max(32, min(window - begin, target_pairs // max(1, len(alive))))
        slice_seeds = round_seeds[begin:begin + span]
        digits = [str(int(sd)).encode() for sd in slice_seeds]

        values = SH.experiment_seeds_gpu(digits, [suffix_of[g] for g in alive], xp)
        _mh, cut = first_two_draws_gpu(values, xp)
        cut = cut.reshape(len(alive), slice_seeds.size)
        thresholds = xp.asarray([cut_p_of(g) for g in alive], dtype=xp.float64)[:, None]
        over = xp.asnumpy(cut > thresholds)

        survivors = []
        for index, guide in enumerate(alive):
            hits = over[index]
            if hits.any():
                fails[guide].extend(int(x) for x in slice_seeds[hits])
                if len(fails[guide]) > max_fail:
                    continue
            survivors.append(guide)
        alive = survivors
        begin += span

    return {g: np.array(sorted(fails[g]), dtype=np.int64) for g in alive}


def verify_gpu(sample: int = 20000, seed: int = 999) -> dict:
    """Check the GPU kernel against the CPU one bit-for-bit.

    The draws feed a threshold comparison, so anything short of exact equality could flip a guide's
    verdict — and GPU floating point is where that would silently happen. The construction is
    integer arithmetic plus one multiply-add in float64, so exactness is expected, not hoped for.
    """
    xp = _gpu()
    if xp is None:
        return {"available": False}
    rng = random.Random(seed)
    values = np.array([0, 1, 2, 2 ** 31, 2 ** 32 - 1]
                      + [rng.randrange(0, 2 ** 32) for _ in range(sample)], dtype=np.uint32)
    cpu_first, cpu_second = first_two_draws(values)
    gpu_first, gpu_second = first_two_draws_gpu(values, xp)
    return {
        "available": True,
        "checked": int(values.size),
        "first_draw_mismatches": int(np.count_nonzero(xp.asnumpy(gpu_first) != cpu_first)),
        "second_draw_mismatches": int(np.count_nonzero(xp.asnumpy(gpu_second) != cpu_second)),
    }


def first_two_draws(seeds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The first two ``random()`` values for each 32-bit seed, vectorised over the batch.

    ``seeds`` is a uint32 array of ``experiment_seed`` values. Returns (mh_draw, cut_draw) as float64
    arrays — draw 1 is the microhomology coin, draw 2 the cut coin, exactly the order stage 3 uses.

    The state is held as a (624, batch) uint32 array: init_by_array's two mixing passes are
    sequential in the *state index* (mt[i] depends on mt[i-1]) but independent across seeds, so each
    step is one array operation over the whole batch.
    """
    seeds = np.ascontiguousarray(seeds, dtype=np.uint32)
    batch = seeds.shape[0]
    mt = np.repeat(_INIT_CONST[:, None], batch, axis=1)
    over = np.errstate(over="ignore")
    over.__enter__()

    # Pass 1: 624 mixes folding in the key. CPython's key array for an int < 2**32 is a single word,
    # so `init_key[j] + j` is just the seed on every iteration (j stays 0).
    i = 1
    for _ in range(N):
        prev = mt[i - 1]
        mt[i] = (mt[i] ^ ((prev ^ (prev >> np.uint32(30))) * _MULT1)) + seeds
        i += 1
        if i >= N:
            mt[0] = mt[N - 1]
            i = 1

    # Pass 2: 623 further mixes, subtracting the index.
    for _ in range(N - 1):
        prev = mt[i - 1]
        mt[i] = (mt[i] ^ ((prev ^ (prev >> np.uint32(30))) * _MULT2)) - np.uint32(i)
        i += 1
        if i >= N:
            mt[0] = mt[N - 1]
            i = 1
    mt[0] = UPPER_MASK          # "MSB is 1; assuring non-zero initial array"

    # init leaves index == N, so the first draw triggers the twist; only words 0..3 are ever read,
    # and each needs mt[kk], mt[kk+1] and mt[kk+M] from the pre-twist state.
    out = np.empty((4, batch), dtype=np.uint32)
    for kk in range(4):
        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK)
        out[kk] = mt[kk + M] ^ (y >> np.uint32(1)) ^ np.where(
            (y & np.uint32(1)).astype(bool), MATRIX_A, np.uint32(0)).astype(np.uint32)

    for kk in range(4):
        y = out[kk]
        y = y ^ (y >> np.uint32(11))
        y = y ^ ((y << np.uint32(7)) & np.uint32(0x9D2C5680))
        y = y ^ ((y << np.uint32(15)) & np.uint32(0xEFC60000))
        out[kk] = y ^ (y >> np.uint32(18))

    # random() consumes two words: (a >> 5) * 2**26 + (b >> 6), scaled by 2**-53.
    first = ((out[0] >> np.uint32(5)).astype(np.float64) * _TWO26
             + (out[1] >> np.uint32(6)).astype(np.float64)) / _TWO53
    second = ((out[2] >> np.uint32(5)).astype(np.float64) * _TWO26
              + (out[3] >> np.uint32(6)).astype(np.float64)) / _TWO53
    over.__exit__(None, None, None)
    return first, second


def experiment_seeds(round_seeds, mutation: str, cas: str, guide: str,
                     start: int, strand: str) -> np.ndarray:
    """``stage3.experiment_seed`` for many round seeds against one design.

    The digest's low 32 bits are what ``int(hexdigest, 16) % 2**32`` selects, i.e. its last four
    bytes big-endian. sha256 is not vectorisable here — the round seed leads the key, so no prefix
    state can be shared — so this stays a Python loop and is the floor on batch cost.
    """
    suffix = f"|{mutation}|{cas}|{guide}|{start}|{strand}".encode()
    return np.fromiter(
        (int.from_bytes(hashlib.sha256(str(sd).encode() + suffix).digest()[-4:], "big")
         for sd in round_seeds),
        dtype=np.uint32, count=len(round_seeds))


def cut_fails(round_seeds: np.ndarray, seed_values: np.ndarray, cut_p: float,
              max_fail: int) -> np.ndarray | None:
    """Round seeds whose cut coin exceeds ``cut_p``, or None if there are more than ``max_fail``."""
    _mh, cut = first_two_draws(seed_values)
    failed = round_seeds[cut > cut_p]
    if failed.size > max_fail:
        return None
    return failed


def screen_guides(guides: list[str], round_seeds: np.ndarray, mutation: str, cas: str,
                  start: int, strand: str, cut_p_of, max_fail: int,
                  chunk_seeds: int = 32) -> dict[str, np.ndarray]:
    """Survivor-screening across a whole guide list: which guides fail at most ``max_fail`` seeds.

    Batching along the *guide* axis is what makes the kernel pay off — a single guide is only 900
    seeds, small enough that array overhead dominates (measured 11 us/pair against 1.3 us/pair at
    450k). So this walks the seed window in slices and, after each slice, drops the guides that have
    already exceeded their budget. That keeps every kernel call wide while preserving the early-out
    that makes a strict search affordable: with cut_p 0.99 and max_fail 0, ~99% of guides are gone
    within the first hundred seeds and are never hashed again.

    ``cut_p_of`` maps a guide to its cut probability (it varies with GC across a target's variants).
    Returns {guide: failed_seeds} for the survivors only.
    """
    alive = list(guides)
    fails: dict[str, list[int]] = {g: [] for g in alive}
    suffix = {g: f"|{mutation}|{cas}|{g}|{start}|{strand}".encode() for g in alive}
    cut_p = np.array([cut_p_of(g) for g in alive], dtype=np.float64)

    for begin in range(0, round_seeds.size, chunk_seeds):
        if not alive:
            break
        slice_seeds = round_seeds[begin:begin + chunk_seeds]
        prefixes = [str(int(sd)).encode() for sd in slice_seeds]
        # one flat batch of (guide x seed-slice) seed values
        values = np.empty(len(alive) * slice_seeds.size, dtype=np.uint32)
        at = 0
        for guide in alive:
            tail = suffix[guide]
            for prefix in prefixes:
                values[at] = int.from_bytes(
                    hashlib.sha256(prefix + tail).digest()[-4:], "big")
                at += 1
        _mh, cut = first_two_draws(values)
        cut = cut.reshape(len(alive), slice_seeds.size)

        survivors = []
        for index, guide in enumerate(alive):
            over_budget = cut[index] > cut_p[index]
            if over_budget.any():
                fails[guide].extend(int(x) for x in slice_seeds[over_budget])
                if len(fails[guide]) > max_fail:
                    continue          # dropped: never hashed for the remaining seeds
            survivors.append(guide)
        alive = survivors
        cut_p = np.array([cut_p_of(g) for g in alive], dtype=np.float64)

    return {g: np.array(sorted(fails[g]), dtype=np.int64) for g in alive}


def verify(sample: int = 20000, seed: int = 12345) -> dict:
    """Check the batched kernel against ``random.Random`` on a sample of 32-bit seeds.

    Compares the exact float bit patterns, not an approximate closeness: the draws feed a threshold
    comparison, so a one-ulp difference could flip a guide's verdict.
    """
    rng = random.Random(seed)
    values = np.array(
        [0, 1, 2, 2 ** 31, 2 ** 32 - 1]
        + [rng.randrange(0, 2 ** 32) for _ in range(sample)], dtype=np.uint32)
    got_first, got_second = first_two_draws(values)

    bad_first = bad_second = 0
    for index, value in enumerate(values):
        reference = random.Random(int(value))
        want_first = reference.random()
        want_second = reference.random()
        if want_first != float(got_first[index]):
            bad_first += 1
        if want_second != float(got_second[index]):
            bad_second += 1
    return {"checked": int(values.size), "first_draw_mismatches": bad_first,
            "second_draw_mismatches": bad_second,
            "exact": bad_first == 0 and bad_second == 0}
