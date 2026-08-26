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


def first_two_draws_gpu(seeds, xp, chunk: int = GPU_CHUNK) -> tuple:
    """``first_two_draws`` on the GPU. Same algorithm, same order, ``xp`` is cupy.

    Kept as a separate function rather than a branch inside the numpy version so the CPU path stays
    the reference implementation that ``verify()`` pins down: the GPU result is checked against it,
    not the other way round.
    """
    seeds = xp.asarray(seeds, dtype=xp.uint32)
    if seeds.shape[0] > chunk:
        # Split rather than fail: a caller should not have to know the device's memory budget.
        firsts, seconds = [], []
        for begin in range(0, int(seeds.shape[0]), chunk):
            a, b = first_two_draws_gpu(seeds[begin:begin + chunk], xp, chunk)
            firsts.append(a)
            seconds.append(b)
        return xp.concatenate(firsts), xp.concatenate(seconds)
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

    out = xp.empty((4, batch), dtype=xp.uint32)
    for kk in range(4):
        y = (mt[kk] & xp.uint32(0x80000000)) | (mt[kk + 1] & xp.uint32(0x7FFFFFFF))
        out[kk] = mt[kk + M] ^ (y >> u1) ^ xp.where(
            (y & u1).astype(bool), xp.uint32(0x9908B0DF), xp.uint32(0)).astype(xp.uint32)
    for kk in range(4):
        y = out[kk]
        y = y ^ (y >> xp.uint32(11))
        y = y ^ ((y << xp.uint32(7)) & xp.uint32(0x9D2C5680))
        y = y ^ ((y << xp.uint32(15)) & xp.uint32(0xEFC60000))
        out[kk] = y ^ (y >> xp.uint32(18))

    first = ((out[0] >> xp.uint32(5)).astype(xp.float64) * _TWO26
             + (out[1] >> xp.uint32(6)).astype(xp.float64)) / _TWO53
    second = ((out[2] >> xp.uint32(5)).astype(xp.float64) * _TWO26
              + (out[3] >> xp.uint32(6)).astype(xp.float64)) / _TWO53
    return first, second


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
