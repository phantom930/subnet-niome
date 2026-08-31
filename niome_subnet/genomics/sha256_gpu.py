"""sha256_gpu.py — device-side sha256 of the stage-3 experiment keys.

Once the Mersenne-Twister draws moved to the GPU, hashing became the binding cost of a seed-agnostic
scan: the key is ``f"{round_seed}|{mutation}|{cas}|{guide}|{start}|{strand}"`` and one is needed per
(guide, seed) pair, millions per build. This computes them on the device.

The trick that makes it cheap is not the hash itself but what gets uploaded. A naive port would send
every message — O(guides x seeds) bytes. Instead the key splits at the first separator: the round
seed varies down one axis, and everything after it is constant per guide. So only a per-seed digit
table and a per-guide suffix table cross the bus (O(guides + seeds)), and each thread concatenates
its own message in registers before hashing.

Only the low 32 bits of the digest are returned, because that is all ``experiment_seed`` keeps
(``int(hexdigest, 16) % 2**32`` is the last four bytes big-endian, i.e. the final state word).
"""

from __future__ import annotations

import hashlib

import numpy as np

# Longest message the kernel will assemble. The real keys run ~55-80 bytes; 119 is the most that
# still pads into two 64-byte blocks, and the host asserts against it rather than truncating.
MAX_MSG = 119
_BUF = 128

_SOURCE = r"""
#define ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))

__constant__ unsigned int K256[64] = {
  0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
  0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
  0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
  0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
  0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
  0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
  0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
  0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
};

extern "C" __global__
void sha256_low32(
    const unsigned char* __restrict__ seed_bytes, const int* __restrict__ seed_lens, int seed_stride,
    const unsigned char* __restrict__ suf_bytes,  const int* __restrict__ suf_lens,  int suf_stride,
    int n_seeds, long long total, unsigned int* __restrict__ out)
{
    long long tid = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (tid >= total) return;
    int g = (int)(tid / n_seeds);
    int s = (int)(tid % n_seeds);

    unsigned char msg[128];
    int len = 0;
    int sl = seed_lens[s];
    const unsigned char* sp = seed_bytes + (long long)s * seed_stride;
    for (int i = 0; i < sl; ++i) msg[len++] = sp[i];
    int ul = suf_lens[g];
    const unsigned char* up = suf_bytes + (long long)g * suf_stride;
    for (int i = 0; i < ul; ++i) msg[len++] = up[i];

    /* sha256 padding: 0x80, zeros, then the bit length big-endian in the last 8 bytes. */
    int nblocks = (len + 9 + 63) / 64;
    int padded = nblocks * 64;
    msg[len] = 0x80;
    for (int i = len + 1; i < padded; ++i) msg[i] = 0;
    unsigned long long bits = (unsigned long long)len * 8ULL;
    for (int i = 0; i < 8; ++i) msg[padded - 1 - i] = (unsigned char)((bits >> (8 * i)) & 0xffULL);

    unsigned int h0=0x6a09e667u,h1=0xbb67ae85u,h2=0x3c6ef372u,h3=0xa54ff53au,
                 h4=0x510e527fu,h5=0x9b05688cu,h6=0x1f83d9abu,h7=0x5be0cd19u;

    for (int blk = 0; blk < nblocks; ++blk) {
        unsigned int w[64];
        const unsigned char* p = msg + blk * 64;
        for (int i = 0; i < 16; ++i)
            w[i] = ((unsigned int)p[i*4] << 24) | ((unsigned int)p[i*4+1] << 16)
                 | ((unsigned int)p[i*4+2] << 8) | (unsigned int)p[i*4+3];
        for (int i = 16; i < 64; ++i) {
            unsigned int s0 = ROTR(w[i-15],7) ^ ROTR(w[i-15],18) ^ (w[i-15] >> 3);
            unsigned int s1 = ROTR(w[i-2],17) ^ ROTR(w[i-2],19) ^ (w[i-2] >> 10);
            w[i] = w[i-16] + s0 + w[i-7] + s1;
        }
        unsigned int a=h0,b=h1,c=h2,d=h3,e=h4,f=h5,g2=h6,hh=h7;
        for (int i = 0; i < 64; ++i) {
            unsigned int S1 = ROTR(e,6) ^ ROTR(e,11) ^ ROTR(e,25);
            unsigned int ch = (e & f) ^ ((~e) & g2);
            unsigned int t1 = hh + S1 + ch + K256[i] + w[i];
            unsigned int S0 = ROTR(a,2) ^ ROTR(a,13) ^ ROTR(a,22);
            unsigned int mj = (a & b) ^ (a & c) ^ (b & c);
            unsigned int t2 = S0 + mj;
            hh=g2; g2=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
        }
        h0+=a; h1+=b; h2+=c; h3+=d; h4+=e; h5+=f; h6+=g2; h7+=hh;
    }
    /* experiment_seed keeps only the low 32 bits, i.e. the final digest word. */
    out[tid] = h7;
}
"""

_KERNEL = None


def _kernel(xp):
    global _KERNEL
    if _KERNEL is None:
        _KERNEL = xp.RawKernel(_SOURCE, "sha256_low32")
    return _KERNEL


def pack(strings: list[bytes]) -> tuple[np.ndarray, np.ndarray]:
    """Byte matrix + lengths for a ragged list, so the device can index it flat."""
    if not strings:
        return np.zeros((0, 1), dtype=np.uint8), np.zeros(0, dtype=np.int32)
    stride = max(len(s) for s in strings)
    out = np.zeros((len(strings), stride), dtype=np.uint8)
    for row, value in enumerate(strings):
        out[row, :len(value)] = np.frombuffer(value, dtype=np.uint8)
    return out, np.array([len(s) for s in strings], dtype=np.int32)


def experiment_seeds_gpu_packed(seed_digits: list[bytes], d_sufb, d_ulen, xp,
                                block: int = 256):
    """``experiment_seeds_gpu`` with the suffix table already packed and resident on the device.

    The suffix table is constant for a target while the seed slice moves, so packing and
    re-uploading it per slice is pure waste — and it dominated a large screen (37k suffixes x ~23
    slices x 915 targets). The caller packs once and slices the device array to pick survivors,
    which is a device-side gather rather than host work plus a PCIe copy.
    """
    seed_bytes, seed_lens = pack(seed_digits)
    longest = int(seed_bytes.shape[1]) + int(d_sufb.shape[1])
    if longest > MAX_MSG:
        raise ValueError(f"key of {longest} bytes exceeds the kernel's {MAX_MSG}-byte buffer")
    d_seed = xp.asarray(seed_bytes)
    d_slen = xp.asarray(seed_lens)
    n_seeds = len(seed_digits)
    total = int(d_sufb.shape[0]) * n_seeds
    out = xp.empty(total, dtype=xp.uint32)
    grid = (total + block - 1) // block
    _kernel(xp)((grid,), (block,), (
        d_seed, d_slen, np.int32(seed_bytes.shape[1]),
        d_sufb, d_ulen, np.int32(d_sufb.shape[1]),
        np.int32(n_seeds), np.int64(total), out))
    return out


def experiment_seeds_gpu(seed_digits: list[bytes], suffixes: list[bytes], xp,
                         block: int = 256):
    """Low 32 bits of sha256(seed_digits[s] + suffixes[g]) for every (g, s), on the device.

    Returns a uint32 device array of shape (len(suffixes) * len(seed_digits),), guide-major — the
    same order ``screen_guides`` reshapes to (guides, seeds).
    """
    longest = max(len(a) for a in seed_digits) + max(len(b) for b in suffixes)
    if longest > MAX_MSG:
        raise ValueError(f"key of {longest} bytes exceeds the kernel's {MAX_MSG}-byte buffer")

    seed_bytes, seed_lens = pack(seed_digits)
    suf_bytes, suf_lens = pack(suffixes)
    d_seed = xp.asarray(seed_bytes)
    d_sufb = xp.asarray(suf_bytes)
    d_slen = xp.asarray(seed_lens)
    d_ulen = xp.asarray(suf_lens)

    n_seeds = len(seed_digits)
    total = len(suffixes) * n_seeds
    out = xp.empty(total, dtype=xp.uint32)
    grid = (total + block - 1) // block
    _kernel(xp)((grid,), (block,), (
        d_seed, d_slen, np.int32(seed_bytes.shape[1]),
        d_sufb, d_ulen, np.int32(suf_bytes.shape[1]),
        np.int32(n_seeds), np.int64(total), out))
    return out


def verify(xp, guides: int = 40, seeds: int = 60) -> dict:
    """Check the device hash against hashlib on a grid of realistic keys."""
    seed_digits = [str(s).encode() for s in range(100, 100 + seeds)]
    suffixes = [f"|NC_000011.10:g.5226521G>A|Cas12a|{'ACGT' * 5 + 'GCA'}|{5225490 + i}|-".encode()
                for i in range(guides)]
    got = xp.asnumpy(experiment_seeds_gpu(seed_digits, suffixes, xp))

    bad = 0
    for gi, suffix in enumerate(suffixes):
        for si, digits in enumerate(seed_digits):
            want = int.from_bytes(hashlib.sha256(digits + suffix).digest()[-4:], "big")
            if int(got[gi * len(seed_digits) + si]) != want:
                bad += 1
    return {"checked": len(suffixes) * len(seed_digits), "mismatches": bad, "exact": bad == 0}
