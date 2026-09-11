#!/usr/bin/env python3
"""hdr_pool.py — STAGE 2 of the floor rebuild: the HDR-clean pool as a function of band width.

The rebuild's premise, from `floor_bound.py` and `band_floor.py`:

* the off-band floor CANNOT be raised by row distribution -- an oracle regressor on the most
  heterogeneous 250 rows the design space allows reaches cons 0.085-0.111 against our live
  0.095-0.113, because `is_hdr`'s nmae is a structural 0.997 and energy is clamped at 1.0 on every
  erythroid row. So an elevated floor needs `is_cut` PINNED, i.e. a wide cut-clean set.
* a wide cut-clean set alongside an HDR band is blocked by pool freedom: "only 82 of 60,000 bank
  guides are HDR-clean on the band and 80 are needed, so the Cas12a min-union has no freedom."
* `band_floor.py` prices narrowing the band to buy that freedom. Break-even floors: band 10 needs
  0.157 (371 of 900 cut-clean), band 8 needs 0.170 (457), band 7 needs 0.186 (562). Band 10 pays
  +63.5% at floor 0.237 and needs only two thirds of the clean set all-cut already reaches.

Everything in that pricing rides on ONE extrapolated number: a per-seed P(HDR) of 0.602, inferred
from the single (82 of 60,000 at band 13) datapoint, giving pool = 60000 * 0.602**B. That
extrapolation assumes seeds are independent ACROSS guides, which they are not -- a high-energy
guide is HDR-clean more often on every seed -- so the true curve should sit ABOVE it. This measures
it instead of assuming it.

No GPU and no build: `AH.build_bank` already stores each guide's HDR-failed seeds over the window
in `r["fails"]`, so the pool of any candidate band S is just

    pool(S) = #{r : fails(r) INTERSECT S = empty}

and the whole curve is set arithmetic over banks the live fleet has already built.

Two bands per width are reported, because they answer different questions:
  GREEDY   the widest-pool band of that width (add the seed costing the fewest survivors) -- the
           band a rebuild would deliberately choose, and the right input to stage 3's min-union.
  RANDOM   a uniformly sampled band of that width -- what an arbitrary window hands you, and the
           quantity the 0.602**B model actually describes.
"""
import glob
import json
import os
import sys

os.environ.setdefault("NIOME_INSTANCE", "hdr_pool")
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

import logging
logging.basicConfig(level=logging.ERROR)

import numpy as np

from niome_subnet.genomics import all_cut as AC

WIDTHS = [int(x) for x in os.environ.get("HP_WIDTHS", "4,6,7,8,10,12,13,16,20").split(",")]
NRAND = int(os.environ.get("HP_NRAND", "40"))
GROUP = int(os.environ.get("HP_GROUP", "42"))
NBANKS = int(os.environ.get("HP_NBANKS", "6"))


def load_banks(n):
    """The most recently built all-HDR Cas12a banks the fleet has on disk."""
    from niome_subnet.genomics import all_hdr as AH
    paths = sorted(glob.glob(os.path.join(AH.HDR_BANK_DIR, "cas12a-*.npz")),
                   key=os.path.getmtime, reverse=True)
    out = []
    for p in paths:
        if len(out) >= n:
            break
        try:
            rec = AC.load_bank(p)
        except Exception as e:
            print("  skip %s (%s)" % (os.path.basename(p), type(e).__name__))
            continue
        if not rec:
            continue
        out.append((os.path.basename(p), rec))
    return out


def screened_seeds(fail_sets):
    """The seeds the bank was actually screened over.

    Do NOT use range(min, max): the fleet's windows are JOINED (window_plan.py stitches rotated
    slices), so that range includes seeds no guide was ever tested on. A greedy then picks those
    unscreened gaps at zero cost and the pool never drops -- which is exactly what the first run of
    this script reported (pool flat at the full bank size for every width).

    Any seed that WAS screened is failed by ~30% of a 14k-guide bank, so it appears in some fail
    set with probability 1 - 0.7**14000. The union of the fail sets is therefore the screened set."""
    u = set()
    for f in fail_sets:
        u |= f
    return sorted(u)


def pool_of(fail_sets, S):
    if not S:
        return len(fail_sets)
    s = set(S)
    return sum(1 for f in fail_sets if not (f & s))


def greedy_band(fail_sets, seeds, width):
    """Add the seed that kills the fewest surviving guides. Returns (band, pool_after_each_add)."""
    alive = list(range(len(fail_sets)))
    band, curve = [], []
    remaining = set(seeds)
    for _ in range(width):
        best, best_keep = None, -1
        for sd in remaining:
            keep = sum(1 for i in alive if sd not in fail_sets[i])
            if keep > best_keep:
                best, best_keep = sd, keep
        if best is None:
            break
        band.append(best); remaining.discard(best)
        alive = [i for i in alive if best not in fail_sets[i]]
        curve.append(len(alive))
    return band, curve


def main():
    banks = load_banks(NBANKS)
    if not banks:
        print("no all-HDR banks on disk"); return
    print("=== STAGE 2: HDR-clean pool vs band width, over %d live banks ===\n" % len(banks))
    rng = np.random.default_rng(20260910)
    agg = {w: {"greedy": [], "rand": [], "q": []} for w in WIDTHS}
    for name, rec in banks:
        fail_sets = [set(int(x) for x in (r.get("fails") if r.get("fails") is not None else []))
                     for r in rec]
        seeds = screened_seeds(fail_sets)
        if not seeds:
            print("%s: no fails recorded, skipping" % name); continue
        W = len(seeds)
        contig = (seeds[-1] - seeds[0] + 1) == W
        rate = 1.0 - np.mean([len(f) for f in fail_sets]) / W
        print("%s  guides %6d  screened %d seeds in %d-%d (%s)  per-seed HDR-clean rate %.4f"
              % (name, len(rec), W, seeds[0], seeds[-1],
                 "contiguous" if contig else "JOINED", rate))
        gb, curve = greedy_band(fail_sets, seeds, max(WIDTHS))
        for w in WIDTHS:
            if w <= len(curve):
                agg[w]["greedy"].append(curve[w - 1])
            rs = []
            for _ in range(NRAND):
                S = rng.choice(seeds, size=min(w, W), replace=False).tolist()
                rs.append(pool_of(fail_sets, S))
            agg[w]["rand"].append(float(np.mean(rs)))
            agg[w]["q"].append(len(rec) * rate ** w)
        print("   greedy band pool by width: %s"
              % "  ".join("%d:%d" % (w, curve[w - 1]) for w in WIDTHS if w <= len(curve)))

    print("\n=== pooled across banks ===")
    print("  the extrapolation this replaces: pool = 60000 * 0.602**B")
    print("\n  width   GREEDY pool (median)   RANDOM pool (median)   independence model   greedy freedom @g%d" % GROUP)
    for w in WIDTHS:
        g = agg[w]["greedy"]; r = agg[w]["rand"]; q = agg[w]["q"]
        if not g:
            continue
        gm = float(np.median(g))
        print("   %3d    %14.0f       %14.0f       %14.0f       %10.1fx"
              % (w, gm, float(np.median(r)), float(np.median(q)), gm / GROUP))
    print("\n  break-even floors from band_floor.py, and the clean set each needs:")
    for w, f, need in ((10, 0.157, 371), (8, 0.170, 457), (7, 0.186, 562), (6, 0.192, 602)):
        g = agg.get(w, {}).get("greedy") or []
        if g:
            print("   band %2d: pool %5.0f (%5.1fx group %d) -> must reach floor %.3f = %d of 900 clean"
                  % (w, float(np.median(g)), float(np.median(g)) / GROUP, GROUP, f, need))
    json.dump({str(w): {k: agg[w][k] for k in agg[w]} for w in WIDTHS},
              open("hdr_pool.json", "w"), indent=1)
    print("\n  wrote hdr_pool.json")


if __name__ == "__main__":
    main()
