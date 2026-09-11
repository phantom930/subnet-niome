#!/usr/bin/env python3
"""band_floor.py -- the indifference curve between band WIDTH and off-band FLOOR.

Why this is the go/no-go for a rebuild. Three things are now measured:

 1. `floor_bound.py` -- the floor cannot be raised by row DISTRIBUTION. An oracle regressor on the
    most heterogeneous 250 rows the design space allows reaches cons 0.085-0.111 against our live
    0.095-0.113, because avg_nmae rises faster than avg_r2 (`is_hdr`'s nmae is 0.997, a structural
    fair coin). So an elevated floor requires PINNING `is_cut`, i.e. a wide cut-clean set.
 2. A wide cut-clean set alongside an HDR band is blocked by pool freedom, and the blocker is
    parametric: only 82 of 60,000 bank guides are HDR-clean on a 13-seed band (implied per-seed
    P(HDR) = 0.602) against the 80 a group needs. The pool goes as 60000*0.602^B, so it is 1036 at
    band 8 and 1721 at band 7 -- and all-cut min-unions group 42 from a pool of ~1749 to reach a
    cut-fail union of 341/900. Comparable freedom appears at band 6-7.
 3. Narrowing the band costs spike frequency: P(k>=1) = 1-(1-B/900)^3 falls 0.0427 -> 0.0264 ->
    0.0199 at B = 13 -> 8 -> 6.

So the rebuild is a single trade -- give up band width, buy the floor -- and this prices it against
the real post-shift fields BEFORE any GPU time. What it produces is the floor a given band width
would have to reach to beat the shipped (band 12, floor 0.101) fleet. If band 7 needs a floor the
conjunction cannot deliver, the rebuild is not worth starting.

Reuses floor_price.simulate, so the field handling, the per-cell weighted x fidelity and the
correlated-floor/independent-band modelling are identical.
"""
import numpy as np

import floor_price as FP


def main():
    fs = FP.fields()
    print("fields: %d post-%s erythroid rounds, our hotkeys excluded by ss58 address"
          % (len(fs), FP.SINCE))
    print("per-cell weighted x fidelity: %s"
          % "  ".join("%s %.1f" % kv for kv in sorted(FP.OUR_WXF_CELL.items())))
    bands = [13, 12, 10, 8, 7, 6, 4]
    floors = [0.101, 0.15, 0.198, 0.237]
    trials = 3000
    print("%d trials per round per cell\n" % trials)

    def price(f, band):
        vals = []
        for _tid, cell, _d, fld in fs:
            rng = np.random.default_rng(20260910)
            vals.append(FP.simulate(fld, f, FP.OUR_WXF_CELL[cell], band, rng, trials))
        return float(np.mean(vals))

    base = price(0.101, 12)
    print("shipped: band 12, floor 0.101 -> E[share] %.5f\n" % base)
    print("=== E[share] over the (band, floor) grid; %+ is against shipped ===")
    print("  band   pool@band " + "".join("  f=%-16s" % f for f in floors))
    import math
    q = 0.6021
    for b in bands:
        cells = []
        for f in floors:
            v = price(f, b)
            cells.append("%.5f (%+6.1f%%)" % (v, 100 * (v / base - 1)))
        print("  %4d   %8.0f  %s" % (b, 60000 * q ** b, "  ".join("%-18s" % c for c in cells)))

    print("\n=== break-even floor for each band width (linear interpolation) ===")
    print("  the floor a rebuild at that band must reach just to MATCH the shipped fleet")
    for b in bands:
        vs = [(f, price(f, b)) for f in floors]
        hit = None
        for (f0, v0), (f1, v1) in zip(vs, vs[1:]):
            if (v0 - base) * (v1 - base) <= 0 and v1 != v0:
                hit = f0 + (f1 - f0) * (base - v0) / (v1 - v0)
                break
        if vs[0][1] >= base:
            print("  band %2d: already at or above shipped at floor 0.101" % b)
        elif hit is None:
            print("  band %2d: does not reach shipped even at floor %.3f (%.5f vs %.5f)"
                  % (b, floors[-1], vs[-1][1], base))
        else:
            print("  band %2d: break-even floor %.3f   %s" % (b, hit,
                  "REACHABLE-ish, all-cut's pinned value is 0.237" if hit <= 0.237
                  else "ABOVE all-cut's 0.237 -- not reachable"))


if __name__ == "__main__":
    main()
