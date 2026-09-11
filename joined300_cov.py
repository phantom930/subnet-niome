#!/usr/bin/env python3
"""joined300_cov.py — what the joined-300 layout costs in fleet coverage, per cell type.

Reads `joined300.json` (the per-cell band measurements from joined300.py) and turns the per-hotkey
band sizes into the only fleet quantity that matters: `|union of the eleven bands| / 900`, and
`1 - (1 - union/900)**3`, the chance a round's three seeds touch some sibling's band.

Two layouts at the same eleven hotkeys:

  joined 300   h0 on the whole joined space (band measured at width 300) + h1-h10 on rotated
               width-225 slices at stride 30 (band measured at width 225). Bands are independent
               draws inside their own windows -- verified in band_overlap.py, mean pairwise overlap
               0.56 against 0.58 expected -- so the union is computed exactly under independence
               rather than assumed disjoint, which is the whole point: eleven 12-seed bands cannot
               sit distinctly in 300 seeds.
  spread 900   the same eleven hotkeys on DISJOINT windows tiling 100-999 (width ~82). Disjoint
               windows give disjoint bands, so the union is the plain sum. The band at width 82 is
               taken from the width-100 measurement, which is conservative: CLAUDE.md's
               `narrow_width` row measures the band GROWING toward the window as the window
               narrows, so a width-82 window returns at least the width-100 band.

Coverage is not the whole answer -- a wider window also moves `weighted x fidelity`, and placement
is a threshold, not a linear function of frequency. That half is priced by fleet_price.py with
FP_BAND / FP_JOINED / FP_AH_WXF. This script answers only "how often does some hotkey spike".
"""
import json
import os
import sys

N_SLICE = 10          # h1-h10
SUB_WIDTH = 225       # JOINED_SUB_WIDTH
STRIDE = 30           # JOINED_STRIDE
JOINED = 300
N_HK = 11
NSEED = 900


def slice_membership():
    """How many of the ten rotated slices contain each index of the joined space."""
    counts = []
    for i in range(JOINED):
        n = sum(1 for h in range(N_SLICE) if (i - h * STRIDE) % JOINED < SUB_WIDTH)
        counts.append(n)
    return counts


def joined_union(b300, b225):
    """Expected |union| of h0's width-300 band and ten width-225 slice bands, under independence."""
    total = 0.0
    for k in slice_membership():
        p_miss = (1 - b300 / JOINED) * (1 - b225 / SUB_WIDTH) ** k
        total += 1 - p_miss
    return total


def p_hit(union):
    return 1 - (1 - union / NSEED) ** 3


def main(path="joined300.json"):
    rows = [r for r in json.load(open(path)) if "band" in r]
    by = {}
    for r in rows:
        by.setdefault(r["cell"], {})[r["arm"]] = r
    print("Fleet coverage from the measured bands — 11 hotkeys, two layouts\n")
    print("%-12s %-26s %-28s %s" % ("cell", "joined 300 (live)", "spread over 900", "cost"))
    print("%-12s %-26s %-28s %s" % ("", "b300 b225 -> union  P(hit)", "b100 x11 -> union  P(hit)", ""))
    for cell, arms in by.items():
        if not {"joined300", "joined225", "contig100"} <= set(arms):
            print(f"{cell:<12} incomplete ({sorted(arms)})")
            continue
        b300 = arms["joined300"]["band"]
        b225 = arms["joined225"]["band"]
        b100 = arms["contig100"]["band"]
        uj = joined_union(b300, b225)
        us = min(NSEED, b300 + N_SLICE * b100)      # h0 keeps its width-300 window in both layouts
        print("%-12s %4d %4d -> %5.1f  %5.1f%%   %4d x10+%d -> %5.1f  %5.1f%%   %+.1f%% P(hit)"
              % (cell, b300, b225, uj, 100 * p_hit(uj), b100, b300, us, 100 * p_hit(us),
                 100 * (p_hit(uj) / p_hit(us) - 1)))
    print("\nP(hit) is the chance at least one of the round's three seeds lands in some hotkey's "
          "band.\nIt is frequency only: what a hit is WORTH also moves with the window "
          "(see wxfid per arm).")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "joined300.json")
