#!/usr/bin/env python3
"""joined300_price.py — price ONE all-HDR window layout against a cell type's real fields.

fleet_price.py sweeps all-cut/all-HDR compositions at a fixed layout. This asks the other question
with the same machinery and the same real fields: hold the composition at eleven all-HDR hotkeys
(what the fleet ships) and vary the LAYOUT — the band space each hotkey searches (`FP_JOINED`), the
band that space returns (`FP_BAND`, measured by joined300.py) and the `weighted x fidelity` that
comes with it (`FP_AH_WXF`, measured on the same build).

Every parameter is read from the environment by fleet_price at import, so one process prices one
layout:

    FP_CELL=K562 FP_JOINED=300 FP_BAND=12 FP_AH_WXF=216.0 python joined300_price.py
"""
import os
import sys

import numpy as np

import fleet_price as fp


def main():
    fs = fp.fields()
    rng = np.random.default_rng(fp.RNG)
    vals = [fp.simulate(f, (0, 11), rng, fp.TRIALS) for _tid, f in fs]
    print("%-12s space %-4d band %-3d wxfid %6.1f   fields %-3d  E[share] %.5f"
          % (fp.CELL, fp.JOINED, fp.AH_BAND, fp.AH_WXF, len(fs), float(np.mean(vals))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
