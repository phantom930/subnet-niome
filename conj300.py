#!/usr/bin/env python3
"""conj300.py — all-cut over 300 seeds, HDR band inside a nested 100, and no-cut rows to pay for it.

The construction, at group 80 (Cas12a 80 / Cas9 170):

    cut window    100-399 (300 seeds), cas12a_max_fail scaled 100 -> 33
    band          k seeds chosen inside 100-199, a sub-window of the cut window
    Cas9 rows     HDR on every band seed, cut over the clean set  -> is_cut True at the band
    Cas12a rows   NO-cut on every band seed, cut elsewhere        -> is_cut False at the band

Two things make this different from conj_test.py, and both come from the same observation: the
binding constraint on every band so far has been the *conditional Cas9 fill*, which dies as
P(HDR)**k.

1. **No-cut Cas12a rows are free.** The band lies inside the cut window, so a row cannot both cut
   everywhere and skip the band — but the all-cut bank already stores each guide's FAIL set, and for
   rule ``cut`` a fail IS a no-cut seed. Guides whose fail set contains the band are exactly the
   no-cut rows, at zero extra screening cost.
2. **They substitute for Cas9.** Shifting n rows from Cas9 to Cas12a leaves 170-n Cas9 rows needing
   HDR on the band, which is the requirement that ran out at 176-of-208. That is the whole idea
   under test: buy band width by needing fewer Cas9 candidates.

The cost is that the band seed is no longer a pure pin. Outcome is now a function of cas system —
Cas9 cuts and repairs, Cas12a does not cut — which stage 4 can only recover through gc-over-length
and energy, measured at **0.9483** rather than 1.0000. A mixture NOT aligned to cas measured
0.31-0.41, so the alignment must stay clean: every Cas12a row no-cut, every Cas9 row HDR.

    python conj300.py <task_id> [k1,k2,...] [n1,n2,...]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else "83f430e9-1c9a-48b4-b594-e17477a046ea"
KS = [int(x) for x in (_ARGV[2].split(",") if len(_ARGV) > 2 else ["2", "4", "8"])]
NS = [int(x) for x in (_ARGV[3].split(",") if len(_ARGV) > 3 else ["0", "40", "80"])]
os.environ["NIOME_INSTANCE"] = "conj300"

import dataclasses as _dc  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import random  # noqa: E402
import statistics as st  # noqa: E402
import time  # noqa: E402
from collections import Counter  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import fastgreedy as FG  # noqa: E402
from niome_subnet.genomics import mt19937 as MT  # noqa: E402
from conj_test import hdr_compliance  # noqa: E402
from sd_task import score, task_content  # noqa: E402

CUT_LO, CUT_HI = 100, 399          # all-cut window (width 300)
BAND_LO, BAND_HI = 100, 199        # the band is chosen inside this nested width-100 window
GROUP = 80                         # Cas12a rows before substitution; Cas9 takes 250 - GROUP
ROWS = 250
N_SAMPLE = int(os.getenv("C3_SAMPLE", "8"))
OUT = os.getenv("C3_OUT", "conj300.json")


def band_by_nocut(records, k):
    """The k seeds inside the band window that the MOST bank guides fail to cut on.

    Bank fails are no-cut seeds, so the seeds most often failed are the cheapest place to put a
    no-cut band. Greedy on the intersection, because the rows must skip every band seed, not any.
    """
    alive = set(range(len(records)))
    fails = [set(int(x) for x in r["fails"]) for r in records]
    band = []
    for _ in range(k):
        best, keep = None, None
        for seed in range(BAND_LO, BAND_HI + 1):
            if seed in band:
                continue
            got = {i for i in alive if seed in fails[i]}
            if keep is None or len(got) > len(keep):
                best, keep = seed, got
        if best is None:
            break
        band.append(best)
        alive = keep
    return band, alive


def regimes(rows, contract, reference, cell_types, band, clean, dirty, rng):
    """Consistency in each of the four regimes a round's seed can land in."""
    out = {}
    out["band"] = st.mean([score(rows, contract, reference, cell_types, seed=s)["consistency"]
                           for s in band[:min(3, len(band))]]) if band else None
    pools = {
        "clean": [s for s in clean if s not in set(band)],
        "dirty": sorted(dirty),
        "off": [s for s in range(100, 1000) if not (CUT_LO <= s <= CUT_HI)],
    }
    for name, pool in pools.items():
        if not pool:
            out[name] = None
            continue
        out[name] = st.mean([score(rows, contract, reference, cell_types, seed=s)["consistency"]
                             for s in rng.sample(pool, min(N_SAMPLE, len(pool)))])
    return out


def main():
    task, contract, reference = task_content(TASK)
    cell = contract.get("cell_type")
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AC.config_for(cell) or AC.AllCutConfig()
    span = CUT_HI - CUT_LO + 1
    cfg = _dc.replace(base, start_seed=CUT_LO, end_seed=CUT_HI, group_size=GROUP,
                      cas12a_max_fail=max(1, round(base.cas12a_max_fail * span / 900)))
    print(f"task {task['id'][:8]}  {cell}  cut {CUT_LO}-{CUT_HI} (span {span}, mf "
          f"{cfg.cas12a_max_fail})  band from {BAND_LO}-{BAND_HI}  group {GROUP}\n")

    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        print("building the width-300 all-cut bank")
        bank = AC.build_bank(contract, reference, cell_types, G.build_context(contract, reference,
                             cell_types), G.enumerate_sites(G.build_context(contract, reference,
                             cell_types), 3000, (20, 23)), cfg)
        if not bank:
            raise SystemExit("bank scan produced nothing")
        AC.save_bank(path, bank)
    records = AC.load_bank(path)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    fails_of = [set(int(x) for x in r["fails"]) for r in records]
    print(f"bank: {len(records)} Cas12a guides, mean fails "
          f"{st.mean(len(f) for f in fails_of):.1f} of {span}\n")

    rng = random.Random(5)
    out = {}
    print(f"  {'k':>3} {'n':>4} {'nocut pool':>10} {'hdr9 pool':>9} {'c12':>4} {'c9':>4} "
          f"{'rows':>5} {'cells':>5} {'weighted':>9} {'band':>7} {'clean':>7} {'dirty':>7} "
          f"{'off':>7} {'E[cons]':>8} {'E[final]':>9}")
    for k in KS:
        band, nocut_pool = band_by_nocut(records, k)
        if len(band) < k:
            print(f"  {k:>3}  band of {k} unreachable inside {BAND_LO}-{BAND_HI}")
            continue
        # Cas9 rows must repair by HDR on every band seed. Scanned over the clean set below, then
        # filtered — this is the pool the substitution is meant to economise on.
        for n in NS:
            n_c12, n_c9 = GROUP + n, ROWS - GROUP - n
            pool = [records[i] for i in sorted(nocut_pool)]
            if len(pool) < n_c12:
                print(f"  {k:>3} {n:>4} {len(pool):>10}  no-cut pool below {n_c12} Cas12a rows")
                out[f"k{k}n{n}"] = {"reason": f"nocut pool {len(pool)}"}
                continue
            sel = FG.FastGreedy(pool, window_lo=CUT_LO, window_hi=CUT_HI)
            idx, _u = sel.best(n_c12, restarts=8)
            group = [pool[i] for i in idx]
            bad = set()
            for rec in group:
                bad.update(int(x) for x in rec["fails"])
            bad -= set(band)                 # band failures are intentional, not coverage loss
            clean = sorted(set(range(CUT_LO, CUT_HI + 1)) - bad - set(band))
            cas9 = AC.scan_cas9(np.array(clean, dtype=np.int64), contract, cell_types, ctx, sites,
                                cfg, n_c9)
            if cas9:
                ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
                cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
            if len(cas9) < n_c9:
                print(f"  {k:>3} {n:>4} {len(pool):>10} {len(cas9):>9}  Cas9 half short of {n_c9}")
                out[f"k{k}n{n}"] = {"reason": f"cas9 {len(cas9)} of {n_c9}",
                                    "nocut_pool": len(pool)}
                continue
            rows = AC.assemble(group, cas9, contract, ctx, _dc.replace(cfg, group_size=n_c12),
                               ROWS)
            cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
            mix = Counter(r["cas_system"] for r in rows)
            reg = regimes(rows, contract, reference, cell_types, band, clean, bad, rng)
            n_off = 900 - span
            e_cons = ((len(band) * (reg["band"] or 0.0)
                       + len(clean) * (reg["clean"] or 0.0)
                       + len(bad) * (reg["dirty"] or 0.0)
                       + n_off * (reg["off"] or 0.0)) / 900.0)
            base_s = score(rows, contract, reference, cell_types, seed=band[0])
            e_final = base_s["weighted"] * e_cons * base_s["fidelity"]
            out[f"k{k}n{n}"] = {"band": band, "nocut_pool": len(pool), "cas9_pool": len(cas9),
                                "cas12a_rows": mix.get("Cas12a", 0), "cas9_rows": mix.get("Cas9", 0),
                                "clean": len(clean), "dirty": len(bad), "cells": cells,
                                "weighted": base_s["weighted"], "fidelity": base_s["fidelity"],
                                **{f"cons_{x}": reg[x] for x in reg},
                                "e_cons": e_cons, "e_final": e_final}
            print(f"  {k:>3} {n:>4} {len(pool):>10} {len(cas9):>9} {mix.get('Cas12a',0):>4} "
                  f"{mix.get('Cas9',0):>4} {len(rows):>5} {cells:>5} {base_s['weighted']:>9.1f} "
                  f"{(reg['band'] or float('nan')):>7.4f} {(reg['clean'] or float('nan')):>7.4f} "
                  f"{(reg['dirty'] or float('nan')):>7.4f} {(reg['off'] or float('nan')):>7.4f} "
                  f"{e_cons:>8.4f} {e_final:>9.2f}")
            MT.free_gpu_memory()

    print(f"\n  references on this cell type: pure all-cut over 100-999 E[cons] 0.2163 / "
          f"E[final] 57.24; all-HDR E[cons] ~0.15")
    json.dump({"task": TASK, "cell": cell, "cut_window": [CUT_LO, CUT_HI],
               "band_window": [BAND_LO, BAND_HI], "group": GROUP, "bank": len(records),
               "arms": out}, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
