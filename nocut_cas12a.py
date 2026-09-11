#!/usr/bin/env python3
"""nocut_cas12a.py — the no-cut band from Cas12a guides driven to LOW ENERGY, width 75, K562.

`nocut_mf.py` measured the no-cut band at the shipped Cas12a config and got 2 (buildable) / 3
(Cas12a-only, unbuildable because the Cas9 half cannot no-cut). Both were limited by the same
quantity: P(no cut) = 1 - cut_p, and at the shipped band (gc 0.40-0.60, d<=400) energy sits at
0.98 so cut_p is 0.957 and P(no cut) only 0.043.

`cut_p = min(0.99, 0.78 + 0.18*energy)` bottoms out at **0.78** when energy reaches 0, which is
P(no cut) = 0.220 -- 5x the shipped rate. Energy is
`clamp(accessibility * (1.8*gc + 0.6*exp(-d/1500) + region_offset), 0, 1)`, so on K562
(accessibility 0.77) reaching it needs BOTH low GC and large distance; at d <= 400 the distance
term alone contributes 0.46 and pins energy above 0.35 whatever the GC. Measured reachable regimes:

    gc band      max_d   mean gc   mean energy   mean cut_p   P(nocut)   best cut_p
    0.00-0.25      400     0.222         0.712       0.9082     0.0918       0.8641
    0.00-0.25     3000     0.221         0.480       0.8664     0.1336       0.7966   <- this module
    0.40-0.60      400     0.451         0.984       0.9572     0.0428       0.9437   <- shipped

so gc 0.00-0.25 at d 3000 triples P(no cut) and puts the best guides at cut_p 0.797, essentially
the 0.78 floor. The mismatch budget is 3, which is what stops GC going lower.

The band is then a pure bank-size question -- a bank of N holds about `N * p**k` guides that all
skip the same k seeds -- so GROUP SIZE is a real variable here and not a free parameter: every row
of a Cas12a-only submission must come from the group, so filling 250 rows needs 250 coincident
guides, against the 42-80 the mixed constructions need.

**What this costs, and it is not small.** A submission of only Cas12a rows leaves stage 5's
cas-coverage entropy at zero and its joint coverage at 0.667, and `geometric_mean` floors zeros at
1e-9 rather than returning 0 -- so `distribution_fidelity_score` lands near **0.029** against the
~0.89 a mixed build gets. The band measured here is real; a submission built only from it is worth
roughly a thirtieth of its weighted score.

    python nocut_cas12a.py
    NCC_MF=60,65 NCC_GROUPS=80,250 python nocut_cas12a.py
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "ncc")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank   # noqa: E402
from allhdr_mf import pick                                  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402

CELL = os.getenv("NCC_CELL", "K562")
WINDOW = tuple(int(x) for x in os.getenv("NCC_WINDOW", "100-174").split("-"))
MFS = [int(x) for x in os.getenv("NCC_MF", "50,55,60,63,65,68").split(",")]
GROUPS = [int(x) for x in os.getenv("NCC_GROUPS", "42,80,125,170,250").split(",")]
GC = tuple(float(x) for x in os.getenv("NCC_GC", "0.00,0.25").split(","))
MAXD = int(os.getenv("NCC_MAXD", "3000"))
VARIANTS = int(os.getenv("NCC_VARIANTS", "4000"))
LOAD = int(os.getenv("NCC_LOAD", "300000"))
OUT = os.getenv("NCC_OUT", "nocut_cas12a.json")


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = WINDOW
    span = hi - lo + 1
    base = AC.config_for(CELL) or AC.AllCutConfig()
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  window {lo}-{hi} (span {span})"
          f"  rule nocut  gc {GC}  max_distance {MAXD}  variants {VARIANTS}\n")
    print(f"{'mf':>4}{'needs':>7}{'raw bank':>10}{'loaded':>8}{'nocut/guide':>14}"
          + "".join(f"{('g'+str(g)):>7}" for g in GROUPS) + f"{'s':>7}")
    out = []
    for mf in MFS:
        t0 = time.monotonic()
        cfg = _dc.replace(base, start_seed=lo, end_seed=hi, cas12a_max_fail=mf, rule="nocut",
                          cas12a_gc=GC, max_distance=MAXD, variants=VARIANTS)
        rec = {"cell": CELL, "window": [lo, hi], "span": span, "mf": mf, "needs": span - mf,
               "gc": list(GC), "max_distance": MAXD, "variants": VARIANTS,
               "task": (task.get("task_id") or task["id"])[:8]}
        path = os.path.join(AC.BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
            MT.free_gpu_memory()
            rec["raw_bank"] = len(bank)
            if not bank:
                rec["empty"] = True
                print(f"{mf:>4}{span-mf:>7}{0:>10}    nothing qualifies"
                      f"{time.monotonic()-t0:>7.0f}")
                out.append(rec); json.dump(out, open(OUT, "w"), indent=1); continue
            save_bank(path, bank)
        records = load_bank(path, limit=LOAD)
        rec["loaded_bank"] = len(records)
        nocuts = np.asarray([span - len(r["fails"]) for r in records])
        rec["nocut_per_guide"] = {"max": int(nocuts.max()), "mean": round(float(nocuts.mean()), 2)}
        line = (f"{mf:>4}{span-mf:>7}{rec.get('raw_bank', 'cached'):>10}{len(records):>8}"
                f"{rec['nocut_per_guide']['mean']:>8.1f}/{rec['nocut_per_guide']['max']:<5}")
        rec["bands"] = {}
        for g in GROUPS:
            if len(records) < g:
                rec["bands"][g] = None
                line += f"{'-':>7}"
                continue
            sel = FG.FastGreedy(records, window_lo=lo, window_hi=hi)
            idx, _u = sel.best(g, restarts=12)
            bad = set()
            for i in idx:
                bad.update(int(x) for x in records[i]["fails"])
            clean = sorted(set(range(lo, hi + 1)) - bad)
            rec["bands"][g] = {"union": len(bad), "clean": len(clean), "seeds": clean}
            line += f"{len(clean):>7}"
        rec["elapsed_s"] = round(time.monotonic() - t0, 1)
        print(line + f"{rec['elapsed_s']:>7.0f}")
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}   (columns under g<N> are the CLEAN BAND for that group size)")


if __name__ == "__main__":
    main()
