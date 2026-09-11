#!/usr/bin/env python3
"""nocut_mf.py — the NO-CUT clean band against `main_max_fail`, at width 75 on K562.

The "nocut" rule (mt19937.RULE_SPECS) is the mirror of "cut": every row must FAIL to cut. Like
"hdr" it pins all three of stage 4's targets at once (is_cut False, is_hdr False, indel_length 0),
so a seed where every row complies scores consistency 1.000. What makes it hard is the probability:
`cut_p = min(0.99, 0.78 + 0.18*energy)` for Cas12a, and on K562 (accessibility 0.77) energy sits at
0.96-1.00 across the Cas12a GC band, so

    P(no cut) = 1 - cut_p  ~ 0.040 - 0.047 per (guide, seed)

against the ~0.95-0.99 that the "cut" rule asks for. A guide therefore no-cuts on only ~3 of 75
window seeds, so `main_max_fail` has to sit just under the window width for anything to bank at all
-- the opposite end of the range from every other rule.

The clean band is then the seeds where EVERY group guide no-cuts, i.e. the complement of the
min-union of failed seeds. With each guide failing ~72 of 75, the group can only hold a band where
its guides' no-cut seeds COINCIDE, and the arithmetic is a bank-size question: a bank of N holds
about N*p**k guides sharing any given k-seed set, so k=2 needs ~N*0.002 and k=3 needs ~N*9e-5 to
exceed the group size.

`nocut_band.py` already measured the ceiling a different way -- greedy seed choice over the whole
2.25M-guide space, no window and no screen -- and got a common band of **3** for every row count
from 20 to 170 (its greedy: 1 seed -> 95,816 guides, 2 -> 4,213, 3 -> 212, 4 -> 15, 5 -> 1). This
asks what the shipped bank/min-union path reaches at width 75, and at which mf.

Reports the min-union band at both group sizes in use (all-cut's 42, all-HDR's 80), since a smaller
group needs fewer coincident guides and the two can differ here.

    python nocut_mf.py
    NC_MF=70,72,74 NC_WINDOW=100-174 python nocut_mf.py
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "ncmf")

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

CELL = os.getenv("NC_CELL", "K562")
WINDOW = tuple(int(x) for x in os.getenv("NC_WINDOW", "100-174").split("-"))
MFS = [int(x) for x in os.getenv("NC_MF", "60,65,68,70,72,73,74").split(",")]
GROUPS = [int(x) for x in os.getenv("NC_GROUPS", "42,80").split(",")]
LOAD_LIMIT = int(os.getenv("NC_LOAD", "300000"))
OUT = os.getenv("NC_OUT", "nocut_mf.json")


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
    acc = cell_types.get(CELL, {}).get("accessibility")
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  window {lo}-{hi} "
          f"(span {span})  rule nocut  accessibility {acc}  "
          f"cas12a gc {base.cas12a_gc}\n")
    print(f"{'mf':>4}{'needs':>7}{'raw bank':>10}{'loaded':>8}{'nocut/guide':>13}"
          + "".join(f"{('g'+str(g)+' union'):>10}{('band'):>6}" for g in GROUPS)
          + f"{'s':>6}")
    out = []
    for mf in MFS:
        t0 = time.monotonic()
        cfg = _dc.replace(base, start_seed=lo, end_seed=hi, cas12a_max_fail=mf, rule="nocut")
        rec = {"cell": CELL, "window": [lo, hi], "span": span, "mf": mf,
               "needs_nocut": span - mf, "rule": "nocut",
               "task": (task.get("task_id") or task["id"])[:8]}
        path = os.path.join(AC.BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
            MT.free_gpu_memory()
            rec["raw_bank"] = len(bank)
            if not bank:
                rec["empty"] = True
                print(f"{mf:>4}{span-mf:>7}{0:>10}    no guide no-cuts on {span-mf} of {span} seeds"
                      f"{time.monotonic()-t0:>6.0f}")
                out.append(rec); json.dump(out, open(OUT, "w"), indent=1); continue
            save_bank(path, bank)
        records = load_bank(path, limit=LOAD_LIMIT)
        rec["loaded_bank"] = len(records)
        # how many of the window's seeds each banked guide actually skips
        nocuts = np.asarray([span - len(r["fails"]) for r in records])
        rec["nocut_per_guide"] = {"max": int(nocuts.max()), "mean": round(float(nocuts.mean()), 2)}
        line = (f"{mf:>4}{span-mf:>7}{rec.get('raw_bank', 'cached'):>10}{len(records):>8}"
                f"{rec['nocut_per_guide']['mean']:>7.1f}/{rec['nocut_per_guide']['max']:<5}")
        for g in GROUPS:
            if len(records) < g:
                rec[f"g{g}"] = {"error": f"bank {len(records)} < group {g}"}
                line += f"{'-':>10}{'-':>6}"
                continue
            sel = FG.FastGreedy(records, window_lo=lo, window_hi=hi)
            idx, _u = sel.best(g, restarts=12)
            bad = set()
            for i in idx:
                bad.update(int(x) for x in records[i]["fails"])
            clean = sorted(set(range(lo, hi + 1)) - bad)
            rec[f"g{g}"] = {"union": len(bad), "clean": len(clean), "clean_seeds": clean}
            line += f"{len(bad):>10}{len(clean):>6}"
        rec["elapsed_s"] = round(time.monotonic() - t0, 1)
        print(line + f"{rec['elapsed_s']:>6.0f}")
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
