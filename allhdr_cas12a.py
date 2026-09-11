#!/usr/bin/env python3
"""allhdr_cas12a.py — the all-HDR band from a Cas12a-ONLY group, width 30, K562.

The shipped all-HDR construction min-unions a Cas12a group of 80 on the `hdr` rule and then fills
the remaining ~170 rows with Cas9 guides required to reach HDR on the resulting clean band. That
second half is the binding constraint everywhere it has been measured -- it decays as P(HDR) per
extra band seed, and CLAUDE.md's "the band cannot be widened" entry is entirely about it.

Dropping Cas9 removes that constraint and replaces it with a harder one: every row now comes from
the min-union group, so filling 250 rows needs 250 guides that all reach HDR on the SAME band,
against the 80 the mixed build needs. Group size is therefore the real variable, exactly as in
`nocut_cas12a.py`.

Unlike the no-cut rule, HDR wants HIGH energy: `hdr_w = 0.24 + 0.35*energy` for Cas12a and the row
must cut first, so P(HDR) climbs with energy and the shipped `cas12a_gc` is the wide 0.40-0.95
rather than all-cut's 0.40-0.60. This sweep therefore runs the shipped config and varies only
`main_max_fail` and group size.

For scale: P(HDR) ~ 0.49 on K562, against the no-cut rule's ~0.13 at its best. A bank of N holds
about `N * p**k` guides sharing any k-seed band, so the same 300,000-guide bank that gave a no-cut
band of 5-6 should reach roughly `ln(250/300000)/ln(0.49)` ~ 10 here.

    python allhdr_cas12a.py
    AHC_MF=12,14,16 AHC_GROUPS=80,250 python allhdr_cas12a.py
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "ahc")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank   # noqa: E402
from allhdr_mf import pick                                  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402

CELL = os.getenv("AHC_CELL", "K562")
WINDOW = tuple(int(x) for x in os.getenv("AHC_WINDOW", "100-129").split("-"))
MFS = [int(x) for x in os.getenv("AHC_MF", "8,10,12,14,16,18,20").split(",")]
GROUPS = [int(x) for x in os.getenv("AHC_GROUPS", "42,80,125,170,250").split(",")]
LOAD = int(os.getenv("AHC_LOAD", "300000"))
OUT = os.getenv("AHC_OUT", "allhdr_cas12a.json")


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = WINDOW
    span = hi - lo + 1
    base = AH.config_for(CELL)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    linear = AH._scaled_max_fail(CELL, base.main_max_fail, span)
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  window {lo}-{hi} (span {span})"
          f"  rule hdr  gc {base.cas12a_gc}  max_distance {base.max_distance}"
          f"  variants {base.variants}")
    print(f"shipped mf at this span (linear from {base.main_max_fail}/100) = {linear}\n")
    print(f"{'mf':>4}{'needs':>7}{'raw bank':>10}{'loaded':>8}{'hdr/guide':>12}"
          + "".join(f"{('g'+str(g)):>7}" for g in GROUPS) + f"{'s':>7}")
    out = []
    for mf in MFS:
        t0 = time.monotonic()
        cfg = _dc.replace(base, hdr_range=(lo, hi), main_max_fail=mf)
        rec = {"cell": CELL, "window": [lo, hi], "span": span, "mf": mf, "needs": span - mf,
               "shipped_mf": linear, "rule": "hdr", "cas12a_gc": list(base.cas12a_gc),
               "max_distance": base.max_distance,
               "task": (task.get("task_id") or task["id"])[:8]}
        path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg)
            MT.free_gpu_memory()
            rec["raw_bank"] = len(bank)
            if not bank:
                rec["empty"] = True
                print(f"{mf:>4}{span-mf:>7}{0:>10}    nothing qualifies{time.monotonic()-t0:>7.0f}")
                out.append(rec); json.dump(out, open(OUT, "w"), indent=1); continue
            save_bank(path, bank)
        records = load_bank(path, limit=LOAD)
        rec["loaded_bank"] = len(records)
        hdrs = np.asarray([span - len(r["fails"]) for r in records])
        rec["hdr_per_guide"] = {"max": int(hdrs.max()), "mean": round(float(hdrs.mean()), 2)}
        line = (f"{mf:>4}{span-mf:>7}{rec.get('raw_bank', 'cached'):>10}{len(records):>8}"
                f"{rec['hdr_per_guide']['mean']:>7.1f}/{rec['hdr_per_guide']['max']:<4}")
        rec["bands"] = {}
        for g in GROUPS:
            if len(records) < g:
                rec["bands"][g] = None
                line += f"{'-':>7}"
                continue
            sel = FG.FastGreedy(records, window_lo=lo, window_hi=hi,
                                per_cell_min=cfg.per_cell_min)
            idx, _u = sel.best(g, restarts=cfg.restarts)
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
