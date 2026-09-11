#!/usr/bin/env python3
"""joined300_mf.py — is the band collapse at span 900 the WINDOW or the screen?

joined300.py measures K562's band at 12 over a 300-seed space and **4** over the full 900. Before
that is read as "all-HDR does not work at 900", the screen has to be cleared: `main_max_fail` is
calibrated against a 100-seed band and `_scaled_max_fail` falls back to LINEAR scaling for every
cell but HEK293. Linear holds the fail *rate* and not its z-score, and the binomial narrows as the
span grows, so the same rate is a far deeper tail cut at 900 than at 100:

    span   mean fails   sd     mf (linear)   z
     100         51.0   5.0            45   -1.20
     300        153.0   8.7           135   -2.08
     900        459.0  15.0           405   -3.60

That is the same mechanism CLAUDE.md records for HEK293 at span 300 (mf 48 at z -3.32 scaling to
z -5.66, every band hotkey declining) — only milder, so it thins the bank instead of emptying it.

This builds the 900-seed arm at several `main_max_fail` values, holding everything else at the
shipped config. If the band recovers at a z-matched screen, the collapse is the screen; if it stays
at 4, it is the window and the 900-wide arm is a real measurement.

    J3M_CELL=K562 J3M_MF=405,441,480 python joined300_mf.py
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "j300mf")

import dataclasses as _dc           # noqa: E402
import json                         # noqa: E402
import logging                      # noqa: E402
import time                         # noqa: E402
import urllib.request               # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                     # noqa: E402
from niome_subnet.genomics import all_hdr as AH        # noqa: E402
from niome_subnet.genomics import mt19937 as MT        # noqa: E402
from joined300 import fetch_tasks, true_band           # noqa: E402
from sd_task import score                              # noqa: E402

CELL = os.getenv("J3M_CELL", "K562")
MFS = [int(x) for x in os.getenv("J3M_MF", "405,441,480").split(",")]
SPAN = int(os.getenv("J3M_SPAN", "900"))
OUT = os.getenv("J3M_OUT", "joined300_mf.json")


def main():
    items = fetch_tasks()
    task = next(t for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True)
                if ((t.get("content") or {}).get("contract") or {}).get("cell_type") == CELL
                and str(((t.get("content") or {}).get("contract") or {}).get("seed", "0")).strip()
                not in ("0", ""))
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = 100, 100 + SPAN - 1
    base = AH.config_for(CELL)
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  span {SPAN} ({lo}-{hi})  "
          f"group {base.group_size}  linear mf "
          f"{AH._scaled_max_fail(CELL, base.main_max_fail, SPAN)}\n")
    out = []
    for mf in MFS:
        cfg = _dc.replace(base, hdr_range=(lo, hi), main_max_fail=mf,
                          variants=min(base.variants, AH.WIDE_WINDOW_VARIANTS))
        t0 = time.monotonic()
        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=2400)
        el = round(time.monotonic() - t0, 1)
        rec = {"cell": CELL, "span": SPAN, "mf": mf, "build_s": el,
               "meta": {k: v for k, v in meta.items() if k != "bank_path"}}
        if not rows:
            rec["declined"] = meta.get("reason")
            print(f"  mf {mf:<4} DECLINED after {el:.0f}s — {rec['declined']}")
        else:
            band, _n = true_band(rows, contract, reference, cell_types)
            hit = score(rows, contract, reference, cell_types,
                        seed=band[len(band) // 2] if band else lo)
            rec.update(bank=meta.get("bank"), band=len(band), cells=meta.get("cells"),
                       weighted=hit["weighted"], fidelity=hit["fidelity"],
                       wxfid=hit["weighted"] * hit["fidelity"], spike_cons=hit["consistency"])
            print(f"  mf {mf:<4} bank {rec['bank']:<6} band {rec['band']:<3} "
                  f"w {rec['weighted']:7.1f}  fid {rec['fidelity']:.4f}  "
                  f"wxfid {rec['wxfid']:7.1f}  spike {rec['spike_cons']:.4f}  "
                  f"cells {rec['cells']}/8  {el:.0f}s")
        MT.free_gpu_memory()
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
