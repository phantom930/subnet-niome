#!/usr/bin/env python3
"""hek_cutspan.py — how wide should HEK293's CUT window be, now that the band space is the full 900?

Three options are live, and they differ only for `REST_HK` (h6-h9), whose band is drawn from the 600
seeds the plan did NOT predict:

    narrow   cut = the predicted 300              band ⊄ cut for group B  (pre-2026-09-18 shape)
    union    cut = predicted 300 ∪ own band 150   band ⊆ cut, 450 seeds   (the middle option)
    wide     cut = the full 900                   band ⊆ cut              (SHIPS today)

For `BAND_HK` (h0-h5) the band already sits inside the predicted 300, so `union` IS `narrow` and only
two arms exist — group A is carried anyway as the reference point, since a per-cell `WIDE_CUT_CELLS`
decision applies to both groups at once and `union` would have to be expressed per HOTKEY.

**Priced by `widecut_price.evaluate`, i.e. `hud_resolve.py`'s method, and that is the whole point.**
An earlier read of this question compared `meta["clean"]` across arms -- 22-25 of 300 against 17-21
of 900 -- and that is biased by construction, exactly as `widecut_price.py`'s docstring warns: the
count is taken WITHIN each arm's own cut space, so the two numbers do not denote the same thing. The
imported `evaluate` enumerates band and cut-clean from the SHIPPED ROWS over all 900 for every arm,
scores each regime at its measured consistency, and prices the result against the one field that
played that contract. Identical work for all three arms whatever space built them.

`band_k` is NOT swept: `conj_wall.py` and the 2026-09-18 margin run both put the wide wall at 8 and
the narrow wall at 9, so the live 8 is feasible on every arm here and holding it fixed keeps this one
question. The surviving-pool margin is reported per arm, since that is where the wall difference
shows up as a number (narrow 2.6x = one full step of 1/P(HDR), wide 1.4-1.8x).

    HC_N=4 python hek_cutspan.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hekcut")

import dataclasses                                       # noqa: E402
import json                                              # noqa: E402
import logging                                           # noqa: E402
import time                                              # noqa: E402
from collections import defaultdict                      # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402
from widecut_price import SPACE, evaluate, own_fields    # noqa: E402

CELL = os.getenv("HC_CELL", "HEK293")
N_CONTRACTS = int(os.getenv("HC_N", "4"))
BUDGET = float(os.getenv("HC_BUDGET", "1800"))
OUT = os.getenv("HC_JSON", "hek_cutspan.json")
GROUPS = [("niome_hotkey", "A/predicted"), ("niome_hotkey7", "B/complement")]


def arm_cfg(base, cut_seeds, cands, k=None):
    return dataclasses.replace(base, seed_list=tuple(sorted(cut_seeds)),
                               band_candidates=tuple(sorted(cands)),
                               **({"band_k": k} if k else {}))


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    base = CJ.config_for(CELL)
    ftasks = own_fields(CELL)

    # The fixed fallback triple, not the live plan: deterministic and reproducible, and band
    # position is free under a uniform generator so the choice costs nothing.
    predicted = sorted({s for a, b in JW.fallback_for(CELL, "niome_hotkey")
                        for s in range(a, b + 1)})

    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        s = str(c.get("seed", "") or "")
        if c.get("cell_type") != CELL or len(
                [x for x in s.split(",") if x.strip().isdigit()]) != 3:
            continue
        if (t.get("task_id") or t["id"]) not in ftasks:
            continue
        tasks.append(t)
        if len(tasks) >= N_CONTRACTS:
            break

    print(f"=== {CELL}  live k={base.band_k} group={base.group_size} width={base.band_width} "
          f"light={base.light_cell_rows} ca={base.band_cell_aware} | predicted {len(predicted)} "
          f"| {len(tasks)} contracts with own fields ===", flush=True)

    recs, agg = [], defaultdict(list)
    for t in tasks:
        tid = (t.get("task_id") or t["id"])
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        field = ftasks[tid]
        print(f"\n  {tid[:8]}  own cut10 {field[9]:.1f}  top {field[0]:.1f}", flush=True)
        ctx = G.build_context(contract, reference, cell_types)
        for hk, gname in GROUPS:
            bspace = JW.band_space(hk, predicted)
            cands = CJ.sub_window(bspace, base.band_width, JW.band_offset_frac(hk) or 0.0)
            arms = [("narrow", predicted), ("wide", SPACE)]
            if not set(cands) <= set(predicted):
                arms.insert(1, ("union", sorted(set(predicted) | set(cands))))
            print(f"   {gname:13} band cands {min(cands)}-{max(cands)} ({len(cands)})", flush=True)
            for arm, cut in arms:
                cfg = arm_cfg(base, cut, cands)
                t0 = time.monotonic()
                rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                                 budget_s=BUDGET)
                pool, grp = meta.get("pool"), meta.get("group_size")
                key = f"{gname[0]}/{arm}"
                if not rows:
                    print(f"    {arm:6s} cut {len(cut):3d} DECLINED — {meta.get('reason')} "
                          f"({time.monotonic()-t0:.0f}s)", flush=True)
                    agg[key].append(0.0)
                    recs.append({"task": tid[:8], "group": gname, "arm": arm, "cut": len(cut),
                                 "built": False, "reason": meta.get("reason"), "share": 0.0})
                    continue
                r = evaluate(rows, contract, reference, cell_types, ctx, field,
                             f"{arm} cut {len(cut)}", t0)
                agg[key].append(r["share"] if r else 0.0)
                recs.append(dict(r or {}, task=tid[:8], group=gname, arm=arm, cut=len(cut),
                                 built=True, pool=pool, margin=(pool / grp) if pool else None,
                                 cas9=meta.get("cas9_pool"), bank=meta.get("bank"),
                                 clean_in_cut=meta.get("clean"),
                                 band_seeds=sorted(int(x) for x in meta["band_seeds"]),
                                 elapsed=round(time.monotonic() - t0, 1), cut10=field[9]))
                print(f"           pool {pool}/{grp} = {pool/grp:.2f}x  cas9 {meta.get('cas9_pool')}"
                      f"  bank {meta.get('bank')}", flush=True)

    with open(OUT, "w") as fh:
        json.dump(recs, fh, indent=1)
    print(f"\n=== mean E[own-field share], {len(tasks)} contracts ===", flush=True)
    for key in sorted(agg):
        v = agg[key]
        print(f"  {key:12} {np.mean(v):.6f}   per contract "
              f"{' '.join(f'{x:.6f}' for x in v)}", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
