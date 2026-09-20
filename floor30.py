#!/usr/bin/env python3
"""floor30.py — does the width-30 tiling's frequency gain survive the FLOOR?

`width30_test.py` settled band formation (100/100) and coverage (union 104 vs 70 of 900, +43% on
band-hit frequency) but stopped before the clean set. That gap matters: the off-band floor comes
from `FastGreedy`'s min-union over the group, and at width 30 the group is drawn from a surviving
pool of 163-206 against 199-232 at width 150. Less freedom there can mean a larger failed-seed union
and a thinner floor, which is the term the +43% would have to be paid out of.

**Per-hotkey E[own-field share] isolates the floor exactly, and that is why it is the right unit
here.** The two arms build the SAME band size (k is unchanged), so per-hotkey `P(band hit)` is
identical between them and cancels. Whatever difference the pricing shows is therefore the clean set
and term 1, nothing else. The coverage half is already measured and is a property of the tiling, so
the two halves multiply rather than needing to be measured together:

    fleet value  ~  (band union, measured: 104 vs 70)  x  (per-hotkey value given the band, here)

Two hotkeys per contract rather than one, because at width 30 the windows are disjoint and a single
hotkey's slice is only 30 of the 300 — h0 (100-129) and h5 (250-279) sit at opposite ends, so a
per-slice artefact would show up as a disagreement between them.

Priced with `widecut_price.evaluate`, which enumerates band and cut-clean from the SHIPPED ROWS over
all 900 and scores each regime at its measured consistency — the only method in this repo that is
unbiased across arms whose seed spaces differ.

    F30_N=3 python floor30.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "floor30")

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
from widecut_price import evaluate, own_fields           # noqa: E402

CELL = os.getenv("F30_CELL", "CD34+_HSPC")
N = int(os.getenv("F30_N", "3"))
HKS = [h for h in os.getenv("F30_HK", "niome_hotkey,niome_hotkey5").split(",") if h]
ARMS = [("w30", 30), ("w150", 150)]
OUT = os.getenv("F30_JSON", "floor30.json")


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    plan = json.load(open("data/window_plan.json"))
    base = CJ.config_for(CELL)
    ftasks = own_fields(CELL)
    raw = plan["assignments"][CELL][HKS[0]]
    pairs = raw if isinstance(raw[0], (list, tuple)) else [raw]
    predicted = sorted({s for a, b in pairs for s in range(a, b + 1)})

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
        if len(tasks) >= N:
            break

    print(f"=== {CELL}  k={base.band_k} group={base.group_size} light={base.light_cell_rows} "
          f"| window {pairs} | {len(tasks)} contracts x {len(HKS)} hotkeys x 2 arms ===",
          flush=True)

    recs, agg = [], defaultdict(list)
    for t in tasks:
        tid = (t.get("task_id") or t["id"])
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        field = ftasks[tid]
        print(f"\n  {tid[:8]}  own cut10 {field[9]:.1f}  top {field[0]:.1f}", flush=True)
        ctx = G.build_context(contract, reference, cell_types)
        for hk in HKS:
            for arm, width in ARMS:
                cands = CJ.sub_window(predicted, width, JW.band_offset_frac(hk) or 0.0)
                cut = JW.conjunction_cut_seeds(hk, CELL, predicted=predicted,
                                               band_candidates=cands)
                cfg = dataclasses.replace(base, seed_list=tuple(sorted(cut or predicted)),
                                          band_candidates=tuple(sorted(cands)),
                                          band_width=width)
                t0 = time.monotonic()
                # Build on seed 0: the miner prefetches before the seeds are stamped.
                rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                                 cfg=cfg, budget_s=1800)
                key = f"{hk[-2:] or 'h0'}/{arm}"
                if not rows:
                    print(f"    {key:9} DECLINED — {meta.get('reason')}", flush=True)
                    agg[arm].append(0.0)
                    recs.append({"task": tid[:8], "hk": hk, "arm": arm, "built": False,
                                 "reason": meta.get("reason"), "share": 0.0})
                    continue
                r = evaluate(rows, contract, reference, cell_types, ctx, field, key, t0)
                agg[arm].append(r["share"] if r else 0.0)
                pool = meta.get("pool")
                recs.append(dict(r or {}, task=tid[:8], hk=hk, arm=arm, built=True,
                                 pool=pool, margin=pool / meta["group_size"],
                                 cas9=meta.get("cas9_pool"),
                                 elapsed=round(time.monotonic() - t0, 1)))
                print(f"              pool {pool}/{meta['group_size']} = "
                      f"{pool/meta['group_size']:.2f}x  cas9 {meta.get('cas9_pool')}", flush=True)

    with open(OUT, "w") as fh:
        json.dump(recs, fh, indent=1)

    print(f"\n=== per-hotkey means, {len(tasks)} contracts x {len(HKS)} hotkeys ===", flush=True)
    for arm, _w in ARMS:
        v = agg[arm]
        b = [r for r in recs if r["arm"] == arm and r.get("built")]
        print(f"  {arm:5} E[own] {np.mean(v):.6f}  built {len(b)}/{len(v)}  "
              f"clean {np.mean([r['cutonly'] + r['band'] for r in b]):.1f}/900  "
              f"w*fid {np.mean([r['wxf'] for r in b]):.1f}  "
              f"margin {np.mean([r['margin'] for r in b]):.2f}x", flush=True)

    a = [r for r in recs if r["arm"] == "w30"]
    b = {(r["task"], r["hk"]): r for r in recs if r["arm"] == "w150"}
    d = [(x["share"], b[(x["task"], x["hk"])]["share"]) for x in a
         if (x["task"], x["hk"]) in b]
    w = sum(1 for p, q in d if p > q)
    l = sum(1 for p, q in d if p < q)
    ma, mb = np.mean([p for p, _ in d]), np.mean([q for _, q in d])
    print(f"\n=== paired within (contract, hotkey) ===", flush=True)
    print(f"  w30 {ma:.6f} vs w150 {mb:.6f} = {(ma/mb if mb else float('nan')):.3f}x  "
          f"{w}W/{l}L/{len(d)-w-l}T", flush=True)
    print(f"\n  FLEET, folding in the measured union (104 vs 70 of 900):", flush=True)
    print(f"    w30  {ma/mb if mb else float('nan'):.3f}x per hotkey  x  "
          f"{104/70:.3f}x coverage  =  {(ma/mb if mb else 0)*104/70:.3f}x", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
