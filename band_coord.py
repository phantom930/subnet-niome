#!/usr/bin/env python3
"""band_coord.py — what per-class band coordination costs in depth.

Under ORACLE class knowledge the three round seeds sit one per 100-seed class, so a hotkey holding
`(k1,k2,k3)` band seeds across those classes doubles with probability `(k1k2+k1k3+k2k3)/100**2` --
zero if the band is confined to one class. The shipped `sub_window` takes a CONTIGUOUS 150-seed
slice of a 300-seed joined space, so every band spans at most 1.5 classes: measured over 10 real
fleet bands x 10 tasks, **0 doubles in 100 hotkey-rounds** at a mean pair-sum of 14.1 of 21.

The fix is to draw the band from a block inside EVERY class (`class_blocks`) under a per-class quota,
which also makes the fleet's rectangles disjoint so the union equals the sum. **What that costs is
depth**, and this measures it: `choose_band` earns its k by picking whichever seeds the most guides
comply on, and `hdr_pool.py` prices that freedom at an 11x pool edge over an arbitrary band. A small
block spends it.

Per (cell, contract, block width) it descends k until a build succeeds, so the output is the
reachable depth at each level of coordination, plus the band's actual class split and pair-sum.

    BC_N=2 python band_coord.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "bandcoord")

import dataclasses
import json
import logging
import math
import time
from collections import Counter, defaultdict

logging.basicConfig(level=logging.ERROR)

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402

N = int(os.getenv("BC_N", "2"))
CELLS = os.getenv("BC_CELLS", "K562,HEK293").split(",")
BLOCKS = [int(x) for x in os.getenv("BC_BLOCKS", "10,20,30,100").split(",")]
KMAX = {"HEK293": 9}
HK = int(os.getenv("BC_HK", "0"))


def pairsum(split):
    a, b, c = (split + [0, 0, 0])[:3]
    return a * b + a * c + b * c


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    out = []
    for cell in CELLS:
        base = CJ.config_for(cell)
        spans = JW.fallback_for(cell, "niome_hotkey")
        joined = sorted({s for a, b in spans for s in range(a, b + 1)})
        classes = sorted({s // 100 for s in joined})
        kmax = KMAX.get(cell, 12)
        cands = []
        for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
            c = (t.get("content") or {}).get("contract") or {}
            s = str(c.get("seed", "") or "")
            if c.get("cell_type") != cell or len(
                    [x for x in s.split(",") if x.strip().isdigit()]) != 3:
                continue
            cands.append(t)
            if len(cands) >= N:
                break
        print(f"########## {cell}  joined {spans}  group={base.group_size} "
              f"k_start={kmax} ##########", flush=True)
        for t in cands:
            tid = (t.get("task_id") or t["id"])[:8]
            contract = dict(t["content"]["contract"])
            reference = t["content"]["hbb_reference"]
            print(f"=== {tid} ===", flush=True)
            for block in BLOCKS + [0]:          # 0 = shipped sub_window, no quota
                cset = tuple(CJ.class_blocks(joined, block, HK)) if block else ()
                got = None
                for k in range(kmax, 3, -1):
                    cap = math.ceil(k / max(1, len(classes)))
                    quota = tuple((c * 100, c * 100 + 99, cap) for c in classes) if block else ()
                    cfg = dataclasses.replace(base, seed_list=tuple(joined), band_k=k,
                                              band_candidates=cset, band_quota=quota)
                    t0 = time.monotonic()
                    rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                                     budget_s=1200.0)
                    if rows:
                        band = [int(x) for x in meta["band_seeds"]]
                        cnt = Counter(s // 100 for s in band)
                        split = [cnt.get(c, 0) for c in classes]
                        got = {"cell": cell, "task": tid, "block": block, "k": k,
                               "split": split, "pairsum": pairsum(split),
                               "cas9": meta.get("cas9_pool"), "pool": meta.get("pool"),
                               "secs": round(time.monotonic() - t0)}
                        print(f"  block {block if block else 'shipped':>7}: k={k:2d} "
                              f"split {'/'.join(map(str, split))} pairsum {got['pairsum']:2d} "
                              f"cas9 {meta.get('cas9_pool')} ({got['secs']}s)", flush=True)
                        break
                if got is None:
                    print(f"  block {block if block else 'shipped':>7}: no build at k>=4",
                          flush=True)
                    got = {"cell": cell, "task": tid, "block": block, "k": 0,
                           "split": [], "pairsum": 0}
                out.append(got)
        print(flush=True)
    with open("band_coord.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print("=== summary: reachable depth and double-hit pair-sum vs coordination ===")
    print(f"{'cell':9s} {'block':>7s} {'mean k':>7s} {'mean pairsum':>13s} {'P(double)/hk':>13s} "
          f"{'fleet x10':>10s}")
    agg = defaultdict(list)
    for r in out:
        agg[(r["cell"], r["block"])].append(r)
    for (cell, block), rs in sorted(agg.items(), key=lambda x: (x[0][0], x[0][1])):
        mk = sum(r["k"] for r in rs) / len(rs)
        mp = sum(r["pairsum"] for r in rs) / len(rs)
        print(f"{cell:9s} {block if block else 'shipped':>7} {mk:7.1f} {mp:13.1f} "
              f"{mp/1e4:13.5f} {min(1.0, 10*mp/1e4):10.4f}")


if __name__ == "__main__":
    main()
