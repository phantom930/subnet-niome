#!/usr/bin/env python3
"""width30_test.py — can k=11/k=8 form a band from a 30-seed candidate window?

The 2026-09-19 layout puts all ten hotkeys on the predicted 300 at `BAND_SUB_WIDTH` 30 and
`BAND_STRIDE` 30, so the ten candidate windows are DISJOINT and the fleet's bands are provably
disjoint too -- a first for this layout, which until now only reached ~86% distinct. What it costs
is greedy freedom: `choose_band` must pick k seeds out of 30 instead of out of 150, and CLAUDE.md's
width sweep records the band GROWING toward the window while the surviving pool collapses. So the
question this settles is availability, not score: does the band still reach k, and does the pool
still clear `group_size`?

`w150` is carried as the control at the same k, same stride, same hotkeys and the same prepare --
only the width differs -- because a decline at width 30 means nothing unless width 150 builds on the
same contract.

**Cheap by construction.** `choose_band` and `cas9_cell_probe` restrict to their `candidates`
argument by column lookup, so compliance computed ONCE over the whole 300-seed window is valid for
every sub-window of it. One bank scan and one `hdr_compliance` pass per task, then 20 greedy runs.

    W30_N=10 python width30_test.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "w30test")

import dataclasses                                       # noqa: E402
import json                                              # noqa: E402
import logging                                           # noqa: E402
import time                                              # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from niome_subnet.genomics import mt19937 as MT          # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds  # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402

N_TASKS = int(os.getenv("W30_N", "10"))
HKS = JW.BAND_HK + JW.REST_HK
ARMS = [("w30", 30), ("w150", 150)]
OUT = os.getenv("W30_JSON", "width30_test.json")


def prepare(contract, reference, cell_types, cell, predicted):
    base = CJ.config_for(cell)
    cut = JW.conjunction_cut_seeds(HKS[0], cell, predicted=predicted,
                                   band_candidates=predicted) or predicted
    cfg = dataclasses.replace(base, seed_list=tuple(sorted(cut)),
                              band_candidates=tuple(predicted))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    path = os.path.join(CJ.BANK_DIR, f"cas12a-{CJ.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        os.makedirs(CJ.BANK_DIR, exist_ok=True)
        bank = CJ.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return None
        CJ.save_bank(path, bank)
    records = CJ.load_bank(path, limit=cfg.bank_keep)
    ok = CJ.hdr_compliance(records, contract, cell_types, ctx, predicted, cfg.band_rule)
    MT.free_gpu_memory()
    cell_ok = None
    if cfg.band_cell_aware:
        cell_ok = CJ.cas9_cell_probe(contract, cell_types, ctx, sites, cfg, predicted,
                                     cfg.band_rule)
        MT.free_gpu_memory()
    return {"cfg": cfg, "ok": ok, "cell_ok": cell_ok, "bank": len(records), "cut": len(cut)}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    plan = json.load(open("data/window_plan.json"))

    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if len(_parse_seeds(c.get("seed"))) != 3:
            continue
        tasks.append(t)
        if len(tasks) >= N_TASKS:
            break

    recs = []
    for t in tasks:
        tid = (t.get("task_id") or t["id"])
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        cell = contract["cell_type"]
        seeds = _parse_seeds(contract["seed"])
        raw = plan["assignments"][cell][HKS[0]]
        pairs = raw if isinstance(raw[0], (list, tuple)) else [raw]
        predicted = sorted({s for a, b in pairs for s in range(a, b + 1)})
        t0 = time.monotonic()
        prep = prepare(dict(contract, seed=0), reference, cell_types, cell, predicted)
        if prep is None:
            print(f"{tid[:8]} {cell}: bank empty, skipped", flush=True)
            continue
        cfg, ok, cell_ok = prep["cfg"], prep["ok"], prep["cell_ok"]
        in_win = [s for s in seeds if s in set(predicted)]
        print(f"\n=== {tid[:8]} {cell:12} seeds {seeds} | k={cfg.band_k} group={cfg.group_size} "
              f"| bank {prep['bank']} cut {prep['cut']} | prep {time.monotonic()-t0:.0f}s ===",
              flush=True)
        for arm, width in ARMS:
            built, bands, pools, union = 0, [], [], set()
            for hk in HKS:
                cands = CJ.sub_window(predicted, width, JW.band_offset_frac(hk) or 0.0)
                band, alive = CJ.choose_band(ok, cfg.band_k, cands, cfg.group_size, cell_ok,
                                             {}, cfg.band_quota)
                okbuild = len(band) >= cfg.band_k and len(alive) >= cfg.group_size
                built += okbuild
                bands.append(len(band))
                pools.append(len(alive))
                if okbuild:
                    union |= set(int(x) for x in band)
            hits = [s for s in seeds if s in union]
            print(f"    {arm:5} built {built:>2}/10  band {min(bands)}-{max(bands)} "
                  f"(need {cfg.band_k})  pool {min(pools)}-{max(pools)} "
                  f"(need {cfg.group_size}, margin {min(pools)/cfg.group_size:.2f}x)  "
                  f"union {len(union):>3}/900" + (f"   SEED HIT {hits}" if hits else ""),
                  flush=True)
            recs.append({"task": tid[:8], "cell": cell, "arm": arm, "width": width,
                         "k": cfg.band_k, "group": cfg.group_size, "built": built,
                         "band_min": min(bands), "band_max": max(bands),
                         "pool_min": min(pools), "pool_max": max(pools),
                         "union": len(union), "hits": hits, "in_window_seeds": in_win})

    with open(OUT, "w") as fh:
        json.dump(recs, fh, indent=1)
    print("\n=== summary over %d tasks ===" % len(tasks), flush=True)
    for arm, _w in ARMS:
        a = [r for r in recs if r["arm"] == arm]
        if not a:
            continue
        print(f"  {arm:5} built {sum(r['built'] for r in a)}/{10*len(a)} hotkey-builds  "
              f"({sum(1 for r in a if r['built']==10)}/{len(a)} tasks all-10)  "
              f"band {min(r['band_min'] for r in a)}-{max(r['band_max'] for r in a)}  "
              f"pool {min(r['pool_min'] for r in a)}-{max(r['pool_max'] for r in a)}  "
              f"mean union {np.mean([r['union'] for r in a]):.0f}/900  "
              f"seed hits {sum(1 for r in a if r['hits'])}/{len(a)} tasks", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
