#!/usr/bin/env python3
"""conj_mhany.py — the conjunction with its band pinned on `mh_any` instead of `hdr`.

The shipped conjunction min-unions a Cas12a group on CUT over the joined 300-seed space, keeping
only guides that also repair by HDR on every seed of a k-seed band. This swaps that band rule for
`mh_any` ("HDR if the microhomology coin came up, else BLUNT_NHEJ") and narrows the band sub-window
to 100.

What the swap trades, both halves measured on HEK293 at width 100 standalone (`rule_window.py`):

    rule      per-row comply   band (standalone)   band-seed cons
    hdr       0.4471           7                   1.0000
    mh_any    0.3807           8                   0.7778

so `mh_any` is band-ADVANTAGED on HEK293 (unlike K562, where it is 9 against 12) and pays with a
~0.22 lower band value. Inside the conjunction the trade differs again: k is chosen explicitly and
the Cas9 conditional fill decays by P(rule) per band seed -- 0.381 against 0.447 -- so a LOWER k may
be feasible from the same pool, and the arm DECLINES rather than shipping a shallower band.

Three arms per contract, so rule and width are separable:

    A  hdr     width 300   the shipped HEK293 arm
    B  mh_any  width 100   the requested construction
    C  hdr     width 100   isolates the width change from the rule change

Each arm is priced against the ONE field that played its contract, never a pooled field, and the
band regime is scored at its MEASURED value rather than an assumed 1.0 -- `conj_stageb.price`
hardcodes 1.0, which is right for `hdr` and wrong for any rule that leaves a target merely
predictable instead of constant.

    CM_N=3 CM_NS=10 python conj_mhany.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "cmhany")

import dataclasses
import json
import logging
import random
import time
from collections import defaultdict
from math import factorial

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from conj_stageb import DIST, NSEED, OURS, API           # noqa: E402
from sd_task import score, fetch                         # noqa: E402

CELL = os.getenv("CM_CELL", "HEK293")
N_CONTRACTS = int(os.getenv("CM_N", "3"))
NS = int(os.getenv("CM_NS", "10"))

ARMS = [
    ("A hdr    w300", "hdr", 300),
    ("B mh_any w100", "mh_any", 100),
    ("C hdr    w100", "hdr", 100),
]


def price(v_band, band, clean_excl_band, v_clean, v_rest, wxf, fs):
    """E[curve share] with the band regime worth `v_band` rather than an assumed 1.0."""
    pb, pc = band / NSEED, clean_excl_band / NSEED
    pr = max(0.0, 1.0 - pb - pc)
    tot = 0.0
    for nb in range(4):
        for nc in range(4 - nb):
            nr = 3 - nb - nc
            w = (factorial(3) / (factorial(nb) * factorial(nc) * factorial(nr))
                 * pb ** nb * pc ** nc * pr ** nr)
            if w <= 0:
                continue
            final = wxf * (nb * v_band + nc * v_clean + nr * v_rest) / 3.0
            share = 0.0
            for f in fs:
                rank = 1 + sum(1 for x in f if x > final)
                if rank <= 10:
                    share += DIST[rank - 1]
            tot += w * share / max(1, len(fs))
    return tot


def own_fields():
    sc = fetch(f"{API}/miners/scores?limit=40000", cache="sd_task_scores.json")
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    tk = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    tk = tk if isinstance(tk, list) else (tk.get("items") or tk.get("data") or [])
    meta = {}
    for t in tk:
        c = (t.get("content") or {}).get("contract", {})
        s = str(c.get("seed", "") or "")
        meta[t.get("task_id") or t.get("id")] = (
            c.get("cell_type"), len([x for x in s.split(",") if x.strip().isdigit()]))
    best = defaultdict(dict)
    for r in sc:
        c, n = meta.get(r.get("task_id"), (None, 0))
        if c != CELL or n != 3 or r.get("miner_hotkey") in OURS:
            continue
        hk, f = r["miner_hotkey"], float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    return {tid: sorted(v.values(), reverse=True)
            for tid, v in best.items() if len(v) >= 10 and max(v.values()) > 0}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ftasks = own_fields()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    cands = []
    for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        s = str(c.get("seed", "") or "")
        if c.get("cell_type") != CELL or len([x for x in s.split(",") if x.strip().isdigit()]) != 3:
            continue
        tid = t.get("task_id") or t["id"]
        if tid in ftasks:
            cands.append(t)
        if len(cands) >= N_CONTRACTS:
            break

    spans = JW.fallback_for(CELL, "niome_hotkey")
    joined = sorted({s for a, b in spans for s in range(a, b + 1)})
    base = CJ.config_for(CELL)
    print(f"{CELL}: {len(cands)} contracts, joined space {spans} ({len(joined)} seeds), "
          f"k={base.band_k} group={base.group_size} light={base.light_cell_rows}\n")

    agg = defaultdict(list)
    for t in cands:
        tid = t.get("task_id") or t["id"]
        contract = dict(t["content"]["contract"])
        reference = t["content"]["hbb_reference"]
        fs = [ftasks[tid]]
        print(f"=== {tid[:8]}  own cut10 {ftasks[tid][9]:.1f}  top {ftasks[tid][0]:.1f} ===",
              flush=True)
        for label, rule, width in ARMS:
            cfg = dataclasses.replace(base, seed_list=tuple(joined),
                                      band_rule=rule, band_width=width)
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=900.0)
            if not rows:
                print(f"  {label}: DECLINED — {meta.get('reason')}", flush=True)
                continue
            band = [int(x) for x in meta["band_seeds"]]
            clean = set(int(x) for x in meta.get("clean_seeds", [])) if meta.get("clean_seeds") \
                else None
            rng = random.Random(0)
            sb = score(rows, contract, reference, cell_types, seed=band[0])
            wxf = sb["weighted"] * sb["fidelity"]
            # clean-off-band and dirty, sampled from the joined space
            from conj_stageb import hdr_clean_seeds  # noqa
            cutclean = [s for s in joined if s not in set(band)]
            # the clean set the meta reports is over the joined space
            nclean = int(meta.get("clean", 0))
            pool_clean = cutclean[:]  # sampled below; meta['clean'] gives the count
            vc = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(pool_clean, min(NS, len(pool_clean)))]
            offs = [s for s in range(100, 1000) if s not in set(joined)]
            vr = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(offs, min(NS, len(offs)))]
            v_band, v_clean, v_rest = sb["consistency"], float(np.mean(vc)), float(np.mean(vr))
            e = price(v_band, len(band), max(0, nclean - len(band)), v_clean, v_rest, wxf, fs)
            agg[label].append(e)
            print(f"  {label}: band {len(band)} clean {nclean}/{len(joined)} "
                  f"pool {meta.get('pool')} cas9 {meta.get('cas9_pool')} "
                  f"wxf {wxf:.1f} | v_band {v_band:.4f} v_clean {v_clean:.4f} "
                  f"v_rest {v_rest:.4f} | E[own] {e:.6f}  ({time.monotonic()-t0:.0f}s)",
                  flush=True)
        print(flush=True)

    print("=== aggregate E[own-field share] ===")
    for label, _r, _w in ARMS:
        v = agg.get(label) or [0.0]
        print(f"  {label}: mean {float(np.mean(v)):.6f}  over {len(v)} contracts  {v}")


if __name__ == "__main__":
    main()
