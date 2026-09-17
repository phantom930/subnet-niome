#!/usr/bin/env python3
"""conj_nomh.py — the conjunction with its band pinned on `not_mhnhej` instead of `hdr`.

**RESULT (2026-09-16, 3 HEK293 contracts; 2026-09-17, 2 K562 contracts): dead on both cells.**
`not_mhnhej` only pins `is_cut` inside the conjunction too, so the "band" is not a third regime —
its measured value lands ON TOP of the ordinary off-band clean value, not above it (HEK293
0.23-0.29 vs clean 0.17-0.27; K562 0.19-0.28 vs clean 0.24-0.30 — on K562's k=25 arm the clean
value was even slightly ABOVE the band). It buys a much wider clean set via a far bigger Cas12a
pool (HEK293 14,000+ against `hdr`'s ~100; K562 ~26,000 against ~230-240) — HEK293 65-67 of 300
against 13-15; K562 219-220 of 300 against 95-98 — but every regime tops out under ~0.30, so the
best reachable round still scores below every field's cutoff. `E[own-field share]` measured
**exactly 0.000000 on 3/3 HEK293 contracts and 2/2 K562 contracts** at every `not_mhnhej` depth
tried, against `hdr`'s real positive value on every one of those same contracts (HEK293
~0.00009-0.00010; K562 0.000133/0.000491). The Cas12a band-formation wall differs sharply by cell
— HEK293 26 seeds, K562 ~39-40 -- and was tested AT the wall on K562 too (k=39/40, the deepest
either contract reaches): band value rises to 0.21-0.36 there but the clean set collapses in
lockstep (clean value 0.16-0.18, down from 0.24-0.30 at k=11-25), and `E[own-field share]` is
still exactly 0.000000 on both contracts. Depth is not a lever that ever rescues this rule, from
k=11 to its build limit, on either cell. Do not re-run; see the module docstring of
`niome_subnet/genomics/conjunction.py` and CLAUDE.md's falsified table for the same finding
recorded in-line where other builders would look for it.

Supports any cell via `CN_CELL` (default HEK293); the baseline/isolate arms (A/B) use that cell's
OWN shipped `band_k` from `conjunction.CELL_CONFIG`, not a literal — HEK293 and the erythroid cells
no longer share one value (8 vs 11 as of 2026-09-17).

The shipped conjunction (`niome_subnet/genomics/conjunction.py`) min-unions a Cas12a group on CUT
over the joined 300-seed space, keeping only guides that ALSO repair by HDR on every seed of a
k-seed band. `not_mhnhej` ({"any": ("HDR", "BLUNT_NHEJ")}) was already measured as a STANDALONE
band rule (CLAUDE.md, "no-MH_NHEJ (`not_mhnhej`) for a wider band") and found dead: it only pins
`is_cut` (a no-cut row still fails it, but both admitted outcomes leave `is_hdr` and `indel_length`
free), so a band seed scores ~0.123 against `hdr`'s exact 1.0 — even though the band itself reaches
~45 seeds against `hdr`'s ~12-13 at the same conditions.

That measurement was never run INSIDE the conjunction, where every seed is already cut-clean before
the band is even considered — so the question here is different: does `not_mhnhej`'s much higher
per-row compliance (fewer guides get excluded per band seed than `hdr` excludes) buy a bigger
Cas12a pool for `FastGreedy`, hence a WIDER cut-clean set (the term this file's own docstring says
is what actually pays), even though the band itself is close to worthless? Or does the extra filter
just shrink the pool for nothing, since a `not_mhnhej`-compliant seed is worth roughly what an
ordinary cut-clean seed is worth already?

Arms, at the shipped HEK293 config (group 80, width 150, light_cell_rows 12, band_cell_aware on)
otherwise unchanged — only `band_rule` and `band_k` vary:

    A  hdr         k=9   the shipped arm
    B  not_mhnhej  k=9   isolates the rule swap alone, band_k held fixed
    C  not_mhnhej  k=20
    D  not_mhnhej  k=35
    E  not_mhnhej  k=50

Each arm is priced against the ONE field that played its contract (never pooled — CLAUDE.md,
"Pricing a construction"), and the band regime is scored at its MEASURED consistency rather than an
assumed 1.0, exactly as `conj_mhany.py` does for `mh_any`.

No network: `sd_task_listing.json` / `sd_task_scores.json` are read directly (both cached earlier
this session) rather than through `sd_task.fetch()`, which would block on the unreachable backend
once its 900s freshness window has passed. Cell-type accessibility is read from
`test_hek/cell_types.json` rather than `genExp.fetch_cell_types()` for the same reason — that
function degrades to accessibility 1.0 on a failed fetch, which would silently change every energy
and cut_p computation relative to a real validator run.

    CN_N=2 CN_NS=10 python conj_nomh.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "cnomh")

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
from sd_task import score, OURS                          # noqa: E402

CELL = os.getenv("CN_CELL", "HEK293")
N_CONTRACTS = int(os.getenv("CN_N", "2"))
NS = int(os.getenv("CN_NS", "10"))
NSEED = 900
DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]

# k for arms A/B is the CELL's OWN shipped band_k (resolved in main() from CJ.config_for), not a
# literal -- HEK293 and the three erythroid cells no longer share one value (8 vs 11 as of
# 2026-09-17), and a HEK293-tuned baseline would not isolate the rule change on another cell.
ARMS_TEMPLATE = [
    ("A hdr        k{k}", "hdr", None),
    ("B not_mhnhej k{k}", "not_mhnhej", None),
    ("C not_mhnhej k25", "not_mhnhej", 25),
    ("D not_mhnhej k39", "not_mhnhej", 39),
    ("E not_mhnhej k40", "not_mhnhej", 40),
]
BAND_WIDTH = 150


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def cell_types_table():
    """Real accessibility table (HEK293 0.35), not the network-dependent live fetch."""
    return load_json("test_hek/cell_types.json")


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


def own_fields(items, scores):
    meta = {}
    for t in items:
        c = (t.get("content") or {}).get("contract", {})
        s = str(c.get("seed", "") or "")
        meta[t.get("task_id") or t.get("id")] = (
            c.get("cell_type"), len([x for x in s.split(",") if x.strip().isdigit()]))
    best = defaultdict(dict)
    for r in scores:
        c, n = meta.get(r.get("task_id"), (None, 0))
        if c != CELL or n != 3 or r.get("miner_hotkey") in OURS:
            continue
        hk, f = r["miner_hotkey"], float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    return {tid: sorted(v.values(), reverse=True)
            for tid, v in best.items() if len(v) >= 10 and max(v.values()) > 0}


def main():
    cell_types = cell_types_table()
    G.load_sequence()
    items = load_json("sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    scores = load_json("sd_task_scores.json")
    scores = scores if isinstance(scores, list) else (scores.get("data") or scores.get("items") or [])

    ftasks = own_fields(items, scores)
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
          f"shipped k={base.band_k} group={base.group_size} light={base.light_cell_rows}\n",
          flush=True)
    ARMS = [(label.format(k=base.band_k), rule, k if k is not None else base.band_k)
            for label, rule, k in ARMS_TEMPLATE]

    agg = defaultdict(list)
    for t in cands:
        tid = t.get("task_id") or t["id"]
        contract = dict(t["content"]["contract"])
        reference = t["content"]["hbb_reference"]
        fs = [ftasks[tid]]
        print(f"=== {tid[:8]}  own cut10 {ftasks[tid][9]:.1f}  top {ftasks[tid][0]:.1f} "
              f"seed {contract.get('seed')} ===", flush=True)
        for label, rule, k in ARMS:
            cfg = dataclasses.replace(base, seed_list=tuple(joined),
                                      band_rule=rule, band_width=BAND_WIDTH, band_k=k)
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=900.0)
            if not rows:
                print(f"  {label}: DECLINED — {meta.get('reason')} "
                      f"({time.monotonic()-t0:.0f}s)", flush=True)
                continue
            band = [int(x) for x in meta["band_seeds"]]
            clean_seeds = set(int(x) for x in meta.get("clean_seeds", []))
            off_band_clean = sorted(clean_seeds - set(band))
            rng = random.Random(0)
            sb = score(rows, contract, reference, cell_types, seed=band[0])
            wxf = sb["weighted"] * sb["fidelity"]
            vc = ([score(rows, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(off_band_clean, min(NS, len(off_band_clean)))]
                  if off_band_clean else [0.0])
            offs = [s for s in range(100, 1000) if s not in set(joined)]
            vr = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(offs, min(NS, len(offs)))]
            v_band, v_clean, v_rest = sb["consistency"], float(np.mean(vc)), float(np.mean(vr))
            e = price(v_band, len(band), len(off_band_clean), v_clean, v_rest, wxf, fs)
            agg[label].append(e)
            print(f"  {label}: band {len(band)} clean {len(clean_seeds)}/{len(joined)} "
                  f"pool {meta.get('pool')} cas9 {meta.get('cas9_pool')} "
                  f"wxf {wxf:.1f} | v_band {v_band:.4f} v_clean {v_clean:.4f} "
                  f"v_rest {v_rest:.4f} | E[own] {e:.6f}  ({time.monotonic()-t0:.0f}s)",
                  flush=True)
        print(flush=True)

    print("=== aggregate E[own-field share] ===")
    for label, _r, _k in ARMS:
        v = agg.get(label) or [0.0]
        print(f"  {label}: mean {float(np.mean(v)):.6f}  over {len(v)} contracts  {v}")


if __name__ == "__main__":
    main()
