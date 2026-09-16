#!/usr/bin/env python3
"""hud_resolve.py — does the conjunction actually lose to all-HDR on HUDEP-2?

CLAUDE.md excludes HUDEP-2 from the conjunction on a 12-contract replication measuring **0.89x**
matched all-HDR (3/12 wins, 0/6 on the fresh half). That exclusion is no longer live: the operator
added a HUDEP-2 entry to `CELL_CONFIG` and `Miner.CONJUNCTION_CELL_TYPES` now names all four cells.

**The 0.89x was measured at an arm that no longer exists** — k=8 / group 80 / width 225 / light 12.
The live arm is k=8 / group 100 / width 150 / light 6, and the band-depth sweep puts HUDEP-2's best
at **k=12 + the band-scaled Cas9 cell floor**, worth 2.06x the live k=8 on own-field E[share]. So
the open question is not "was 0.89x right" but "does the conjunction at its BEST arm beat all-HDR",
which nothing has measured.

Three arms per contract, each priced against the ONE field that played that contract:

    A  all-HDR, matched          `AH.build_for_cell` on the same joined window
    B  conjunction k=8           the live arm
    C  conjunction k=12 +floor   the best measured arm

Both constructions are priced the SAME way and neither is given a modelled regime: for each build
the band set is enumerated over all 900 seeds (every row satisfying `hdr`) and so is the cut-clean
set (every row satisfying `cut`), then each regime is scored at its MEASURED consistency. That
matters because all-HDR carries a handful of cut-clean seeds beyond its band — small, but assuming
them away would stack the comparison.

    HR_N=6 python hud_resolve.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hudres")

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
from niome_subnet.genomics import all_hdr as AH          # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from niome_subnet.genomics import mt19937 as MT          # noqa: E402
from conj_stageb import DIST, OURS, API                  # noqa: E402
from sd_task import score, fetch                         # noqa: E402

CELL = "HUDEP-2"
N = int(os.getenv("HR_N", "6"))
NS = int(os.getenv("HR_NS", "5"))
SPACE = list(range(100, 1000))


def records_of(rows):
    return [{"guide": r["guideRNA"], "start": int(r["target_alignment_start"]),
             "length": int(r["target_alignment_end"]) - int(r["target_alignment_start"]),
             "strand": r["strand"], "mutation": r["mutation"],
             "cas_system": r["cas_system"]} for r in rows]


def compliant(rows, contract, cell_types, ctx, rule):
    """Every seed of 100-999 on which EVERY row satisfies `rule`."""
    ok = CJ.hdr_compliance(records_of(rows), contract, cell_types, ctx, SPACE, rule)
    MT.free_gpu_memory()
    if not ok:
        return []
    keep = set(SPACE)
    for i in sorted(ok):
        keep &= ok[i]
        if not keep:
            break
    return sorted(keep)


def price(nb, nc, v_band, v_clean, v_rest, wxf, field):
    pb, pc = nb / 900.0, max(0, nc) / 900.0
    pr = max(0.0, 1.0 - pb - pc)
    tot = 0.0
    for j in range(4):
        for i in range(4 - j):
            r = 3 - j - i
            w = (factorial(3) / (factorial(j) * factorial(i) * factorial(r))
                 * pb ** j * pc ** i * pr ** r)
            if w <= 0:
                continue
            final = wxf * (j * v_band + i * v_clean + r * v_rest) / 3.0
            rank = 1 + sum(1 for x in field if x > final)
            if rank <= 10:
                tot += w * DIST[rank - 1]
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


def evaluate(rows, contract, reference, cell_types, ctx, joined, field, label, t0):
    band = compliant(rows, contract, cell_types, ctx, "hdr")
    cut = compliant(rows, contract, cell_types, ctx, "cut")
    if not band:
        print(f"  {label}: no hdr-clean seed over 900 — unpriceable", flush=True)
        return None
    cutonly = [s for s in cut if s not in set(band)]
    sb = score(rows, contract, reference, cell_types, seed=band[0])
    wxf = sb["weighted"] * sb["fidelity"]
    rng = random.Random(0)
    vc = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(cutonly, min(NS, len(cutonly)))] if cutonly else []
    rest = [s for s in SPACE if s not in set(cut)]
    vr = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(rest, min(NS, len(rest)))]
    v_band = sb["consistency"]
    v_clean = float(np.mean(vc)) if vc else float(np.mean(vr))
    v_rest = float(np.mean(vr))
    e = price(len(band), len(cutonly), v_band, v_clean, v_rest, wxf, field)
    print(f"  {label}: band {len(band):3d} cut-clean {len(cut):3d} (+{len(cutonly)}) "
          f"wxf {wxf:6.1f} | v_band {v_band:.4f} v_clean {v_clean:.4f} v_rest {v_rest:.4f} "
          f"| E[own] {e:.6f}  ({time.monotonic()-t0:.0f}s)", flush=True)
    return {"label": label, "band": len(band), "cut": len(cut), "wxf": wxf,
            "v_band": v_band, "v_clean": v_clean, "v_rest": v_rest, "share": e}


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
        if len(cands) >= N:
            break

    spans = JW.fallback_for(CELL, "niome_hotkey")
    joined = sorted({s for a, b in spans for s in range(a, b + 1)})
    base = CJ.config_for(CELL)
    print(f"{CELL}: {len(cands)} contracts, joined {spans} ({len(joined)} of 900), "
          f"live conjunction k={base.band_k} group={base.group_size} "
          f"width={base.band_width} light={base.light_cell_rows}\n", flush=True)

    out = defaultdict(list)
    recs = []
    for t in cands:
        tid = t.get("task_id") or t["id"]
        contract = dict(t["content"]["contract"])
        reference = t["content"]["hbb_reference"]
        field = ftasks[tid]
        print(f"=== {tid[:8]}  cut10 {field[9]:.1f}  top {field[0]:.1f} ===", flush=True)
        ctx = G.build_context(contract, reference, cell_types)

        t0 = time.monotonic()
        rows, meta = AH.build_for_cell(contract, reference, cell_types, seed_list=tuple(joined))
        if rows:
            r = evaluate(rows, contract, reference, cell_types, ctx, joined, field,
                         "A all-HDR         ", t0)
            if r:
                out["A all-HDR         "].append(r["share"]); recs.append(dict(r, task=tid[:8]))
        else:
            print(f"  A all-HDR         : DECLINED — {meta.get('reason')}", flush=True)
            out["A all-HDR         "].append(0.0)

        for label, k, ca in (("B conj k=8 live   ", 8, False), ("C conj k=12 +floor", 12, True)):
            cfg = dataclasses.replace(base, seed_list=tuple(joined), band_k=k,
                                      band_cell_aware=ca)
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=1200.0)
            if not rows:
                print(f"  {label}: DECLINED — {meta.get('reason')}", flush=True)
                out[label].append(0.0)
                continue
            r = evaluate(rows, contract, reference, cell_types, ctx, joined, field, label, t0)
            out[label].append(r["share"] if r else 0.0)
            if r:
                recs.append(dict(r, task=tid[:8]))
        print(flush=True)

    print("=== aggregate E[own-field share] ===")
    ref = np.mean(out.get("A all-HDR         ") or [0.0])
    for label in ("A all-HDR         ", "B conj k=8 live   ", "C conj k=12 +floor"):
        v = out.get(label) or [0.0]
        ratio = (np.mean(v) / ref) if ref > 0 else float("nan")
        print(f"  {label}: mean {np.mean(v):.6f}  vs all-HDR {ratio:5.2f}x  "
              f"wins {sum(1 for a, b in zip(v, out['A all-HDR         ']) if a > b)}/{len(v)}")
    with open("hud_resolve.json", "w") as fh:
        json.dump(recs, fh, indent=1)


if __name__ == "__main__":
    main()
