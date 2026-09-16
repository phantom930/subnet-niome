#!/usr/bin/env python3
"""band_hit.py — band HIT ACCURACY of the conjunction, `mh_any` band against `hdr` band.

The question this answers is whether a band rule can be *more often right* rather than merely
wider. Under a uniform seed generator those are the same quantity -- a round draws 3 seeds from
100-999 and the band is a fixed set, so

    P(exactly j of 3 seeds land in the band) = C(3,j) * p**j * (1-p)**(3-j),   p = |B| / 900

and nothing about the RULE enters except through |B|. So the only way `mh_any` can beat `hdr` on
hit accuracy is by holding a larger |B|, and the only way it can beat it on payout is by holding a
larger |B| * (what a band seed is worth).

**|B| is not the `band_k` the build pinned.** `choose_band` pins k seeds inside a width-100/300
sub-window of the joined space; the shipped rows may ALSO comply on seeds nobody searched -- the
band "leaks" (h3's 12-seed band held 5 seeds below its own window). So this measures |B| the only
honest way: replay the 250 shipped rows through the rule screen at every one of the 900 seeds and
intersect. That is the set a round seed can actually hit.

Reported per arm per contract:

    band_k      what the build pinned (6 on HEK293)
    |B900|      the full compliant set over 100-999 -- the real hit surface
    v_band      the MEASURED consistency of a band seed (1.0000 for hdr, ~0.78 for mh_any)
    hit         how many of the contract's 3 real stamped seeds landed in B900
    P(j)        exact hit distribution from |B900|
    E[cons]     expected round consistency, and P(place) against the contract's own field

    BH_N=5 BH_CELL=HEK293 python band_hit.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "bandhit")

import dataclasses
import json
import logging
import random
import time
from collections import defaultdict

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from conj_stageb import DIST, OURS, API                   # noqa: E402
from sd_task import score, fetch                          # noqa: E402

CELL = os.getenv("BH_CELL", "HEK293")
N_CONTRACTS = int(os.getenv("BH_N", "5"))
NS = int(os.getenv("BH_NS", "8"))
GROUP = os.getenv("BH_GROUP", "")
SPACE = list(range(100, 1000))
REC: list = []
OUT = os.getenv("BH_JSON", "")

ARMS = [
    ("hdr    w300 k6   ", "hdr", 300, 6, False),
    ("mh_any w100 k6   ", "mh_any", 100, 6, False),
]
EXTRA = os.getenv("BH_ARMS", "")
if EXTRA:
    ARMS = []
    for spec in EXTRA.split(","):
        parts = spec.split(":")
        rule, width, k = parts[0], parts[1], parts[2]
        ca = len(parts) > 3 and parts[3] == "ca"
        ARMS.append((f"{rule:6s} w{width} k{k}{' CA' if ca else '   '}",
                     rule, int(width), int(k), ca))


def records_of(rows):
    """The shipped submission rows, back in the shape `hdr_compliance` screens."""
    return [{"guide": r["guideRNA"], "start": int(r["target_alignment_start"]),
             "length": int(r["target_alignment_end"]) - int(r["target_alignment_start"]),
             "strand": r["strand"], "mutation": r["mutation"],
             "cas_system": r["cas_system"]} for r in rows]


def full_band(rows, contract, cell_types, ctx, rule):
    """Every seed of 100-999 on which EVERY shipped row satisfies `rule`."""
    ok = CJ.hdr_compliance(records_of(rows), contract, cell_types, ctx, SPACE, rule)
    from niome_subnet.genomics import mt19937 as MT
    MT.free_gpu_memory()
    if not ok:
        return []
    keep = set(SPACE)
    for i in sorted(ok):
        keep &= ok[i]
        if not keep:
            break
    return sorted(keep)


def hit_dist(nb, nc, ntot=900):
    """P(j band seeds, i clean-off-band seeds) over a 3-seed draw, as a flat dict."""
    pb, pc = nb / ntot, max(0, nc) / ntot
    pr = max(0.0, 1.0 - pb - pc)
    from math import factorial
    out = {}
    for j in range(4):
        for i in range(4 - j):
            r = 3 - j - i
            w = (factorial(3) / (factorial(j) * factorial(i) * factorial(r))
                 * pb ** j * pc ** i * pr ** r)
            if w > 0:
                out[(j, i)] = w
    return out


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
    print(f"{CELL}: {len(cands)} contracts, joined {spans} ({len(joined)} of 900), "
          f"group={GROUP or base.group_size} light={base.light_cell_rows}\n", flush=True)

    agg = defaultdict(list)
    for t in cands:
        tid = t.get("task_id") or t["id"]
        contract = dict(t["content"]["contract"])
        reference = t["content"]["hbb_reference"]
        seeds3 = [int(x) for x in str(contract["seed"]).split(",") if x.strip().isdigit()]
        field = ftasks[tid]
        cut10 = field[9]
        print(f"=== {tid[:8]}  seeds {seeds3}  own cut10 {cut10:.1f}  top {field[0]:.1f} ===",
              flush=True)
        ctx = G.build_context(contract, reference, cell_types)
        for label, rule, width, k, ca in ARMS:
            cfg = dataclasses.replace(base, seed_list=tuple(joined), band_rule=rule,
                                      band_width=width, band_k=k, band_cell_aware=ca)
            if GROUP:
                cfg = dataclasses.replace(cfg, group_size=int(GROUP))
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=1200.0)
            if not rows:
                print(f"  {label}: DECLINED — {meta.get('reason')}  "
                      f"({time.monotonic()-t0:.0f}s)", flush=True)
                agg[label].append(None)
                REC.append({"task": tid[:8], "arm": label, "rule": rule, "width": width,
                            "k": k, "ca": ca, "built": False, "reason": meta.get("reason"),
                            "share": 0.0, "cut10": cut10})
                continue
            pinned = [int(x) for x in meta["band_seeds"]]
            B = full_band(rows, contract, cell_types, ctx, rule)
            leak = sorted(set(B) - set(pinned))
            sb = score(rows, contract, reference, cell_types, seed=B[0])
            wxf = sb["weighted"] * sb["fidelity"]
            v_band = sb["consistency"]
            # clean-off-band and dirty values, sampled
            rng = random.Random(0)
            cl = [s for s in joined if s not in set(B)]
            vc = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(cl, min(NS, len(cl)))]
            off = [s for s in SPACE if s not in set(joined)]
            vr = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(off, min(NS, len(off)))]
            v_clean, v_rest = float(np.mean(vc)), float(np.mean(vr))
            nclean = int(meta.get("clean", 0))
            nclean_excl = max(0, nclean - len([s for s in B if s in set(joined)]))

            hits = [s for s in seeds3 if s in set(B)]
            d = hit_dist(len(B), nclean_excl)
            pj = [sum(w for (j, _i), w in d.items() if j == jj) for jj in range(4)]
            econs = sum(w * (j * v_band + i * v_clean + (3 - j - i) * v_rest) / 3.0
                        for (j, i), w in d.items())
            pplace = share = 0.0
            for (j, i), w in d.items():
                final = wxf * (j * v_band + i * v_clean + (3 - j - i) * v_rest) / 3.0
                rank = 1 + sum(1 for x in field if x > final)
                if rank <= 10:
                    pplace += w
                    share += w * DIST[rank - 1]
            agg[label].append((len(B), v_band, wxf, pj, econs, pplace, share, len(hits)))
            REC.append({"task": tid[:8], "arm": label, "rule": rule, "width": width,
                        "k": k, "ca": ca, "built": True, "B900": len(B), "v_band": v_band,
                        "wxf": wxf, "cas9": meta.get("cas9_pool"), "clean": nclean,
                        "v_clean": v_clean, "v_rest": v_rest, "econs": econs,
                        "pplace": pplace, "share": share, "cut10": cut10,
                        "real_hits": len(hits)})
            print(f"  {label}: pinned {len(pinned)} |B900| {len(B):3d} leak {len(leak)} "
                  f"cas9 {meta.get('cas9_pool')} steer {meta.get('band_steered')}/"
                  f"{meta.get('band_fell_back')} wxf {wxf:.1f} v_band {v_band:.4f} "
                  f"v_clean {v_clean:.4f} v_rest {v_rest:.4f}", flush=True)
            print(f"      P(j hits) 0:{pj[0]:.5f} 1:{pj[1]:.6f} 2:{pj[2]:.3e} 3:{pj[3]:.3e} | "
                  f"E[cons] {econs:.4f} P(place) {pplace:.5f} E[share] {share:.6f} | "
                  f"real seeds hit {len(hits)}  ({time.monotonic()-t0:.0f}s)", flush=True)
        print(flush=True)

    print("=== aggregate ===")
    for label, _r, _w, _k, _c in ARMS:
        v = [x for x in agg.get(label, []) if x]
        if not v:
            print(f"  {label}: no builds")
            continue
        B = float(np.mean([x[0] for x in v]))
        vb = float(np.mean([x[1] for x in v]))
        p1 = float(np.mean([x[3][1] for x in v]))
        p2 = float(np.mean([x[3][2] for x in v]))
        p3 = float(np.mean([x[3][3] for x in v]))
        print(f"  {label}: |B900| {B:5.1f}  v_band {vb:.4f}  P(>=1 hit) {p1+p2+p3:.5f}  "
              f"P(2) {p2:.3e}  P(3) {p3:.3e}  E[cons] {np.mean([x[4] for x in v]):.4f}  "
              f"P(place) {np.mean([x[5] for x in v]):.5f}  E[share] {np.mean([x[6] for x in v]):.6f}"
              f"  builds {len(v)}/{len(agg[label])}  real hits {sum(x[7] for x in v)}")


if __name__ == "__main__":
    try:
        main()
    finally:
        if OUT:
            with open(OUT, "w") as fh:
                json.dump(REC, fh, indent=1)
            print(f"\nwrote {len(REC)} rows to {OUT}")
