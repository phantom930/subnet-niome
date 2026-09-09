#!/usr/bin/env python3
"""cmp_k562.py — all-cut vs all-HDR on K562, each priced against its OWN contract's field.

This is the measurement CLAUDE.md's "Pricing a construction" section asks for. The live config has
h0 on all-cut for K562 on the strength of E[share] 0.0023 against a POOL of 19 fields — but that
score was measured on one contract (83f430e9, weighted 304.6) and the pool contains softer
contracts' fields, so it prices field softness as well as the construction. Against its own field
the same build earns 0.0000.

Fixing it needs both arms on the SAME contracts:

    all-cut   whole window 100-999, group 42     clean set of ~549-570 at cons ~0.26, else ~0.10
    all-HDR   band inside hdr_range, group 80    band of ~13 at cons 1.000, else the ~0.126 floor

Each arm is then a regime lottery over the round's three seeds, and every possible round score is
placed in the one real field that played that contract. Per-contract measurements are cached to
cmp_k562_<task8>.json, so a re-run only builds what is missing.

    python cmp_k562.py <task_id> [task_id ...]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "cmpk562")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import random                        # noqa: E402
import statistics as st              # noqa: E402
import time                          # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                   # noqa: E402
from conj_price import lottery, share   # noqa: E402
from niome_subnet.genomics import all_cut as AC    # noqa: E402
from niome_subnet.genomics import all_hdr as AH    # noqa: E402
from niome_subnet.genomics import fastgreedy as FG  # noqa: E402
from niome_subnet.genomics import mt19937 as MT    # noqa: E402
from sd_task import score, task_content            # noqa: E402

N_SAMPLE = int(os.getenv("CMP_SAMPLE", "12"))
BAND_HUNT = int(os.getenv("CMP_BAND_HUNT", "15"))


def own_field(task_id):
    from sd_task import OURS
    rows = json.load(open("sd_task_scores.json"))
    rows = rows if isinstance(rows, list) else (rows.get("data") or rows.get("items") or [])
    best = {}
    for r in rows:
        if r.get("task_id") == task_id and r.get("miner_hotkey") not in OURS:
            h = r["miner_hotkey"]
            best[h] = max(float(r.get("final_score") or 0), best.get(h, 0.0))
    return sorted(best.values(), reverse=True)


def measure_all_cut(contract, reference, cell_types, rng):
    """Whole-window all-cut at the shipped config: clean set, and cons clean vs dirty."""
    cell = contract.get("cell_type")
    cfg = AC.config_for(cell) or AC.AllCutConfig()
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    t0 = time.monotonic()
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        if not bank:
            return {"reason": "bank scan produced nothing"}
        AC.save_bank(path, bank)
    records = AC.load_bank(path)
    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed)
    idx, _u = sel.best(cfg.group_size, restarts=12)
    group = [records[i] for i in idx]
    bad = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    clean = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
    if not clean:
        return {"reason": "no clean seed"}
    cas9 = AC.scan_cas9(np.array(clean, dtype=np.int64), contract, cell_types, ctx, sites, cfg,
                        n_rows - cfg.group_size)
    if len(cas9) < n_rows - cfg.group_size:
        return {"reason": f"cas9 pool {len(cas9)} short of {n_rows - cfg.group_size}"}
    rows = AC.assemble(group, cas9, contract, ctx, cfg, n_rows)
    base = score(rows, contract, reference, cell_types, seed=clean[0])
    cs = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(clean, min(N_SAMPLE, len(clean)))]
    ds = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(sorted(bad), min(N_SAMPLE, len(bad)))] or [0.10]
    MT.free_gpu_memory()
    return {"method": "all-cut", "group": cfg.group_size, "hit": len(clean),
            "union": len(bad), "cells": len(Counter(
                (r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
            "weighted": base["weighted"], "fidelity": base["fidelity"],
            "hit_cons": st.mean(cs), "floor_cons": st.mean(ds),
            "build_s": round(time.monotonic() - t0, 1)}


def measure_all_hdr(contract, reference, cell_types, rng):
    """all-HDR at the shipped config: band size, the 1.000 spike, and the off-band floor."""
    t0 = time.monotonic()
    rows, meta = AH.build_for_cell(contract, reference, cell_types, budget_s=900)
    if not rows:
        return {"reason": meta.get("reason", "declined"), "meta": meta}
    cfg = AH.config_for(contract.get("cell_type"))
    lo, hi = cfg.hdr_range
    # The band is the subset of hdr_range where EVERY row repairs by HDR; hunt one to verify the
    # spike is really 1.000 rather than assuming it.
    band_cons, band_seed = None, None
    for s in rng.sample(range(lo, hi + 1), min(BAND_HUNT, hi - lo + 1)):
        got = score(rows, contract, reference, cell_types, seed=s)
        if got["consistency"] > 0.99:
            band_cons, band_seed, base = got["consistency"], s, got
            break
    off = [s for s in range(100, 1000) if not (lo <= s <= hi)]
    fs = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(off, min(N_SAMPLE, len(off)))]
    if band_cons is None:                       # band too narrow to find by sampling
        base = score(rows, contract, reference, cell_types, seed=lo)
    MT.free_gpu_memory()
    return {"method": "all-HDR", "group": cfg.group_size, "band_range": [lo, hi],
            "hit": meta.get("clean"), "union": meta.get("union"),
            "cells": meta.get("cells"), "cas9_pool": meta.get("cas9_pool"),
            "weighted": base["weighted"], "fidelity": base["fidelity"],
            "hit_cons": band_cons if band_cons else 1.0,
            "hit_cons_verified": band_cons is not None, "band_seed": band_seed,
            "floor_cons": st.mean(fs), "build_s": round(time.monotonic() - t0, 1)}


def price(arm, field):
    """Expected curve share of this arm in the one field that played its contract."""
    hit, tot = arm["hit"], 900
    lot = lottery((hit, tot - hit), (arm["hit_cons"], arm["floor_cons"]))
    wtd, fid = arm["weighted"], arm["fidelity"]
    e_final = sum(p * wtd * c * fid for p, c in lot)
    e_share = sum(p * share(wtd * c * fid, field) for p, c in lot)
    p_place = sum(p for p, c in lot if share(wtd * c * fid, field) > 0)
    return {"e_final": e_final, "e_share": e_share, "p_place": p_place,
            "best_case": max(wtd * c * fid for _p, c in lot)}


def main(tasks):
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    rng = random.Random(29)
    out = []
    for tid in tasks:
        cache = f"cmp_k562_{tid[:8]}.json"
        if os.path.exists(cache):
            rec = json.load(open(cache))
            print(f"{tid[:8]}  cached")
        else:
            task, contract, reference = task_content(tid)
            print(f"{tid[:8]}  {contract.get('cell_type')}  seeds {contract.get('seed')}")
            rec = {"task": tid, "seeds": contract.get("seed"),
                   "all_cut": measure_all_cut(contract, reference, cell_types, rng),
                   "all_hdr": measure_all_hdr(contract, reference, cell_types, rng)}
            json.dump(rec, open(cache, "w"), indent=1)
        rec["field"] = own_field(tid)
        out.append(rec)
        for k in ("all_cut", "all_hdr"):
            a = rec[k]
            if "reason" in a:
                print(f"    {k:<8} declined: {a['reason']}")
            else:
                print(f"    {a['method']:<8} g{a['group']:<3} hit {a['hit']:>4} of 900  "
                      f"weighted {a['weighted']:>6.1f}  fid {a['fidelity']:.3f}  "
                      f"hit cons {a['hit_cons']:.4f}  floor {a['floor_cons']:.4f}  "
                      f"{a['build_s']:>6.1f}s")

    print(f"\n{'':<10} {'cut10':>7}  {'all-cut':>28}   {'all-HDR':>28}")
    print(f"  {'task':<8} {'':>7}  {'E[final]':>8} {'P(pl)':>6} {'E[shr]':>7} {'best':>5}   "
          f"{'E[final]':>8} {'P(pl)':>6} {'E[shr]':>7} {'best':>5}   winner")
    tot = {"all-cut": 0.0, "all-HDR": 0.0}
    n = 0
    for rec in out:
        f = rec.get("field") or []
        if len(f) < 50 or "reason" in rec["all_cut"] or "reason" in rec["all_hdr"]:
            continue
        pc, ph = price(rec["all_cut"], f), price(rec["all_hdr"], f)
        tot["all-cut"] += pc["e_share"]
        tot["all-HDR"] += ph["e_share"]
        n += 1
        win = "all-cut" if pc["e_share"] > ph["e_share"] else (
            "all-HDR" if ph["e_share"] > pc["e_share"] else "tie(0)")
        print(f"  {rec['task'][:8]:<8} {f[9] if len(f) >= 10 else 0:>7.1f}  "
              f"{pc['e_final']:>8.2f} {pc['p_place']:>5.1%} {pc['e_share']:>7.4f} "
              f"{pc['best_case']:>5.0f}   "
              f"{ph['e_final']:>8.2f} {ph['p_place']:>5.1%} {ph['e_share']:>7.4f} "
              f"{ph['best_case']:>5.0f}   {win}")
    if n:
        print(f"\n  mean over {n} matched contracts:  all-cut {tot['all-cut']/n:.4f}   "
              f"all-HDR {tot['all-HDR']/n:.4f}")
        a, h = tot["all-cut"] / n, tot["all-HDR"] / n
        if max(a, h) > 0:
            better = "all-cut" if a > h else "all-HDR"
            fac = f"{max(a,h)/min(a,h):.2f}x" if min(a, h) > 0 else "infinitely (other is 0)"
            print(f"  -> {better} better by {fac} on expected curve share")
    json.dump([{k: v for k, v in r.items() if k != "field"} for r in out],
              open("cmp_k562.json", "w"), indent=1)
    print("\nwrote cmp_k562.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(_ARGV[1:] or ["83f430e9-1c9a-48b4-b594-e17477a046ea"]))
