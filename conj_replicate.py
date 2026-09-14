#!/usr/bin/env python3
"""conj_replicate.py — does the conjunction's stage-B verdict hold across SIX contracts?

Stage B built each arm on ONE contract per cell. `total_weighted_score` moves up to 54% with the
contract and the field moves with it, so a single contract cannot separate a construction from a
soft field -- which is exactly how `price_cell.py` put h0 on all-cut and had to be withdrawn. This
rebuilds the leading arms on six contracts per cell and prices each against **its own field**, the
one that actually played that contract, which is `cmp_k562.py`'s method and the standard CLAUDE.md
sets for this question.

Both figures are reported per contract:

    own     the arm's E[curve share] in the single field that played its contract   <- primary
    pooled  the same score placed in every field of that cell type                  <- stage B's
                                                                                       basis, kept
                                                                                       for comparison

The matched all-HDR baseline is rebuilt on each contract too, on the same joined window, so every
ratio is within-contract.

    CR_CELL=HEK293 python conj_replicate.py
    CR_CELL=K562 CR_N=3 CR_NS=8 python conj_replicate.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import random                        # noqa: E402
import statistics as st              # noqa: E402
import time                          # noqa: E402
from collections import defaultdict  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

CELL = os.getenv("CR_CELL", "HEK293")
os.environ.setdefault("NIOME_INSTANCE", "crep_" + CELL.replace("+", "").replace("-", "_"))

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from conj_test import hdr_compliance                        # noqa: E402
from conj_joined import JoinedCut                           # noqa: E402
from conj_feas import greedy_prefix                         # noqa: E402
from conj_grid import sub_window                            # noqa: E402
from conj_stageb import price, hdr_clean_seeds, DIST, NSEED, OURS, API  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402
from sd_task import score, fetch                            # noqa: E402

CLASSES = [(100, 199), (400, 499), (700, 799)]
N_CONTRACTS = int(os.getenv("CR_N", "6"))
NS = int(os.getenv("CR_NS", "12"))
POOL_TARGET = int(os.getenv("CR_POOL_TARGET", "500"))

# The leading arm per cell from stage B, plus one near-best, as (k, group, width, light).
ARMS = {
    "K562":       [(9, 80, 75, 6), (8, 80, 150, 6)],
    "HUDEP-2":    [(8, 80, 225, 12), (8, 80, 225, 25)],
    "CD34+_HSPC": [(8, 80, 100, 6), (8, 80, 100, 12)],
    "HEK293":     [(6, 80, 300, 12), (6, 100, 75, 6)],
}


def fields_by_task(cell):
    """{task_id: sorted finals} for every current-regime 3-seed field of this cell, ours removed."""
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
        if c != cell or n != 3 or r.get("miner_hotkey") in OURS:
            continue
        hk, f = r["miner_hotkey"], float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    return {tid: sorted(v.values(), reverse=True)
            for tid, v in best.items() if len(v) >= 10 and max(v.values()) > 0}


def main():
    tag = CELL.replace("+", "").replace("-", "_")
    ftasks = fields_by_task(CELL)
    all_fields = list(ftasks.values())
    items = fetch_tasks(limit=500)
    # Newest stamped tasks of this cell that also have a scored field of their own.
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
    print(f"{CELL}: replicating {len(ARMS[CELL])} arms over {len(cands)} contracts, "
          f"{NS} samples per regime")
    print(f"  own-field pricing is primary; pooled is over all {len(all_fields)} fields\n")

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    rng = random.Random(20260913)
    jset = set(joined)
    offw = [s for s in range(100, 1000) if s not in jset]
    out = []
    for t in cands:
        tid = (t.get("task_id") or t["id"])
        short = tid[:8]
        own = [ftasks[tid]]
        content = t["content"]
        contract, reference = content["contract"], content["hbb_reference"]
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
        base = AC.config_for(CELL) or AC.AllCutConfig()
        mf = max(1, round(base.cas12a_max_fail * len(joined) / 900))
        cfg0 = _dc.replace(JoinedCut(**{f.name: getattr(base, f.name)
                                        for f in _dc.fields(AC.AllCutConfig)}),
                           seed_list=tuple(joined), start_seed=min(joined),
                           end_seed=max(joined), cas12a_max_fail=mf)
        path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg0)}.npz")
        t0 = time.monotonic()
        if not os.path.exists(path):
            bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg0)
            MT.free_gpu_memory()
            if not bank:
                print(f"  {short}: cut bank produced nothing, skipped"); continue
            AC.save_bank(path, bank)
        records = AC.load_bank(path, limit=300_000)
        ok = hdr_compliance(records, contract, cell_types, ctx, joined)
        print(f"  {short}  own field cut10 {own[0][9]:>6.1f}  bank {len(records)}  "
              f"prep {time.monotonic()-t0:.0f}s")

        row = {"cell": CELL, "task": short, "own_cut10": own[0][9], "arms": []}
        for (k, g, width, light) in ARMS[CELL]:
            cand = sub_window(joined, 700, width)
            okw = {i: (v & set(cand)) for i, v in ok.items()}
            steps = greedy_prefix(okw, k, cand, g)
            if k >= len(steps):
                print(f"    k={k} g={g} w={width} l={light}: band unreachable"); continue
            band, alive = steps[k]
            pool = [records[i] for i in alive]
            if len(pool) < g:
                print(f"    k={k} g={g} w={width} l={light}: pool {len(pool)} < group"); continue
            sel = FG.FastGreedy(pool, window_lo=min(joined), window_hi=max(joined),
                                seeds=np.asarray(joined, dtype=np.int64))
            idx, _u = sel.best(g, restarts=12)
            grp = [pool[i] for i in idx]
            bad = set()
            for rec in grp:
                bad.update(int(x) for x in rec["fails"])
            clean = sorted(jset - bad)
            cfg_scan = _dc.replace(cfg0, group_size=g, pool_target=POOL_TARGET)
            cas9 = AC.scan_cas9(np.array(clean, dtype=np.int64), contract, cell_types, ctx,
                                sites, cfg_scan, n_rows - g)
            if band and cas9:
                ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
                cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
            if len(cas9) < n_rows - g:
                print(f"    k={k} g={g} w={width} l={light}: cas9 {len(cas9)} short"); continue
            cfg = _dc.replace(cfg0, group_size=g, light_cell_rows=light)
            rw = AC.assemble(grp, cas9, contract, ctx, cfg, n_rows)
            pc = [s for s in clean if s not in set(band)]
            rest = sorted(bad) + offw
            vc = [score(rw, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(pc, min(NS, len(pc)))] or [0.0]
            vr = [score(rw, contract, reference, cell_types, seed=s)["consistency"]
                  for s in rng.sample(rest, min(NS, len(rest)))]
            bs = score(rw, contract, reference, cell_types, seed=int(band[0]))
            wxf = bs["weighted"] * bs["fidelity"]
            e_own, p_own = price(len(band), len(pc), st.mean(vc), st.mean(vr), wxf, own)
            e_pool, _ = price(len(band), len(pc), st.mean(vc), st.mean(vr), wxf, all_fields)
            row["arms"].append({"k": k, "group": g, "width": width, "light": light,
                                "band": len(band), "clean": len(clean), "wxfid": wxf,
                                "v_clean": st.mean(vc), "v_rest": st.mean(vr),
                                "e_own": e_own, "e_pool": e_pool})
            print(f"    conj k={k} g={g} w={width} l={light}  band {len(band):>2} clean "
                  f"{len(clean):>3}  wxf {wxf:>6.1f}  own {e_own:.5f}  pooled {e_pool:.5f}")
            MT.free_gpu_memory()

        from niome_subnet.genomics import all_hdr as AH
        rows_ah, meta = AH.build_for_cell(contract, reference, cell_types, seed_list=joined)
        if rows_ah is None:
            print(f"    all-HDR declined: {meta.get('reason')}")
        else:
            score(rows_ah, contract, reference, cell_types, seed=500)
            aband = hdr_clean_seeds()
            if aband:
                aoff = [s for s in range(100, 1000) if s not in set(aband)]
                afl = [score(rows_ah, contract, reference, cell_types, seed=s)["consistency"]
                       for s in rng.sample(aoff, min(NS, len(aoff)))]
                abs_ = score(rows_ah, contract, reference, cell_types, seed=int(aband[0]))
                awxf = abs_["weighted"] * abs_["fidelity"]
                ae_own, _ = price(len(aband), 0, 0.0, st.mean(afl), awxf, own)
                ae_pool, _ = price(len(aband), 0, 0.0, st.mean(afl), awxf, all_fields)
                row["all_hdr"] = {"band": len(aband), "wxfid": awxf, "floor": st.mean(afl),
                                  "e_own": ae_own, "e_pool": ae_pool}
                print(f"    all-HDR    band {len(aband):>2}            wxf {awxf:>6.1f}  "
                      f"own {ae_own:.5f}  pooled {ae_pool:.5f}")
            MT.free_gpu_memory()
        out.append(row)
        json.dump(out, open(f"conj_replicate_{tag}.json", "w"), indent=1)
        print()

    print(f"===== {CELL}: {len(out)} contracts =====")
    for keyn, lbl in (("e_own", "OWN FIELD"), ("e_pool", "pooled")):
        print(f"  {lbl}")
        for i, arm in enumerate(ARMS[CELL]):
            rs = [(r, a) for r in out for a in r["arms"]
                  if (a["k"], a["group"], a["width"], a["light"]) == arm and "all_hdr" in r]
            if not rs:
                continue
            ratios = [a[keyn] / r["all_hdr"][keyn] if r["all_hdr"][keyn] > 0 else float("inf")
                      for r, a in rs]
            fin = [x for x in ratios if x != float("inf")]
            wins = sum(1 for r, a in rs if a[keyn] > r["all_hdr"][keyn])
            cs = sum(a[keyn] for r, a in rs); bs = sum(r["all_hdr"][keyn] for r, a in rs)
            print(f"    k={arm[0]} g={arm[1]} w={arm[2]} l={arm[3]}: conj {cs:.5f} vs "
                  f"all-HDR {bs:.5f} -> {cs/bs:.2f}x aggregate, wins {wins}/{len(rs)}"
                  + (f", per-contract median {st.median(fin):.2f}x" if fin else ""))
    print(f"\nwrote conj_replicate_{tag}.json")


if __name__ == "__main__":
    main()
