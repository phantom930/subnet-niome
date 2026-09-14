#!/usr/bin/env python3
"""conj_feas.py — which (k, group) the conjunction can actually FILL, per cell type.

The conjunction filters the Cas9 half to guides that repair by HDR on every one of the k band
seeds, so the surviving pool decays as P(HDR)**k — ~0.57 on the erythroid cells, ~0.37 on HEK293.
`assemble` needs `n_rows - group` Cas9 rows, so there is a hard k above which no config fills,
and it differs per cell. This measures that ceiling before a scoring grid is worth running.

Cheap by construction: ONE `scan_cas9` per cell (at k=0, the widest clean set, so the survivor
counts below are an UPPER bound), then ONE GPU HDR screen of those candidates over the band
candidate seeds. `choose_band` is greedy and therefore prefix-nested in k, so every k comes from
the same screen with no extra work.

    python conj_feas.py
    CF_CELLS=K562 CF_GROUPS=80 python conj_feas.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjfeas")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from conj_test import hdr_compliance                        # noqa: E402
from conj_joined import JoinedCut                           # noqa: E402
from joined300 import fetch_tasks, pick_tasks               # noqa: E402

CELLS = os.getenv("CF_CELLS", "K562,HUDEP-2,CD34+_HSPC,HEK293").split(",")
CLASSES = [tuple(int(y) for y in x.split("-"))
           for x in os.getenv("CF_JOINED", "100-199,400-499,700-799").split(",")]
GROUPS = [int(x) for x in os.getenv("CF_GROUPS", "42,50,80,100,125").split(",")]
KS = [int(x) for x in os.getenv("CF_KS", "0,2,4,6,8,10").split(",")]
POOL_TARGET = int(os.getenv("CF_POOL_TARGET", "200"))
OUT = os.getenv("CF_OUT", "conj_feas.json")


def greedy_prefix(ok, kmax, candidates, need):
    """Greedy band to depth kmax, returning (band, alive) at EVERY prefix length."""
    col = {s: j for j, s in enumerate(candidates)}
    idx = sorted(ok)
    M = np.zeros((len(idx), len(candidates)), dtype=bool)
    for r, i in enumerate(idx):
        for s in ok[i]:
            if s in col:
                M[r, col[s]] = True
    alive = np.ones(len(idx), dtype=bool)
    band, taken, steps = [], set(), [([], [idx[r] for r in np.flatnonzero(alive)])]
    for _ in range(kmax):
        counts = M[alive].sum(axis=0)
        for j in taken:
            counts[j] = -1
        j = int(np.argmax(counts))
        if counts[j] < need:
            break
        taken.add(j)
        band.append(candidates[j])
        alive &= M[:, j]
        steps.append((list(band), [idx[r] for r in np.flatnonzero(alive)]))
    return steps


def main():
    tasks = pick_tasks(fetch_tasks())
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    out = []
    for cell in CELLS:
        task = tasks.get(cell)
        if task is None:
            print(f"{cell}: no stamped task in the listing, skipped\n")
            continue
        content = task["content"]
        contract, reference = content["contract"], content["hbb_reference"]
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
        base = AC.config_for(cell) or AC.AllCutConfig()
        mf = max(1, round(base.cas12a_max_fail * len(joined) / 900))
        acc = cell_types.get(cell, {}).get("accessibility")
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
                print(f"{cell}: cut bank scan produced nothing\n"); continue
            AC.save_bank(path, bank)
            print(f"{cell}: built the joined cut bank, {len(bank)} guides in "
                  f"{time.monotonic()-t0:.0f}s")
        records = AC.load_bank(path, limit=300_000)
        tid = (task.get("task_id") or task["id"])[:8]
        print(f"=== {cell}  task {tid}  acc {acc}  mf {mf}  rows {n_rows}  "
              f"bank {len(records)} Cas12a guides ===")

        t0 = time.monotonic()
        ok = hdr_compliance(records, contract, cell_types, ctx, joined)
        p_hdr = sum(len(v) for v in ok.values()) / max(1, len(ok) * len(joined))
        print(f"  HDR screen over {len(joined)} joined seeds in {time.monotonic()-t0:.0f}s "
              f"| per-seed P(HDR) on bank guides = {p_hdr:.3f}")
        steps = greedy_prefix(ok, max(KS), joined, min(GROUPS))
        print(f"  greedy band reaches k={len(steps)-1} of {max(KS)} "
              f"(Cas12a pool at each k: {[len(s[1]) for s in steps]})")

        for g in GROUPS:
            kmax_pool = max(k for k in range(len(steps)) if len(steps[k][1]) >= g)
            band0, alive0 = steps[0]
            pool0 = [records[i] for i in alive0]
            sel = FG.FastGreedy(pool0, window_lo=min(joined), window_hi=max(joined),
                                seeds=np.asarray(joined, dtype=np.int64))
            idx, _u = sel.best(g, restarts=12)
            grp = [pool0[i] for i in idx]
            bad = set()
            for rec in grp:
                bad.update(int(x) for x in rec["fails"])
            clean = np.array(sorted(set(joined) - bad), dtype=np.int64)
            want = n_rows - g
            cfg_scan = _dc.replace(cfg0, group_size=g, pool_target=POOL_TARGET)
            t0 = time.monotonic()
            cas9 = AC.scan_cas9(clean, contract, cell_types, ctx, sites, cfg_scan, want)
            ts = time.monotonic() - t0
            ok9 = hdr_compliance(cas9, contract, cell_types, ctx, joined) if cas9 else {}
            print(f"  group {g:>3}  clean {clean.size:>3}/{len(joined)}  want {want:>3}  "
                  f"cas9(k=0) {len(cas9):>6} in {ts:>4.0f}s  pool-limited kmax {kmax_pool}")
            for k in KS:
                if k >= len(steps):
                    print(f"     k={k:<3} band unreachable (Cas12a pool)"); continue
                band = steps[k][0]
                surv = sum(1 for i in ok9 if all(s in ok9[i] for s in band)) if band else len(cas9)
                fills = surv >= want
                print(f"     k={k:<3} band {str(band[:4]) + ('...' if len(band) > 4 else ''):<26}"
                      f" cas9 survivors {surv:>6}  need {want:>3}  "
                      f"{'FILLS' if fills else 'SHORT'}")
                out.append({"cell": cell, "task": tid, "group": g, "k": k,
                            "clean_k0": int(clean.size), "want": want, "cas9_k0": len(cas9),
                            "survivors": surv, "fills": bool(fills),
                            "cas12a_pool": len(steps[k][1]), "p_hdr": p_hdr,
                            "band": band, "scan_s": round(ts, 1)})
                json.dump(out, open(OUT, "w"), indent=1)
            MT.free_gpu_memory()
        print()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
