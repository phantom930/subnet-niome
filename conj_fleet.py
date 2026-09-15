#!/usr/bin/env python3
"""conj_fleet.py — a 10-hotkey conjunction fleet under ORACLE window knowledge, scored for real.

Assumes what no live miner knows: the 100-seed classes the round's three seeds fall in. The joined
cut space is exactly those classes -- 300 seeds when the three seeds land in three distinct classes,
200 when two share one, 100 when all three do -- every hotkey min-unions on cut over the WHOLE of
it, and hotkeys differ only in where their `SUB_WIDTH` HDR sub-window sits, rotated at `STRIDE`.

Unlike `conj_oracle.py` this builds the Cas9 half and ASSEMBLES, so each hotkey produces real rows
that are scored through the pipeline at the round's three actual seeds and placed in the field that
actually played that contract. That is the expensive part: the Cas9 scan runs at `pool_target` 500
because the on-band HDR filter keeps only ~P(HDR)**k of the candidates.

**Where the fleet has no decorrelation to give, this says so rather than hiding it.** When the
joined space is not wider than `SUB_WIDTH` every hotkey searches the identical candidate set and
builds the identical rows; the build is memoised on (band, clean) so that costs one build instead of
ten, and `distinct_builds` reports it.

    CF_TASKS=1 python conj_fleet.py        # pilot
    CF_TASKS=10 python conj_fleet.py
"""
import dataclasses as _dc
import json
import logging
import os
import statistics as st
import sys
import time
from collections import defaultdict

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjfleet")
logging.basicConfig(level=logging.ERROR)

import numpy as np                                          # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import conjunction as CJ         # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from conj_oracle import greedy_prefix, slice_at  # noqa: E402
from conj_replicate import fields_by_task                   # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402
from sd_task import score                                   # noqa: E402

NTASKS = int(os.getenv("CF_TASKS", "10"))
SUB_WIDTH = int(os.getenv("CF_SUB_WIDTH", "150"))
STRIDE = int(os.getenv("CF_STRIDE", "30"))
NHK = int(os.getenv("CF_NHK", "10"))
POOL_TARGET = int(os.getenv("CF_POOL_TARGET", "500"))
OUT = os.getenv("CF_OUT", "conj_fleet.json")

N_WINDOWS = int(os.getenv("CF_WINDOWS", "3"))


def windows_of(seeds, n_target=N_WINDOWS):
    """The classes the seeds occupy, PADDED out to `n_target` distinct 100-seed classes.

    Two of three seeds share a class on ~31% of rounds and all three share one on ~1.2% (chance,
    for 9 classes), which would otherwise leave a 200- or 100-seed joined space. That breaks the
    fleet in a way that has nothing to do with the construction: at a 150-wide sub-window over 100
    seeds every hotkey searches the identical candidate set and builds identical rows, so the
    10-hotkey layout silently collapses to one. Padding keeps the space at 300 so `STRIDE` 30 x 10
    hotkeys tiles it exactly once on every task.

    Padding classes contain NO round seed by construction, so they can only dilute: a hotkey whose
    sub-window lands wholly inside one wastes its band. That cost is real and is the price of a
    uniform 10-hotkey layout -- it is not a modelling artefact to be explained away.

    The filler is chosen farthest-point (maximise the minimum class-index distance to what is
    already chosen, ties to the higher index), which is deterministic and reproducible. Any other
    choice would do as well: band position is measured free under a uniform generator, and the
    filler holds no seeds either way. Note this does NOT reproduce the illustrative examples
    exactly -- for a lone 700-799 it picks 100-199 and 400-499 rather than 100-199 and 900-999.
    """
    used = sorted({(int(s) - 100) // 100 for s in seeds})
    chosen, pool = list(used), [c for c in range(9) if c not in used]
    while len(chosen) < n_target and pool:
        best = max(pool, key=lambda c: (min(abs(c - x) for x in chosen), c))
        chosen.append(best)
        pool.remove(best)
    chosen.sort()
    return sorted(s for c in chosen for s in range(c * 100 + 100, c * 100 + 200)), len(used), chosen


# (group, k, light) per cell type, as specified.
ARM = {
    "CD34+_HSPC": (125, 10, 6),
    "K562": (125, 10, 6),
    "HUDEP-2": (125, 10, 6),
    "HEK293": (80, 8, 12),
}


def build_one(band, clean, group, contract, cell_types, ctx, sites, cfg, n_rows):
    """Cas9 half + assemble for one (band, clean, group) -> (rows, None) or (None, reason)."""
    cas9 = AC.scan_cas9(np.array(sorted(clean), dtype=np.int64), contract, cell_types, ctx, sites,
                        _dc.replace(cfg, pool_target=POOL_TARGET), n_rows - cfg.group_size)
    MT.free_gpu_memory()
    n_cut = len(cas9)
    if band and cas9:
        ok9 = CJ.hdr_compliance(cas9, contract, cell_types, ctx, sorted(band))
        cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
        MT.free_gpu_memory()
    cells = len({(r["mutation"], r["strand"]) for r in cas9})
    if len(cas9) < n_rows - cfg.group_size or cells < 4:
        return None, f"cas9 {len(cas9)}/{n_cut} over {cells} cells short of {n_rows-cfg.group_size}"
    rows = AC.assemble(group, cas9, contract, ctx, cfg, n_rows)
    if len(rows) < n_rows:
        return None, f"assembled {len(rows)} rows"
    return rows, None


def run_task(task, cell_types, fields):
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell = contract.get("cell_type")
    tid = (task.get("task_id") or task["id"])
    seeds = sorted(int(x) for x in str(contract.get("seed", "")).split(",") if x.strip().isdigit())
    joined, nwin, classes = windows_of(seeds)
    group, k, light = ARM[cell]
    base = CJ.config_for(cell) or CJ.ConjunctionConfig()
    mf = max(1, round(base.cas12a_max_fail * len(joined) / 900))
    cfg = _dc.replace(base, seed_list=tuple(joined), start_seed=joined[0], end_seed=joined[-1],
                      cas12a_max_fail=mf, group_size=group, light_cell_rows=light, band_k=k)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    os.makedirs(AC.BANK_DIR, exist_ok=True)
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    t0 = time.monotonic()
    if not os.path.exists(path):
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        MT.free_gpu_memory()
        if not bank:
            return {"task": tid[:8], "cell": cell, "reason": f"cut bank empty at mf {mf}"}
        AC.save_bank(path, bank)
    records = AC.load_bank(path, limit=cfg.bank_keep)
    ok = CJ.hdr_compliance(records, contract, cell_types, ctx, joined)
    MT.free_gpu_memory()
    t_prep = time.monotonic() - t0

    jset, sset = set(joined), set(seeds)
    width = min(SUB_WIDTH, len(joined))
    memo, hks = {}, []
    t0 = time.monotonic()
    for h in range(NHK):
        cand = set(slice_at(joined, (h * STRIDE) % len(joined), width))
        steps = greedy_prefix({i: (v & cand) for i, v in ok.items()}, k, sorted(cand), group)
        if k >= len(steps):
            hks.append({"hk": h, "reason": f"band reached {len(steps)-1} of {k}"}); continue
        band, alive = steps[k]
        if len(alive) < group:
            hks.append({"hk": h, "reason": f"pool {len(alive)} < group {group}"}); continue
        pool = [records[i] for i in alive]
        sel = FG.FastGreedy(pool, window_lo=joined[0], window_hi=joined[-1],
                            seeds=np.asarray(joined, dtype=np.int64))
        idx, _u = sel.best(group, restarts=cfg.restarts)
        gidx = [alive[i] for i in idx]
        bad = set()
        for r in (records[i] for i in gidx):
            bad.update(int(x) for x in r["fails"])
        clean = jset - bad
        key = (tuple(sorted(band)), tuple(sorted(clean)))
        if key not in memo:
            rows, why = build_one(band, clean, [records[i] for i in gidx], contract, cell_types,
                                  ctx, sites, cfg, n_rows)
            if rows is None:
                memo[key] = {"reason": why}
            else:
                s = [score(rows, contract, reference, cell_types, seed=int(x)) for x in seeds]
                cons = st.mean(x["consistency"] for x in s)
                memo[key] = {"weighted": s[0]["weighted"], "fidelity": s[0]["fidelity"],
                             "cons_per_seed": [round(x["consistency"], 4) for x in s],
                             "cons": cons,
                             "final": s[0]["weighted"] * cons * s[0]["fidelity"]}
                MT.free_gpu_memory()
        hks.append({"hk": h, "band": sorted(band), "clean": sorted(clean),
                    "clean_n": len(clean),
                    "band_hits": sorted(set(band) & sset),
                    "clean_hits": sorted(clean & sset), **memo[key]})
    t_build = time.monotonic() - t0

    ub = set().union(*[set(x.get("band", [])) for x in hks]) if hks else set()
    uc = set().union(*[set(x.get("clean", [])) for x in hks]) if hks else set()
    fld = fields.get(tid)
    built = [x for x in hks if "final" in x]
    out = {"task": tid[:8], "cell": cell, "at": task.get("created_at", "")[:16], "seeds": seeds,
           "n_windows": nwin, "classes": [f"{c*100+100}-{c*100+199}" for c in classes],
           "padded": N_WINDOWS - nwin, "span": len(joined), "mf": mf, "group": group, "k": k,
           "light": light, "bank": len(records), "sub_width": width,
           "prep_s": round(t_prep, 1), "build_s": round(t_build, 1),
           "distinct_builds": len(memo), "hotkeys": hks,
           "band_total": sorted(ub), "band_hits": sorted(ub & sset),
           "clean_total_n": len(uc), "clean_hits": sorted(uc & sset),
           "field_top10": fld[:10] if fld else None,
           "best_final": max((x["final"] for x in built), default=None),
           "n_built": len(built)}
    return out


def main():
    items = fetch_tasks(limit=500)
    cands = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        sd = [x for x in str(c.get("seed", "") or "").split(",") if x.strip().isdigit()]
        if len(sd) == 3 and c.get("cell_type") in ARM:
            cands.append(t)
        if len(cands) >= NTASKS:
            break
    fields = {}
    for cell in ARM:
        fields.update(fields_by_task(cell))
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    print(f"{len(cands)} tasks | {NHK} hotkeys, sub-window {SUB_WIDTH} stride {STRIDE}\n")
    out = []
    for n, t in enumerate(cands, 1):
        r = run_task(t, cell_types, fields)
        out.append(r)
        json.dump(out, open(OUT, "w"), indent=1)
        if "reason" in r:
            print(f"  [{n}] {r['task']} {r.get('cell')}: {r['reason']}")
            continue
        f10 = r["field_top10"]
        print(f"  [{n}/{len(cands)}] {r['task']} {r['cell']:<11} seeds {r['seeds']} win "
              f"{r['n_windows']} | built {r['n_built']}/{NHK} distinct {r['distinct_builds']} "
              f"| band {len(r['band_total'])} hit {len(r['band_hits'])} "
              f"clean {r['clean_total_n']} hit {len(r['clean_hits'])} "
              f"| best {r['best_final'] or 0:.1f} vs cut10 "
              f"{f10[9] if f10 and len(f10) >= 10 else float('nan'):.1f} "
              f"top1 {f10[0] if f10 else float('nan'):.1f} "
              f"| prep {r['prep_s']:.0f}s build {r['build_s']:.0f}s")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
