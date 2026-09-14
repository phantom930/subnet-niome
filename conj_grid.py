#!/usr/bin/env python3
"""conj_grid.py — stage A of the conjunction config sweep: the DETERMINISTIC quantities.

Sweeps k x group_size x band sub-window width x light_cell_rows for the conjunction on the joined
300 window, one cell type per run, and records only what can be measured without sampling seeds:

    band reached, Cas12a pool, clean-set size, Cas9 survivors, rows, stage-5 cells,
    total_weighted_score, distribution_fidelity_score

`total_weighted_score` and `distribution_fidelity_score` are both SEED-INDEPENDENT -- stage 12 is
seed-independent by construction and `distribution_fidelity_score` reads only stage-12 design
fields (verified on a live submission: per-seed 315.1750..315.1750 and 0.8864..0.8864 across three
seeds). So one `score()` call per config gives both exactly, against the ~10 that a per-config
regime measurement costs. The regime values (`v_clean`, `v_rest`) are NOT measured here: at the 3
samples the earlier sweeps used they swing more than the configs differ (`conj_bandwidth.json` is
flat at w x fid 231-234 while its `e_final` ranges 33.8-45.8 on 3 samples), so they are deferred to
stage B on the shortlist, at a sample size that can actually separate arms.

Three savings make the full grid affordable:

  * ONE HDR screen over all 300 joined seeds. `hdr_compliance` returns, per bank guide, the set of
    candidate seeds it repairs by HDR on, so every sub-window is a subset of one screen.
  * `choose_band` is greedy, so the k=KMAX run passes through every shorter band on the way. One
    greedy per width yields every k.
  * `light_cell_rows` only governs how `assemble` apportions rows, so the Cas9 scan is hoisted
    above it (already true in conj_joined.py).

    python conj_grid.py
    CG_CELL=HEK293 CG_KS=4,6 CG_WIDTHS=30,300 python conj_grid.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

CELL = os.getenv("CG_CELL", "K562")
os.environ.setdefault("NIOME_INSTANCE", "conjgrid_" + CELL.replace("+", "").replace("-", "_"))

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from conj_test import hdr_compliance                        # noqa: E402
from conj_joined import JoinedCut                           # noqa: E402
from conj_feas import greedy_prefix                         # noqa: E402
from joined300 import fetch_tasks, pick_tasks               # noqa: E402
from sd_task import score                                   # noqa: E402

CLASSES = [tuple(int(y) for y in x.split("-"))
           for x in os.getenv("CG_JOINED", "100-199,400-499,700-799").split(",")]
KS = [int(x) for x in os.getenv("CG_KS", "6,8,10").split(",")]
GROUPS = [int(x) for x in os.getenv("CG_GROUPS", "42,50,80,100,125").split(",")]
WIDTHS = [int(x) for x in os.getenv("CG_WIDTHS", "30,75,100,150,225,300").split(",")]
LIGHTS = [None if x.lower() in ("none", "") else int(x)
          for x in os.getenv("CG_LIGHT", "none,6,12,25").split(",")]
POOL_TARGET = int(os.getenv("CG_POOL_TARGET", "500"))
MF = os.getenv("CG_MF", "")
START = int(os.getenv("CG_START", "700"))
OUT = os.getenv("CG_OUT", "")


def sub_window(joined, start, width):
    """`width` seeds of the joined space starting at `start`, wrapping circularly."""
    i = joined.index(start)
    return [joined[(i + j) % len(joined)] for j in range(min(width, len(joined)))]


def main():
    out_path = OUT or f"conj_grid_{CELL.replace('+', '').replace('-', '_')}.json"
    task = pick_tasks(fetch_tasks()).get(CELL)
    if task is None:
        raise SystemExit(f"no stamped task for {CELL}")
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    base = AC.config_for(CELL) or AC.AllCutConfig()
    mf = int(MF) if MF else max(1, round(base.cas12a_max_fail * len(joined) / 900))
    cfg0 = _dc.replace(JoinedCut(**{f.name: getattr(base, f.name)
                                    for f in _dc.fields(AC.AllCutConfig)}),
                       seed_list=tuple(joined), start_seed=min(joined),
                       end_seed=max(joined), cas12a_max_fail=mf)
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg0)}.npz")
    if not os.path.exists(path):
        t0 = time.monotonic()
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg0)
        MT.free_gpu_memory()
        if not bank:
            raise SystemExit("joined cut bank scan produced nothing")
        AC.save_bank(path, bank)
        print(f"built the joined cut bank: {len(bank)} guides in {time.monotonic()-t0:.0f}s")
    records = AC.load_bank(path, limit=300_000)
    tid = (task.get("task_id") or task["id"])[:8]
    print(f"{CELL}  task {tid}  acc {cell_types.get(CELL, {}).get('accessibility')}  "
          f"joined {','.join(f'{a}-{b}' for a, b in CLASSES)}  mf {mf}  rows {n_rows}  "
          f"bank {len(records)}\n")

    t0 = time.monotonic()
    ok = hdr_compliance(records, contract, cell_types, ctx, joined)
    p_hdr = sum(len(v) for v in ok.values()) / max(1, len(ok) * len(joined))
    print(f"one HDR screen over {len(joined)} joined seeds in {time.monotonic()-t0:.0f}s "
          f"| P(HDR) {p_hdr:.3f}\n")

    print(f"{'wid':>4}{'k':>3}{'grp':>5}{'light':>7}{'pool':>8}{'clean':>7}{'cas9':>8}"
          f"{'rows':>6}{'cells':>6}{'weighted':>10}{'fidelity':>10}{'w x fid':>10}{'s':>6}")
    out = []
    for width in WIDTHS:
        cand = sub_window(joined, START, width)
        okw = {i: (v & set(cand)) for i, v in ok.items()}
        steps = greedy_prefix(okw, max(KS), cand, min(GROUPS))
        for k in KS:
            if k >= len(steps):
                out.append({"cell": CELL, "task": tid, "width": width, "k": k,
                            "reason": "band unreachable"}); continue
            band, alive = steps[k]
            pool = [records[i] for i in alive]
            for g in GROUPS:
                if len(pool) < g:
                    out.append({"cell": CELL, "task": tid, "width": width, "k": k, "group": g,
                                "pool": len(pool), "reason": "pool below group"}); continue
                sel = FG.FastGreedy(pool, window_lo=min(joined), window_hi=max(joined),
                                    seeds=np.asarray(joined, dtype=np.int64))
                idx, _u = sel.best(g, restarts=12)
                grp = [pool[i] for i in idx]
                bad = set()
                for rec in grp:
                    bad.update(int(x) for x in rec["fails"])
                clean = np.array(sorted(set(joined) - bad), dtype=np.int64)
                if clean.size == 0:
                    out.append({"cell": CELL, "task": tid, "width": width, "k": k, "group": g,
                                "pool": len(pool), "reason": "clean empty"}); continue
                want = n_rows - g
                cfg_scan = _dc.replace(cfg0, group_size=g, pool_target=POOL_TARGET)
                cas9 = AC.scan_cas9(clean, contract, cell_types, ctx, sites, cfg_scan, want)
                if band and cas9:
                    ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
                    cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
                if len(cas9) < want:
                    print(f"{width:>4}{k:>3}{g:>5}{'-':>7}{len(pool):>8}{clean.size:>7}"
                          f"{len(cas9):>8}   Cas9 half short of {want}")
                    out.append({"cell": CELL, "task": tid, "width": width, "k": k, "group": g,
                                "pool": len(pool), "clean": int(clean.size), "cas9": len(cas9),
                                "want": want, "reason": "cas9 short"}); continue
                for light in LIGHTS:
                    t0 = time.monotonic()
                    cfg = _dc.replace(cfg0, group_size=g, light_cell_rows=light)
                    rows = AC.assemble(grp, cas9, contract, ctx, cfg, n_rows)
                    cells = len(Counter((r["mutation"], r["cas_system"], r["strand"])
                                        for r in rows))
                    s = score(rows, contract, reference, cell_types,
                              seed=int(band[0]) if band else int(clean[0]))
                    wf = s["weighted"] * s["fidelity"]
                    dt = time.monotonic() - t0
                    print(f"{width:>4}{k:>3}{g:>5}{str(light):>7}{len(pool):>8}{clean.size:>7}"
                          f"{len(cas9):>8}{len(rows):>6}{cells:>6}{s['weighted']:>10.2f}"
                          f"{s['fidelity']:>10.4f}{wf:>10.2f}{dt:>6.1f}")
                    out.append({"cell": CELL, "task": tid, "width": width, "k": k, "group": g,
                                "light": light, "band": band, "pool": len(pool),
                                "clean": int(clean.size), "cas9": len(cas9), "want": want,
                                "rows": len(rows), "cells": cells, "weighted": s["weighted"],
                                "fidelity": s["fidelity"], "wxfid": wf,
                                "band_cons": s["consistency"], "p_hdr": p_hdr, "mf": mf})
                    json.dump(out, open(out_path, "w"), indent=1)
                MT.free_gpu_memory()
    json.dump(out, open(out_path, "w"), indent=1)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
