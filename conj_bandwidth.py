#!/usr/bin/env python3
"""conj_bandwidth.py — the conjunction's HDR-band SUB-WINDOW width, at group 80 on K562.

`conj_joined.py` fixed the band sub-window at 30 seeds and swept k. This sweeps the sub-window
itself — 30 / 75 / 100 / 150 / 225 / 300 seeds — holding the cut window (joined 300), the group
(80) and the band width (k) fixed.

What the sub-window controls is CHOICE, not the band: the band is always k seeds, but a wider
sub-window gives the greedy more candidate seeds to choose among, so it can pick the k seeds that
the largest number of bank guides already satisfy. A wider window should therefore leave a LARGER
surviving pool and so a better cut-clean floor — the opposite of the all-HDR case, where widening
the window shrinks the band.

The sub-window is a slice of the joined space starting at 700 and wrapping, so width 30 is exactly
`conj_joined.py`'s 700-729 and the results are directly comparable. Widths above 100 are
necessarily non-contiguous, since the joined space is three disjoint width-100 classes.

Two changes from conj_joined.py, both forced by the wider windows:

  * `choose_band` there is O(candidates x pool) in Python — 300 candidates against a 300,000-guide
    pool is 90M set lookups per pick. Replaced with a boolean matrix and vectorised column sums.
  * candidates are every seed of the sub-window rather than a fixed sample of 30, which is what
    makes "width" mean anything here.

    python conj_bandwidth.py
    CB_WIDTHS=30,300 CB_K=6 python conj_bandwidth.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjbw")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import random                        # noqa: E402
import statistics as st              # noqa: E402
import time                          # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from conj_joined import JoinedCut                           # noqa: E402
from conj_test import hdr_compliance                        # noqa: E402
from allhdr_mf import pick                                  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402
from sd_task import score                                   # noqa: E402

CLASSES = [tuple(int(y) for y in x.split("-"))
           for x in os.getenv("CB_JOINED", "100-199,400-499,700-799").split(",")]
WIDTHS = [int(x) for x in os.getenv("CB_WIDTHS", "30,75,100,150,225,300").split(",")]
K = int(os.getenv("CB_K", "6"))
GROUP = int(os.getenv("CB_GROUP", "80"))
LIGHT = None if os.getenv("CB_LIGHT", "6").lower() == "none" else int(os.getenv("CB_LIGHT", "6"))
POOL_TARGET = int(os.getenv("CB_POOL_TARGET", "200"))
N_SAMPLE = int(os.getenv("CB_SAMPLE", "3"))
START = int(os.getenv("CB_START", "700"))
OUT = os.getenv("CB_OUT", "conj_bandwidth.json")


def fast_choose_band(ok, n_records, candidates, k, need):
    """Greedy k seeds maximising the surviving pool, vectorised.

    `ok` maps guide index -> set of compliant seeds. Materialised as a bool matrix so each greedy
    step is a column sum over the alive rows instead of a Python set intersection per candidate.
    """
    col = {s: j for j, s in enumerate(candidates)}
    M = np.zeros((n_records, len(candidates)), dtype=bool)
    for i, seeds in ok.items():
        for s in seeds:
            M[i, col[s]] = True
    alive = np.ones(n_records, dtype=bool)
    band, taken = [], set()
    for _ in range(k):
        counts = M[alive].sum(axis=0)
        for j in taken:
            counts[j] = -1
        j = int(np.argmax(counts))
        if counts[j] < need:
            break
        taken.add(j)
        band.append(candidates[j])
        alive &= M[:, j]
    return band, np.flatnonzero(alive)


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    base = AC.config_for("K562") or AC.AllCutConfig()
    mf = max(1, round(base.cas12a_max_fail * len(joined) / 900))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    rng = random.Random(20260911)
    cfg0 = _dc.replace(JoinedCut(**{f.name: getattr(base, f.name)
                                    for f in _dc.fields(AC.AllCutConfig)}),
                       seed_list=tuple(joined), start_seed=min(joined), end_seed=max(joined),
                       cas12a_max_fail=mf, group_size=GROUP, light_cell_rows=LIGHT)
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg0)}.npz")
    if not os.path.exists(path):
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg0)
        MT.free_gpu_memory()
        AC.save_bank(path, bank)
    records = AC.load_bank(path, limit=300_000)
    start_i = joined.index(START)
    print(f"K562  task {(task.get('task_id') or task['id'])[:8]}  cut window joined "
          f"{len(joined)} seeds  group {GROUP}  k {K}  light {LIGHT}  bank {len(records)}\n")
    print(f"{'bandw':>6}{'cands':>7}{'pool':>8}{'union':>7}{'clean':>7}{'cas9':>7}{'weighted':>10}"
          f"{'fid':>8}{'w x fid':>9}{'bandC':>7}{'cleanC':>8}{'E[cons]':>9}{'E[final]':>10}{'s':>6}")
    jset, out = set(joined), []
    offw = [s for s in range(100, 1000) if s not in jset]
    for width in WIDTHS:
        t0 = time.monotonic()
        cands = [joined[(start_i + j) % len(joined)] for j in range(min(width, len(joined)))]
        cands = sorted(set(cands))
        ok = hdr_compliance(records, contract, cell_types, ctx, cands)
        band, alive = fast_choose_band(ok, len(records), cands, K, GROUP)
        if len(band) < K:
            print(f"{width:>6}{len(cands):>7}   band of {K} unreachable")
            out.append({"band_width": width, "reason": "band unreachable"}); continue
        pool = [records[int(i)] for i in alive]
        sel = FG.FastGreedy(pool, window_lo=min(joined), window_hi=max(joined),
                            seeds=np.asarray(joined, dtype=np.int64))
        idx, _u = sel.best(GROUP, restarts=12)
        group = [pool[i] for i in idx]
        bad = set()
        for rec in group:
            bad.update(int(x) for x in rec["fails"])
        clean = np.array(sorted(jset - bad), dtype=np.int64)
        cfg_scan = _dc.replace(cfg0, pool_target=POOL_TARGET)
        cas9 = AC.scan_cas9(clean, contract, cell_types, ctx, sites, cfg_scan, n_rows - GROUP)
        if cas9:
            ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
            cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
        if len(cas9) < n_rows - GROUP:
            print(f"{width:>6}{len(cands):>7}{len(pool):>8}{len(bad):>7}{clean.size:>7}"
                  f"{len(cas9):>7}   Cas9 half short of {n_rows-GROUP}")
            out.append({"band_width": width, "pool": len(pool), "clean": int(clean.size),
                        "cas9": len(cas9), "reason": "cas9 short"}); continue
        rows = AC.assemble(group, cas9, contract, ctx, cfg0, n_rows)
        cleanset = set(int(x) for x in clean)
        s_band = score(rows, contract, reference, cell_types, seed=band[0])
        pc = [s for s in cleanset if s not in set(band)]
        cs = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(pc, min(N_SAMPLE, len(pc)))] or [0.0]
        ds = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(sorted(bad), min(N_SAMPLE, len(bad)))] or [0.10]
        os_ = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
               for s in rng.sample(offw, min(N_SAMPLE, len(offw)))]
        e_cons = ((len(band) * s_band["consistency"] + (len(cleanset) - len(band)) * st.mean(cs)
                   + len(bad) * st.mean(ds) + len(offw) * st.mean(os_)) / 900.0)
        wf = s_band["weighted"] * s_band["fidelity"]
        el = time.monotonic() - t0
        print(f"{width:>6}{len(cands):>7}{len(pool):>8}{len(bad):>7}{len(cleanset):>7}"
              f"{len(cas9):>7}{s_band['weighted']:>10.1f}{s_band['fidelity']:>8.4f}{wf:>9.2f}"
              f"{s_band['consistency']:>7.3f}{st.mean(cs):>8.3f}{e_cons:>9.4f}"
              f"{wf*e_cons:>10.2f}{el:>6.0f}")
        out.append({"band_width": width, "candidates": len(cands), "k": K, "group": GROUP,
                    "light": LIGHT, "band": band, "pool": len(pool), "union": len(bad),
                    "clean": len(cleanset), "cas9": len(cas9), "rows": len(rows),
                    "cells": len(Counter((r["mutation"], r["cas_system"], r["strand"])
                                         for r in rows)),
                    "weighted": s_band["weighted"], "fidelity": s_band["fidelity"], "wxfid": wf,
                    "band_cons": s_band["consistency"], "clean_cons": st.mean(cs),
                    "dirty_cons": st.mean(ds), "offwindow_cons": st.mean(os_),
                    "e_cons": e_cons, "e_final": wf * e_cons, "build_s": round(el, 1)})
        json.dump(out, open(OUT, "w"), indent=1)
        MT.free_gpu_memory()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
