#!/usr/bin/env python3
"""conj_hdrmf.py — an HDR `max_fail` screen at the BAND stage, instead of a strict k-seed greedy.

Every conjunction arm so far picked k band seeds greedily and kept only guides HDR-compliant on
ALL of them. That caps the pool at ~P(HDR)**k and hit a ceiling at k=10 (`conj_joined_k12.json`).

This replaces that stage with the construction `all_hdr` itself uses: screen the cut bank on HDR
over the sub-window at a `max_fail` budget, then MIN-UNION the survivors' HDR failures. The band is
then the sub-window minus that union — an outcome rather than an input, and the min-union is free
to find seeds the greedy could not, because it optimises the union directly instead of choosing
seeds one at a time.

The trade is the usual one and is what the sweep is for:

    tight mf_band -> few guides qualify, but each is HDR on most of the sub-window -> wide band
    loose mf_band -> a large pool, each guide failing more seeds -> the union covers everything

Both halves of the conjunction are still measured: the HDR min-union sets the BAND (consistency
1.000), and the same group's CUT failures set the floor over the joined 300 (consistency ~0.24-0.29).
The group is selected on HDR, so the cut floor is whatever that selection leaves — the reverse of
`conj_joined.py`, where the group was selected on cut and HDR was the filter.

    python conj_hdrmf.py
    CH_SUB=30 CH_MF=12,15,18 python conj_hdrmf.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjhm")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import random                        # noqa: E402
import statistics as st              # noqa: E402
import time                          # noqa: E402
from collections import Counter, defaultdict   # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.validation import stage3         # noqa: E402
from conj_joined import JoinedCut                           # noqa: E402
from allhdr_mf import pick                                  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402
from sd_task import score                                   # noqa: E402

CLASSES = [tuple(int(y) for y in x.split("-"))
           for x in os.getenv("CH_JOINED", "100-199,400-499,700-799").split(",")]
SUB = int(os.getenv("CH_SUB", "30"))          # sub-window size, sliced from START and wrapping
START = int(os.getenv("CH_START", "700"))
GROUP = int(os.getenv("CH_GROUP", "80"))
LIGHT = None if os.getenv("CH_LIGHT", "6").lower() == "none" else int(os.getenv("CH_LIGHT", "6"))
MFS = [int(x) for x in os.getenv("CH_MF", "").split(",") if x] or None
POOL_TARGET = int(os.getenv("CH_POOL_TARGET", "200"))
FG_CAP = int(os.getenv("CH_FG_CAP", "60000"))  # guides handed to FastGreedy, fewest-fails first
N_SAMPLE = int(os.getenv("CH_SAMPLE", "3"))
OUT = os.getenv("CH_OUT", "conj_hdrmf.json")


def hdr_matrix(records, contract, cell_types, ctx, seeds):
    """Bool matrix [guide, seed] of HDR compliance. A matrix, not sets: at sub-window 300 the
    per-guide sets would hold ~45M ints across the bank."""
    acc = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    arr = np.asarray(seeds, dtype=np.int64)
    col = {int(s): j for j, s in enumerate(arr)}
    M = np.zeros((len(records), len(arr)), dtype=bool)
    groups = defaultdict(list)
    for i, rec in enumerate(records):
        groups[(rec["start"], rec["strand"], rec["cas_system"], rec["length"],
                rec["mutation"])].append(i)
    for (start, strand, cas, length, mutation), idxs in groups.items():
        site = type("S", (), {"start": start, "strand": strand, "cas": cas, "length": length})()
        distance = abs(start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        params_of = AC._params_fn(site, distance, acc, offset)
        guides = [records[i]["guide"] for i in idxs]
        got = MT.screen_guides_rule_gpu(guides, arr, mutation, cas, start, strand,
                                        params_of, "hdr", len(arr))
        for i, guide in zip(idxs, guides):
            fails = got.get(guide)
            row = np.ones(len(arr), dtype=bool)
            if fails is None:
                row[:] = False
            else:
                for f in fails:
                    row[col[int(f)]] = False
            M[i] = row
    MT.free_gpu_memory()
    return M


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    i0 = joined.index(START)
    sub = sorted({joined[(i0 + j) % len(joined)] for j in range(min(SUB, len(joined)))})
    base = AC.config_for("K562") or AC.AllCutConfig()
    mf_cut = max(1, round(base.cas12a_max_fail * len(joined) / 900))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    rng = random.Random(20260911)
    cfg0 = _dc.replace(JoinedCut(**{f.name: getattr(base, f.name)
                                    for f in _dc.fields(AC.AllCutConfig)}),
                       seed_list=tuple(joined), start_seed=min(joined), end_seed=max(joined),
                       cas12a_max_fail=mf_cut, group_size=GROUP, light_cell_rows=LIGHT)
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg0)}.npz")
    if not os.path.exists(path):
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg0)
        MT.free_gpu_memory()
        AC.save_bank(path, bank)
    records = AC.load_bank(path, limit=300_000)
    t0 = time.monotonic()
    M = hdr_matrix(records, contract, cell_types, ctx, sub)
    hdr_ok = M.sum(axis=1)
    print(f"K562 {(task.get('task_id') or task['id'])[:8]}  cut joined {len(joined)}  "
          f"sub-window {len(sub)} seeds from {START}  group {GROUP}  light {LIGHT}  "
          f"bank {len(records)}")
    print(f"HDR screen {time.monotonic()-t0:.0f}s — HDR seeds per guide over the sub-window: "
          f"mean {hdr_ok.mean():.1f}, max {hdr_ok.max()}\n")
    mfs = MFS or sorted({max(1, round(len(sub) * r)) for r in
                         (0.40, 0.50, 0.60, 0.67, 0.73, 0.80)})
    print(f"{'mf_band':>8}{'needs':>7}{'pool':>8}{'BAND':>6}{'clean':>7}{'cas9':>7}"
          f"{'weighted':>10}{'fid':>8}{'w x fid':>9}{'bandC':>7}{'cleanC':>8}"
          f"{'E[cons]':>9}{'E[fin]':>8}{'s':>6}")
    offw = [s for s in range(100, 1000) if s not in set(joined)]
    out = []
    for mfb in mfs:
        t1 = time.monotonic()
        keep = np.flatnonzero(hdr_ok >= len(sub) - mfb)
        if keep.size < GROUP:
            print(f"{mfb:>8}{len(sub)-mfb:>7}{keep.size:>8}   pool below the group")
            out.append({"mf_band": mfb, "pool": int(keep.size), "reason": "pool short"}); continue
        order = keep[np.argsort(-hdr_ok[keep])][:FG_CAP]
        pseudo = [{**records[int(i)],
                   "fails": np.asarray([s for j, s in enumerate(sub) if not M[i, j]],
                                       dtype=np.int16)} for i in order]
        sel = FG.FastGreedy(pseudo, window_lo=min(sub), window_hi=max(sub),
                            seeds=np.asarray(sub, dtype=np.int64))
        idx, _u = sel.best(GROUP, restarts=12)
        gi = [int(order[i]) for i in idx]
        hbad = set()
        for i in idx:
            hbad.update(int(x) for x in pseudo[i]["fails"])
        band = sorted(set(sub) - hbad)
        group = [records[i] for i in gi]
        cbad = set()
        for rec in group:
            cbad.update(int(x) for x in rec["fails"])
        clean = np.array(sorted(set(joined) - cbad), dtype=np.int64)
        if not band or clean.size == 0:
            print(f"{mfb:>8}{len(sub)-mfb:>7}{keep.size:>8}{len(band):>6}{clean.size:>7}"
                  f"   band or clean empty")
            out.append({"mf_band": mfb, "pool": int(keep.size), "band": len(band),
                        "clean": int(clean.size), "reason": "empty"}); continue
        cfg_scan = _dc.replace(cfg0, pool_target=POOL_TARGET)
        cas9 = AC.scan_cas9(clean, contract, cell_types, ctx, sites, cfg_scan, n_rows - GROUP)
        if cas9:
            M9 = hdr_matrix(cas9, contract, cell_types, ctx, band)
            cas9 = [cas9[i] for i in range(len(cas9)) if M9[i].all()]
        if len(cas9) < n_rows - GROUP:
            print(f"{mfb:>8}{len(sub)-mfb:>7}{keep.size:>8}{len(band):>6}{clean.size:>7}"
                  f"{len(cas9):>7}   Cas9 half short of {n_rows-GROUP}")
            out.append({"mf_band": mfb, "pool": int(keep.size), "band": len(band),
                        "clean": int(clean.size), "cas9": len(cas9),
                        "reason": "cas9 short"}); continue
        rows = AC.assemble(group, cas9, contract, ctx, cfg0, n_rows)
        cleanset = set(int(x) for x in clean)
        s_band = score(rows, contract, reference, cell_types, seed=band[0])
        pc = [s for s in cleanset if s not in set(band)]
        cs = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(pc, min(N_SAMPLE, len(pc)))] or [0.0]
        ds = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(sorted(cbad), min(N_SAMPLE, len(cbad)))] or [0.10]
        os_ = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
               for s in rng.sample(offw, min(N_SAMPLE, len(offw)))]
        nb = len([s for s in band if s in cleanset or True])
        e_cons = ((len(band) * s_band["consistency"]
                   + (len(cleanset - set(band))) * st.mean(cs)
                   + len(cbad) * st.mean(ds) + len(offw) * st.mean(os_)) / 900.0)
        wf = s_band["weighted"] * s_band["fidelity"]
        el = time.monotonic() - t1
        print(f"{mfb:>8}{len(sub)-mfb:>7}{keep.size:>8}{len(band):>6}{len(cleanset):>7}"
              f"{len(cas9):>7}{s_band['weighted']:>10.1f}{s_band['fidelity']:>8.4f}{wf:>9.2f}"
              f"{s_band['consistency']:>7.3f}{st.mean(cs):>8.3f}{e_cons:>9.4f}"
              f"{wf*e_cons:>8.2f}{el:>6.0f}")
        out.append({"mf_band": mfb, "sub": len(sub), "group": GROUP, "light": LIGHT,
                    "pool": int(keep.size), "band": len(band), "band_seeds": band,
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
