#!/usr/bin/env python3
"""conj_joined.py — all-cut over a JOINED 300 window AND an HDR band in a width-30 sub-window.

The conjunction CLAUDE.md records as falsified was measured at CONTIGUOUS cut windows of 100 / 300
/ 900 and groups 42 / 50 / 80. This runs it over the joined 300-seed space the fleet actually uses,
with the HDR band confined to a 30-seed sub-window inside it — the width that `allhdr_cas12a.py`
measured as the band optimum — and sweeps group size across 42 / 50 / 80 / 100 / 125.

The construction, in order:

  1. bank Cas12a guides on the `cut` rule over the JOINED 300 seeds
  2. keep only those that also repair by HDR on every seed of a k-seed band, chosen greedily
     inside the width-30 sub-window
  3. min-union the survivors on cut -> the clean set (cut-pinned, cons ~0.24)
  4. fill Cas9 strictly on cut over that clean set AND HDR over the band

A band seed then pins all three stage-4 targets (cons 1.0) while the rest of the clean set pins
`is_cut` alone — the "spike on an elevated floor" that no construction measured here has produced.
The cost is pool shrinkage: step 2 keeps ~P(HDR)**k of the bank, and a smaller pool min-unions to a
LARGER failed-seed union, so the floor itself falls. k=0 is the pure all-cut reference.

`AllCutConfig` has no `seed_list` (only `AllHdrConfig` does), so the joined space is carried by the
subclass below. `bank_key` already folds a `seed_list` in when present, so the joined bank cannot
collide with the contiguous range that spans it.

    python conj_joined.py
    CJ_GROUPS=80 CJ_KS=0,1,2 CJ_LIGHT=6 python conj_joined.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjj")

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
from conj_test import choose_band, hdr_compliance           # noqa: E402
from allhdr_mf import pick                                  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402
from sd_task import score                                   # noqa: E402

# Three width-100 classes, deliberately spaced so the space is genuinely non-contiguous.
CLASSES = [tuple(int(y) for y in x.split("-"))
           for x in os.getenv("CJ_JOINED", "100-199,400-499,700-799").split(",")]
BAND_WIN = tuple(int(x) for x in os.getenv("CJ_BAND", "700-729").split("-"))
GROUPS = [int(x) for x in os.getenv("CJ_GROUPS", "42,50,80,100,125").split(",")]
KS = [int(x) for x in os.getenv("CJ_KS", "0,1").split(",")]
LIGHTS = [None if x.lower() in ("none", "") else int(x)
          for x in os.getenv("CJ_LIGHT", "none,6").split(",")]
MF = os.getenv("CJ_MF", "")            # blank -> linear scaling of the 900-seed default
N_SAMPLE = int(os.getenv("CJ_SAMPLE", "3"))
# `scan_cas9` stops at `pool_target * want` candidates. The band filter then keeps only
# ~P(HDR)**k of them, so a high k needs a deeper scan or it declines on an early exit rather
# than on a real shortage: at k=8 the default 8x leaves ~83 of 21,144 against 170 rows needed.
POOL_TARGET = int(os.getenv("CJ_POOL_TARGET", "8"))
BAND_POOL = int(os.getenv("CJ_BAND_POOL", "30"))
OUT = os.getenv("CJ_OUT", "conj_joined.json")


@_dc.dataclass(frozen=True)
class JoinedCut(AC.AllCutConfig):
    """all-cut over an explicit, possibly non-contiguous seed set."""

    seed_list: tuple = ()

    @property
    def seeds(self):
        if self.seed_list:
            return np.asarray(sorted(set(int(x) for x in self.seed_list)), dtype=np.int64)
        return np.arange(self.start_seed, self.end_seed + 1, dtype=np.int64)


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    span = len(joined)
    base = AC.config_for("K562") or AC.AllCutConfig()
    mf = int(MF) if MF else max(1, round(base.cas12a_max_fail * span / 900))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    rng = random.Random(20260911)
    print(f"K562  task {(task.get('task_id') or task['id'])[:8]}  joined cut window "
          f"{','.join(f'{a}-{b}' for a, b in CLASSES)} ({span} seeds)  "
          f"HDR band sub-window {BAND_WIN[0]}-{BAND_WIN[1]}  cut mf {mf}  rows {n_rows}\n")

    cfg0 = _dc.replace(JoinedCut(**{f.name: getattr(base, f.name)
                                    for f in _dc.fields(AC.AllCutConfig)}),
                       seed_list=tuple(joined), start_seed=min(joined), end_seed=max(joined),
                       cas12a_max_fail=mf)
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
    print(f"bank: {len(records)} Cas12a guides\n")

    band_candidates = sorted(rng.sample(range(BAND_WIN[0], BAND_WIN[1] + 1),
                                        min(BAND_POOL, BAND_WIN[1] - BAND_WIN[0] + 1)))
    t0 = time.monotonic()
    ok = hdr_compliance(records, contract, cell_types, ctx, band_candidates)
    print(f"HDR screen over {len(band_candidates)} candidate band seeds "
          f"({time.monotonic()-t0:.0f}s)\n")

    jset = set(joined)
    offw = [s for s in range(100, 1000) if s not in jset]
    print(f"{'grp':>4}{'k':>3}{'light':>7}{'pool':>8}{'union':>7}{'clean':>7}{'cas9':>7}"
          f"{'cells':>6}{'weighted':>10}{'fid':>8}{'w x fid':>9}{'band':>7}{'clean':>7}"
          f"{'dirty':>7}{'offwin':>7}{'E[cons]':>9}{'E[final]':>10}")
    out = []
    for g in GROUPS:
        for k in KS:
            band, alive = choose_band(ok, k, band_candidates, g)
            if k and len(band) < k:
                print(f"{g:>4}{k:>3}   band of {k} unreachable (pool below the group)")
                out.append({"group": g, "k": k, "reason": "band unreachable"}); continue
            pool = [records[i] for i in sorted(alive)]
            if len(pool) < g:
                print(f"{g:>4}{k:>3}   pool {len(pool)} below group {g}")
                out.append({"group": g, "k": k, "reason": "pool below group"}); continue
            sel = FG.FastGreedy(pool, window_lo=min(joined), window_hi=max(joined),
                                seeds=np.asarray(joined, dtype=np.int64))
            idx, _u = sel.best(g, restarts=12)
            group = [pool[i] for i in idx]
            bad = set()
            for rec in group:
                bad.update(int(x) for x in rec["fails"])
            clean = np.array(sorted(jset - bad), dtype=np.int64)
            if clean.size == 0:
                print(f"{g:>4}{k:>3}   the group's failures cover the whole joined window")
                out.append({"group": g, "k": k, "reason": "clean empty"}); continue
            # The Cas9 scan depends on the clean set and the band, NOT on light_cell_rows, which
            # only governs how `assemble` apportions the rows it returns. Hoisted out of the light
            # loop: re-scanning per light value doubled this sweep for identical candidates.
            cfg_scan = _dc.replace(cfg0, group_size=g, pool_target=POOL_TARGET)
            cas9 = AC.scan_cas9(clean, contract, cell_types, ctx, sites, cfg_scan, n_rows - g)
            if band and cas9:
                ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
                cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
            for light in LIGHTS:
                cfg = _dc.replace(cfg0, group_size=g, light_cell_rows=light)
                if len(cas9) < n_rows - g:
                    print(f"{g:>4}{k:>3}{str(light):>7}{len(pool):>8}{len(bad):>7}"
                          f"{clean.size:>7}{len(cas9):>7}   Cas9 half short of {n_rows-g}")
                    out.append({"group": g, "k": k, "light": light, "pool": len(pool),
                                "clean": int(clean.size), "cas9": len(cas9),
                                "reason": "cas9 short"}); continue
                rows = AC.assemble(group, cas9, contract, ctx, cfg, n_rows)
                cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
                cleanset = set(int(x) for x in clean)
                s_band = score(rows, contract, reference, cell_types, seed=band[0]) if band else None
                pc = [s for s in cleanset if s not in set(band)]
                cs = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                      for s in rng.sample(pc, min(N_SAMPLE, len(pc)))] or [0.0]
                ds = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                      for s in rng.sample(sorted(bad), min(N_SAMPLE, len(bad)))] or [0.10]
                os_ = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
                       for s in rng.sample(offw, min(N_SAMPLE, len(offw)))]
                bse = s_band or score(rows, contract, reference, cell_types, seed=int(clean[0]))
                e_cons = ((len(band) * (s_band["consistency"] if s_band else 0.0)
                           + (len(cleanset) - len(band)) * st.mean(cs)
                           + len(bad) * st.mean(ds)
                           + len(offw) * st.mean(os_)) / 900.0)
                wf = bse["weighted"] * bse["fidelity"]
                print(f"{g:>4}{k:>3}{str(light):>7}{len(pool):>8}{len(bad):>7}{len(cleanset):>7}"
                      f"{len(cas9):>7}{cells:>6}{bse['weighted']:>10.1f}{bse['fidelity']:>8.4f}"
                      f"{wf:>9.2f}{(s_band['consistency'] if s_band else float('nan')):>7.3f}"
                      f"{st.mean(cs):>7.3f}{st.mean(ds):>7.3f}{st.mean(os_):>7.3f}"
                      f"{e_cons:>9.4f}{wf*e_cons:>10.2f}")
                out.append({"group": g, "k": k, "light": light, "band": band,
                            "pool": len(pool), "union": len(bad), "clean": len(cleanset),
                            "cas9": len(cas9), "rows": len(rows), "cells": cells,
                            "weighted": bse["weighted"], "fidelity": bse["fidelity"],
                            "wxfid": wf, "band_cons": s_band["consistency"] if s_band else None,
                            "clean_cons": st.mean(cs), "dirty_cons": st.mean(ds),
                            "offwindow_cons": st.mean(os_), "e_cons": e_cons,
                            "e_final": wf * e_cons})
                json.dump(out, open(OUT, "w"), indent=1)
                MT.free_gpu_memory()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
