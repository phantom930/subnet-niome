#!/usr/bin/env python3
"""cell_grid.py — the full all-HDR mixed config grid for one cell type.

Sweeps group size x window width x light_cell_rows, with `main_max_fail` optimised per width, and
records the two quantities that decide payout: the clean BAND and `weighted x fidelity`.

Why those two are enough. all-HDR mixed has no elevated floor — every seed outside the band sits at
the ~0.10 arithmetic floor (CLAUDE.md's "the ~0.10 floor is arithmetic, not a design failure"), so a
round's consistency is `(n_band + (3 - n_band) * floor) / 3` and expected payout is a function of
`band` and `w x fid` alone. Both are design-only and seed-independent, which is what makes a
400-arm grid affordable:

  * `total_weighted_score` is exactly the sum of stage 12's per-row `stage2.weighted_score`
    (verified to 1e-6), so stage 4's RandomForest is never run.
  * `distribution_fidelity_score` uses only stage 12's valid experiments — `cas_shift` is a
    diagnostic outside the score — so stage 3 is never run either (verified identical).

`main_max_fail` is found ADAPTIVELY rather than swept: the optimum is always the first value whose
raw bank reaches the 300,000 `bank_keep` cap (measured across five widths on K562, and the same
z-score in every case), so the search starts at `_scaled_max_fail` and loosens 15% at a time until
the bank caps. That replaces a 4-value mf sweep per width with 1-2 bank builds.

The Cas9 scan depends on (width, group) but not on light, so it is hoisted above the light loop.

    CG_CELL=K562 python cell_grid.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "cgrid")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402
import urllib.request                # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank    # noqa: E402
from niome_subnet.genomics.validation.stage12 import run_stage12            # noqa: E402
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity  # noqa: E402
from niome_subnet.utils import settings                     # noqa: E402

CELL = os.getenv("CG_CELL", "K562")
GROUPS = [int(x) for x in os.getenv("CG_GROUPS", "42,50,80,100,125").split(",")]
WIDTHS = [int(x) for x in os.getenv("CG_WIDTHS", "30,75,100,150,225").split(",")]
# CG_PAIRS restricts the sweep to explicit "width:group" combinations instead of the full
# cross product — the validation run only needs the few configs the grid short-listed.
PAIRS = [tuple(int(y) for y in x.split(":"))
         for x in os.getenv("CG_PAIRS", "").split(",") if x]
LIGHTS = [None if x.lower() in ("none", "") else int(x)
          for x in os.getenv("CG_LIGHTS", "none,6,12,25").split(",")]
START = int(os.getenv("CG_START", "100"))
TASK = os.getenv("CG_TASK", "")
MAXLOAD = float(os.getenv("CG_MAXLOAD", "4"))
OUT = os.getenv("CG_OUT", "")


def pick_task():
    tk = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks?limit=500"))
    tk = tk if isinstance(tk, list) else (tk.get("items") or tk.get("data"))
    best = None
    for t in sorted(tk, key=lambda t: t.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL:
            continue
        tid = t.get("task_id") or t["id"]
        if TASK:
            if tid.startswith(TASK):
                return t
            continue
        s = str(c.get("seed", "") or "")
        if len([x for x in s.split(",") if x.strip().isdigit()]) == 3 and best is None:
            best = t
    if best is None:
        raise SystemExit(f"no stamped 3-seed task for {CELL}")
    return best


def wait_idle():
    while True:
        with open("/proc/loadavg") as fh:
            if float(fh.read().split()[0]) < MAXLOAD:
                return
        time.sleep(20)


def score_design(rows, contract, reference, cell_types):
    """weighted and fidelity from stage 12 alone — both are design-only, so no stage 3/4."""
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    weighted = sum(e["stage2"]["weighted_score"] for e in valid)
    summary = compute_distribution_fidelity(valid, [], contract, k=12)
    fid = max(0.0, min(1.0, summary.get("distribution_fidelity_score", 0.0)))
    return weighted, fid, len(valid), summary


def bank_for(contract, reference, cell_types, ctx, sites, cfg, span):
    """Build/load the bank, loosening `main_max_fail` until the raw scan reaches the keep cap."""
    tried = []
    mf = cfg.main_max_fail
    for _ in range(3):
        c = _dc.replace(cfg, main_max_fail=mf)
        path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, c)}.npz")
        raw = None
        if not os.path.exists(path):
            bank = AH.build_bank(contract, reference, cell_types, ctx, sites, c)
            MT.free_gpu_memory()
            raw = len(bank)
            if bank:
                save_bank(path, bank)
        recs = load_bank(path, limit=300_000) if os.path.exists(path) else []
        tried.append((mf, raw, len(recs)))
        # capped (or already huge) -> on the plateau; otherwise loosen and retry
        if len(recs) >= 300_000 or (raw is not None and raw >= 300_000):
            return c, recs, tried
        nxt = max(mf + 1, int(round(mf * 1.15)))
        if nxt >= span:
            return c, recs, tried
        mf = nxt
    return c, recs, tried


def main():
    task = pick_task()
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    base = AH.config_for(CELL)
    w = contract.get("mutation_weights", {})
    spread = (max(w.values()) / min(w.values())) if len(w) == 2 else float("nan")
    out = OUT or f"cell_grid_{CELL.replace('+','').replace('-','_')}.json"
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  rows {n_rows}  "
          f"weight spread {spread:.2f}  groups {GROUPS}  widths {WIDTHS}  lights {LIGHTS}")
    print(f"window = contiguous, starting at {START} (band position is free)\n")
    print(f"{'width':>6}{'mf':>5}{'bank':>8}{'grp':>5}{'band':>6}{'cas9':>7}"
          + "".join(f"{('L'+str(l)):>9}" for l in LIGHTS) + f"{'best':>6}{'s':>6}")
    res = []
    widths = sorted({w for w, _g in PAIRS}) if PAIRS else WIDTHS
    for width in widths:
        wait_idle()
        lo, hi = START, START + width - 1
        # `_scaled_max_fail` treats its second argument as the value calibrated at span 100.
        # `base.main_max_fail` is NOT that for every cell: K562's CELL_CONFIG now carries 69, a
        # span-150 value, so scaling it would start the search 50% too loose. Use the span-100
        # reference explicitly — the shared 45, or HEK293's own 48.
        ref100 = 48 if CELL == "HEK293" else 45
        cfg0 = _dc.replace(base, hdr_range=(lo, hi), seed_list=None,
                           main_max_fail=(ref100 if width == 100
                                          else AH._scaled_max_fail(CELL, ref100, width)),
                           variants=(min(base.variants, AH.WIDE_WINDOW_VARIANTS)
                                     if width > 100 else base.variants))
        t0 = time.monotonic()
        cfg0, records, tried = bank_for(contract, reference, cell_types, ctx, sites, cfg0, width)
        if len(records) < min(GROUPS):
            print(f"{width:>6}{cfg0.main_max_fail:>5}{len(records):>8}   bank short of "
                  f"the smallest group; tried {tried}")
            continue
        for g in ([g for w, g in PAIRS if w == width] if PAIRS else GROUPS):
            wait_idle()
            t1 = time.monotonic()
            if len(records) < g:
                print(f"{width:>6}{cfg0.main_max_fail:>5}{len(records):>8}{g:>5}"
                      f"   bank < group"); continue
            cfg = _dc.replace(cfg0, group_size=g)
            sel = FG.FastGreedy(records, window_lo=lo, window_hi=hi,
                                per_cell_min=cfg.per_cell_min)
            idx, _u = sel.best(g, restarts=cfg.restarts)
            group = [records[i] for i in idx]
            bad = set()
            for rec in group:
                bad.update(int(x) for x in rec["fails"])
            clean = np.array(sorted(set(range(lo, hi + 1)) - bad), dtype=np.int64)
            if clean.size == 0:
                print(f"{width:>6}{cfg0.main_max_fail:>5}{len(records):>8}{g:>5}"
                      f"   band empty"); continue
            cas9 = AH.scan_cas9(clean, contract, cell_types, ctx, sites, cfg, n_rows - g)
            if len(cas9) < n_rows - g:
                print(f"{width:>6}{cfg0.main_max_fail:>5}{len(records):>8}{g:>5}"
                      f"{clean.size:>6}{len(cas9):>7}   Cas9 short of {n_rows-g}")
                res.append({"width": width, "mf": cfg0.main_max_fail, "group": g,
                            "band": int(clean.size), "cas9": len(cas9), "reason": "cas9 short"})
                json.dump(res, open(out, "w"), indent=1)
                MT.free_gpu_memory(); continue
            line = (f"{width:>6}{cfg0.main_max_fail:>5}{len(records):>8}{g:>5}"
                    f"{clean.size:>6}{len(cas9):>7}")
            arms, best, bestv = {}, None, -1
            for light in LIGHTS:
                rows = AC.assemble(group, cas9, contract, ctx,
                                   _dc.replace(cfg, light_cell_rows=light), n_rows)
                wt, fid, nval, summ = score_design(rows, contract, reference, cell_types)
                wf = wt * fid
                cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
                arms[str(light)] = {"weighted": wt, "fidelity": fid, "wxfid": wf,
                                    "rows": len(rows), "valid": nval, "cells": cells}
                line += f"{wf:>9.1f}"
                if wf > bestv:
                    best, bestv = str(light), wf
            print(line + f"{best:>6}{time.monotonic()-t1:>6.0f}")
            res.append({"cell": CELL, "task": (task.get("task_id") or task["id"])[:8],
                        "spread": spread, "width": width, "window": [lo, hi],
                        "mf": cfg0.main_max_fail, "mf_tried": tried, "bank": len(records),
                        "group": g, "band": int(clean.size), "band_seeds": clean.tolist(),
                        "cas9": len(cas9), "arms": arms, "best_light": best,
                        "build_s": round(time.monotonic() - t1, 1)})
            json.dump(res, open(out, "w"), indent=1)
            MT.free_gpu_memory()
        print(f"   [width {width} done in {time.monotonic()-t0:.0f}s]")
    json.dump(res, open(out, "w"), indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
