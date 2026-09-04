#!/usr/bin/env python3
"""win20_width.py — where does the all-HDR band actually peak as a function of window width?

Two measured points bracket an optimum that neither is at: band 13 (60k) / 14 (150k) at width 100,
and 15 at width 20. Coverage at a fixed hotkey count is just the absolute band, so the question is
which width maximises it. Fills in 30 / 40 / 50, several windows each, against the same K562
contract the width-20 sweep used.

Two things the earlier runs established and this one is built around:

* ``load_bank``'s 60,000 cap is binding at width 100 and not at width 20, so a width sweep at a
  fixed limit would confound the two. Every window is greedy-solved at 60,000 and 150,000.
* ``main_max_fail`` bounds what reaches disk; the load cap then takes the lowest-fail records. It
  is z-matched to width 20's ``mf6`` (``0.43w - 0.584*sqrt(w)``) so each width banks a comparable
  fraction of the pool. Constant-*rate* scaling tightens with width and would starve the wide banks
  below the load cap, which is the confound above wearing a different hat.

A wider band is not free — each seed added multiplies the conditional Cas9 requirement by
P(HDR) ~ 0.57 — so the best limit per window is carried through ``scan_cas9`` + ``assemble`` and
scored, rather than reported as a band count that might not build.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import dataclasses
import json
import sys
import time
from collections import Counter

import numpy as np

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

TASK = "task-k562.json"
LIMITS = [int(x) for x in
          (os.getenv("WIN20_LIMITS") or "60000,150000").split(",")]

# width -> (max_fail, windows). 20 and 100 re-use the banks win20.py already wrote.
PLAN = [
    (20, 6, [(300, 319)]),
    (30, 10, [(300, 329), (330, 359), (360, 389)]),
    (40, 14, [(300, 339), (340, 379)]),
    (50, 17, [(300, 349), (350, 399)]),
    (100, 45, [(300, 399)]),
]
OUT = os.getenv("WIN20_OUT", "win20_width.json")
# A looser screen re-run: the z-matched max_fail above under-delivers as width grows
# (36k on disk at w50 against 150k at w20), and pool size moves the band, so the wide
# points needed a second pass at a comparable bank. Set WIN20_PLAN to override.
if os.getenv("WIN20_PLAN"):
    PLAN = [(w, mf, [tuple(x) for x in ws])
            for w, mf, ws in json.loads(os.environ["WIN20_PLAN"])]


def bank_for(contract, cell_types, ctx, sites, cfg, reference):
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return None, 0
        save_bank(path, bank)
    n = int(np.load(path, allow_pickle=False)["fails"].shape[0])
    return path, n


def band_at(path, cfg, ctx, contract, limit):
    records = load_bank(path, limit=limit)
    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg))
    index, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
    bad = set()
    for i in index:
        bad.update(int(x) for x in records[i]["fails"])
    band = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
    return records, index, band


def main():
    task = json.load(open(TASK))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    seeds = [int(x) for x in str(contract["seed"]).split(",")]
    probe = seeds[0]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AH.config_for(contract["cell_type"])
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    want = n_rows - base.group_size

    print(f"task {task['id'][:8]}  {contract['cell_type']}  seeds {seeds}  "
          f"(weighted/fid are seed-independent; cons is the off-band probe {probe})\n")
    print(f"{'w':>4} {'mf':>3} {'window':>9} {'disk':>7} {'b@60k':>6} {'b@150k':>7} {'best':>5} "
          f"{'frac':>6} {'cas9':>6} {'rows':>5} {'cells':>5} {'weighted':>9} {'fid':>6} "
          f"{'build':>7}")
    out = []
    for width, max_fail, windows in PLAN:
        for window in windows:
            started = time.monotonic()
            cfg = dataclasses.replace(base, hdr_range=window, main_max_fail=max_fail)
            path, disk = bank_for(contract, cell_types, ctx, sites, cfg, reference)
            if path is None:
                print(f"{width:>4} {max_fail:>3} {window[0]}-{window[1]:<4}  bank scan failed",
                      flush=True)
                continue
            bands = {}
            for limit in LIMITS:
                records, index, band = band_at(path, cfg, ctx, contract, limit)
                bands[limit] = (records, index, band)
            best_limit = max(LIMITS, key=lambda L: (len(bands[L][2]), -L))
            records, index, band = bands[best_limit]
            rec = {"width": width, "max_fail": max_fail, "window": list(window), "disk": disk,
                   "band_60k": len(bands[60_000][2]), "band_150k": len(bands[150_000][2]),
                   "best_limit": best_limit, "band": band, "clean": len(band),
                   "clean_fraction": round(len(band) / width, 4),
                   "k": len([s for s in seeds if s in set(band)])}

            group = [records[i] for i in index]
            cas9 = AH.scan_cas9(np.array(band, dtype=np.int64), contract, cell_types, ctx, sites,
                                cfg, want, None)
            MT.free_gpu_memory()
            rec["cas9_pool"] = len(cas9)
            cells4 = len({(r["mutation"], r["strand"]) for r in cas9})
            if len(cas9) < want or cells4 < 4:
                rec.update(declined=f"Cas9 pool {len(cas9)} over {cells4} cells short of {want}",
                           build_s=round(time.monotonic() - started, 1))
                out.append(rec)
                print(f"{width:>4} {max_fail:>3} {window[0]}-{window[1]:<4} {disk:>7} "
                      f"{rec['band_60k']:>6} {rec['band_150k']:>7} {len(band):>5} "
                      f"{len(band) / width * 100:>5.1f}% {len(cas9):>6}  DECLINED", flush=True)
                json.dump(out, open(OUT, "w"), indent=1)
                continue

            rows = AH.assemble(group, cas9, contract, ctx, cfg, n_rows)
            build_s = time.monotonic() - started
            json.dump(contract, open(settings.CONTRACT_PATH, "w"))
            json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
            json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
            run_stage12(cell_types)
            run_stage3(seed=probe)
            run_stage4(seed=probe)
            res = run_stage5()
            rec.update(rows=len(rows),
                       cells=len(Counter((r["mutation"], r["cas_system"], r["strand"])
                                         for r in rows)),
                       total_weighted_score=res["total_weighted_score"],
                       distribution_fidelity_factor=res["distribution_fidelity_factor"],
                       probe_consistency_factor=res["consistency_factor"],
                       build_s=round(build_s, 1))
            out.append(rec)
            print(f"{width:>4} {max_fail:>3} {window[0]}-{window[1]:<4} {disk:>7} "
                  f"{rec['band_60k']:>6} {rec['band_150k']:>7} {len(band):>5} "
                  f"{len(band) / width * 100:>5.1f}% {len(cas9):>6} {rec['rows']:>5} "
                  f"{rec['cells']:>3}/8 {rec['total_weighted_score']:>9.1f} "
                  f"{rec['distribution_fidelity_factor']:>6.3f} {build_s:>6.1f}s", flush=True)
            json.dump(out, open(OUT, "w"), indent=1)

    print("\nband by width (best limit), mean over windows:")
    for width, _mf, _w in PLAN:
        got = [r for r in out if r["width"] == width and "rows" in r]
        if got:
            mean = sum(r["clean"] for r in got) / len(got)
            print(f"  width {width:>3}: band {mean:>5.1f}  "
                  f"({'/'.join(str(r['clean']) for r in got)})  "
                  f"{len(got)}/{len([r for r in out if r['width'] == width])} built")
    json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
