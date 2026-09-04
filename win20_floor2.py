#!/usr/bin/env python3
"""win20_floor2.py — off-band consistency floors for the widths the first pass did not cover.

win20_floor.py measured widths 100 and 20 for the three erythroid types. Pricing width 16 for those
and width 12 for HEK293 needs floors there too, and HEK293 has none at any width.

Driven from the stored sweep records rather than reconstructing configs: each build is rebuilt with
the exact window and ``main_max_fail`` that produced it, so the cached bank is guaranteed to match
(``bank_key`` folds in the screen, and ``ensure_bank`` may have bumped it during the sweep).

Note on what the floor is: a miss means the drawn seed is anywhere in 100-999 outside the band, and
seeds are independent hashes, so a non-band seed sampled from inside the window is representative
of one outside it. That matters at width 16, where a band of ~14.5 leaves only one or two off-band
seeds in the window itself.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import dataclasses
import json
import statistics as st
import sys
import time
from collections import defaultdict

sys.path.insert(0, "/root/workspace/subnet-niome")

import numpy as np
import urllib.request

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.all_cut import bank_key, load_bank
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

PLAN = {"CD34+_HSPC": [16], "HUDEP-2": [16], "K562": [16], "HEK293": [100, 20, 12]}
SEEDS_PER = 3
SOURCES = ["win20_cells.json", "win20_narrow.json", "win20_fill.json"]
OUT = "win20_floor2.json"


def main():
    recs, seen, uniq = [], set(), []
    for p in SOURCES:
        if os.path.exists(p):
            recs += json.load(open(p))
    for r in recs:
        k = (r["task_id"], r["width"], r["window"][0])
        if k in seen or "declined" in r:
            continue
        seen.add(k)
        uniq.append(r)

    # one record per (cell, width, contract) -- the lowest window offset
    want = {}
    for r in uniq:
        if r["width"] not in PLAN.get(r["cell_type"], []):
            continue
        key = (r["cell_type"], r["width"], r["task_id"])
        if key not in want or r["window"][0] < want[key]["window"][0]:
            want[key] = r

    ids = {k[2] for k in want}
    tl = json.load(urllib.request.urlopen("https://niome-api.genomes.io/api/v3/tasks", timeout=120))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    tasks = {t["id"]: t for t in tl if t["id"] in ids}
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["cell_type"], r["width"], r["task_id"]) for r in out}

    print(f"{len(want)} (cell, width, contract) cells, {len(done)} done\n")
    print(f"{'cell':<11} {'task':<9} {'w':>4} {'mf':>3} {'band':>5} {'floors':>26} {'mean':>7} "
          f"{'build':>7}")
    started, step = time.monotonic(), 0
    for (cell, width, tid), rec in sorted(want.items()):
        step += 1
        if (cell, width, tid) in done:
            continue
        task = tasks.get(tid)
        if task is None:
            print(f"{cell:<11} {tid[:8]:<9} {width:>4}  task no longer listed", flush=True)
            continue
        contract = task["content"]["contract"]
        reference = task["content"]["hbb_reference"]
        base = AH.config_for(cell)
        cfg = dataclasses.replace(base, hdr_range=tuple(rec["window"]),
                                  main_max_fail=rec["max_fail"])
        path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            print(f"{cell:<11} {tid[:8]:<9} {width:>4}  no cached bank", flush=True)
            continue
        t0 = time.monotonic()
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
        records = load_bank(path, limit=rec.get("best_limit", 60_000))
        sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                            per_cell_min=cfg.per_cell_min,
                            caps=AH._group_caps(contract, ctx, cfg))
        index, _ = sel.best(cfg.group_size, restarts=cfg.restarts)
        bad = set()
        for i in index:
            bad.update(int(x) for x in records[i]["fails"])
        band = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
        cas9 = AH.scan_cas9(np.array(band, dtype=np.int64), contract, cell_types, ctx, sites, cfg,
                            n_rows - cfg.group_size, None)
        MT.free_gpu_memory()
        if len(cas9) < n_rows - cfg.group_size:
            print(f"{cell:<11} {tid[:8]:<9} {width:>4}  Cas9 short; skipped", flush=True)
            continue
        rows = AH.assemble([records[i] for i in index], cas9, contract, ctx, cfg, n_rows)
        json.dump(contract, open(settings.CONTRACT_PATH, "w"))
        json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        off = [s for s in range(cfg.start_seed, cfg.end_seed + 1) if s not in set(band)]
        probes = off[:: max(1, len(off) // SEEDS_PER)][:SEEDS_PER] or [cfg.start_seed]
        floors = []
        for seed in probes:
            run_stage3(seed=seed)
            run_stage4(seed=seed)
            floors.append(run_stage5()["consistency_factor"])
        r = {"task_id": tid, "cell_type": cell, "width": width, "max_fail": rec["max_fail"],
             "band": len(band), "probe_seeds": probes, "floors": floors, "floor": st.mean(floors),
             "build_s": round(time.monotonic() - t0, 1)}
        out.append(r)
        json.dump(out, open(OUT, "w"), indent=1)
        eta = (time.monotonic() - started) / step * (len(want) - step) / 60
        print(f"{cell:<11} {tid[:8]:<9} {width:>4} {rec['max_fail']:>3} {len(band):>5} "
              f"{' '.join(f'{x:.4f}' for x in floors):>26} {r['floor']:>7.4f} "
              f"{r['build_s']:>6.1f}s   eta {eta:.0f}m", flush=True)

    agg = defaultdict(list)
    for r in out:
        agg[(r["cell_type"], r["width"])] += r["floors"]
    print("\nmeasured floor:")
    for k in sorted(agg):
        v = agg[k]
        print(f"  {k[0]:<11} w{k[1]:<4} mean {st.mean(v):.4f} sd {st.stdev(v):.4f} n={len(v)}")
    print("\ndone")


if __name__ == "__main__":
    main()
