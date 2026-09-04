#!/usr/bin/env python3
"""win20_floor.py — measure the off-band consistency floor, the input that flips E[pay]'s sign.

The E[pay] pricing assumed 0.10 for a miss, from CLAUDE.md's k-of-3 table. That assumption is
load-bearing: sweeping it over the range these builds actually measure (0.086-0.106) moves the
width-20 verdict from -4.0% to +4.1% and back to -1.5% at 0.12. The reason is structural rather
than numerical — our k=1 score of 86-107 lands in the densest part of the field, against a rank-10
cutoff median of 74-84, so a small shift crosses several rank steps.

So the floor is measured here rather than assumed: every paired contract is rebuilt at both widths
(banks are cached, so this is greedy + Cas9 + assemble) and scored on SEEDS_PER off-band seeds, and
the per-(cell, width) mean floor is written for win20_epay.py to consume. Off-band seeds are drawn
from inside the config's own window but outside its band, which is the case a real miss hits.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import dataclasses
import json
import statistics as st
import sys
import time

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

CELLS = ["CD34+_HSPC", "HUDEP-2", "K562"]
WIDTHS = [100, 20]
MF20 = 6
SEEDS_PER = 3
OUT = "win20_floor.json"
SOURCES = ["win20_cells.json", "win20_narrow.json", "win20_fill.json"]


def paired_contracts():
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
    out = {}
    for cell in CELLS:
        sets = [{r["task_id"] for r in uniq if r["cell_type"] == cell and r["width"] == w}
                for w in WIDTHS]
        out[cell] = sorted(set.intersection(*sets))
    return out


def main():
    want = paired_contracts()
    ids = {t for v in want.values() for t in v}
    tl = json.load(urllib.request.urlopen("https://niome-api.genomes.io/api/v3/tasks", timeout=120))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    tasks = {t["id"]: t for t in tl if t["id"] in ids}
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["task_id"], r["width"]) for r in out}
    total = sum(len(v) for v in want.values()) * len(WIDTHS)
    print(f"{total} (contract, width) cells, {len(done)} done\n")
    print(f"{'cell':<11} {'task':<9} {'w':>4} {'band':>5} {'floor seeds':>28} {'mean':>6} "
          f"{'build':>7}")
    started, step = time.monotonic(), 0
    for cell in CELLS:
        base = AH.config_for(cell)
        for tid in want[cell]:
            task = tasks.get(tid)
            if task is None:
                continue
            contract = task["content"]["contract"]
            reference = task["content"]["hbb_reference"]
            ctx = G.build_context(contract, reference, cell_types)
            sites = G.enumerate_sites(ctx, 3000, (20, 23))
            n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
            for width in WIDTHS:
                step += 1
                if (tid, width) in done:
                    continue
                cfg = base if width == 100 else dataclasses.replace(
                    base, hdr_range=(base.start_seed, base.start_seed + width - 1),
                    main_max_fail=MF20)
                t0 = time.monotonic()
                path = os.path.join(AH.HDR_BANK_DIR,
                                    f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
                if not os.path.exists(path):
                    print(f"{cell:<11} {tid[:8]:<9} {width:>4}  no cached bank; skipped",
                          flush=True)
                    continue
                records = load_bank(path)
                sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                                    per_cell_min=cfg.per_cell_min,
                                    caps=AH._group_caps(contract, ctx, cfg))
                index, _ = sel.best(cfg.group_size, restarts=cfg.restarts)
                bad = set()
                for i in index:
                    bad.update(int(x) for x in records[i]["fails"])
                band = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
                cas9 = AH.scan_cas9(np.array(band, dtype=np.int64), contract, cell_types, ctx,
                                    sites, cfg, n_rows - cfg.group_size, None)
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
                probes = off[:: max(1, len(off) // SEEDS_PER)][:SEEDS_PER]
                floors = []
                for seed in probes:
                    run_stage3(seed=seed)
                    run_stage4(seed=seed)
                    floors.append(run_stage5()["consistency_factor"])
                rec = {"task_id": tid, "cell_type": cell, "width": width, "band": len(band),
                       "probe_seeds": probes, "floors": floors, "floor": st.mean(floors),
                       "build_s": round(time.monotonic() - t0, 1)}
                out.append(rec)
                json.dump(out, open(OUT, "w"), indent=1)
                eta = (time.monotonic() - started) / step * (total - step) / 60
                print(f"{cell:<11} {tid[:8]:<9} {width:>4} {len(band):>5} "
                      f"{' '.join(f'{x:.4f}' for x in floors):>28} {rec['floor']:>6.4f} "
                      f"{rec['build_s']:>6.1f}s   eta {eta:.0f}m", flush=True)

    print("\nmeasured floor by cell type and width:")
    for cell in CELLS:
        for w in WIDTHS:
            g = [r for r in out if r["cell_type"] == cell and r["width"] == w]
            if g:
                allf = [x for r in g for x in r["floors"]]
                print(f"  {cell:<11} w{w:<4} floor {st.mean(allf):.4f} "
                      f"(sd {st.stdev(allf):.4f}, n={len(allf)}, "
                      f"range {min(allf):.4f}-{max(allf):.4f})")
    print("\ndone")


if __name__ == "__main__":
    main()
