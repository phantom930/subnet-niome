#!/usr/bin/env python3
"""win20_fill.py — complete the K562 matrix so its width comparison is actually paired.

``pick_contracts`` takes evenly-spaced tasks from the backend list, and that list grew between the
width-20 run and the sub-20 run. Every other cell type happened to re-select the same five
contracts; K562's stride shifted and it kept only one. Its widths were therefore measured on
disjoint contract sets, which is why K562's weighted appeared to jump 274.9 -> 317.2 between widths
20 and 18 while its band did not move at all — that was contract mix, not width.

Fills every missing (contract, width) cell across the union of both K562 contract sets, giving a
9-contract x 4-width matrix where each width is measured on identical contracts. Also writes the
chosen ids to win20_contracts.json so later runs pin the sample instead of re-drawing it.
"""
import dataclasses
import json
import os
import sys
import time

sys.path.insert(0, "/root/workspace/subnet-niome")

import urllib.request

from win20_cells import AH, G, MT, one_build  # noqa: E402  (sets NIOME_INSTANCE)

CELL = "K562"
WIDTHS = {100: [0], 20: [0, 40], 18: [0, 40], 16: [0, 40]}
MF = {100: None, 20: 6, 18: 6, 16: 6}          # None = the cell type's shipped main_max_fail
SOURCES = ["win20_cells.json", "win20_narrow.json"]
OUT = "win20_fill.json"


def main():
    have, seen = [], set()
    for path in SOURCES:
        if os.path.exists(path):
            have += json.load(open(path))
    if os.path.exists(OUT):
        have += json.load(open(OUT))
    done = {(r["task_id"], r["width"], r["window"][0]) for r in have}
    wanted_tasks = sorted({r["task_id"] for r in have if r["cell_type"] == CELL})

    tl = json.load(urllib.request.urlopen("https://niome-api.genomes.io/api/v3/tasks", timeout=60))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    tasks = {t["id"]: t for t in tl if t["id"] in wanted_tasks}
    json.dump(wanted_tasks, open("win20_contracts.json", "w"), indent=1)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AH.config_for(CELL)
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    todo = [(tid, w, off) for tid in wanted_tasks for w, offs in WIDTHS.items() for off in offs
            if (tid, w, base.start_seed + off) not in done]
    print(f"{CELL}: {len(wanted_tasks)} contracts, {len(todo)} missing builds\n")
    print(f"{'task':<9} {'w':>3} {'mf':>3} {'window':>9} {'disk':>7} {'b@60k':>6} {'b@150k':>7} "
          f"{'band':>5} {'cas9':>6} {'weighted':>9} {'fid':>6} {'build':>7}")
    started = time.monotonic()
    for i, (tid, width, off) in enumerate(todo, 1):
        task = tasks.get(tid)
        if task is None:
            print(f"{tid[:8]:<9} task no longer listed; skipped", flush=True)
            continue
        contract = task["content"]["contract"]
        reference = task["content"]["hbb_reference"]
        seeds = [int(x) for x in str(contract.get("seed") or "0").split(",") if x.strip()]
        probe = next((s for s in seeds if s >= 100), 500)
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        lo = base.start_seed
        cfg = base if width == 100 else dataclasses.replace(
            base, hdr_range=(lo + off, lo + off + width - 1), main_max_fail=MF[width])
        rec = one_build(contract, reference, cell_types, ctx, sites, cfg, width != 100, probe,
                        seeds)
        rec.update(task_id=tid, cell_type=CELL)
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
        MT.free_gpu_memory()
        eta = (time.monotonic() - started) / i * (len(todo) - i) / 60
        tail = (f"{rec.get('cas9_pool', 0):>6} {rec.get('total_weighted_score', 0):>9.1f} "
                f"{rec.get('distribution_fidelity_factor', 0):>6.3f} "
                f"{rec.get('build_s', 0):>6.1f}s")
        if "declined" in rec:
            tail = f"DECLINED: {rec['declined']}"
        print(f"{tid[:8]:<9} {rec['width']:>3} {rec['max_fail']:>3} "
              f"{rec['window'][0]}-{rec['window'][1]:<4} {rec['disk']:>7} "
              f"{rec.get('band_60k', 0):>6} {rec.get('band_150k', 0):>7} "
              f"{rec.get('clean', 0):>5} {tail}   eta {eta:.0f}m", flush=True)
    print("\ndone")


if __name__ == "__main__":
    main()
