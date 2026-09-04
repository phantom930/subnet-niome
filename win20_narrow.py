#!/usr/bin/env python3
"""win20_narrow.py — below width 20: does the band hold, or does the window start capping it?

Width 20 beat the shipped width 100 on all four cell types (+1.6 to +3.2 band seeds, every
interval clear of zero). This asks what happens further down: widths 16 and 18 for the erythroid
types, 12 and 16 for HEK293, against the same 20 contracts so every number pairs with the width-20
and width-100 results in win20_cells.json.

The two regimes are already at different places against their window, which is why the widths
differ by cell type:

* erythroid — band 14.6-14.9 of 20, i.e. 73-75% of the window. Widths 18 and 16 press straight
  against that ceiling.
* HEK293 — band 10.2 of 20 (51%), but its bank already saturates ``bank_keep`` at 300,000, so its
  pool cannot grow further and its band bound should be roughly fixed. Width 12 tests whether the
  window itself becomes the binding constraint.

``main_max_fail`` goes *looser* as the window narrows, which inverts the first width sweep's
mistake. ``load_bank`` keeps the lowest-fail 60,000 and ``bank_keep`` the lowest-fail 300,000, so a
looser screen only admits high-fail guides that truncation discards — it can enlarge the effective
pool but never degrade it. The erythroid types therefore hold ``mf6`` (a bigger bank at width 16
than at width 20), and HEK293 takes ``mf8``/``mf6``, z-matched to the ``mf10`` that already
saturated bank_keep. ``ensure_bank``'s +2 retry still guards the 60,000 floor.
"""
import json
import os
import sys
import time

sys.path.insert(0, "/root/workspace/subnet-niome")

import dataclasses

from win20_cells import (AH, G, MT, one_build, pick_contracts)  # noqa: E402  (sets NIOME_INSTANCE)

CELL_WIDTHS = {
    "CD34+_HSPC": [(18, 6), (16, 6)],
    "HUDEP-2": [(18, 6), (16, 6)],
    "K562": [(18, 6), (16, 6)],
    "HEK293": [(16, 8), (12, 6)],
}
OFFSETS = [0, 40]
OUT = "win20_narrow.json"
BASE = "win20_cells.json"


def main():
    tasks = pick_contracts()
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    done, out = set(), []
    if os.path.exists(OUT):
        out = json.load(open(OUT))
        done = {(r["task_id"], r["width"], r["window"][0]) for r in out}
    total = sum(len(CELL_WIDTHS.get(t["content"]["contract"]["cell_type"], [])) * len(OFFSETS)
                for t in tasks)
    print(f"{len(tasks)} contracts, {total} builds, {len(done)} already done\n")
    print(f"{'cell':<11} {'task':<9} {'w':>3} {'mf':>3} {'window':>9} {'disk':>7} {'b@60k':>6} "
          f"{'b@150k':>7} {'band':>5} {'frac':>6} {'cas9':>6} {'rows':>5} {'weighted':>9} "
          f"{'fid':>6} {'build':>7}")

    started, step = time.monotonic(), 0
    for task in tasks:
        contract = task["content"]["contract"]
        reference = task["content"]["hbb_reference"]
        cell = contract["cell_type"]
        plan = CELL_WIDTHS.get(cell)
        base = AH.config_for(cell)
        if not plan or base is None:
            continue
        seeds = [int(x) for x in str(contract.get("seed") or "0").split(",") if x.strip()]
        probe = next((s for s in seeds if s >= 100), 500)
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        lo = base.start_seed
        for width, mf in plan:
            for off in OFFSETS:
                step += 1
                key = (task["id"], width, lo + off)
                if key in done:
                    continue
                cfg = dataclasses.replace(base, hdr_range=(lo + off, lo + off + width - 1),
                                          main_max_fail=mf)
                rec = one_build(contract, reference, cell_types, ctx, sites, cfg, True, probe,
                                seeds)
                rec.update(task_id=task["id"], cell_type=cell)
                out.append(rec)
                json.dump(out, open(OUT, "w"), indent=1)
                eta = (time.monotonic() - started) / step * (total - step) / 60
                tail = (f"{rec.get('cas9_pool', 0):>6} {rec.get('rows', 0):>5} "
                        f"{rec.get('total_weighted_score', 0):>9.1f} "
                        f"{rec.get('distribution_fidelity_factor', 0):>6.3f} "
                        f"{rec.get('build_s', 0):>6.1f}s")
                if "declined" in rec:
                    tail = f"DECLINED: {rec['declined']}"
                print(f"{cell:<11} {task['id'][:8]:<9} {rec['width']:>3} {rec['max_fail']:>3} "
                      f"{rec['window'][0]}-{rec['window'][1]:<4} {rec['disk']:>7} "
                      f"{rec.get('band_60k', 0):>6} {rec.get('band_150k', 0):>7} "
                      f"{rec.get('clean', 0):>5} {rec.get('clean_fraction', 0) * 100:>5.1f}% "
                      f"{tail}   eta {eta:.0f}m", flush=True)
                MT.free_gpu_memory()

    prev = json.load(open(BASE)) if os.path.exists(BASE) else []
    print("\nband by cell type and width (mean over contracts x windows, b@60k):")
    print(f"{'cell':<11} " + "".join(f"{f'w{w}':>8}" for w in (100, 20, 18, 16, 12)))
    for cell in sorted({r["cell_type"] for r in out}):
        cells = []
        for w in (100, 20, 18, 16, 12):
            g = [r for r in (prev + out)
                 if r["cell_type"] == cell and r["width"] == w and "declined" not in r]
            cells.append(f"{sum(r['band_60k'] for r in g) / len(g):>8.2f}" if g else f"{'-':>8}")
        print(f"{cell:<11} " + "".join(cells))
    dec = [r for r in out if "declined" in r]
    print(f"\ndeclines: {len(dec)}" + ("" if not dec else
          "\n  " + "\n  ".join(f"{r['cell_type']} w{r['width']} {r['window']}: {r['declined']}"
                               for r in dec)))


if __name__ == "__main__":
    main()
