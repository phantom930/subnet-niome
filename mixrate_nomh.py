#!/usr/bin/env python3
"""mixrate_nomh.py — sweep the HDR:BLUNT row ratio in a not_mhnhej submission (UNSTRUCTURED row
order, no KFold-fold alignment), to find the best non-trivial mix rate for consistency_factor.

`kfold_align_nomh.py` found that deliberately SEGREGATING BLUNT rows into one KFold fold is actively
harmful relative to a naive mix at a similar ratio (200/50 aligned: cons 0.0018, against the
unstructured `nomh_standalone.py` build's 0.2743 at roughly a 159/91 natural ratio) -- structuring
the rows by outcome apparently hands the RandomForest a spurious, wrongly-generalising pattern to
overfit to. That result leaves the RATIO question open on its own: holding row order unstructured
(no alignment attempt at all, just HDR and BLUNT rows shuffled together), how does consistency_factor
vary as the BLUNT count goes from 0 (the trivial, already-known cons=1.0000 case, kept here only as
the top anchor of the curve) up toward a natural not_mhnhej mix?

Uses the SAME candidate pool, task and pinned seed as `nomh_standalone.py` and `kfold_align_nomh.py`
(K562 task 2af71117, seed 401) for direct comparability with those results.

    python mixrate_nomh.py [task_id]

RESULT (2026-09-21, K562, task 2af71117, pinned seed 401): the best non-trivial mix is the smallest
possible contamination -- 1 BLUNT row of 250 (249/1) -- giving consistency_factor ~0.4999, roughly
HALF the trivial all-HDR ceiling. The curve falls steeply and monotonically from there: blunt=2
0.4801, blunt=3 0.4577, blunt=4 0.3946, blunt=5 0.3404. Past ~6 rows it stops falling and flattens
into a noisy plateau, roughly 0.22-0.37, all the way out to 50% BLUNT (125/250) -- a coarser sweep
(0,2,5,...,125) found no further systematic decline past that point, just noise (one outlier crash
to 0.2283 at blunt=8, where is_hdr r2 hit -0.7278, an isolated instance of the same
confidently-wrong-prediction behaviour `kfold_align_nomh.py` found more of when BLUNT rows were
deliberately segregated). So essentially all the damage happens in the first handful of
contaminating rows -- consistent with is_hdr's R2 term contributing nothing either way (clipped to 0
by max(avg_r2,0) whenever negative, which it is at every non-trivial point tested), so
consistency_score is driven almost entirely by is_hdr's normalized MAE, which scales with BLUNT
count roughly as expected for a model that still predicts close to "always HDR" once BLUNT is rare
but is not exploitable at all once BLUNT is common. Practical upshot for anyone trying to use
not_mhnhej this way: minimising BLUNT count, not finding a "sweet spot" ratio, is the only lever --
there is no non-trivial mix that beats near-zero contamination, and once contamination exceeds a
handful of rows the ratio stops mattering much at all.
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "mixratenomh")
CELL = os.getenv("KA_CELL", "K562")
os.environ["EH_CELL"] = CELL

import json                                                # noqa: E402
import logging                                             # noqa: E402
import random                                               # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import seed_depend as SD         # noqa: E402
from niome_subnet.genomics.validation import (              # noqa: E402
    run_stage12, run_stage3, run_stage4, _parse_seeds)
from niome_subnet.utils import settings                     # noqa: E402
from sd_task import task_content                            # noqa: E402
from energy_hdr_k562 import pick_task, candidates_by_cell    # noqa: E402

N_ROWS = 250
BLUNT_SWEEP = (0, 1, 2, 3, 4, 5, 6, 7, 8)


def select(hdr_by_cell, blunt_by_cell, cells, hdr_take, blunt_take):
    n_cells = len(cells)
    per_hdr = hdr_take // n_cells
    per_blunt = blunt_take // n_cells
    hdr_pool, blunt_pool = [], []
    for i, cell in enumerate(cells):
        hn = per_hdr + (1 if i < hdr_take - per_hdr * n_cells else 0)
        bn = per_blunt + (1 if i < blunt_take - per_blunt * n_cells else 0)
        hdr_pool.extend(hdr_by_cell[cell][:hn])
        blunt_pool.extend(blunt_by_cell[cell][:bn])
    return hdr_pool[:hdr_take], blunt_pool[:blunt_take]


def build_rows(hdr_pool, blunt_pool, rng):
    combined = list(hdr_pool) + list(blunt_pool)
    rng.shuffle(combined)          # unstructured: no attempt to align with KFold's fold positions
    rows, seen = [], set()
    for i, rec in enumerate(combined):
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": rec["guide"],
                     "target_alignment_start": rec["start"],
                     "target_alignment_end": rec["start"] + rec["length"],
                     "strand": rec["strand"], "mutation": rec["mutation"],
                     "cas_system": rec["cas_system"], "cell_type": CELL})
    return rows


def main():
    task_id = _ARGV[1] if len(_ARGV) > 1 else pick_task()
    task, contract, reference = task_content(task_id)
    seeds = _parse_seeds(contract["seed"])
    seed = seeds[0]
    print(f"task {task_id[:8]}  {CELL}  real seeds {seeds}  pinned to seed {seed}\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    cfg = SD.SeedDependConfig(rule="not_mhnhej", variants_per_site=300)
    ctx, by_cell = candidates_by_cell(contract, reference, cell_types, cfg, seed)
    print(f"cells with candidates: {len(by_cell)}/8", flush=True)

    hdr_by_cell, blunt_by_cell = {}, {}
    for cell, recs in by_cell.items():
        hdr_by_cell[cell] = [r for r in recs if r["record"]["outcome"] == "HDR"]
        blunt_by_cell[cell] = [r for r in recs if r["record"]["outcome"] == "BLUNT_NHEJ"]
    cells = sorted(by_cell)
    min_blunt_per_cell = min(len(blunt_by_cell[c]) for c in cells)
    print(f"min BLUNT candidates in any one cell: {min_blunt_per_cell} "
          f"(caps how high blunt_take can go per-cell before a cell runs short)\n", flush=True)

    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))

    print("=== consistency_factor vs BLUNT row count, unstructured order, n=250 ===")
    results = []
    for blunt_take in BLUNT_SWEEP:
        hdr_take = N_ROWS - blunt_take
        hdr_pool, blunt_pool = select(hdr_by_cell, blunt_by_cell, cells, hdr_take, blunt_take)
        if len(hdr_pool) < hdr_take or len(blunt_pool) < blunt_take:
            print(f"  blunt={blunt_take:3d}: a cell ran short, skipped")
            continue
        rows = build_rows(hdr_pool, blunt_pool, random.Random(42))
        if len(rows) < N_ROWS:
            print(f"  blunt={blunt_take:3d}: deduped to {len(rows)} of {N_ROWS}, skipped")
            continue
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        run_stage3(seed=seed)
        r4 = run_stage4(seed=seed)
        hdr_r2 = r4["model_results"].get("is_hdr", {}).get("r2_mean", float("nan"))
        hdr_mae = r4["model_results"].get("is_hdr", {}).get("mae_mean", float("nan"))
        results.append((blunt_take, r4["consistency_factor"], hdr_r2, hdr_mae,
                        r4["total_weighted_score"]))
        print(f"  blunt={blunt_take:3d} ({100 * blunt_take / N_ROWS:4.1f}%)  "
              f"hdr={hdr_take:3d}  consistency_factor {r4['consistency_factor']:.4f}  "
              f"is_hdr r2 {hdr_r2:+.4f}  mae {hdr_mae:.4f}  weighted {r4['total_weighted_score']:.1f}",
              flush=True)

    print()
    non_trivial = [r for r in results if r[0] > 0]
    if non_trivial:
        best = max(non_trivial, key=lambda r: r[1])
        print(f"best NON-TRIVIAL mix (excluding blunt=0): blunt={best[0]} "
              f"({100 * best[0] / N_ROWS:.1f}%)  consistency_factor {best[1]:.4f}")
    trivial = [r for r in results if r[0] == 0]
    if trivial:
        print(f"trivial anchor: blunt=0  consistency_factor {trivial[0][1]:.4f}")


if __name__ == "__main__":
    main()
