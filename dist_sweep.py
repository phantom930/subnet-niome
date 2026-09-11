#!/usr/bin/env python3
"""dist_sweep.py — does tightening max_distance ALONE buy base structural without costing band?

`total_weighted_score` = sum over rows of `base_structural * offtarget_factor * mutation_weight`,
with `base_structural = 0.625*gc_score + 0.375*dist_score`. Measured on the live build (571f4843,
CD34+, after the cas9_cell_target fix):

    gc_score 0.9337   dist_score 0.8477   base 0.9015   offtarget 1.0000

gc_score is nearly exhausted; **dist_score at 0.848 is where the room is**. Solving the competitor's
row backwards, ranks 8-11 reach weighted 262.8 at fidelity 0.940, which needs base ~1.005 against
our 0.9015 — they have near-perfect structural terms, not a heavier mutation skew. Taking our base
to 1.000 moves weighted 252.6 -> 280.2 and the product 221.6 -> 245.8, against rank 10's 246.8. That
one term is essentially the whole rank-10/rank-11 gap.

**Why this is not the falsified distance/GC row.** That arm tightened `cas12a_gc`/`cas9_gc` to
0.42-0.58 *together with* max_distance 250, and lost 6-18% on frequency x value by costing band on
3 of 4 conditions. Its own notes record that tightening GC cost fidelity monotonically while
tightening DISTANCE raised it. This sweep moves distance only, leaves the GC bounds alone, and
scores band as a first-class outcome rather than a side effect.

**And not the 210-config sweep either.** That priced E[pay] at max_distance {400, 150, 100} and put
400 ahead (0.0242 / 0.0237 / 0.0235 E[eff]) — but 100 had the best score *among builds* (0.0257)
and lost only on a decline, and the three points are non-monotone, which is the shape of noise
rather than a trend. It never tested 200-300, and it reports no band or dist_score, so it cannot say
whether a band-preserving sweet spot exists. That is exactly what this prints.

Same contract, same joined window, group and every other knob fixed. The decision quantity is
`band x weighted x fidelity` — band is spike frequency, the product is spike value.
"""
import os

os.environ.setdefault("NIOME_INSTANCE", "distsweep")

import dataclasses
import json
import logging
import statistics
import sys
import time
import urllib.request

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
logging.basicConfig(level=logging.ERROR)

from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics.validation.stage12 import run_stage12
from niome_subnet.genomics.validation.stage3 import run_stage3
from niome_subnet.genomics.validation.stage4 import run_stage4
from niome_subnet.genomics.validation.stage5 import run_stage5
from niome_subnet.utils import settings

TASK = os.environ.get("DS_TASK", "571f4843-1671-49c3-951c-47d16a9418b0")
DISTANCES = [int(x) for x in os.environ.get("DS_D", "400,300,250,200,150,100").split(",")]
# the live rotated slice h4 runs (width 225), so band numbers are comparable to the fleet's
WINDOW = sorted(set(range(100, 115)) | set(range(190, 300)) | set(range(800, 900)))
REFERENCE = None


def load_case():
    doc = json.load(urllib.request.urlopen("https://niome-api.genomes.io/api/v3/tasks?limit=40"))
    items = doc if isinstance(doc, list) else doc.get("data") or doc.get("items")
    task = next(t for t in items if (t.get("task_id") or t.get("id")) == TASK)
    cells = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/data/cell-types?format=json"))
    return task["content"]["contract"], task["content"]["hbb_reference"], cells


def score(rows, contract, cells):
    os.makedirs(f"data/inst/{os.environ['NIOME_INSTANCE']}", exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(REFERENCE, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cells); run_stage3(seed=500); run_stage4(seed=500)
    s5 = run_stage5()
    v = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    gc = statistics.mean(i["features"]["gc_score"] for i in v)
    ds = statistics.mean(i["features"]["dist_score"] for i in v)
    dist = statistics.mean(i["features"]["distance_to_mutation"] for i in v)
    w = json.load(open(settings.FINAL_REWARD_PATH))["total_weighted_score"]
    return dict(gc=gc, ds=ds, base=0.625 * gc + 0.375 * ds, dist=dist, n=len(v),
                weighted=w, fid=s5["distribution_fidelity_factor"])


def main():
    global REFERENCE
    contract, reference, cells = load_case()
    REFERENCE = reference
    cell = contract["cell_type"]
    base_cfg = AH.config_for(cell)
    span = len(WINDOW)
    cfg0 = dataclasses.replace(
        base_cfg, hdr_range=(WINDOW[0], WINDOW[-1]), seed_list=tuple(WINDOW),
        main_max_fail=AH._scaled_max_fail(cell, base_cfg.main_max_fail, span),
        variants=min(base_cfg.variants, AH.WIDE_WINDOW_VARIANTS))
    print(f"task {TASK[:8]}  {cell}  joined window {span} seeds, group {cfg0.group_size}, "
          f"mf {cfg0.main_max_fail}\n")
    print("%-7s %-6s %-8s %-8s %-8s %-8s %-9s %-8s %-9s %s" % (
        "maxd", "band", "meandist", "dist_sc", "gc_sc", "base", "weighted", "fid", "w x fid",
        "band x w x fid"))
    for d in DISTANCES:
        t0 = time.monotonic()
        cfg = dataclasses.replace(cfg0, max_distance=d)
        rows, meta = AH.build_submission(contract, reference, cells, cfg=cfg, budget_s=900)
        if not rows:
            print("%-7d DECLINED — %s" % (d, meta.get("reason"))); continue
        m = score(rows, contract, cells)
        band = meta.get("clean")
        prod = m["weighted"] * m["fid"]
        print("%-7d %-6s %-8.1f %-8.4f %-8.4f %-8.4f %-9.1f %-8.4f %-9.1f %-8.0f (%.0fs)" % (
            d, band, m["dist"], m["ds"], m["gc"], m["base"], m["weighted"], m["fid"], prod,
            band * prod, time.monotonic() - t0))
    print("\nshipped is max_distance 400. band is spike FREQUENCY, w x fid is spike VALUE;")
    print("the last column is the product that decides E[pay] at fixed group size.")


if __name__ == "__main__":
    main()
