#!/usr/bin/env python3
"""joined300.py — all-HDR measured at the JOINED 300-seed band space, per cell type.

The fleet's band space is the *joined* window: three width-100 classes concatenated (300 seeds),
with h1-h10 taking rotated width-225 slices of it and h0 the whole 300. Everything CLAUDE.md
records about all-HDR's band, weighted and fidelity was measured on a CONTIGUOUS width-100 window,
and the joined numbers exist only in fragments: band 12 at joined 225 / 11 at joined 300 on K562,
7 at both on HEK293, fidelity 0.888/0.877 on one CD34+ contract. There is no per-cell table at the
width the fleet actually ships, and no measurement at all of what a hotkey searching the whole
900-seed space would get.

Four arms per cell type, so window WIDTH and CONTIGUITY are separable and the "total 900" case is
priced rather than assumed:

    contig100   the CELL_CONFIG default — where every tuned number was measured
    joined225   the live rotated slice (h1's assignment in data/window_plan.json)
    joined300   the live full joined space (h0's assignment)
    contig900   one hotkey over the entire 100-999 space
    joined300nc a deliberately non-contiguous 300 (100-199, 400-499, 800-899) — only needed where
                the live plan's three classes happen to be adjacent, as K562's 200-499 are

What is measured per arm:

  * **the TRUE band** — the seeds where all 250 assembled rows repair by HDR, scanned over the full
    100-999 rather than trusting `meta["clean"]`. `meta["clean"]` is the Cas12a group's clean set
    *inside the band space*; the assembled submission can be clean on seeds the search never saw
    (CLAUDE.md records h3's 12-seed band leaking five seeds below its window), and a band seed is
    worth the same wherever it sits. `leak` reports the difference.
  * **weighted x fidelity** at that band, from the validator's own stages, plus the spike
    consistency at a true-band seed (verified, not assumed) and the floor off it.

Fidelity and weighted are design-only quantities and so seed-independent (verified in
fidelity_window.py); the seed matters only for stage 3/4, which is why one clean seed suffices.

    python joined300.py                       # all four cell types, all four arms
    J3_CELLS=K562 J3_ARMS=joined300 python joined300.py
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "j300")

import json                          # noqa: E402
import logging                       # noqa: E402
import random                        # noqa: E402
import statistics as st              # noqa: E402
import time                          # noqa: E402
import urllib.request                # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.validation import stage3         # noqa: E402
from niome_subnet.genomics.validation.stage12 import run_stage12   # noqa: E402
from niome_subnet.utils import settings                     # noqa: E402
from sd_task import score                                   # noqa: E402

API = "https://niome-api.genomes.io/api/v3"
CELLS = [c for c in os.getenv("J3_CELLS", "K562,HUDEP-2,CD34+_HSPC,HEK293").split(",") if c]
ARMS = [a for a in os.getenv("J3_ARMS",
                             "contig100,joined225,joined300,contig900").split(",") if a]
BUDGET = float(os.getenv("J3_BUDGET", "2400"))
FLOOR_N = int(os.getenv("J3_FLOOR_N", "3"))
PLAN = os.getenv("J3_PLAN", "data/window_plan.json")
OUT = os.getenv("J3_OUT", "joined300.json")
# A non-contiguous 300 for the joined300nc arm: window indices 0, 3, 7 of the nine width-100
# classes. Fixed rather than per cell so contiguity is the only thing that differs from joined300.
NC_CLASSES = (0, 3, 7)


def fetch_tasks(limit=90):
    doc = json.load(urllib.request.urlopen(f"{API}/tasks?limit={limit}"))
    return doc if isinstance(doc, list) else (doc.get("items") or doc.get("data") or [])


def pick_tasks(items):
    """The newest STAMPED task per cell type — a round that really played, one per cell."""
    out = {}
    for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        cell = c.get("cell_type")
        if cell not in CELLS or cell in out:
            continue
        if str(c.get("seed", "0")).strip() in ("0", ""):
            continue                      # unstamped: the live round, still in flight
        out[cell] = t
    for cell in CELLS:
        override = os.getenv("J3_TASK_" + cell.replace("+", "").replace("-", "_"))
        if override:
            out[cell] = next(t for t in items
                             if (t.get("task_id") or t.get("id", "")).startswith(override))
    return out


def seeds_of(ranges):
    return sorted(s for lo, hi in ranges for s in range(lo, hi + 1))


def arm_spec(name, cell, plan):
    """(build kwargs, label) for one arm, or None where the arm does not apply."""
    if name == "contig100":
        return {"hdr_range": AH.CELL_CONFIG[cell]["hdr_range"]}
    if name == "contig900":
        return {"hdr_range": (100, 999)}
    if name == "joined225":
        return {"seed_list": seeds_of(plan[cell]["niome_hotkey1"])}
    if name == "joined300":
        return {"seed_list": seeds_of(plan[cell]["niome_hotkey"])}
    if name == "joined300nc":
        return {"seed_list": seeds_of([[w * 100 + 100, w * 100 + 199] for w in NC_CLASSES])}
    if name.startswith("contig"):
        # contig<N>: the first N seeds of the cell's default window. The decomposition below the
        # main table says the whole layout cost is BAND SIZE and not the space the fleet shares,
        # and the band grows as the window narrows (CLAUDE.md's `narrow_width` row measures band
        # 10 of 10 at width 10 on HEK293 against 7 at width 100), so whether a sub-100 window pays
        # more band is the one lever the pricing points at.
        lo = AH.CELL_CONFIG[cell]["hdr_range"][0]
        return {"hdr_range": (lo, lo + int(name[len("contig"):]) - 1)}
    raise SystemExit(f"unknown arm {name}")


def true_band(rows, contract, reference, cell_types, lo=100, hi=999):
    """Seeds in [lo, hi] where EVERY assembled row repairs by HDR.

    Runs stage 12 once and then stage 3's own `simulate` per seed, bailing on the first row that is
    not HDR — which is nearly always the first or second, so the full 900-seed sweep costs seconds.
    """
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    band = []
    for seed in range(lo, hi + 1):
        if all(stage3.simulate(exp, seed)["outcome"] == "HDR" for exp in valid):
            band.append(seed)
    return band, len(valid)


def measure(cell, task, arm, plan, rng, contract, reference, cell_types):
    kw = arm_spec(arm, cell, plan)
    space = (kw["seed_list"] if "seed_list" in kw
             else list(range(kw["hdr_range"][0], kw["hdr_range"][1] + 1)))
    span = len(space)
    base = AH.config_for(cell)
    mf = (base.main_max_fail if span == 100
          else AH._scaled_max_fail(cell, base.main_max_fail, span))
    rec = {"cell": cell, "arm": arm, "task": (task.get("task_id") or task["id"])[:8],
           "span": span, "window": AH._range_label(space), "main_max_fail": mf,
           "group_size": base.group_size}
    t0 = time.monotonic()
    rows, meta = AH.build_for_cell(contract, reference, cell_types, budget_s=BUDGET, **kw)
    rec["build_s"] = round(time.monotonic() - t0, 1)
    rec["meta"] = {k: v for k, v in meta.items() if k != "bank_path"}
    if not rows:
        rec["declined"] = meta.get("reason", "declined")
        print(f"  {arm:<12} DECLINED after {rec['build_s']:.0f}s — {rec['declined']}")
        MT.free_gpu_memory()
        return rec

    t1 = time.monotonic()
    band, n_valid = true_band(rows, contract, reference, cell_types)
    inside = [s for s in band if s in set(space)]
    rec.update(rows=len(rows), valid=n_valid, band=len(band), band_in_window=len(inside),
               leak=len(band) - len(inside), band_seeds=band,
               group_clean=meta.get("clean"), bank=meta.get("bank"),
               cas9_pool=meta.get("cas9_pool"), cells=meta.get("cells"),
               cas_mix=meta.get("cas_mix"), scan_s=round(time.monotonic() - t1, 1))

    # One full 5-stage pass on a true-band seed: weighted and fidelity (both seed-independent)
    # plus the spike itself, verified rather than assumed.
    if band:
        hit = score(rows, contract, reference, cell_types, seed=band[len(band) // 2])
        rec["band_seed"] = band[len(band) // 2]
    else:
        hit = score(rows, contract, reference, cell_types, seed=space[0])
    off = [s for s in range(100, 1000) if s not in set(band)]
    floors = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(off, FLOOR_N)]
    rec.update(weighted=hit["weighted"], fidelity=hit["fidelity"],
               wxfid=hit["weighted"] * hit["fidelity"], spike_cons=hit["consistency"],
               floor_cons=st.mean(floors))
    MT.free_gpu_memory()
    print(f"  {arm:<12} span {span:<4} mf {mf:<4} bank {rec['bank']:<6} "
          f"band {rec['band']:<3} (in-win {rec['band_in_window']}, leak {rec['leak']})  "
          f"w {rec['weighted']:7.1f}  fid {rec['fidelity']:.4f}  wxfid {rec['wxfid']:7.1f}  "
          f"spike {rec['spike_cons']:.4f}  floor {rec['floor_cons']:.4f}  "
          f"cells {rec['cells']}/8  {rec['build_s']:.0f}s")
    return rec


def main():
    plan = json.load(open(PLAN))["assignments"]
    tasks = pick_tasks(fetch_tasks())
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    rng = random.Random(20260910)
    results = []
    for cell in CELLS:
        task = tasks.get(cell)
        if task is None:
            print(f"{cell}: no stamped task in the feed, skipped")
            continue
        content = task["content"]
        contract, reference = content["contract"], content["hbb_reference"]
        print(f"\n{cell}  task {(task.get('task_id') or task['id'])[:8]}  "
              f"{task.get('created_at', '')[:16]}  seeds {contract.get('seed')}  "
              f"joined {AH._range_label(seeds_of(plan[cell]['niome_hotkey']))}")
        for arm in ARMS:
            try:
                results.append(measure(cell, task, arm, plan, rng, contract, reference,
                                       cell_types))
            except Exception as exc:                       # one arm must not lose the run
                print(f"  {arm:<12} ERROR {type(exc).__name__}: {exc}")
                results.append({"cell": cell, "arm": arm, "error": f"{type(exc).__name__}: {exc}"})
            json.dump(results, open(OUT, "w"), indent=1)
    json.dump(results, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
