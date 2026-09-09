#!/usr/bin/env python3
"""cas9_floor.py — give scan_cas9's early break a PER-CELL floor, and price the result.

Measured defect (see the funnel in this session): `all_hdr.scan_cas9` walks Cas9 sites nearest-first
and breaks on

    len(found) >= want * pool_target and len({(mutation, strand)}) == 4

The cell check only requires each (mutation, strand) cell to be NON-EMPTY. On HEK293 a8b9f1bb the
break fires at job 34 of 395 with HEAVY+ holding **12** candidates against **5531** available, so
`assemble` cannot fill its heavy quota (it wants ~79 on that strand) and backfills with light rows.
That is the whole of HEK293's mean mutation_weight 0.78 against a reachable ~1.02, i.e. weighted
164 against the rank-1 miner's 216.

This arm requires every cell to reach ``per_cell_floor`` before breaking. Scanning further is the
cost: the full 395-job scan takes 67s against the break's 5s, and heavy(+) is 17% of jobs, so a
floor of ~80 should land near job 90.

**The result may well be negative.** With availability fixed, `assemble`'s quota finally realises
the 158 heavy / 12 light split that ``light_cell_rows = 6`` encodes -- and all_cut.py records that
capped split as measured dead, because it drives stage 5's mutation-coverage entropy to its floor
and fidelity falls faster than weighted rises. If that reproduces here, the defect is real but
fixing it buys nothing, and the honest conclusion is that HEK293's weighted is capped by stage 5.

    python cas9_floor.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "c9floor")

import dataclasses as _dc     # noqa: E402
import json                   # noqa: E402
import logging                # noqa: E402
import statistics as st       # noqa: E402
import time                   # noqa: E402
from collections import Counter  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np            # noqa: E402

import genExp as G            # noqa: E402
from niome_subnet.genomics import all_hdr as AH      # noqa: E402
from niome_subnet.genomics import mt19937 as MT      # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3  # noqa: E402
from niome_subnet.utils import settings              # noqa: E402
from sd_task import score, task_content              # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "a8b9f1bb-f25d-42af-92a1-e0aff16bee3a"
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "300-399"
FLOOR = 0.10
FLOORS = [int(x) for x in (os.getenv("C9_FLOORS") or "42,80").split(",")]
_REAL_SCAN = AH.scan_cas9


def scan_floored(clean, contract, cell_types, ctx, sites, cfg, want, deadline=None,
                 per_cell_floor=0):
    """all_hdr.scan_cas9 with a per-cell floor on the early break. Body kept in step with the
    original deliberately -- only the break condition differs."""
    cell = contract.get("cell_type")
    accessibility = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    jobs = [(i, m, abs(sites[i].start - ctx.mutation_map[m]))
            for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(sites[i].start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda job: job[2])
    found = []
    per = Counter()
    for site_index, mutation, distance in jobs:
        if deadline is not None and time.monotonic() > deadline:
            break
        site = sites[site_index]
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = AH._params_fn(site, distance, accessibility, offset)
        for guide in MT.screen_guides_rule_gpu(guides, clean, mutation, "Cas9", site.start,
                                               site.strand, params_of, "hdr", 0):
            gc, _e, _c = params_of(guide)
            found.append({"guide": guide, "mutation": mutation, "cas_system": "Cas9",
                          "strand": site.strand, "start": site.start, "length": site.length,
                          "gc": gc, "distance": distance})
            per[(mutation, site.strand)] += 1
        enough_cells = len(per) == 4 and min(per.values()) >= per_cell_floor
        if len(found) >= want * cfg.pool_target and enough_cells:
            break
    return found


def report(rows, meta, contract, cell_types, reference, cfg, group_n, dt, label, lo):
    if not rows:
        return {"arm": label, "reason": meta.get("reason"), "build_s": round(dt, 1)}
    s = score(rows, contract, reference, cell_types, seed=lo)
    feats = [d["features"] for d in json.load(open(settings.VALID_EXPERIMENTS_PATH))]
    mw = contract.get("mutation_weights", {})
    heavy = max(mw, key=lambda m: mw.get(m, 1.0))
    band = meta["clean"]
    per = Counter(("HEAVY" if r["mutation"] == heavy else "light", r["strand"],
                   r["cas_system"]) for r in rows)
    k1 = s["weighted"] * ((1 + 2 * FLOOR) / 3) * s["fidelity"]
    cov = 1 - (1 - 7 * band / 900.0) ** 3
    return {"arm": label, "band": band, "cas9": meta.get("cas9_pool"),
            "heavy_rows": sum(1 for r in rows if r["mutation"] == heavy),
            "mean_weight": st.mean(mw.get(r["mutation"], 1.0) for r in rows),
            "weighted": s["weighted"], "fidelity": s["fidelity"],
            "structural": st.mean(0.625 * f["gc_score"] + 0.375 * f["dist_score"] for f in feats),
            "cells": len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
            "k1_final": k1, "fleet_cov": cov, "freq_x_val": cov * k1,
            "split": {f"{a}{b}": n for (a, b, c), n in sorted(per.items()) if c == "Cas9"},
            "build_s": round(dt, 1)}


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = (int(x) for x in WINDOW.split("-"))
    span = hi - lo + 1
    base = AH.config_for(cell)
    cfg = _dc.replace(base, hdr_range=(lo, hi),
                      main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span),
                      variants=(min(base.variants, AH.WIDE_WINDOW_VARIANTS)
                                if span > 100 else base.variants))
    print(f"task {TASK[:8]}  {cell}  window {lo}-{hi}  group {cfg.group_size}  "
          f"mutation_weights {contract.get('mutation_weights')}\n")
    print(f"  {'arm':<22} {'band':>5} {'cas9':>6} {'heavy':>6} {'mean_w':>7} {'struct':>7} "
          f"{'weighted':>9} {'fid':>7} {'k=1':>7} {'cov':>6} {'freq*val':>9} {'s':>5}  Cas9 split")
    out = []
    for label, floor in [("baseline (non-empty)", None)] + [(f"per-cell floor {f}", f)
                                                            for f in FLOORS]:
        if floor is None:
            AH.scan_cas9 = _REAL_SCAN
        else:
            AH.scan_cas9 = (lambda *a, _f=floor, **k:
                            scan_floored(*a, **k, per_cell_floor=_f))
        t0 = time.monotonic()
        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=800)
        dt = time.monotonic() - t0
        rec = report(rows, meta, contract, cell_types, reference, cfg, cfg.group_size, dt,
                     label, lo)
        MT.free_gpu_memory()
        out.append(rec)
        if "reason" in rec:
            print(f"  {label:<22}  declined: {rec['reason']} ({rec['build_s']}s)")
            continue
        print(f"  {label:<22} {rec['band']:>5} {rec['cas9']:>6} {rec['heavy_rows']:>6} "
              f"{rec['mean_weight']:>7.4f} {rec['structural']:>7.4f} {rec['weighted']:>9.2f} "
              f"{rec['fidelity']:>7.4f} {rec['k1_final']:>7.2f} {rec['fleet_cov']:>5.1%} "
              f"{rec['freq_x_val']:>9.2f} {rec['build_s']:>5.0f}  {rec['split']}")
    AH.scan_cas9 = _REAL_SCAN
    good = [r for r in out if r.get("band")]
    if len(good) > 1:
        b = good[0]
        print(f"\n  against baseline (weighted {b['weighted']:.1f}, fid {b['fidelity']:.4f}, "
              f"freq*val {b['freq_x_val']:.2f}):")
        for r in good[1:]:
            print(f"    {r['arm']:<22} weighted {(r['weighted']/b['weighted']-1)*100:+6.1f}%  "
                  f"fidelity {(r['fidelity']/b['fidelity']-1)*100:+6.1f}%  "
                  f"band {r['band']-b['band']:+d}  "
                  f"k=1 {(r['k1_final']/b['k1_final']-1)*100:+6.1f}%  "
                  f"freq*val {(r['freq_x_val']/b['freq_x_val']-1)*100:+6.1f}%  "
                  f"(+{r['build_s']-b['build_s']:.0f}s)")
    json.dump({"task": TASK, "cell": cell, "window": [lo, hi], "arms": out},
              open(os.getenv("C9_OUT", "cas9_floor.json"), "w"), indent=1)
    print(f"\nwrote {os.getenv('C9_OUT', 'cas9_floor.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
