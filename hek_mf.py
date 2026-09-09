#!/usr/bin/env python3
"""hek_mf.py — HEK293's all-HDR `main_max_fail` curve at a width-300 band window.

The width-300 fleet layout breaks all-HDR on HEK293: linear span scaling gives `mf 144` and the
bank comes back with 53 qualifying guides against a group of 80, so all seven band hotkeys decline
and the fleet ships eight flat builds on ~a quarter of rounds.

The mechanism is that HEK293's P(HDR) is ~0.37 against the erythroid types' ~0.57, so `mf` sits far
out in the tail rather than near the mean. Linear scaling preserves the THRESHOLD RATE (52% of
seeds must repair by HDR at both widths) but the binomial narrows as the window grows, so the same
rate is a much deeper tail cut at 300 than at 100. Matching the z-score instead:

    mf(span) = (1-p)*span - z*sqrt(span*p*(1-p))

Raising `mf` costs selectivity though, and selectivity is what the min-union turns into a band — so
the question this answers is whether any `mf` both fills the group AND keeps a band worth having,
or whether HEK293 simply needs to stay at width 100.

    python hek_mf.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hekmf")

import dataclasses as _dc      # noqa: E402
import json                    # noqa: E402
import logging                 # noqa: E402
import math                    # noqa: E402
import time                    # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np             # noqa: E402

import genExp as G             # noqa: E402
from niome_subnet.genomics import all_cut as AC       # noqa: E402
from niome_subnet.genomics import all_hdr as AH       # noqa: E402
from niome_subnet.genomics import mt19937 as MT       # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3   # noqa: E402
from sd_task import task_content                      # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "414dab89-1bf9-4f16-ba6f-93c43186b8cf"
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "100-399"


def median_fails(contract, reference, cell_types, ctx, sites, cfg, n_targets=6, n_var=2000):
    """The true fail-count distribution at this window, with no cap filtering it."""
    cell = contract["cell_type"]
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    seeds = np.arange(cfg.start_seed, cfg.end_seed + 1, dtype=np.int64)
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == "Cas12a" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance][:n_targets]
    counts = []
    for si, mut in jobs:
        site = sites[si]
        d = abs(site.start - ctx.mutation_map[mut])
        off = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mut), 0.0)
        gs = SA.enumerate_variants(site, ctx, cfg.cas12a_gc[0], cfg.cas12a_gc[1],
                                   ctx.max_mismatches, True, n_var)
        if not gs:
            continue
        pf = AH._params_fn(site, d, acc, off)
        got = MT.screen_guides_rule_gpu(gs, seeds, mut, site.cas, site.start, site.strand,
                                        pf, "hdr", len(seeds))
        counts.extend(len(f) for f in got.values())
    MT.free_gpu_memory()
    return np.asarray(counts)


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = (int(x) for x in WINDOW.split("-"))
    span = hi - lo + 1
    base = AH.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    print(f"task {TASK[:8]}  {cell}  accessibility "
          f"{cell_types.get(cell, {}).get('accessibility')}  window {lo}-{hi} (span {span})  "
          f"group {base.group_size}\n")

    # Ground the z arithmetic in the real distribution rather than the nominal P(HDR).
    for probe_span, probe_w in ((100, (lo, lo + 99)), (span, (lo, hi))):
        c = _dc.replace(base, hdr_range=probe_w)
        f = median_fails(contract, reference, cell_types, ctx, sites, c)
        p_hdr = 1.0 - f.mean() / probe_span
        sd = f.std()
        mf_ref = base.main_max_fail if probe_span == 100 else None
        z = (mf_ref - f.mean()) / sd if mf_ref and sd else None
        print(f"  span {probe_span:>3}: mean fails {f.mean():6.1f}  sd {sd:5.1f}  "
              f"-> P(HDR) {p_hdr:.3f}"
              + (f"   mf {mf_ref} sits at z {z:+.2f}" if z is not None else ""))
        if probe_span == 100:
            z_target, mean100, sd100 = z, f.mean(), sd
        else:
            mean300, sd300 = f.mean(), sd
    mf_z = max(1, round(mean300 + z_target * sd300))
    print(f"\n  z-matched mf at span {span}: {mf_z}  (linear scaling gives "
          f"{round(base.main_max_fail * span / 100.0)})\n")

    mfs = sorted({round(base.main_max_fail * span / 100.0), mf_z, mf_z + 20, mf_z + 40,
                  int(mean300)})
    print(f"  {'mf':>5} {'z':>6} {'bank':>7} {'band':>5} {'union':>6} {'cas9':>6} {'cells':>5} "
          f"{'rows':>5} {'build s':>8}  outcome")
    out = []
    for mf in mfs:
        cfg = _dc.replace(base, hdr_range=(lo, hi), main_max_fail=mf,
                          variants=min(base.variants, AH.WIDE_WINDOW_VARIANTS))
        p = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
        if os.path.exists(p):
            os.remove(p)
        t0 = time.monotonic()
        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=600)
        dt = time.monotonic() - t0
        MT.free_gpu_memory()
        z = (mf - mean300) / sd300 if sd300 else float("nan")
        rec = {"mf": mf, "z": round(z, 2), "bank": meta.get("bank"),
               "band": meta.get("clean"), "union": meta.get("union"),
               "cas9": meta.get("cas9_pool"), "cells": meta.get("cells"),
               "rows": len(rows) if rows else 0, "build_s": round(dt, 1),
               "reason": None if rows else meta.get("reason")}
        out.append(rec)
        print(f"  {mf:>5} {z:>+6.2f} {str(meta.get('bank') or 0):>7} "
              f"{str(meta.get('clean') or '-'):>5} {str(meta.get('union') or '-'):>6} "
              f"{str(meta.get('cas9_pool') or '-'):>6} "
              f"{(str(meta.get('cells'))+'/8') if meta.get('cells') else '-':>5} "
              f"{rec['rows']:>5} {dt:>8.0f}  "
              + ("built" if rows else f"declined: {meta.get('reason')}"))
    good = [r for r in out if r["rows"]]
    print()
    if good:
        best = max(good, key=lambda r: r["band"])
        print(f"  best band at width {span}: mf {best['mf']} -> band {best['band']}  "
              f"(7 hotkeys = {7*best['band']} seeds, "
              f"{1-(1-7*best['band']/900)**3:.1%} spike rate)")
        print(f"  HEK293 at width 100 for comparison: band 7 -> 6 hotkeys = 42 seeds, "
              f"{1-(1-42/900)**3:.1%}")
    else:
        print(f"  NO mf builds at width {span} -> HEK293 needs a per-cell window width.")
    json.dump({"task": TASK, "cell": cell, "window": [lo, hi], "z_target": z_target,
               "mean_fails_300": mean300, "sd_300": sd300, "arms": out},
              open("hek_mf.json", "w"), indent=1)
    print("\nwrote hek_mf.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
