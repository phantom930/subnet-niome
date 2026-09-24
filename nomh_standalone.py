#!/usr/bin/env python3
"""nomh_standalone.py — the not_mhnhej rule as its OWN construction, no all-cut stage at all,
and a direct test of the REVERSED energy hypothesis (low energy = HDR, high energy = BLUNT).

Every earlier not_mhnhej test in this repo combined it with something else: `conj_nomh.py` used it
as `conjunction.py`'s BAND rule on top of an outer all-cut min-union (falsified in CLAUDE.md --
"dead on both cells"), and `energy_hdr_k562.py`/`energy_combo_k562.py` sampled the raw
`SD.enumerate_candidates` pool by hand (nearest-first / spread / uniform-high) rather than running
the actual `weighted_score x fidelity` allocator. This script runs `seed_depend`'s real allocator
(`SD.allocate`, the same one `SD.build()` calls) with `cfg.rule="not_mhnhej"` and NOTHING else --
no cut-min-union, no band, no conjunction -- so it is the actual standalone construction a hotkey
would ship if it ran not_mhnhej on its own, scored through the real five-stage pipeline.

**Direction under test.** `stage3.repair_mode` gives `hdr_w = hdr_base + 0.35*energy` against a
FIXED `blunt_w = 0.35`, so the formula's OWN prediction is P(HDR) rising WITH energy (high energy =
more HDR-favouring, low energy = more BLUNT-favouring). The hypothesis this script tests is the
OPPOSITE sign -- low energy = HDR, high energy = BLUNT -- prompted by the one negative correlation
already seen in this file's sibling (`energy_hdr_k562.py`'s K562 `spread` arm: Pearson -0.0362,
wrong-signed against the formula, though inside noise at n=240). This run reports which sign the
ACTUAL allocated construction's rows show, rather than assuming either direction, and scores that
construction through the real pipeline at all three real stamped seeds.

    NOMH_CELL=K562 python nomh_standalone.py [task_id]

RESULT (2026-09-21, K562, task 2af71117, real seeds 401/403/249, built pinned to seed 401): the
reversed hypothesis cannot even be tested. The real allocator's own weighted_score ranking (nearest,
highest-GC guides first, same mechanism as energy_hdr_k562.py's `default`/`uniform_high` arms) drives
energy to mean 1.0000 / std 0.0000 across all 250 allocated rows -- zero variance, so
Pearson(energy, is_hdr) = +0.0000 by necessity and neither "low energy = HDR" nor the formula's own
"high energy = HDR" has anything to ride on. Tercile P(HDR) (0.5542/0.7349/0.6190) is noise around
one saturated energy value, not a trend. Separately, and regardless of the energy question: this
standalone construction reaches only consistency_factor 0.2743 at its OWN build seed (not_mhnhej
pins is_cut alone, never is_hdr, so it cannot reach 1.0000 the way hdr/mh do) and averages 0.1558
over the contract's 3 real seeds -- close to the ordinary dirty-seed floor despite every row being
individually rule-compliant. See CLAUDE.md's falsified table for the full writeup and the sibling
findings this extends.
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "nomhstandalone")
CELL = os.getenv("NOMH_CELL", "K562")
os.environ["EH_CELL"] = CELL  # so energy_hdr_k562.pick_task() targets the same cell

import json                                             # noqa: E402
import logging                                          # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
from niome_subnet.genomics import seed_depend as SD      # noqa: E402
from niome_subnet.genomics.validation import (           # noqa: E402
    run_stage12, run_stage3, run_stage4, _parse_seeds)
from niome_subnet.utils import settings                  # noqa: E402
from sd_task import task_content                         # noqa: E402
from energy_hdr_k562 import pick_task, pearson           # noqa: E402


def main():
    task_id = _ARGV[1] if len(_ARGV) > 1 else pick_task()
    task, contract, reference = task_content(task_id)
    seeds = _parse_seeds(contract["seed"])
    seed = seeds[0]
    print(f"task {task_id[:8]}  {CELL}  real seeds {seeds}  building pinned to seed {seed}\n",
          flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    # variants_per_site 300, not the shipped 12000 -- see energy_hdr_k562.py's docstring: not_mhnhej
    # admits far more draws than hdr/mh, so enumerate_candidates' per-survivor dict (entry + kmers +
    # record) reaches tens of GB at the shipped value on a box whose live fleet already holds ~26 GB
    # of 49. 300 is ample for allocate()'s per-cell take.
    cfg = SD.SeedDependConfig(rule="not_mhnhej", variants_per_site=300)

    pinned = dict(contract)
    pinned["seed"] = seed
    ctx = G.build_context(pinned, reference, cell_types)
    sites = G.enumerate_sites(ctx, cfg.flank, cfg.lengths)
    n_rows = pinned["rules"].get("max_experiments") or ctx.max_experiments

    by_cell = SD.enumerate_candidates(ctx, sites, cfg, None)
    print(f"cells with candidates: {len(by_cell)}/8, sizes "
          f"{sorted(len(v) for v in by_cell.values())} "
          f"(total {sum(len(v) for v in by_cell.values())})", flush=True)
    if len(by_cell) < 8:
        raise SystemExit(f"only {len(by_cell)}/8 cells -- allocate() would decline, stop here")

    chosen, alloc = SD.allocate(by_cell, pinned, ctx, n_rows, cfg)
    if not chosen:
        raise SystemExit(f"allocate() declined: {alloc}")
    print(f"allocate(): {len(chosen)} rows, weighted {alloc['weighted']:.1f}  "
          f"fidelity {alloc['fidelity']:.4f}  heavy {alloc['heavy']}  cells {alloc['cells']}/8\n",
          flush=True)

    energies = [r["record"]["energy"] for r in chosen]
    is_hdr = [1.0 if r["record"]["outcome"] == "HDR" else 0.0 for r in chosen]
    corr = pearson(energies, is_hdr)
    print(f"=== the REAL not_mhnhej-only allocation: {len(chosen)} rows ===")
    print(f"  energy: mean {np.mean(energies):.4f} std {np.std(energies):.4f} "
          f"min {min(energies):.4f} max {max(energies):.4f}")
    print(f"  raw P(HDR) overall: {np.mean(is_hdr):.4f}")
    order = np.argsort(energies)
    n = len(order)
    for name, idx in (("low third", order[:n // 3]), ("mid third", order[n // 3:2 * n // 3]),
                      ("high third", order[2 * n // 3:])):
        e = [energies[i] for i in idx]
        h = [is_hdr[i] for i in idx]
        print(f"    {name}: n={len(idx):3d} mean energy {np.mean(e):.4f} "
              f"observed P(HDR) {np.mean(h):.4f}")
    print(f"  Pearson(energy, is_hdr) = {corr:+.4f}")
    if corr < -0.02:
        print("  -> NEGATIVE: matches the 'low energy = HDR' hypothesis under test")
    elif corr > 0.02:
        print("  -> POSITIVE: matches stage3.repair_mode's own formula, OPPOSITE of the hypothesis")
    else:
        print("  -> ~ZERO: neither direction distinguishable from noise")
    print(flush=True)

    rows, seen = [], set()
    for i, rec in enumerate(chosen):
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": rec["guide"],
                     "target_alignment_start": rec["start"],
                     "target_alignment_end": rec["start"] + rec["length"],
                     "strand": rec["strand"], "mutation": rec["mutation"],
                     "cas_system": rec["cas_system"], "cell_type": CELL})
    print(f"deduped to {len(rows)} of {len(chosen)} rows\n", flush=True)

    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)

    print("=== scored through the real five-stage pipeline, per real stamped seed ===")
    cons_list, weighted_list = [], []
    for s in seeds:
        run_stage3(seed=s)
        r4 = run_stage4(seed=s)
        cons_list.append(r4["consistency_factor"])
        weighted_list.append(r4["total_weighted_score"])
        print(f"  seed {s}{'  (build seed)' if s == seed else ''}: "
              f"consistency_factor {r4['consistency_factor']:.4f}  "
              f"weighted {r4['total_weighted_score']:.1f}")
    print(f"\n  averaged over the 3 real seeds: consistency_factor {np.mean(cons_list):.4f}  "
          f"weighted {np.mean(weighted_list):.1f}", flush=True)


if __name__ == "__main__":
    main()
