#!/usr/bin/env python3
"""cas12a_only_score.py — the REAL weighted x fidelity of a 250-row Cas12a-only submission.

`allhdr_cas12a.py` and `nocut_cas12a.py` measured how wide a band a Cas12a-only min-union group can
hold (13 on the `hdr` rule at width 30, 5 on `nocut`). Neither scored one: the fidelity cost was
quoted as ~0.0295, computed from stage 5's `geometric_mean` assuming balanced mutation/strand
coverage, with the cas-coverage term floored at eps=1e-9 because only one cas system is present.

This assembles the group into an actual submission -- `assemble(group, [], ...)` with a 250-guide
group leaves `need = 0`, so the rows ARE the group -- and runs the validator's own stages over it,
on a band seed and off it. What comes back is the real `total_weighted_score`, the real
`distribution_fidelity_score` with its six terms, and the real `consistency_factor`, against the
shipped mixed build on the same contract.

    python cas12a_only_score.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "c12only")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.all_cut import bank_key, load_bank   # noqa: E402
from niome_subnet.genomics.validation import stage3          # noqa: E402
from niome_subnet.genomics.validation.stage12 import run_stage12   # noqa: E402
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity  # noqa: E402
from niome_subnet.utils import settings                      # noqa: E402
from allhdr_mf import pick                                   # noqa: E402
from joined300 import fetch_tasks                            # noqa: E402
from sd_task import score                                    # noqa: E402

WINDOW = (100, 129)
GROUP = int(os.getenv("C12_GROUP", "250"))
TERMS = ["mutation_coverage_entropy_ratio", "cas_system_coverage_entropy_ratio",
         "strand_coverage_entropy_ratio", "joint_coverage_entropy_ratio",
         "kmer_diversity_entropy_ratio", "distinct_guide_ratio"]
# (label, module, config factory, mf, the stage-3 outcome test a "clean" seed must satisfy)
ARMS = [
    ("all-HDR  Cas12a-only", AH, "hdr", 12, lambda o: o == "HDR"),
    ("no-cut   Cas12a-only", AC, "nocut", 24, lambda o: o == "no_cut"),
]


def cfg_for(mod, rule, mf, base_cell="K562"):
    if mod is AH:
        base = AH.config_for(base_cell)
        return _dc.replace(base, hdr_range=WINDOW, main_max_fail=mf)
    base = AC.config_for(base_cell) or AC.AllCutConfig()
    return _dc.replace(base, start_seed=WINDOW[0], end_seed=WINDOW[1], cas12a_max_fail=mf,
                       rule="nocut", cas12a_gc=(0.00, 0.25), max_distance=3000, variants=4000)


def clean_seeds(rows, contract, reference, cell_types, test, lo=100, hi=999):
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    return [s for s in range(lo, hi + 1)
            if all(test(stage3.simulate(e, s)["outcome"]) for e in valid)], valid


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)
    n_rows = contract["rules"].get("max_experiments") or 250
    print(f"K562  task {(task.get('task_id') or task['id'])[:8]}  window {WINDOW[0]}-{WINDOW[1]}"
          f"  group {GROUP}  target rows {n_rows}\n")
    out = []
    for label, mod, rule, mf, test in ARMS:
        cfg = cfg_for(mod, rule, mf)
        bank_dir = AH.HDR_BANK_DIR if mod is AH else AC.BANK_DIR
        path = os.path.join(bank_dir, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            print(f"{label}: bank not cached at mf {mf} — run the sweep first"); continue
        records = load_bank(path, limit=300_000)
        sel = FG.FastGreedy(records, window_lo=WINDOW[0], window_hi=WINDOW[1],
                            per_cell_min=getattr(cfg, "per_cell_min", 2))
        idx, _u = sel.best(GROUP, restarts=getattr(cfg, "restarts", 12))
        group = [records[i] for i in idx]
        bad = set()
        for rec in group:
            bad.update(int(x) for x in rec["fails"])
        band_in_window = sorted(set(range(WINDOW[0], WINDOW[1] + 1)) - bad)

        rows = AC.assemble(group, [], contract, ctx, cfg, n_rows)
        band, valid = clean_seeds(rows, contract, reference, cell_types, test)
        cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
        cas_mix = dict(Counter(r["cas_system"] for r in rows))

        hit = score(rows, contract, reference, cell_types, seed=band[0]) if band else None
        off = [s for s in range(100, 1000) if s not in set(band)]
        dirty = score(rows, contract, reference, cell_types, seed=off[0])
        summary = compute_distribution_fidelity(
            json.load(open(settings.VALID_EXPERIMENTS_PATH)),
            json.load(open(settings.STAGE3_DATASET)), contract, k=12)
        w = dirty["weighted"]
        fid = dirty["fidelity"]
        print(f"{label}   mf {mf}")
        print(f"  rows {len(rows)}  valid {len(valid)}  cells {cells}/8  cas mix {cas_mix}")
        print(f"  band in window {len(band_in_window)}   band over 100-999 {len(band)}")
        print(f"  WEIGHTED {w:.1f}   FIDELITY {fid:.4f}   w x fid {w*fid:.2f}")
        print("  fidelity terms: " + "  ".join(
            f"{t.split('_')[0][:5]} {summary.get(t, float('nan')):.3f}" for t in TERMS))
        if hit:
            print(f"  consistency  on-band {hit['consistency']:.4f} -> final {hit['final']:.2f}"
                  f"   |  off-band {dirty['consistency']:.4f} -> final {dirty['final']:.2f}")
        out.append({"arm": label, "mf": mf, "group": GROUP, "rows": len(rows), "cells": cells,
                    "cas_mix": cas_mix, "band_in_window": len(band_in_window), "band": len(band),
                    "band_seeds": band, "weighted": w, "fidelity": fid, "wxfid": w * fid,
                    "terms": {t: summary.get(t) for t in TERMS},
                    "cons_on_band": hit["consistency"] if hit else None,
                    "final_on_band": hit["final"] if hit else None,
                    "cons_off_band": dirty["consistency"], "final_off_band": dirty["final"]})
        MT.free_gpu_memory()
        print()
    json.dump(out, open("cas12a_only_score.json", "w"), indent=1)
    print("wrote cas12a_only_score.json")


if __name__ == "__main__":
    main()
