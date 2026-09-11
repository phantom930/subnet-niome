#!/usr/bin/env python3
"""intersect_score.py — the REAL weighted x fidelity of a two-group all-HDR submission.

`allhdr_cas9.py` and `allhdr_cas12a.py` measured what each cas system reaches ALONE at width 30,
and `cas12a_only_score.py` measured what a single-cas submission is worth: fidelity 0.0287, because
stage 5's cas-coverage entropy is floored at 1e-9 when only one system is present and joint coverage
caps at 4 of 8 cells. The band was wide and the submission was worth 5.16.

This is the construction that fixes that term: min-union the two banks INDEPENDENTLY and ship both
groups, so the submission carries both cas systems and stage 5 sees a full 8 cells. Every row must
still be HDR on a band seed, so the band is the INTERSECTION of the two groups' bands -- measured at
10 seeds for all three 250-row splits (125/125, 80/170, 170/80), which is what independence predicts
(obs/exp 1.13-1.19 over the window).

Note this is NOT the shipped construction. all_hdr requires the Cas9 half to comply CONDITIONALLY
on the Cas12a group's clean band, which preserves that band entirely rather than intersecting two
independent ones -- strictly better on band size, at the cost of a Cas9 pool that decays as P(HDR)
per band seed. This asks what the unconditional two-group version is worth.

    python intersect_score.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "isect")

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
MF9 = int(os.getenv("IS_MF9", "10"))       # each arm's own plateau entry point
MF12 = int(os.getenv("IS_MF12", "12"))
SPLITS = [tuple(int(x) for x in s.split("/"))
          for s in os.getenv("IS_SPLITS", "125/125,80/170,170/80").split(",")]
TERMS = ["mutation_coverage_entropy_ratio", "cas_system_coverage_entropy_ratio",
         "strand_coverage_entropy_ratio", "joint_coverage_entropy_ratio",
         "kmer_diversity_entropy_ratio", "distinct_guide_ratio"]


def group_from(path, size, cfg, cas):
    # `load_bank` hardcodes cas_system="Cas12a" (all_cut.py:458): save_bank never stored the field,
    # because the only bank production builds is a Cas12a one. A Cas9 bank reloads mislabelled, and
    # the rows then fail stage 12 as malformed Cas12a -- 83 of 250 valid in the first run of this
    # script. The `fails` arrays are unaffected, so the band measurements stand; only the assembled
    # rows need the true label put back.
    records = load_bank(path, limit=300_000)
    for rec in records:
        rec["cas_system"] = cas
    sel = FG.FastGreedy(records, window_lo=WINDOW[0], window_hi=WINDOW[1],
                        per_cell_min=cfg.per_cell_min)
    idx, _u = sel.best(size, restarts=cfg.restarts)
    group = [records[i] for i in idx]
    bad = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    return group, sorted(set(range(WINDOW[0], WINDOW[1] + 1)) - bad)


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)
    n_rows = contract["rules"].get("max_experiments") or 250
    base = AH.config_for("K562")
    cfg9 = _dc.replace(base, hdr_range=WINDOW, main_max_fail=MF9)
    cfg12 = _dc.replace(base, hdr_range=WINDOW, main_max_fail=MF12)
    p9 = os.path.join(AH.HDR_BANK_DIR, f"cas9-{bank_key(contract, cell_types, cfg9)}.npz")
    p12 = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg12)}.npz")
    for p in (p9, p12):
        if not os.path.exists(p):
            raise SystemExit(f"bank missing: {p} — run the sweeps first")
    print(f"K562  task {(task.get('task_id') or task['id'])[:8]}  window {WINDOW[0]}-{WINDOW[1]}"
          f"  Cas9 mf {MF9} / Cas12a mf {MF12}  target rows {n_rows}\n")
    out = []
    for g9, g12 in SPLITS:
        grp9, band9 = group_from(p9, g9, cfg9, "Cas9")
        grp12, band12 = group_from(p12, g12, cfg12, "Cas12a")
        shared = sorted(set(band9) & set(band12))
        rows = AC.assemble(grp9 + grp12, [], contract, ctx, cfg12, n_rows)

        os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
        json.dump(contract, open(settings.CONTRACT_PATH, "w"))
        json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
        band = [s for s in range(100, 1000)
                if all(stage3.simulate(e, s)["outcome"] == "HDR" for e in valid)]

        hit = score(rows, contract, reference, cell_types, seed=band[0]) if band else None
        off = [s for s in range(100, 1000) if s not in set(band)]
        dirty = score(rows, contract, reference, cell_types, seed=off[0])
        summary = compute_distribution_fidelity(
            json.load(open(settings.VALID_EXPERIMENTS_PATH)),
            json.load(open(settings.STAGE3_DATASET)), contract, k=12)
        w, fid = dirty["weighted"], dirty["fidelity"]
        cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
        mix = dict(Counter(r["cas_system"] for r in rows))
        print(f"Cas9 {g9} + Cas12a {g12}")
        print(f"  rows {len(rows)}  valid {len(valid)}  cells {cells}/8  cas mix {mix}")
        print(f"  band: Cas9 {len(band9)}  Cas12a {len(band12)}  shared {len(shared)}  "
              f"-> TRUE band over 100-999: {len(band)}  {band}")
        print(f"  WEIGHTED {w:.1f}   FIDELITY {fid:.4f}   w x fid {w*fid:.2f}")
        print("  terms: " + "  ".join(
            f"{t.split('_')[0][:5]} {summary.get(t, float('nan')):.3f}" for t in TERMS))
        if hit:
            print(f"  cons on-band {hit['consistency']:.4f} -> final {hit['final']:.2f}"
                  f"   |  off-band {dirty['consistency']:.4f} -> final {dirty['final']:.2f}")
        out.append({"split": f"{g9}/{g12}", "cas9_group": g9, "cas12a_group": g12,
                    "rows": len(rows), "cells": cells, "cas_mix": mix,
                    "band_cas9": len(band9), "band_cas12a": len(band12),
                    "band_shared": len(shared), "band_true": len(band), "band_seeds": band,
                    "weighted": w, "fidelity": fid, "wxfid": w * fid,
                    "terms": {t: summary.get(t) for t in TERMS},
                    "cons_on_band": hit["consistency"] if hit else None,
                    "final_on_band": hit["final"] if hit else None,
                    "cons_off_band": dirty["consistency"], "final_off_band": dirty["final"]})
        MT.free_gpu_memory()
        print()
    json.dump(out, open("intersect_score.json", "w"), indent=1)
    print("wrote intersect_score.json")


if __name__ == "__main__":
    main()
