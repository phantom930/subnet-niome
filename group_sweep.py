#!/usr/bin/env python3
"""group_sweep.py — sweep all-cut's group_size to close the fidelity gap.

Diagnosis (gap_probe.py, round b9051bc7): our all-cut's fidelity 0.8938 is carried almost entirely
by one term, cas_system_coverage_entropy_ratio = 0.6531. The cas mix is 42 Cas12a / 208 Cas9, and
binary entropy at p = 0.168 is exactly 0.653. Lifting that ratio alone to 1.0 would give 0.9596,
past the 0.9473 that 95 miners in that round's field share.

group_size sets the Cas12a count directly, so it is the lever. It is not free and pulls three ways:

  more Cas12a  ->  cas entropy UP        (fidelity)
               ->  min-union grows       (fewer clean seeds -> consistency DOWN)
               ->  Cas12a rows score below the Cas9 rows they displace (weighted DOWN)

AllCutConfig's group_size 42 comes from a table measured on a *different bank* (it records clean
14.6%; today's mf100/d400 bank gives 62.9%), so that choice is stale and has to be re-made here.

bank_key excludes group_size, so every config below reuses one Cas12a scan — the sweep costs the
Cas9 rescan and the scoring, not five bank builds.
"""
import os
os.environ["NIOME_INSTANCE"] = "hedge"

import dataclasses
import json
import math
import statistics as st
import sys
import time
from collections import Counter

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_cut as AC
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity
from niome_subnet.utils import settings

TASK_FILE = os.getenv("GS_TASK", "task-b9051bc7.json")
GROUPS = [int(x) for x in os.getenv("GS_GROUPS", "42,60,80,100,125").split(",")]
# Seeds sampled uniformly from the window: their clean/dirty split is the true one, so the mean
# consistency over them estimates E[cons] without needing the clean set itself.
N_SEEDS = int(os.getenv("GS_SEEDS", "15"))
OUT = os.getenv("GS_OUT", "group_sweep.json")


def sample_seeds(n):
    import random
    rng = random.Random(12345)          # fixed, so every config is scored on the same seeds
    return sorted(rng.sample(range(100, 1000), n))


def main():
    task = json.load(open(TASK_FILE))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell = contract["cell_type"]
    real = [int(x) for x in str(contract["seed"]).split(",")]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AC.config_for(cell)
    seeds = sample_seeds(N_SEEDS)
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))

    print(f"{task['id'][:8]}  {cell}  round seeds {real}  |  {N_SEEDS} sampled seeds for E[cons]")
    print(f"groups {GROUPS}  (bank is shared: group_size is not in bank_key)\n")
    print(f"{'grp':>4} {'cas12a':>7} {'clean':>6} {'weighted':>9} {'casR':>6} {'jointR':>7} "
          f"{'fid':>7} {'E[cons]':>8} {'E[final]':>9} {'build':>7}")
    out = []
    for g in GROUPS:
        t0 = time.monotonic()
        cfg = dataclasses.replace(base, group_size=g)
        try:
            rows, meta = AC.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=1800.0)
        except Exception as exc:
            print(f"{g:>4}  FAILED: {exc}", flush=True)
            continue
        if not rows:
            print(f"{g:>4}  DECLINED: {meta.get('reason')}", flush=True)
            out.append({"group": g, "declined": meta.get("reason")})
            continue
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
        # stage 3/4 per sampled seed -> E[cons]; weighted and fidelity are seed-independent
        cons, weighted = [], None
        fid = None
        for s in seeds:
            run_stage3(seed=s)
            r4 = run_stage4(seed=s)
            cons.append(r4["consistency_factor"])
            if weighted is None:
                weighted = r4["total_weighted_score"]
                fid = compute_distribution_fidelity(
                    valid, json.load(open(settings.STAGE3_DATASET)), contract, k=12)
        ec = st.mean(cons)
        f = fid["distribution_fidelity_score"]
        rec = {"group": g, "cas12a": sum(1 for e in valid
                                         if e["experiment"]["cas_system"] == "Cas12a"),
               "clean_pct": meta.get("clean_pct"), "clean": meta.get("clean"),
               "weighted": weighted, "fidelity": f,
               "cas_ratio": fid["cas_system_coverage_entropy_ratio"],
               "joint_ratio": fid["joint_coverage_entropy_ratio"],
               "mutation_ratio": fid["mutation_coverage_entropy_ratio"],
               "kmer_ratio": fid["kmer_diversity_entropy_ratio"],
               "E_cons": ec, "cons_sd": st.pstdev(cons),
               "E_final": weighted * ec * f,
               "per_seed_cons": dict(zip(map(str, seeds), [round(c, 4) for c in cons])),
               "cells": len(Counter((e["experiment"]["mutation"], e["experiment"]["cas_system"],
                                     e["experiment"].get("strand")) for e in valid)),
               "build_s": round(time.monotonic() - t0, 1)}
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"{g:>4} {rec['cas12a']:>7} {str(meta.get('clean')):>6} {weighted:>9.1f} "
              f"{rec['cas_ratio']:>6.3f} {rec['joint_ratio']:>7.3f} {f:>7.4f} {ec:>8.4f} "
              f"{rec['E_final']:>9.2f} {rec['build_s']:>6.1f}s", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
    ok = [r for r in out if "E_final" in r]
    if ok:
        b = max(ok, key=lambda r: r["E_final"])
        cur = next((r for r in ok if r["group"] == 42), None)
        print(f"\nbest: group {b['group']}  E[final] {b['E_final']:.2f}"
              + (f"  vs group 42's {cur['E_final']:.2f}  ({b['E_final'] - cur['E_final']:+.2f}, "
                 f"{b['E_final'] / cur['E_final'] - 1:+.1%})" if cur else ""))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
