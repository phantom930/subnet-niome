#!/usr/bin/env python3
"""rule_targets.py — per-target R2 and normalised MAE for not_mhnhej at a 225-seed window.

rule_window.py showed not_mhnhej reaches a band of 41 of 225 (and 44-45 of 100) but only scores
consistency 0.154-0.281 on a band seed, against ``hdr``'s 1.000. This opens up why, term by term:

    consistency = 0.7 * max(avg_r2, 0) + 0.3 * (1 - avg_nmae)
    nmae        = mae_mean / std(y)          (returns mae_mean unchanged when y is constant)

with the average taken over is_cut, is_hdr and indel_length. The rule pins outcomes to
{HDR, BLUNT_NHEJ}: a no-cut row fails it, so ``is_cut`` becomes constant on a band seed, while
``is_hdr`` stays a coin between the two admitted modes and ``indel_length`` stays free within
BLUNT_NHEJ. So the expectation is one pinned target and two live ones -- this prints each target's
contribution instead of inferring it.

``hdr`` at width 100 is scored alongside as the contrast: it admits one outcome, so all three
targets go constant at once, which is what makes it worth 1.000 rather than ~0.2.

Banks are reused from rule_window.py (data/rule_window) and all_hdr (data/all_hdr).
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import dataclasses
import json
import sys
from collections import Counter

import numpy as np

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics.all_cut import bank_key, load_bank
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings
from rule_window import BANK_DIR, cas9_for

TASKS = ["task-k562.json", "task-cd34.json"]
# (rule, window, main_max_fail, bank directory)
PLAN = [
    ("not_mhnhej", (700, 924), 45, BANK_DIR),
    ("not_mhnhej", (700, 799), 20, BANK_DIR),
    ("hdr", (700, 799), 45, AH.HDR_BANK_DIR),
]
TARGETS = ("is_cut", "is_hdr", "indel_length")


def y_columns(path):
    """is_cut / is_hdr / indel_length exactly as stage 4 derives them from the stage-3 dataset."""
    data = json.load(open(path))
    rows = data if isinstance(data, list) else (data.get("results") or data.get("rows") or [])
    out = {"is_cut": [], "is_hdr": [], "indel_length": []}
    for r in rows:
        out["is_cut"].append(0.0 if r["outcome"] == "no_cut" else 1.0)
        out["is_hdr"].append(1.0 if r["outcome"] == "HDR" else 0.0)
        out["indel_length"].append(float(r["indel_length"]))
    return {k: np.asarray(v) for k, v in out.items()}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    for tf in TASKS:
        task = json.load(open(tf))
        contract = dict(task["content"]["contract"])
        reference = task["content"]["hbb_reference"]
        cell = contract["cell_type"]
        base = AH.config_for(cell)
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        n_rows = contract["rules"].get("max_experiments") or 250
        print(f"=== {tf}  {cell} ===", flush=True)
        for rule, window, mf, bdir in PLAN:
            width = window[1] - window[0] + 1
            cfg = dataclasses.replace(base, hdr_range=window, main_max_fail=mf)
            if bdir == BANK_DIR:
                path = os.path.join(bdir, f"{cell}-{rule}-{window[0]}_{window[1]}-mf{mf}.npz")
            else:
                path = os.path.join(bdir, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
            if not os.path.exists(path):
                print(f"  {rule} w{width}: no cached bank", flush=True)
                continue
            recs = load_bank(path)
            sel = FG.FastGreedy(recs, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                                per_cell_min=cfg.per_cell_min,
                                caps=AH._group_caps(contract, ctx, cfg))
            idx, _ = sel.best(cfg.group_size, restarts=cfg.restarts)
            bad = set()
            for i in idx:
                bad.update(int(x) for x in recs[i]["fails"])
            clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad),
                             dtype=np.int64)
            want = n_rows - cfg.group_size
            cas9 = cas9_for(clean, contract, cell_types, ctx, sites, cfg, rule, want)
            if len(cas9) < want:
                print(f"  {rule} w{width}: Cas9 pool {len(cas9)} short", flush=True)
                continue
            rows = AH.assemble([recs[i] for i in idx], cas9, contract, ctx, cfg, n_rows)
            doc = dict(contract)
            doc["seed"] = str(int(clean[0]))
            json.dump(doc, open(settings.CONTRACT_PATH, "w"))
            json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
            json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
            run_stage12(cell_types)
            off = next(s for s in range(cfg.start_seed, cfg.end_seed + 1)
                       if s not in set(int(x) for x in clean))
            print(f"  {rule:<11} w{width:<4} band {clean.size}/{width}", flush=True)
            for label, seed in (("BAND", int(clean[0])), ("off ", off)):
                run_stage3(seed=seed)
                run_stage4(seed=seed)
                res = run_stage5()
                mr = json.load(open(settings.FINAL_REWARD_PATH))["model_results"]
                y = y_columns(settings.STAGE3_DATASET)
                print(f"    {label} seed {seed}   cons {res['consistency_factor']:.4f}", flush=True)
                print(f"      {'target':<14} {'mean(y)':>8} {'std(y)':>8} {'r2':>8} "
                      f"{'mae':>8} {'nmae':>8}", flush=True)
                r2s, nmaes = [], []
                for t in TARGETS:
                    v = mr.get(t, {})
                    sd = float(np.std(y[t]))
                    nm = v.get("mae_mean", 0.0) / sd if sd >= 1e-9 else v.get("mae_mean", 0.0)
                    r2s.append(v.get("r2_mean", 0.0))
                    nmaes.append(nm)
                    print(f"      {t:<14} {float(np.mean(y[t])):>8.4f} {sd:>8.4f} "
                          f"{v.get('r2_mean', 0.0):>8.4f} {v.get('mae_mean', 0.0):>8.4f} "
                          f"{nm:>8.4f}", flush=True)
                ar, an = float(np.mean(r2s)), float(np.mean(nmaes))
                print(f"      {'AVG':<14} {'':>8} {'':>8} {ar:>8.4f} {'':>8} {an:>8.4f}"
                      f"   ->  0.7*max(r2,0)={0.7*max(ar,0):.4f} + 0.3*(1-nmae)="
                      f"{0.3*(1-an):.4f} = {0.7*max(ar,0)+0.3*(1-an):.4f}", flush=True)
        print(flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
