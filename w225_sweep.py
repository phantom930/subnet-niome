#!/usr/bin/env python3
"""w225_sweep.py — is a 225-wide all-cut worth it, and how large can its group go?

clean225.py established the enabling fact: over a 225-seed window every guide in the min-union
group is STRICT (union 0, clean 225/225), and 246-277 strict guides exist per window. Strict guides
add nothing to the union, so unlike the whole-window build, group_size no longer trades against
clean seeds — the cas mix is free, and with it the fidelity term that group 80 could not buy.

What it costs is coverage: the clean set is 225 seeds, not 557, so P(a drawn seed is clean) falls
61.9% -> 25%. Seeds are therefore sampled from the FULL 100-999 range, exactly as a round draws
them; sampling inside the window would measure a case that never happens.

Baseline to beat (whole window, group 42, 120 seeds): weighted 246.6, fid 0.8938, E[cons] 0.2092,
E[final] 46.11.
"""
import os
os.environ["NIOME_INSTANCE"] = "c225"

import dataclasses, json, random, statistics as st, sys, time
sys.argv = ["x"]; sys.path.insert(0, "/root/workspace/subnet-niome")
import logging; logging.basicConfig(level=logging.ERROR)

from collections import Counter
import genExp as G
from niome_subnet.genomics import all_cut as AC
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity
from niome_subnet.utils import settings

TASK = os.getenv("W_TASK", "task-b9051bc7.json")
WINDOW = tuple(int(x) for x in os.getenv("W_WINDOW", "100-324").split("-"))
GROUPS = [int(x) for x in os.getenv("W_GROUPS", "42,80,125").split(",")]
N_SEEDS = int(os.getenv("W_SEEDS", "120"))
MF = int(os.getenv("W_MF", "25"))
RULE = os.getenv("W_RULE", "cut")
OUT = os.getenv("W_OUT", "w225_sweep.json")
BASE = {"weighted": 246.6, "fid": 0.8938, "E_cons": 0.2092, "E_final": 46.11}


def main():
    task = json.load(open(TASK))
    contract, reference = task["content"]["contract"], task["content"]["hbb_reference"]
    cell_types = G.fetch_cell_types(); G.load_sequence()
    base = AC.config_for(task["content"]["contract"]["cell_type"])
    # full-range sample: 75% of these fall outside the window, which is the point
    seeds = sorted(random.Random(12345).sample(range(100, 1000), N_SEEDS))
    inw = sum(1 for s in seeds if WINDOW[0] <= s <= WINDOW[1])
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    print(f"{task['id'][:8]}  window {WINDOW[0]}-{WINDOW[1]} (mf {MF}, rule {RULE})  groups {GROUPS}")
    print(f"{N_SEEDS} seeds sampled from the full 100-999; {inw} ({inw/N_SEEDS:.0%}) land in window")
    print(f"baseline whole-window g42: weighted {BASE['weighted']}, fid {BASE['fid']}, "
          f"E[cons] {BASE['E_cons']}, E[final] {BASE['E_final']}\n")
    print(f"{'grp':>4} {'cas12a':>7} {'clean':>6} {'union':>6} {'weighted':>9} {'casR':>6} "
          f"{'fid':>7} {'E[cons]':>8} {'E[final]':>9} {'vs base':>8} {'build':>7}")
    out = []
    for g in GROUPS:
        t0 = time.monotonic()
        cfg = dataclasses.replace(base, start_seed=WINDOW[0], end_seed=WINDOW[1],
                                  cas12a_max_fail=MF, group_size=g, rule=RULE)
        try:
            rows, meta = AC.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=1800.0)
        except Exception as exc:
            print(f"{g:>4}  FAILED: {exc}", flush=True); continue
        if not rows:
            print(f"{g:>4}  DECLINED: {meta.get('reason')}", flush=True)
            out.append({"group": g, "declined": meta.get("reason")}); continue
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
        cons, weighted, fid = [], None, None
        for s in seeds:
            run_stage3(seed=s); r4 = run_stage4(seed=s)
            cons.append(r4["consistency_factor"])
            if weighted is None:
                weighted = r4["total_weighted_score"]
                fid = compute_distribution_fidelity(
                    valid, json.load(open(settings.STAGE3_DATASET)), contract, k=12)
        ec = st.mean(cons); f = fid["distribution_fidelity_score"]
        ef = weighted * ec * f
        rec = {"group": g, "window": list(WINDOW), "max_fail": MF,
               "cas12a": sum(1 for e in valid if e["experiment"]["cas_system"] == "Cas12a"),
               "clean": meta.get("clean"), "union": meta.get("union"),
               "weighted": weighted, "fidelity": f,
               "cas_ratio": fid["cas_system_coverage_entropy_ratio"],
               "joint_ratio": fid["joint_coverage_entropy_ratio"],
               "E_cons": ec, "E_final": ef,
               "cons_in_window": st.mean([c for c, s in zip(cons, seeds)
                                          if WINDOW[0] <= s <= WINDOW[1]] or [0.0]),
               "cons_outside": st.mean([c for c, s in zip(cons, seeds)
                                        if not WINDOW[0] <= s <= WINDOW[1]] or [0.0]),
               "clean_set": [int(x) for x in (meta.get("clean_seeds") or [])][:250],
               "per_seed_cons": dict(zip(map(str, seeds), [round(c, 4) for c in cons])),
               "cells": len(Counter((e["experiment"]["mutation"], e["experiment"]["cas_system"],
                                     e["experiment"].get("strand")) for e in valid)),
               "build_s": round(time.monotonic() - t0, 1)}
        out.append(rec); json.dump(out, open(OUT, "w"), indent=1)
        print(f"{g:>4} {rec['cas12a']:>7} {str(meta.get('clean')):>6} {str(meta.get('union')):>6} "
              f"{weighted:>9.1f} {rec['cas_ratio']:>6.3f} {f:>7.4f} {ec:>8.4f} {ef:>9.2f} "
              f"{ef-BASE['E_final']:>+8.2f} {rec['build_s']:>6.1f}s", flush=True)
        print(f"      cons in-window {rec['cons_in_window']:.4f} vs outside "
              f"{rec['cons_outside']:.4f} | cells {rec['cells']}/8", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
    ok = [r for r in out if "E_final" in r]
    if ok:
        b = max(ok, key=lambda r: r["E_final"])
        print(f"\nbest at this window: group {b['group']} E[final] {b['E_final']:.2f} "
              f"({b['E_final']/BASE['E_final']-1:+.1%} vs whole-window baseline)")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
