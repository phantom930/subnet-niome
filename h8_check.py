#!/usr/bin/env python3
"""h8_check.py — why did the all-cut hedge score the floor on e9a67b4a?

It built all-cut (424.6s, clean 566/900 = 62.9%) and shipped 250 rows, yet scored consistency
0.1066 — indistinguishable from the eight all-HDR hotkeys sitting at their floor. With 62.9% clean,
all three seeds missing the clean set is a ~5% event. Two candidate explanations:

  (a) bad luck — 362/713/178 all fell outside this build's clean set
  (b) a clean seed on THIS contract is not worth what it was on b9051bc7 (0.2935)

Scores the archived submission against its own three seeds AND a uniform sample, so the round's
draw is separated from the build's expectation.
"""
import os
os.environ["NIOME_INSTANCE"] = "h8chk"

import json, random, statistics as st, sys
sys.argv = ["x"]; sys.path.insert(0, "/root/workspace/subnet-niome")
import logging; logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity
from niome_subnet.utils import settings

SRC = "data/inst/niome_hotkey8/result/2026-09-04T16:44:59"
REAL = [362, 713, 178]
N = int(os.getenv("H8_N", "60"))


def main():
    contract = json.load(open(f"{SRC}/contract.json"))
    reference = json.load(open(f"{SRC}/hbb_reference.json"))
    rows = json.load(open(f"{SRC}/submission.json"))
    cell_types = G.fetch_cell_types(); G.load_sequence()
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    print(f"h8's shipped e9a67b4a submission: {len(rows)} rows, {len(valid)} valid")
    print(f"contract seed field: {contract.get('seed')!r}\n")

    sample = sorted(random.Random(4242).sample(range(100, 1000), N))
    seeds = REAL + [s for s in sample if s not in REAL]
    out = {}
    fid = None
    for s in seeds:
        run_stage3(seed=s); r4 = run_stage4(seed=s)
        out[s] = r4["consistency_factor"]
        if fid is None:
            fid = compute_distribution_fidelity(
                valid, json.load(open(settings.STAGE3_DATASET)), contract, k=12)
            w = r4["total_weighted_score"]
    print(f"the round's own three seeds:")
    for s in REAL:
        print(f"  seed {s:>3}: consistency {out[s]:.4f}")
    print(f"  mean {st.mean(out[s] for s in REAL):.4f}   (validator reported 0.1066)")
    rest = [out[s] for s in seeds if s not in REAL]
    hi = [c for c in rest if c > 0.15]
    print(f"\nuniform sample of {len(rest)} other seeds:")
    print(f"  mean {st.mean(rest):.4f}  min {min(rest):.4f}  max {max(rest):.4f}")
    print(f"  seeds scoring > 0.15 (i.e. clean): {len(hi)}/{len(rest)} = {len(hi)/len(rest):.1%}")
    if hi:
        print(f"  their mean consistency: {st.mean(hi):.4f}")
    lo = [c for c in rest if c <= 0.15]
    if lo:
        print(f"  floor seeds mean: {st.mean(lo):.4f}")
    print(f"\nweighted {w:.1f}  fidelity {fid['distribution_fidelity_score']:.4f}")
    print(f"E[cons] over the sample: {st.mean(rest):.4f} -> E[final] "
          f"{w*st.mean(rest)*fid['distribution_fidelity_score']:.2f}")
    json.dump({"real": {str(s): out[s] for s in REAL},
               "sample": {str(s): out[s] for s in seeds if s not in REAL},
               "weighted": w, "fidelity": fid["distribution_fidelity_score"]},
              open("h8_check.json", "w"), indent=1)


if __name__ == "__main__":
    main()
