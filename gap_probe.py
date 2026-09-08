#!/usr/bin/env python3
"""gap_probe.py — locate the structural deficit in our all-cut build.

`final = total_weighted_score x consistency_factor x distribution_fidelity_factor`. On round
b9051bc7 our all-cut scored weighted 246.6 / fidelity 0.8938 where 95 miners sharing one row
structure scored 258-262 / exactly 0.9473 — an 11% product gap that, unlike consistency, trades
against nothing. This dumps both halves per component so the fix targets the term that is actually
short rather than the one that is easiest to change.
"""
import os
os.environ["NIOME_INSTANCE"] = "hedge"

import json
import statistics as st
import sys
from collections import Counter

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity
from niome_subnet.utils import settings

TASK_FILE = "task-b9051bc7.json"
ROWS_CACHE = "hedge_rows.json"


def main():
    task = json.load(open(TASK_FILE))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    seeds = [int(x) for x in str(contract["seed"]).split(",")]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    rows = json.load(open(ROWS_CACHE))["rows"]["all-cut"]

    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    run_stage3(seed=seeds[0]); run_stage4(seed=seeds[0])
    s3 = json.load(open(settings.STAGE3_DATASET))
    fid = compute_distribution_fidelity(valid, s3, contract, k=12)

    print(f"all-cut on {task['id'][:8]}: {len(valid)} valid of {len(rows)} rows\n")

    print("=== stage 5: the six ratios (geometric mean = fidelity) ===")
    keys = ["mutation_coverage_entropy_ratio", "cas_system_coverage_entropy_ratio",
            "strand_coverage_entropy_ratio", "joint_coverage_entropy_ratio",
            "kmer_diversity_entropy_ratio", "distinct_guide_ratio"]
    for k in keys:
        v = fid.get(k)
        bar = "#" * int(round((v or 0) * 40))
        print(f"  {k:<38} {v:.4f} {bar}")
    print(f"  {'-> geometric mean':<38} {fid['distribution_fidelity_score']:.4f}")
    # what the geometric mean would be with each ratio lifted to 1.0 in turn
    print("\n  headroom per component (fidelity if that one ratio were 1.0):")
    import math
    vals = [fid[k] for k in keys]
    for i, k in enumerate(keys):
        alt = vals[:i] + [1.0] + vals[i + 1:]
        gm = math.exp(sum(math.log(max(v, 1e-12)) for v in alt) / len(alt))
        print(f"    {k:<38} {gm:.4f}  ({gm - fid['distribution_fidelity_score']:+.4f})")

    print("\n=== stage 1-2: the structural score, per row ===")
    def stat(name, xs):
        print(f"  {name:<22} mean {st.mean(xs):.4f}  median {st.median(xs):.4f}  "
              f"min {min(xs):.4f}  max {max(xs):.4f}")
    sample = valid[0]
    print(f"  fields available: {sorted(sample.keys())}")
    if "scores" in sample:
        print(f"  score fields:     {sorted(sample['scores'].keys())}")
    for field in ("gc_score", "dist_score", "offtarget_factor", "mutation_weight",
                  "weighted_score", "structural_score"):
        xs = [e.get(field, (e.get("scores") or {}).get(field)) for e in valid]
        xs = [x for x in xs if isinstance(x, (int, float))]
        if xs:
            stat(field, xs)

    print("\n=== composition ===")
    print("  by cas:      ", dict(Counter(e["experiment"]["cas_system"] for e in valid)))
    print("  by strand:   ", dict(Counter(e["experiment"].get("strand") for e in valid)))
    print("  by mutation: ", dict(Counter(e["experiment"]["mutation"] for e in valid)))
    jc = Counter((e["experiment"]["mutation"], e["experiment"]["cas_system"],
                  e["experiment"].get("strand")) for e in valid)
    print(f"  joint cells: {len(jc)} occupied, counts {sorted(jc.values())}")
    guides = [e["experiment"]["guideRNA"] for e in valid]
    print(f"  distinct guides: {len(set(guides))} of {len(guides)}")
    json.dump({"fidelity": fid, "n_valid": len(valid)}, open("gap_probe.json", "w"), indent=1)


if __name__ == "__main__":
    main()
