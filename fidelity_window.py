#!/usr/bin/env python3
"""fidelity_window.py — does the JOINED window cost distribution_fidelity_factor?

Measured on the live fleet, our all-HDR fidelity sits at 0.877-0.912 on every round and every
hotkey, while CLAUDE.md's group-size table records 0.96 for group 80 on the erythroid types. The
field's top-10 median runs 0.92-0.94. That 0.07 shortfall is systematic, not round noise, and the
obvious suspect is the window: 0.96 was measured on a CONTIGUOUS 100-seed band, and everything we
ship now comes from a joined, non-contiguous 225- or 300-seed window.

`distribution_fidelity_score` is the geometric mean of six ratios -- mutation / cas / strand / joint
coverage entropy, k-mer diversity, distinct-guide -- ALL computed from stage 12's valid experiments,
i.e. the design fields alone. `cas_shift` (repair-mode JS divergence, indel Wasserstein) is a
diagnostic and does NOT enter the score. So fidelity is seed-independent and one stage12+stage5 pass
per arm settles it; verify_seed_independence() checks that claim rather than trusting it.

Four arms separate the two confounds -- window WIDTH and non-CONTIGUITY:

    A  contiguous 100   the CELL_CONFIG default, where 0.96 was measured
    B  contiguous 225   width without joining
    C  joined 225       the live rotated slice
    D  joined 300       h0's full joined space

If B ~ A and C < B, non-contiguity is the cost. If B < A, it is width. If all four match, the 0.96
figure was optimistic and the window is exonerated.
"""
import os

os.environ.setdefault("NIOME_INSTANCE", "fidwin")

import json
import logging
import sys
import time
import urllib.request

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
logging.basicConfig(level=logging.ERROR)

from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics.validation.stage12 import run_stage12
from niome_subnet.genomics.validation.stage3 import run_stage3
from niome_subnet.genomics.validation.stage4 import run_stage4
from niome_subnet.genomics.validation.stage5 import run_stage5
from niome_subnet.utils import settings

TASK = os.environ.get("FW_TASK", "571f4843-1671-49c3-951c-47d16a9418b0")
DST = f"data/inst/{os.environ['NIOME_INSTANCE']}"

REFERENCE = None

TERMS = ["mutation_coverage_entropy_ratio", "cas_system_coverage_entropy_ratio",
         "strand_coverage_entropy_ratio", "joint_coverage_entropy_ratio",
         "kmer_diversity_entropy_ratio", "distinct_guide_ratio"]


def load_case():
    doc = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks?limit=40"))
    items = doc if isinstance(doc, list) else doc.get("data") or doc.get("items")
    task = next(t for t in items if (t.get("task_id") or t.get("id")) == TASK)
    content = task["content"]
    cells = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/data/cell-types?format=json"))
    return content["contract"], content["hbb_reference"], cells


def score(rows, contract, cells, seed):
    os.makedirs(DST, exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(REFERENCE, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cells)
    run_stage3(seed=seed)
    run_stage4(seed=seed)
    return run_stage5()


def main():
    global REFERENCE
    contract, reference, cells = load_case()
    REFERENCE = reference
    cell = contract["cell_type"]
    print(f"task {TASK[:8]}  {cell}  (fidelity is design-only; seed shown for stage 3/4 only)\n")

    arms = [
        ("A contiguous 100", dict(hdr_range=(500, 599))),
        ("B contiguous 225", dict(hdr_range=(500, 724))),
        ("C joined 225", dict(seed_list=sorted(
            set(range(100, 115)) | set(range(190, 300)) | set(range(800, 900))))),
        ("D joined 300", dict(seed_list=sorted(
            set(range(100, 300)) | set(range(800, 900))))),
    ]

    print("%-18s %-7s %-9s %-9s | %s" % ("arm", "band", "weighted", "FIDELITY",
                                         "  ".join(t.split("_")[0][:5] for t in TERMS)))
    results = {}
    for label, kw in arms:
        t0 = time.monotonic()
        rows, meta = AH.build_for_cell(contract, reference, cells, budget_s=900, **kw)
        if not rows:
            print("%-18s DECLINED (%s)" % (label, meta.get("reason"))); continue
        s5 = score(rows, contract, cells, 500)
        fid = s5["distribution_fidelity_factor"]
        d = json.load(open(settings.DISTRIBUTION_FIDELITY_PATH))
        w = json.load(open(settings.FINAL_REWARD_PATH))["total_weighted_score"]
        print("%-18s %-7s %-9.1f %-9.4f | %s   (%.0fs)" % (
            label, meta.get("clean"), w, fid,
            "  ".join("%.3f" % d[t] for t in TERMS), time.monotonic() - t0))
        results[label] = (rows, fid)

    # the claim this script rests on
    if results:
        label, (rows, fid) = next(iter(results.items()))
        alt = score(rows, contract, cells, 731)["distribution_fidelity_factor"]
        print(f"\nseed-independence check on {label}: seed 500 -> {fid:.6f}, "
              f"seed 731 -> {alt:.6f}  {'OK' if abs(fid-alt) < 1e-9 else 'DIFFERS — rerun per seed'}")
    print("\nCLAUDE.md records 0.96 for group 80 (erythroid); the field's top-10 median is 0.92-0.94.")


if __name__ == "__main__":
    main()
