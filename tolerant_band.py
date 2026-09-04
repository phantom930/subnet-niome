#!/usr/bin/env python3
"""tolerant_band.py — how much wider is the band if a few rows are allowed to fail?

partial_probe measured what a near-miss is worth: 1 non-compliant row of 250 scores consistency
0.4946-0.5171 (six of six contract/seed combinations), 2-3 rows usually 0.40-0.51, and it decays to
the floor past five. That is exactly the 0.45-0.65 regime the leaders occupy and we never enter --
our band is defined by *zero* failures, which is why our consistency distribution has 1.000, the
floor, and nothing between 0.42 and 1.0.

So the question is purely how many seeds sit at 1, 2 or 3 failures. This takes an existing all-HDR
submission and counts, for every seed in 100-999, how many of its 250 rows break the ``hdr`` rule.
No search and no rebuild: the row set is fixed, and the failure count per seed is deterministic.

  band(0)   is today's band -- the seeds worth consistency 1.0
  band(<=k) adds the seeds worth ~0.50

The trade is frequency against per-seed value, and the round arithmetic decides it: a k=1 hit at
1.0 gives round consistency 0.40, while a k=1 hit at 0.50 gives only 0.23 -- below the rank-10
cutoff -- so the wider band only pays if it is *much* wider.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import json
import sys
import time
from collections import Counter

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH

TASKS = ["task-k562.json", "task-cd34.json"]
SEEDS = range(100, 1000)
OUT = "tolerant_band.json"


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = []
    for tf in TASKS:
        task = json.load(open(tf))
        contract = dict(task["content"]["contract"])
        reference = task["content"]["hbb_reference"]
        cell = contract["cell_type"]
        t0 = time.monotonic()
        rows, meta = AH.build_for_cell(dict(task["content"]["contract"]), reference, cell_types,
                                       budget_s=900.0)
        if not rows:
            print(f"{tf}: all-HDR declined ({meta.get('reason')})", flush=True)
            continue
        print(f"=== {tf}  {cell}  band {meta['band']} clean {meta['clean']} "
              f"({time.monotonic() - t0:.0f}s build) ===", flush=True)

        # Failure count per seed for the fixed row set. ctx carries the seed, so rebuild per seed
        # and simulate each row: deterministic, no search.
        counts = {}
        t0 = time.monotonic()
        for seed in SEEDS:
            c = dict(contract)
            c["seed"] = seed
            ctx = G.build_context(c, reference, cell_types)
            bad = 0
            for r in rows:
                exp = {"experiment_id": r["experiment_id"], "guideRNA": r["guideRNA"],
                       "target_alignment_start": r["target_alignment_start"],
                       "target_alignment_end": r["target_alignment_end"],
                       "strand": r["strand"], "mutation": r["mutation"],
                       "cas_system": r["cas_system"], "cell_type": r["cell_type"]}
                entry = G.build_valid_entry(exp, ctx)
                if entry is None or G.simulate(entry, ctx)["outcome"] != "HDR":
                    bad += 1
            counts[seed] = bad
            if seed % 200 == 100:
                print(f"  ... seed {seed}: {bad} failing rows "
                      f"({time.monotonic() - t0:.0f}s)", flush=True)
        hist = Counter(counts.values())
        cum = {}
        for k in (0, 1, 2, 3, 5, 8, 12, 25):
            cum[k] = sum(1 for v in counts.values() if v <= k)
        print(f"  seeds by failing-row count: "
              + ", ".join(f"{k}:{hist[k]}" for k in sorted(hist)[:12]), flush=True)
        print(f"  band(<=k) of 900: " + "  ".join(f"{k}:{cum[k]}" for k in cum), flush=True)
        print(f"  min failures on a non-band seed: "
              f"{min((v for v in counts.values() if v > 0), default=None)}", flush=True)
        out.append({"task": task["id"], "cell_type": cell, "band0": meta["clean"],
                    "cumulative": cum, "histogram": dict(hist),
                    "counts": {str(k): v for k, v in counts.items()}})
        json.dump(out, open(OUT, "w"), indent=1)
        print(flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
