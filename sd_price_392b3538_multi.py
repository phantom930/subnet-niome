#!/usr/bin/env python3
"""sd_price_392b3538_multi.py — 2-seed and 3-seed pins on 392b3538 (K562, seed 181,178,127).

Extends sd_price_392b3538.py's single-seed pricing: `seed_depend.build_multi` pins a submission to
SEVERAL known seeds at once (every row satisfies the rule at every pinned seed, not just one), so
it reaches consistency exactly 1.0 on all of them simultaneously rather than 1.0 on one and the
~0.10-0.13 floor on the rest. The single-seed results from the earlier run are reused rather than
rebuilt (135.74 / 131.24 / 133.35 for pins 181 / 178 / 127) since they are unchanged by this file.

Four new builds: the three 2-of-3 combinations, plus the one 3-of-3 (all real seeds at once, which
IS achievable here only because the round is already stamped and scored -- unreachable on a live,
unstamped round). Each is scored the same way as before: the validator's own stages, averaged over
the real three-seed round exactly as `validation/__init__.py` averages a round.

    python sd_price_392b3538_multi.py
"""
import os
import sys
from itertools import combinations

_TASK = "392b3538-c4df-48c6-be18-009f928fc665"
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "sdprice392bm")

import logging  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                     # noqa: E402
from niome_subnet.genomics import seed_depend as SD    # noqa: E402
from sd_task import task_content                        # noqa: E402
from sd_price_392b3538 import DIST, own_field, round_score  # noqa: E402

REAL_SEEDS = (181, 178, 127)
# From the earlier single-seed run -- not rebuilt here, unchanged by anything in this file.
SINGLE_SEED_FINAL = {(181,): 135.74, (178,): 131.24, (127,): 133.35}


def main():
    task, contract, reference = task_content(_TASK)
    cell = contract.get("cell_type")
    print(f"task {_TASK[:8]}  cell {cell}  seed field {contract.get('seed')!r}\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()

    plans = list(combinations(REAL_SEEDS, 2)) + [REAL_SEEDS]
    results = dict(SINGLE_SEED_FINAL)
    for pins in plans:
        print(f"=== building seed-depend pinned to ALL of {pins} ===", flush=True)
        cfg = SD.SeedDependConfig()
        rows, meta = SD.build_multi(contract, reference, cell_types, seeds=pins, cfg=cfg,
                                    budget_s=1200.0)
        if not rows:
            print(f"  DECLINED — {meta.get('reason')}  "
                  f"(candidates {meta.get('candidates')}, cells {meta.get('cells_with_candidates')}/8)",
                  flush=True)
            continue
        print(f"  built {meta['rows']} rows in {meta['elapsed_s']}s  "
              f"(candidates {meta['candidates']}, cells {meta['cells_with_candidates']}/8)",
              flush=True)
        rs = round_score(rows, contract, reference, cell_types, REAL_SEEDS)
        results[pins] = rs["final"]
        print(f"  round score: weighted {rs['weighted']:.1f}  consistency {rs['consistency']:.4f}  "
              f"fidelity {rs['fidelity']:.4f}  final {rs['final']:.2f}", flush=True)
        for s in REAL_SEEDS:
            p = rs["per_seed"][s]
            tag = " <-- pinned" if s in pins else ""
            print(f"    seed {s}: cons {p['consistency']:.4f} final {p['final']:.2f}{tag}",
                  flush=True)
        print(flush=True)

    print("=== summary: round final by pin set ===")
    for pins, final in sorted(results.items(), key=lambda kv: -kv[1]):
        label = "+".join(str(s) for s in pins)
        print(f"  {len(pins)}-seed pin {{{label}}}: final {final:.2f}")

    field_scores = own_field(_TASK)
    print(f"\nreal field for this task: {len(field_scores)} other miners scored")
    print("rank  final_score")
    for i, f in enumerate(field_scores[:12], start=1):
        print(f"  {i:>3}  {f:.2f}")
    cutoff10 = field_scores[9] if len(field_scores) >= 10 else (min(field_scores) if field_scores else 0)

    print(f"\nreal rank-10 cutoff on this task: {cutoff10:.2f}\n")
    for pins, final in sorted(results.items(), key=lambda kv: -kv[1]):
        rank = 1 + sum(1 for f in field_scores if f > final)
        share = DIST[rank - 1] if rank <= len(DIST) else 0.0
        label = "+".join(str(s) for s in pins)
        print(f"  {len(pins)}-seed {{{label}}}: final {final:.2f} -> rank {rank} of "
              f"{len(field_scores) + 1}, curve share {share:.0%}")


if __name__ == "__main__":
    main()
