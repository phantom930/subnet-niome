#!/usr/bin/env python3
"""sd_price_392b3538.py — what a seed-depend submission would score on 392b3538 (K562, seed 181,178,127).

This round IS stamped (seed field "181,178,127", not "0"), so seed-depend is not what the ladder
would actually reach for here -- the question is a pricing one: had we somehow known one of the
three real seeds in advance and pinned a seed-depend build to it, what would that submission's
ROUND score (averaged over all three real seeds, exactly as the validator scores it) have been, and
where would that land against the real field that was actually scored on this task?

Three candidate builds, one per real seed (181, 178, 127) -- seed-depend pins to exactly one seed,
so there is no single build that is "pinned to all three"; "best submission" here means whichever
of the three pins scores highest once averaged over the real three-seed round. Each is built with
`seed_depend.build` (pinned generation) and then scored with the VALIDATOR's own stages at all
three real seeds via `sd_task.score`, averaged the same way `validation/__init__.py` averages a
round's breakdown fields.

The real field is read via `sd_task.field`, with our own hotkeys (h0-h5, not just sd_task.OURS'
four) excluded -- h4/h5 registered after that constant was written, and both built real submissions
for this exact task earlier today.

    python sd_price_392b3538.py
"""
import os
import sys

_TASK = "392b3538-c4df-48c6-be18-009f928fc665"
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "sdprice392b")

import logging        # noqa: E402
import statistics as st  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                     # noqa: E402
from niome_subnet.genomics import seed_depend as SD    # noqa: E402
from sd_task import field, score, task_content         # noqa: E402
import sd_task as ST                                    # noqa: E402

REAL_SEEDS = (181, 178, 127)
# sd_task.OURS is missing h4/h5 (added to the fleet after that constant was written).
OURS = {
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW",   # niome_hotkey  (h0)
    "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",   # niome_hotkey1 (h1)
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU",   # niome_hotkey2 (h2)
    "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",   # niome_hotkey3 (h3)
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2",   # niome_hotkey4 (h4)
    "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",   # niome_hotkey5 (h5)
}
DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]


def round_score(rows, contract, reference, cell_types, seeds):
    parts = [score(rows, contract, reference, cell_types, seed=s) for s in seeds]
    return {
        "weighted": st.mean(p["weighted"] for p in parts),
        "consistency": st.mean(p["consistency"] for p in parts),
        "fidelity": st.mean(p["fidelity"] for p in parts),
        "final": st.mean(p["final"] for p in parts),
        "per_seed": {s: p for s, p in zip(seeds, parts)},
    }


def own_field(task_id):
    # ST.field() forces a fresh fetch through ST.fetch() when the on-disk cache is stale -- this
    # task was created and scored today, so the pre-session cache would not contain it.
    rows = ST.field(task_id, exclude_ours=False)
    best = {}
    for final, _cons, _wtd, _fid, hk in rows:
        if hk in OURS:
            continue
        if hk not in best or final > best[hk]:
            best[hk] = final
    return sorted(best.values(), reverse=True)


def main():
    task, contract, reference = task_content(_TASK)
    cell = contract.get("cell_type")
    print(f"task {_TASK[:8]}  cell {cell}  seed field {contract.get('seed')!r}  "
          f"created {task.get('created_at', '?')[:19]}\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()

    results = {}
    for pin in REAL_SEEDS:
        print(f"=== building seed-depend pinned to {pin} ===", flush=True)
        cfg = SD.SeedDependConfig()
        rows, meta = SD.build(contract, reference, cell_types, seed=pin, cfg=cfg, budget_s=900.0)
        if not rows:
            print(f"  DECLINED — {meta.get('reason')}", flush=True)
            continue
        print(f"  built {meta['rows']} rows in {meta['elapsed_s']}s", flush=True)
        rs = round_score(rows, contract, reference, cell_types, REAL_SEEDS)
        results[pin] = rs
        print(f"  round score (avg over real seeds {REAL_SEEDS}): "
              f"weighted {rs['weighted']:.1f}  consistency {rs['consistency']:.4f}  "
              f"fidelity {rs['fidelity']:.4f}  final {rs['final']:.2f}", flush=True)
        for s in REAL_SEEDS:
            p = rs["per_seed"][s]
            tag = " <-- pinned seed" if s == pin else ""
            print(f"    seed {s}: cons {p['consistency']:.4f} final {p['final']:.2f}{tag}",
                  flush=True)
        print(flush=True)

    if not results:
        print("every pin declined; nothing to compare")
        return

    best_pin = max(results, key=lambda p: results[p]["final"])
    best = results[best_pin]
    print(f"=== best of the three pins: seed {best_pin}, round final {best['final']:.2f} ===\n")

    field_scores = own_field(_TASK)
    print(f"real field for this task: {len(field_scores)} other miners scored\n")
    print("rank  final_score")
    for i, f in enumerate(field_scores[:12], start=1):
        print(f"  {i:>3}  {f:.2f}")

    rank = 1 + sum(1 for f in field_scores if f > best["final"])
    share = DIST[rank - 1] if rank <= len(DIST) else 0.0
    cutoff10 = field_scores[9] if len(field_scores) >= 10 else (min(field_scores) if field_scores else 0)
    print(f"\nour best hypothetical submission (seed-depend pinned to {best_pin}): "
          f"final {best['final']:.2f} -> rank {rank} of {len(field_scores) + 1}, "
          f"curve share {share:.0%}")
    print(f"real rank-10 cutoff on this task: {cutoff10:.2f}")


if __name__ == "__main__":
    main()
