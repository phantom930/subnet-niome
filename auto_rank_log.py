#!/usr/bin/env python3
"""auto_rank_log.py — what auto_rank would have PICKED, task by task, and whether it hit.

A strict walk-forward: for target task t the selector is handed only the tasks strictly before t,
so the ranking that chooses the strategy and the strategy's own prediction both see history alone.
Anything else measures the selector against data it used to pick itself.

**`model` is replayed from `walk_record.json`, not from `next_prediction.json`.** `next_windows`
reads the latter, which is ONE current snapshot: asked for a task from last week it returns today's
prediction, which is lookahead. `walk_record.json` holds the model's real one-step-ahead pick per
task, which is the only honest source for a historical row. Where a task is missing from it the row
is marked and `rank_table` has already dropped `model` from that task's pool anyway.

`matched` counts SEEDS covered, not distinct windows -- the same `_seeds_covered` the ranking is
built on, so a round drawing two seeds in one predicted window scores 2. Chance is 1.00 of 3.

    ARL_SIZES=10 ARL_ROWS=20 python -u auto_rank_log.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import seed_model.data as D
import strategy_rank as SR

SIZES = tuple(int(x) for x in os.getenv("ARL_SIZES", "10").split(",") if x.strip())
N_ROWS = int(os.getenv("ARL_ROWS", "20"))


def main():
    from seed_window_model import load_tasks
    rows = sorted(load_tasks(), key=lambda t: t["at"])
    rec = SR._walk_record()
    targets = rows[-N_ROWS:]

    print(f"auto_rank walk-forward   width {D.WIDTH}   ranking window(s) "
          f"{'/'.join(str(s) for s in SIZES)}   latest {len(targets)} tasks\n")
    print(f"  {'No':>2} {'Created':<16} {'Cell':<11} {'Hist':>4} {'Picked':<12} "
          f"{'Predicted':<20} {'Actual':<20} {'Match':>5}")
    tot = 0
    hist_rows = []
    for i, t in enumerate(targets, 1):
        before = [r for r in rows if r["at"] < t["at"]]
        cell = t["cell"]
        hist = sum(1 for r in before if r["cell"] == cell)
        tbl = SR.rank_table(before, cell, SIZES)
        win = tbl.get("winner")
        if not win:
            print(f"  {i:>2} {t['at'][:16]:<16} {cell:<11} {hist:>4} {'-':<12} "
                  f"{'(no ranking)':<20} {'':<20} {'-':>5}")
            hist_rows.append(None)
            continue
        if win == "model":
            pred = rec.get(str(t["at"])[:19])
            note = "" if pred else " [no walk_record]"
        else:
            pred, _n = SR.next_windows(before, cell, win)
            note = ""
        actual_v = SR._drawn(t["seeds"])
        actual = sorted(int(c) for c in np.nonzero(actual_v)[0]
                        for _ in range(int(actual_v[c])))
        m = SR._seeds_covered(pred, actual_v / max(actual_v.sum(), 1.0)) if pred else 0
        tot += m
        hist_rows.append(m)
        ps = "/".join(D.class_name(c).split("-")[0] for c in sorted(pred)) if pred else "-"
        as_ = "/".join(D.class_name(c).split("-")[0] for c in actual)
        print(f"  {i:>2} {t['at'][:16]:<16} {cell:<11} {hist:>4} {win:<12} "
              f"{ps:<20} {as_:<20} {m:>5}{note}")
    n = len([x for x in hist_rows if x is not None])
    print(f"\n  latest {len(targets)}: {tot} seeds matched of {3*n} "
          f"= {tot/max(1,n):.3f}/task   (chance 1.000)")
    last10 = [x for x in hist_rows[-10:] if x is not None]
    if last10:
        print(f"  latest 10: {sum(last10)} of {3*len(last10)} "
              f"= {sum(last10)/len(last10):.3f}/task")
    return 0


if __name__ == "__main__":
    sys.exit(main())
