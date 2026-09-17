#!/usr/bin/env python3
"""auto_rank_alt.py — the operator's selector spec against the shipped `strategy_rank.auto_rank`.

Same walk-forward footing for both: for target task t each selector is handed only the tasks
strictly before t, so the ranking that chooses the method and the method's own prediction are both
history-only. Identical scoring chain, identical fractional-rank rule, identical tie-break
(mean at the longest qualifying horizon, then name).

**Three differences between the spec and what ships, and only one of them can bind here:**

  1. CANDIDATE POOL. The spec ranks 5 candidates -- the baselines. Shipped ranks 6, adding `model`
     (its walk-forward picks from walk_record.json). Removing a candidate shifts the ranks of the
     others, and `model` does not sit at the same rank in every horizon, so this can reorder the
     average even on tasks where `model` would never have been picked. **This is the live one.**
  2. PARTIAL HORIZONS. Shipped drops horizon k unless it can be filled to k exactly. The spec drops
     it only if fewer than MIN_HISTORY tasks exist, otherwise truncating. Inert wherever every cell
     already holds >= 30 prior tasks.
  3. UNIFORM FALLBACK. The spec falls back to `uniform` and counts it; shipped returns no winner and
     `window_plan` drops to its `rank_freq` fallback. Only reachable when 2 leaves no horizon.

2 and 3 are reported as counts so the claim that they are inert is checked rather than assumed.

    ARA_ROWS=20 ARA_MIN_HISTORY=10 python -u auto_rank_alt.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import seed_model.data as D
import strategy_rank as SR

N_ROWS = int(os.getenv("ARA_ROWS", "20"))
MIN_HISTORY = int(os.getenv("ARA_MIN_HISTORY", "10"))
HORIZONS = (10, 20, 30)
CANDIDATES = tuple(s for s in SR.STRATEGIES if s != "model")     # the 5 baselines


def _frac_rank(means, pool):
    order = sorted(pool, key=lambda s: -means[s])
    ranks, i = {}, 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and means[order[j + 1]] == means[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            ranks[order[t]] = shared
        i = j + 1
    return ranks


def alt_select(before, cell):
    """The spec: 5 baselines, horizons truncated to MIN_HISTORY, uniform fallback."""
    scored, _note = SR.per_task_scores(before, cell, max(HORIZONS))
    avail = len(scored)
    used, totals, last_means = [], {s: 0.0 for s in CANDIDATES}, None
    partial = False
    for k in HORIZONS:
        take = min(k, avail)
        if take < MIN_HISTORY:
            continue
        if take < k:
            partial = True
        win = scored[-take:]
        means = {s: float(np.mean([r[s] for r in win])) for s in CANDIDATES}
        ranks = _frac_rank(means, CANDIDATES)
        for s in CANDIDATES:
            totals[s] += ranks[s]
        used.append(k)
        last_means = means
    if not used:
        return "uniform", True, partial, avail
    avg = {s: totals[s] / len(used) for s in CANDIDATES}
    winner = min(CANDIDATES, key=lambda s: (avg[s], -last_means[s], s))
    return winner, False, partial, avail


def main():
    from seed_window_model import load_tasks
    rows = sorted(load_tasks(), key=lambda t: t["at"])
    rec = SR._walk_record()
    targets = rows[-N_ROWS:]

    def predict(before, cell, method, at):
        if method == "model":
            return rec.get(str(at)[:19])
        w, _n = SR.next_windows(before, cell, method)
        return w

    print(f"width {D.WIDTH}   horizons {HORIZONS}   MIN_HISTORY {MIN_HISTORY}   "
          f"latest {len(targets)} tasks")
    print(f"  shipped pool: {len(SR.STRATEGIES)} ({', '.join(SR.STRATEGIES)})")
    print(f"  spec    pool: {len(CANDIDATES)} ({', '.join(CANDIDATES)})\n")
    print(f"  {'No':>2} {'Created':<16} {'Cell':<11} {'Hist':>4} | "
          f"{'SHIPPED':<12} {'pred':<13} {'M':>2} | {'SPEC':<12} {'pred':<13} {'M':>2} | diff")
    a_tot = b_tot = 0
    n_diff = n_fb = n_part = 0
    for i, t in enumerate(targets, 1):
        before = [r for r in rows if r["at"] < t["at"]]
        cell = t["cell"]
        hist = sum(1 for r in before if r["cell"] == cell)
        y = SR._drawn(t["seeds"])
        ynorm = y / max(y.sum(), 1.0)

        tbl = SR.rank_table(before, cell, HORIZONS)
        a_win = tbl.get("winner")
        a_pred = predict(before, cell, a_win, t["at"]) if a_win else None
        a_m = SR._seeds_covered(a_pred, ynorm) if a_pred else 0

        b_win, fb, part, avail = alt_select(before, cell)
        b_pred = predict(before, cell, b_win, t["at"])
        b_m = SR._seeds_covered(b_pred, ynorm) if b_pred else 0

        a_tot += a_m; b_tot += b_m
        n_fb += fb; n_part += part
        same = (a_win == b_win)
        n_diff += (not same)
        fmt = lambda p: "/".join(D.class_name(c).split("-")[0] for c in sorted(p)) if p else "-"
        print(f"  {i:>2} {t['at'][:16]:<16} {cell:<11} {hist:>4} | "
              f"{str(a_win):<12} {fmt(a_pred):<13} {a_m:>2} | "
              f"{b_win:<12} {fmt(b_pred):<13} {b_m:>2} | "
              f"{'' if same else 'PICK'}{' FB' if fb else ''}{' PART' if part else ''}")
    n = len(targets)
    print(f"\n  SHIPPED (6 candidates): {a_tot} of {3*n} = {a_tot/n:.3f}/task")
    print(f"  SPEC    (5 candidates): {b_tot} of {3*n} = {b_tot/n:.3f}/task")
    print(f"  different pick on {n_diff} of {n} tasks; "
          f"partial horizons used on {n_part}; uniform fallbacks {n_fb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
