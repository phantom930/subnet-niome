#!/usr/bin/env python3
"""strategy_pool.py — the strategy rank table POOLED over all cell types.

`strategy_rank.rank_table` ranks within one cell type, because that is what the selector does: the
plan picks a scheme per cell. This asks the other question — over the fleet's latest N tasks
whatever cell they were, which strategy covered the most seeds?

**The per-task score is identical to the per-cell version and is not recomputed differently here.**
Each task is scored by the same `E.BASELINES[name](prior_y, x, mask)` -> `E.top3_from_probs` chain
`strategy_rank.per_task_scores` uses, so a strategy means exactly what it means there. Only the
RANKING WINDOW changes: the latest N tasks in global time order instead of the latest N of one cell.
Note `prior_y` is already global in that chain (it slices the pooled, time-ordered sample list), so
nothing about the prediction depends on this choice; only the set being averaged over does.

Pooling is not strictly more data for the same question — it answers a different one. A strategy
that is strong on two cells and weak on two can pool to chance, and the per-cell tables show exactly
that (`hot_hand` is #1 on HUDEP-2 and #6 on CD34+/K562). Read this as "what would one scheme applied
fleet-wide have scored", not as a tie-break on the per-cell result.

    SP_SIZES=10,20 python -u strategy_pool.py
"""
import os
import sys
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import seed_model.data as D
import seed_model.evaluate as E
import strategy_rank as SR

SIZES = tuple(int(x) for x in os.getenv("SP_SIZES", "10,20").split(",") if x.strip())


def pooled_scores(rows):
    """Every task in global time order, scored by every strategy. Newest last."""
    samples = SR._samples(rows)
    rec = SR._walk_record()
    out, missing = [], 0
    for i, s in enumerate(samples):
        prior_y = (np.stack([p["y"] for p in samples[:i]]) if i
                   else np.zeros((0, D.N_CLASSES)))
        row = {"at": s["at"], "cell": s["cell"]}
        for name in SR.STRATEGIES[1:]:
            p = E.BASELINES[name](prior_y, s["x"], s["mask"])
            row[name] = SR._seeds_covered(E.top3_from_probs(p), s["y"])
        wins = rec.get(str(s["at"])[:19])
        row["model"] = None if wins is None else SR._seeds_covered(wins, s["y"])
        if wins is None:
            missing += 1
        out.append(row)
    return out, missing


def ranked(scored, k, pool):
    """Fractional ranking by mean seeds covered over the latest k — same rule as rank_table."""
    win = scored[-k:]
    means = {s: float(np.mean([r[s] for r in win])) for s in pool}
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
    return means, ranks


def main():
    from seed_window_model import load_tasks
    rows = load_tasks()
    scored, missing = pooled_scores(rows)
    biggest = max(SIZES)
    tail = scored[-biggest:]
    model_ok = all(r["model"] is not None for r in tail)
    pool = [s for s in SR.STRATEGIES if s != "model" or model_ok]
    note = ("model included" if model_ok
            else f"model EXCLUDED: missing from walk_record.json within the latest {biggest}")

    chance = D.SEEDS_PER_TASK * (D.SEEDS_PER_TASK / D.N_CLASSES)
    sd = math.sqrt(D.SEEDS_PER_TASK * (D.SEEDS_PER_TASK / D.N_CLASSES)
                   * (1 - D.SEEDS_PER_TASK / D.N_CLASSES))
    print(f"POOLED over all cell types   width {D.WIDTH} -> {D.N_CLASSES} classes, "
          f"{D.SEEDS_PER_TASK} seeds/task")
    print(f"{len(scored)} scored tasks total; {note}")
    print(f"chance = {chance:.3f} seeds covered (per-task sd {sd:.3f})\n")

    from collections import Counter
    for k in SIZES:
        mix = Counter(r["cell"] for r in scored[-k:])
        print(f"  latest {k}: " + ", ".join(f"{c} {n}" for c, n in sorted(mix.items())))
    print()

    tables = {k: ranked(scored, k, pool) for k in SIZES if len(scored) >= k}
    avg = {s: float(np.mean([tables[k][1][s] for k in tables])) for s in pool}
    hdr = "".join(f"{'n=' + str(k):>22}" for k in SIZES)
    print(f"  {'strategy':<12}{hdr}{'avg rank':>10}{'x chance':>10}")
    for s in sorted(pool, key=lambda s: avg[s]):
        cols = ""
        for k in SIZES:
            if k not in tables:
                cols += f"{'—':>22}"
                continue
            m, r = tables[k]
            z = (m[s] - chance) / (sd / math.sqrt(k))
            cols += f"{m[s]:>9.3f} {'#' + format(r[s], 'g'):>5} {z:>+6.2f}"
        big = tables[max(tables)][0][s]
        star = "  <- best avg rank" if avg[s] == min(avg.values()) else ""
        print(f"  {s:<12}{cols}{avg[s]:>10.2f}{big / chance:>9.2f}x{star}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
