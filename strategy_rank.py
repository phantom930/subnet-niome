#!/usr/bin/env python3
"""strategy_rank.py -- pick a seed-window prediction strategy per cell type, by average rank.

The rule, as specified: for the cell type of the task about to be built, score every strategy on
that cell's latest 10, 20 and 30 tasks, rank them within each of those three windows by mean seeds
covered, and take the strategy with the best (lowest) AVERAGE rank. Worked example on HEK293:

    repeat_last  2 + 1 + 1 = 4  -> 1.33   <- selected
    model        1 + 2 + 2 = 5  -> 1.67
    cold_hand    3 + 3 + 4 = 10 -> 3.33

**Read this before trusting a selection.** The three windows are NESTED (10 subset 20 subset 30),
so they are not three independent votes -- a task in the latest 10 is counted three times and its
rank contribution is correlated across all three columns. And the quantity being ranked separates
by far less than its own standard error: over the 160-task cold walk-forward every strategy sits
between 0.93x and 1.06x chance with |z| <= 1.1, and the model's paired difference against `uniform`
is +0.006 seeds (54W/51L/55T, p 0.497). So this selector is choosing the maximum of six correlated
noisy estimates, which is a winner's-curse setup: expect it to behave like an arbitrary pick.

It is nonetheless SAFE to run, for the same reason `JOINED_SOURCE` was safe to switch before it was
proven: band position is free under a uniform generator and the generator is measured uniform (seed
min-gap median 94 observed against 93 simulated). Choosing differently costs nothing. What it must
not be read as is evidence that any strategy has an edge.

Everything here is torch-free and reads the LIVE task feed, so `.venv-ml`, the retraining cron and
`seed_refresh_guard.py` stay off the plan's critical path -- the same property `_repeat_last_windows`
was written to preserve. The one exception is the `model` strategy, whose per-task history cannot be
recomputed without torch: it is read from `seed_model/walk_record.json` and, when that file is
missing or does not cover the ranking window, `model` is dropped from the pool with a note rather
than silently scored as zero.
"""
import json
import os
from pathlib import Path

import numpy as np

import seed_model.data as D
import seed_model.evaluate as E

ROOT = Path(__file__).resolve().parent
WALK_RECORD = ROOT / "seed_model" / "walk_record.json"
NEXT_PREDICTION = ROOT / "seed_model" / "next_prediction.json"

SIZES = (10, 20, 30)
STRATEGIES = ("model", "uniform", "marginal", "hot_hand", "cold_hand", "repeat_last")


def _window_of(seed):
    return (int(seed) - 100) // D.WIDTH


def _drawn(seeds):
    v = np.zeros(D.N_CLASSES, dtype=np.float32)
    for s in seeds:
        w = _window_of(s)
        if 0 <= w < D.N_CLASSES:
            v[w] += 1.0
    return v


def _samples(rows):
    """Live task feed -> the minimal per-cell sample arrays the BASELINES actually read.

    `evaluate`'s baselines touch exactly two feature blocks -- `DRAWN` (the round's own class
    multiset, /3) and `SHARE` (that cell's cumulative class share) -- so the rest of FEATURE_DIM is
    left zero. The functions themselves are imported rather than reimplemented, because the pairing
    of a baseline with `top3_from_probs`'s stable argsort IS the strategy: a different tie-break is
    a different strategy and would not be the one that was measured.

    Returns a list in global time order of dicts: cell, at, y, x, mask.
    """
    by_cell = {}
    for r in rows:
        by_cell.setdefault(r["cell"], []).append(r)
    out = []
    for cell, mine in by_cell.items():
        mine = sorted(mine, key=lambda t: t["at"])
        feats, cum = [], np.zeros(D.N_CLASSES, dtype=np.float32)
        for r in mine:
            drawn = _drawn(r["seeds"])
            cum = cum + drawn
            f = np.zeros(D.FEATURE_DIM, dtype=np.float32)
            f[D.DRAWN] = drawn / D.SEEDS_PER_TASK
            total = cum.sum()
            f[D.SHARE] = (cum / total if total > 0
                          else np.full(D.N_CLASSES, 1.0 / D.N_CLASSES, dtype=np.float32))
            feats.append(f)
        # sample t predicts round t+1 from history[..t], exactly as build_samples does
        for t in range(len(mine) - 1):
            block, mask = D._history(feats, t, D.CONTEXT)
            nxt = mine[t + 1]
            y = _drawn(nxt["seeds"])
            out.append({"cell": cell, "at": nxt["at"], "id": nxt.get("id"),
                        "y": y / max(y.sum(), 1.0), "x": block, "mask": mask})
    out.sort(key=lambda s: s["at"])
    return out


def _seeds_covered(pred_classes, y):
    counts = np.rint(np.asarray(y) * D.SEEDS_PER_TASK).astype(int)
    return int(sum(counts[c] for c in set(int(x) for x in pred_classes)))


def _walk_record():
    """{created_at -> [class index, ...]} of the model's one-step-ahead walk-forward picks."""
    try:
        doc = json.loads(WALK_RECORD.read_text())
    except Exception:
        return {}
    out = {}
    for e in doc.get("tasks", []):
        at = str(e.get("created_at", ""))[:19]
        wins = []
        for name in e.get("predicted_windows") or []:
            try:
                wins.append((int(str(name).split("-")[0]) - 100) // D.WIDTH)
            except (ValueError, IndexError):
                wins = []
                break
        if at and len(wins) == D.SEEDS_PER_TASK:
            out[at] = wins
    return out


def per_task_scores(rows, cell, need):
    """-> (list of per-task score dicts for `cell`, newest last; note about the model)."""
    samples = _samples(rows)
    mine = [s for s in samples if s["cell"] == cell]
    want = mine[-need:] if need else mine
    if not want:
        return [], "no scored tasks for this cell type"
    idx = {id(s): i for i, s in enumerate(samples)}
    rec = _walk_record()
    model_ok, missing = True, 0
    scored = []
    for s in want:
        i = idx[id(s)]
        prior_y = np.stack([p["y"] for p in samples[:i]]) if i else np.zeros((0, D.N_CLASSES))
        row = {"at": s["at"], "y": s["y"]}
        for name in STRATEGIES[1:]:
            p = E.BASELINES[name](prior_y, s["x"], s["mask"])
            row[name] = _seeds_covered(E.top3_from_probs(p), s["y"])
        wins = rec.get(str(s["at"])[:19])
        if wins is None:
            model_ok, missing = False, missing + 1
            row["model"] = None
        else:
            row["model"] = _seeds_covered(wins, s["y"])
        scored.append(row)
    note = ("model included" if model_ok else
            f"model EXCLUDED: {missing} of {len(want)} tasks missing from walk_record.json")
    return scored, note


def rank_table(rows, cell, sizes=SIZES):
    """-> dict: per size a {strategy: (mean_seeds, rank)}, plus average rank per strategy.

    Ties share the average of the ranks they span (fractional ranking) rather than being broken by
    list order. An arbitrary tie-break is exactly the failure mode CLAUDE.md records for
    `seed_window_model.predict`, whose "rank 1" at beta 0 is a stable-sort tie-break to the lowest
    index and not a prediction at all.
    """
    scored, note = per_task_scores(rows, cell, max(sizes))
    pool = [s for s in STRATEGIES if s != "model" or all(r["model"] is not None for r in scored)]
    out = {"cell": cell, "note": note, "pool": pool, "available": len(scored), "sizes": {}}
    totals = {s: 0.0 for s in pool}
    used = []
    for k in sizes:
        if len(scored) < k:
            out["sizes"][k] = None
            continue
        win = scored[-k:]
        means = {s: float(np.mean([r[s] for r in win])) for s in pool}
        order = sorted(pool, key=lambda s: -means[s])
        ranks, i = {}, 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and means[order[j + 1]] == means[order[i]]:
                j += 1
            shared = (i + j) / 2.0 + 1.0          # average of the tied positions, 1-based
            for t in range(i, j + 1):
                ranks[order[t]] = shared
            i = j + 1
        out["sizes"][k] = {s: {"mean_seeds": means[s], "rank": ranks[s]} for s in pool}
        for s in pool:
            totals[s] += ranks[s]
        used.append(k)
    if not used:
        out["avg_rank"] = {}
        out["winner"] = None
        return out
    out["avg_rank"] = {s: totals[s] / len(used) for s in pool}
    out["used_sizes"] = used
    # tie on average rank -> the strategy that is best on the LARGEST window, then name order, so
    # the choice is at least deterministic and reproducible rather than dict-order dependent.
    big = out["sizes"][used[-1]]
    out["winner"] = min(pool, key=lambda s: (out["avg_rank"][s], -big[s]["mean_seeds"], s))
    return out


def next_windows(rows, cell, strategy):
    """The three class indices `strategy` bets on for this cell's NEXT round -> (wins, note)."""
    if strategy == "model":
        try:
            doc = json.loads(NEXT_PREDICTION.read_text())
        except Exception as exc:
            return None, f"next_prediction.json unreadable ({exc})"
        entry = next((e for e in doc.get("predictions", []) if e.get("cell_type") == cell), None)
        if not entry:
            return None, "no model entry for this cell type"
        wins = []
        for name in entry.get("predicted_classes") or []:
            try:
                wins.append((int(str(name).split("-")[0]) - 100) // D.WIDTH)
            except (ValueError, IndexError):
                return None, f"malformed class {name!r}"
        if len(set(wins)) != D.SEEDS_PER_TASK:
            return None, f"expected 3 distinct classes, got {entry.get('predicted_classes')}"
        seen = (entry.get("last_round") or {}).get("created_at") or ""
        newest = max((t["at"] for t in rows if t.get("cell") == cell), default="")
        if newest and seen and newest > seen:
            return None, f"STALE: predicts the round after {seen[:16]}, {newest[:16]} has drawn"
        return sorted(wins), f"model prediction current as of {seen[:16]}"

    samples = _samples(rows)
    mine = [s for s in samples if s["cell"] == cell]
    if not mine:
        return None, "no history for this cell type"
    # The unlabelled "what comes next" sample: history through this cell's LAST round.
    by_cell = sorted((r for r in rows if r["cell"] == cell), key=lambda t: t["at"])
    feats, cum = [], np.zeros(D.N_CLASSES, dtype=np.float32)
    for r in by_cell:
        drawn = _drawn(r["seeds"])
        cum = cum + drawn
        f = np.zeros(D.FEATURE_DIM, dtype=np.float32)
        f[D.DRAWN] = drawn / D.SEEDS_PER_TASK
        total = cum.sum()
        f[D.SHARE] = (cum / total if total > 0
                      else np.full(D.N_CLASSES, 1.0 / D.N_CLASSES, dtype=np.float32))
        feats.append(f)
    block, mask = D._history(feats, len(by_cell) - 1, D.CONTEXT)
    prior_y = np.stack([s["y"] for s in samples]) if samples else np.zeros((0, D.N_CLASSES))
    p = E.BASELINES[strategy](prior_y, block, mask)
    wins = sorted(int(c) for c in E.top3_from_probs(p))
    return wins, f"{strategy} over {len(by_cell)} rounds of this cell"


def select(rows, cell, sizes=SIZES):
    """-> (windows, note). The picker contract window_plan.py's JOINED_SOURCE dispatch expects."""
    tbl = rank_table(rows, cell, sizes)
    if not tbl.get("winner"):
        return None, f"no ranking possible ({tbl['note']}, {tbl['available']} tasks)"
    win = tbl["winner"]
    wins, note = next_windows(rows, cell, win)
    if wins is None:
        return None, f"selected {win} but it has no prediction: {note}"
    avg = tbl["avg_rank"][win]
    ranks = "/".join(f"{tbl['sizes'][k][win]['rank']:g}" for k in tbl["used_sizes"])
    return wins, (f"auto_rank picked {win} (ranks {ranks} over "
                  f"{'/'.join(str(k) for k in tbl['used_sizes'])}, avg {avg:.2f}; "
                  f"{tbl['note']})")


def _report(rows, sizes=SIZES):
    for cell in sorted({r["cell"] for r in rows}):
        tbl = rank_table(rows, cell, sizes)
        print(f"\n{'=' * 78}\n{cell}   ({tbl['available']} scored tasks; {tbl['note']})")
        hdr = "".join(f"{'n=' + str(k):>18}" for k in sizes)
        print(f"  {'strategy':<12}{hdr}{'avg rank':>11}")
        if not tbl.get("avg_rank"):
            print("  (not enough tasks)")
            continue
        for s in sorted(tbl["pool"], key=lambda s: tbl["avg_rank"][s]):
            cols = ""
            for k in sizes:
                d = tbl["sizes"].get(k)
                cols += (f"{d[s]['mean_seeds']:>11.3f} {'#' + format(d[s]['rank'], 'g'):>6}"
                         if d else f"{'—':>18}")
            mark = "  <- selected" if s == tbl["winner"] else ""
            print(f"  {s:<12}{cols}{tbl['avg_rank'][s]:>11.2f}{mark}")
        wins, note = select(rows, cell, sizes)
        print(f"  -> {[D.class_name(w) for w in wins] if wins else None}   {note}")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(ROOT))
    from seed_window_model import load_tasks
    _report(load_tasks())
