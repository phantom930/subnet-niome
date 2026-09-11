"""Score the model's predictions on the most recent N tasks, honestly.

Every prediction here is made by a model that has never seen the task it is
predicting, nor any later one: the walk steps forward one task at a time, in
task creation order, refitting on strictly-earlier rounds before each step
(warm-started from the previous step, which is also how update.py runs in
production). That is slower than scoring one model against everything, and it
is the only version of these numbers that means anything.

Four fields per task, as asked. The first three count *seeds covered*, not
distinct windows: a predicted window covers every seed that falls in it, so a
correctly predicted duplicate is worth two.

  1_window    at least 1 of the task's 3 seeds sits in a predicted window
  2_windows   at least 2 do
  3_windows   all 3 do
  duplication the task drew one window twice (e.g. 523, 201, 525 -> 500-599
              twice) and the prediction included that window

That makes the duplicate case behave the way it should: naming the duplicated
window alone is a 2-seed success, and naming it plus the task's remaining
window is a 3-seed success even though only two distinct windows were named.
`windows_matched` is still reported alongside for anyone who wants the
distinct-window count.

Each field is shown against what the same field scores for the five reference
strategies and for blind chance, computed exactly over all C(9,3) picks -
conditioned on each task's own multiset, since a duplicate changes the odds.
"""

import argparse
import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import torch

from .data import (CELL_TYPES, CONTEXT, N_CLASSES, SEEDS_PER_TASK, WIDTH,
                   build_samples, class_name, load_rounds, refresh_seeds)
from .evaluate import BASELINES, top3_from_probs
from .model import SeedFormer
from .train import fit, predict_probs

ROOT = Path(__file__).resolve().parents[1]


def score_one(pred_classes, y):
    """-> dict of the four fields for one task."""
    counts = np.rint(y * SEEDS_PER_TASK).astype(int)
    actual = {i for i, c in enumerate(counts) if c > 0}
    dup = {i for i, c in enumerate(counts) if c > 1}
    pred = set(int(c) for c in pred_classes)
    matched = pred & actual
    return {
        "windows_matched": len(matched),
        "seeds_covered": int(sum(counts[c] for c in matched)),
        "has_duplicate": bool(dup),
        "duplicate_hit": bool(dup and dup <= pred),
        "actual_windows": sorted(class_name(i) for i in actual),
        "actual_counts": {class_name(i): int(counts[i]) for i in sorted(actual)},
        "predicted_windows": [class_name(int(c)) for c in pred_classes],
    }


def aggregate(scores):
    n = len(scores)
    dups = [s for s in scores if s["has_duplicate"]]
    if n == 0:
        return {}
    ll = [s["logloss"] for s in scores if "logloss" in s]
    return {
        "n": n,
        **({"logloss": round(float(np.mean(ll)), 4),
            "uniform_logloss": round(float(np.log(N_CLASSES)), 4),
            "max_prob": round(float(np.mean([s["max_prob"] for s in scores])), 4),
            "max_deviation_from_uniform":
                round(float(np.max([s["deviation"] for s in scores])), 4)} if ll else {}),
        # thresholds on seeds covered, so a predicted duplicate counts twice
        "1_window": sum(1 for s in scores if s["seeds_covered"] >= 1),
        "2_windows": sum(1 for s in scores if s["seeds_covered"] >= 2),
        "3_windows": sum(1 for s in scores if s["seeds_covered"] >= 3),
        "duplication": {
            "tasks_with_duplicate": len(dups),
            "hit": sum(1 for s in dups if s["duplicate_hit"]),
        },
        "mean_windows_matched": round(float(np.mean(
            [s["windows_matched"] for s in scores])), 4),
        "mean_seeds_covered": round(float(np.mean(
            [s["seeds_covered"] for s in scores])), 4),
    }


def chance_rates(scores):
    """Exact probabilities for a blind pick of 3 distinct windows out of N_CLASSES.

    Enumerates all C(N_CLASSES,3) picks against each task's own multiset and
    averages, so the duplicate cases carry their true odds: naming the one
    duplicated window already covers two seeds, which makes a 2-seed success
    *more* likely on a duplicate task (1/3) than on a task with three distinct
    windows (15/84), and a 3-seed success far more likely (7/84 vs 1/84).
    """
    all_picks = list(combinations(range(N_CLASSES), SEEDS_PER_TASK))
    acc = {"1_window": [], "2_windows": [], "3_windows": []}
    for s in scores:
        counts = np.zeros(N_CLASSES, dtype=int)
        for name, c in s["actual_counts"].items():
            counts[(int(name.split("-")[0]) - 100) // WIDTH] = c
        covered = np.array([counts[list(p)].sum() for p in all_picks])
        acc["1_window"].append(float((covered >= 1).mean()))
        acc["2_windows"].append(float((covered >= 2).mean()))
        acc["3_windows"].append(float((covered >= 3).mean()))
    out = {k: (round(float(np.mean(v)), 4) if v else None)
           for k, v in acc.items() if k != "duplication"}
    # the duplicated window is one specific class, and a blind pick names 3 of
    # the N_CLASSES, so it is caught with probability 3/N regardless of the task
    out["duplication"] = round(SEEDS_PER_TASK / N_CLASSES, 4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20, help="most recent tasks to score")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--cold", action="store_true",
                    help="do not warm-start from the pooled checkpoint. Required for an honest "
                         "number whenever that checkpoint's training window covers the tasks "
                         "being scored, which it does for the live one -- see below.")
    ap.add_argument("--start", default="2026-08-27T00:00:00")
    ap.add_argument("--end", default="now")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--context", type=int, default=CONTEXT)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    tag = "" if WIDTH == 100 else f"-w{WIDTH}"
    if args.out is None:
        args.out = str(ROOT / "seed_model" / f"backtest{tag}.json")

    rounds, meta = (refresh_seeds(args.start, args.end) if args.refresh
                    else load_rounds())
    sample = build_samples(rounds, args.context)
    order = np.argsort(np.asarray(sample["created"]), kind="stable")
    if len(order) <= args.n:
        raise SystemExit(f"only {len(order)} samples; need more than --n")
    targets = order[-args.n:]

    cfg = dict(dim=64, depth=3, heads=4, dropout=0.2)
    model = SeedFormer(**cfg)
    # WARM-STARTING FROM THE LIVE CHECKPOINT LEAKS. update.py retrains pooled.pt on everything up
    # to "now", so its window covers every task this walk then predicts: step 1 begins from weights
    # that have already seen its own target. The refit is on strictly-earlier rounds, but it starts
    # from contaminated weights, so the docstring's "never seen the task it is predicting" only
    # holds with --cold. This was found by a width sweep in which width 100 was the ONLY arm with a
    # checkpoint to load and the ONLY arm to beat chance (lift 1.35x, log-loss below uniform) while
    # the seven cold arms all sat at or below it. Use --cold for any cross-width comparison.
    pooled = ROOT / "seed_model" / "checkpoints" / f"pooled{tag}.pt"
    if args.cold:
        print(f"  cold start (--cold); width {WIDTH}, {N_CLASSES} classes")
    elif pooled.exists():
        meta = torch.load(pooled, weights_only=True)
        model.load_state_dict(meta["state"])
        end = (meta.get("window") or {}).get("end", "?")
        latest = max(sample["created"][i] for i in targets)
        print(f"  warm-started from {pooled.name} (trained through {end})")
        if end == "?" or str(end) >= str(latest):
            print(f"  !! LEAKAGE: that window covers the scored tasks (latest {latest[:19]}). "
                  f"Re-run with --cold for a number that means anything.")
    else:
        print(f"  cold start (no {pooled.name}); width {WIDTH}, {N_CLASSES} classes")

    rows, model_scores = [], []
    ref_scores = {name: [] for name in BASELINES}
    created = np.asarray(sample["created"])

    for step, i in enumerate(targets, start=1):
        prior = order[created[order] < created[i]]          # strictly earlier
        best = fit(model, sample, prior, epochs=args.epochs, lr=args.lr,
                   weight_decay=0.05, batch_size=16, seed=step, val_frac=0.2)
        probs = predict_probs(best["state"], sample, np.asarray([i]), cfg)[0]
        s = score_one(top3_from_probs(probs), sample["y"][i])
        # soft log-loss against the round's own class multiset (y sums to 1, so this is the
        # cross-entropy of the true distribution under the prediction). Compared against
        # log(N_CLASSES), which is what an exactly-uniform head scores.
        s["logloss"] = float(-(sample["y"][i] * np.log(np.clip(probs, 1e-12, None))).sum())
        s["max_prob"] = float(probs.max())
        s["deviation"] = float(np.abs(probs - 1.0 / N_CLASSES).max())
        model_scores.append(s)
        for name, fn in BASELINES.items():
            p = fn(sample["y"][prior], sample["x"][i], sample["mask"][i])
            ref_scores[name].append(score_one(top3_from_probs(p), sample["y"][i]))
        rows.append({
            "created_at": sample["created"][i][:19],
            "cell_type": CELL_TYPES[int(sample["cell"][i])],
            "round": int(sample["round"][i]),
            "n_train": int(len(prior)),
            **s,
        })
        model.load_state_dict(best["state"])                # carry forward
        print(f"  {step:>2}/{args.n}  {rows[-1]['created_at']}  "
              f"{rows[-1]['cell_type']:<11} round {rows[-1]['round']:>2}  "
              f"pred {'/'.join(x.split('-')[0] for x in s['predicted_windows'])}  "
              f"actual {'/'.join(x.split('-')[0] for x in s['actual_windows'])}"
              f"{' [dup]' if s['has_duplicate'] else ''}  "
              f"-> {s['seeds_covered']} seed(s) "
              f"({s['windows_matched']} window(s))", flush=True)

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "window": meta["window"],
        "protocol": ("one-step-ahead walk-forward: each task is predicted by a "
                     "model refit on strictly earlier rounds only, warm-started "
                     "from the previous step"),
        "n_tasks": args.n,
        "model": aggregate(model_scores),
        "baselines": {k: aggregate(v) for k, v in ref_scores.items()},
        "chance": chance_rates(model_scores),
        "tasks": rows,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n")

    def line(label, agg, n):
        d = agg["duplication"]
        dup = (f"{d['hit']}/{d['tasks_with_duplicate']}"
               if d["tasks_with_duplicate"] else "—")
        return (f"  {label:<14} {agg['1_window']:>2}/{n}  {agg['2_windows']:>2}/{n}  "
                f"{agg['3_windows']:>2}/{n}   {dup:>5}   "
                f"{agg['mean_windows_matched']:.2f}  {agg['mean_seeds_covered']:.2f}")

    n = args.n
    print(f"\n{'=' * 74}")
    print(f"latest {n} tasks, one-step-ahead walk-forward")
    print("  fields count seeds covered: a predicted duplicate window is worth 2")
    print(f"  {'strategy':<14} {'≥1':>5} {'≥2':>6} {'=3':>6}   {'dup':>5}   "
          f"{'win':>4}  {'seed':>4}")
    print(line("model", payload["model"], n))
    for name in sorted(payload["baselines"],
                       key=lambda k: -payload["baselines"][k]["mean_windows_matched"]):
        print(line(name, payload["baselines"][name], n))
    c = payload["chance"]
    print(f"  {'chance':<14} {c['1_window'] * n:>4.1f} {c['2_windows'] * n:>6.1f} "
          f"{c['3_windows'] * n:>6.1f}   {c['duplication'] * 100:>4.0f}%")
    print(f"\n  -> {args.out}")


if __name__ == "__main__":
    main()
