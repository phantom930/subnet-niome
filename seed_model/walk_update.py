"""Keep `seed_model/walk_record.json` current: one refit per NEW task, not a whole re-walk.

`strategy_rank.py` ranks the SeedFormer model against the five torch-free baselines, which means it
needs the model's one-step-ahead prediction for every task in the ranking window. Those cannot be
recomputed without torch, and re-running `seed_model.backtest` over the whole history hourly is not
affordable on a box whose miner is CPU-bound -- the 160-task walk takes ~11 minutes.

So the record is append-only. Tasks arrive about every 2h24m and the plan cron runs hourly, so a
normal invocation finds zero or one new task and costs one refit (~4s). The seed batch came from
`seed_model.backtest --n 160 --cold`.

**The warm-start chain is preserved across invocations.** backtest.py carries the previous step's
weights forward into the next fit, so a cold fit per new task would not be the same protocol that
produced the seeded rows. `walk_state.pt` holds that carried state; it is written after every
appended task and loaded at the start of the next run. Delete it to fall back to a cold fit, which
is honest but no longer matches the seed batch.

There is no leakage: each fit trains only on samples created strictly before the task it predicts,
and the state carried in was itself only ever fit on earlier tasks.

  python -m seed_model.walk_update                 # append whatever is new
  python -m seed_model.walk_update --refresh       # refresh seeds.json from the API first
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch

from .backtest import score_one
from .data import CELL_TYPES, CONTEXT, N_CLASSES, build_samples, class_name, load_rounds, refresh_seeds
from .evaluate import top3_from_probs
from .model import SeedFormer
from .train import fit, predict_probs

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "seed_model" / "walk_record.json"
STATE = ROOT / "seed_model" / "walk_state.pt"
CFG = dict(dim=64, depth=3, heads=4, dropout=0.2)


def load_record():
    if RECORD.exists():
        return json.loads(RECORD.read_text())
    return {"source": "seed_model.walk_update (no seed batch)", "width": 100, "tasks": []}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--start", default="2026-08-27T00:00:00")
    ap.add_argument("--end", default="now")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--horizon", type=int, default=200,
                    help="only backfill tasks inside the newest N samples; older gaps are left "
                         "alone, since the ranking window never reaches them")
    ap.add_argument("--max-new", type=int, default=8,
                    help="refuse to refit more than this in one run. A larger gap means the "
                         "record fell badly behind and should be rebuilt with backtest.py "
                         "instead of walked forward one task at a time.")
    ap.add_argument("--min-prior", type=int, default=8,
                    help="skip a task with fewer than this many strictly-earlier samples. The "
                         "oldest sample has NO prior at all, and `fit` on an empty training set "
                         "returns state=None after evaluating an empty validation batch to NaN "
                         "-- which surfaces much later as a confusing load_state_dict TypeError. "
                         "Such tasks are permanently unscoreable and sit far outside the ranking "
                         "window, so skipping them is the correct answer, not a workaround.")
    args = ap.parse_args()

    rounds, _meta = (refresh_seeds(args.start, args.end) if args.refresh else load_rounds())
    sample = build_samples(rounds, CONTEXT)
    order = np.argsort(np.asarray(sample["created"]), kind="stable")
    created = np.asarray(sample["created"])

    rec = load_record()
    have = {str(t["created_at"])[:19] for t in rec["tasks"]}
    horizon = set(order[-args.horizon:].tolist())
    todo = [i for i in order if i in horizon and sample["created"][i][:19] not in have]
    n_prior = {int(i): int((created[order] < created[i]).sum()) for i in todo}
    skipped = [i for i in todo if n_prior[int(i)] < args.min_prior]
    todo = [i for i in todo if n_prior[int(i)] >= args.min_prior]
    if skipped:
        print(f"  skipping {len(skipped)} task(s) with fewer than {args.min_prior} earlier "
              f"samples (oldest: {sample['created'][skipped[0]][:19]})")
    if not todo:
        print(f"walk_record: up to date ({len(rec['tasks'])} tasks, "
              f"newest {rec['tasks'][-1]['created_at'][:19] if rec['tasks'] else '—'})")
        return 0
    if len(todo) > args.max_new:
        print(f"walk_record: {len(todo)} tasks missing, over --max-new {args.max_new}. "
              f"Rebuild with `seed_model.backtest --n {min(len(order), 200)} --cold` and reseed "
              f"instead of walking forward.")
        return 1

    model = SeedFormer(**CFG)
    if STATE.exists():
        model.load_state_dict(torch.load(STATE, weights_only=True)["state"])
        print(f"  warm-started from {STATE.name} (the walk's carried state)")
    else:
        print("  no walk_state.pt; this batch is a cold fit and does not match the seed batch")

    for step, i in enumerate(todo, start=1):
        prior = order[created[order] < created[i]]
        best = fit(model, sample, prior, epochs=args.epochs, lr=args.lr,
                   weight_decay=0.05, batch_size=16, seed=len(have) + step, val_frac=0.2)
        probs = predict_probs(best["state"], sample, np.asarray([i]), CFG)[0]
        pred = top3_from_probs(probs)
        s = score_one(pred, sample["y"][i])
        rec["tasks"].append({
            "created_at": sample["created"][i][:19],
            "cell_type": CELL_TYPES[int(sample["cell"][i])],
            "round": int(sample["round"][i]),
            "predicted_windows": [class_name(int(c)) for c in pred],
        })
        model.load_state_dict(best["state"])                 # carry forward
        print(f"  +{step}/{len(todo)}  {sample['created'][i][:19]}  "
              f"{CELL_TYPES[int(sample['cell'][i])]:<11} "
              f"pred {'/'.join(x.split('-')[0] for x in rec['tasks'][-1]['predicted_windows'])}  "
              f"-> {s['seeds_covered']} seed(s)")

    rec["tasks"].sort(key=lambda t: t["created_at"])
    rec["updated_at"] = datetime.now(timezone.utc).isoformat()
    RECORD.write_text(json.dumps(rec, indent=2) + "\n")
    torch.save({"state": model.state_dict()}, STATE)
    print(f"  -> {RECORD} now holds {len(rec['tasks'])} tasks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
