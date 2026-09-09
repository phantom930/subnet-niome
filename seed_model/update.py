"""Continuous train / finetune loop.

One invocation does the whole cycle: refresh seeds.json from the API, retrain
and re-backtest, finetune a checkpoint per cell type, then write the next-round
prediction. State is kept in state.json so a run can tell what is new.

A checkpoint is only promoted when the walk-forward result does not regress
against the previous run by more than --tolerance hits per round. That is the
guard that keeps a bad day's data from quietly replacing a better model.

  update.py --end now                 # one cycle
  update.py --end now --loop 3600     # keep cycling, hourly
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ML_PYTHON = ROOT / ".venv-ml" / "bin" / "python"
CKPT_DIR = ROOT / "seed_model" / "checkpoints"
STATE = ROOT / "seed_model" / "state.json"
REPORT = ROOT / "seed_model" / "report.json"


def load_state():
    if STATE.exists():
        return json.loads(STATE.read_text())
    return {"cycles": 0, "history": []}


def run(cmd):
    print("  $ " + " ".join(str(c) for c in cmd[-8:]), flush=True)
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"step failed: {cmd[1:3]}")
    return r.stdout


def cycle(args):
    state = load_state()
    prev = state["history"][-1] if state["history"] else None
    python = str(ML_PYTHON if ML_PYTHON.exists() else sys.executable)

    backup = None
    if CKPT_DIR.exists() and any(CKPT_DIR.glob("*.pt")):
        backup = CKPT_DIR.with_name("checkpoints_prev")
        shutil.rmtree(backup, ignore_errors=True)
        shutil.copytree(CKPT_DIR, backup)

    out = run([python, "-m", "seed_model.train", "--refresh",
               "--start", args.start, "--end", args.end,
               "--epochs", str(args.epochs), "--seeds", str(args.seeds),
               "--folds", str(args.folds)]
              + (["--finetune"] if prev and not args.cold else []))
    print(out[-1400:])

    report = json.loads(REPORT.read_text())
    v = report["verdict"]
    entry = {
        "at": datetime.now(timezone.utc).isoformat(),
        "window": report["args"]["end"],
        "n_samples": sum(f["n_test"] for f in report["folds"]),
        "model_hits": v["model_mean_hits"],
        "best_reference": v["best_reference"],
        "best_reference_hits": v["reference_mean_hits"][v["best_reference"]],
        "p_value": v["vs_best_reference"]["p_value"],
        "promoted": True,
    }

    if prev and entry["model_hits"] < prev["model_hits"] - args.tolerance:
        entry["promoted"] = False
        entry["note"] = (f"regressed {prev['model_hits']:.3f} -> "
                         f"{entry['model_hits']:.3f}; kept previous checkpoints")
        if backup:
            shutil.rmtree(CKPT_DIR, ignore_errors=True)
            shutil.copytree(backup, CKPT_DIR)
        print(f"  NOT promoted: {entry['note']}")
    else:
        print(f"  promoted: {entry['model_hits']:.3f} hits/round "
              f"vs {entry['best_reference']} {entry['best_reference_hits']:.3f}")

    print(run([python, "-m", "seed_model.predict",
               "--start", args.start, "--end", args.end])[-900:])

    state["cycles"] += 1
    state["history"].append(entry)
    state["history"] = state["history"][-50:]
    STATE.write_text(json.dumps(state, indent=2) + "\n")
    return entry


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2026-08-27T00:00:00")
    ap.add_argument("--end", default="now")
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--tolerance", type=float, default=0.10,
                    help="hits/round of regression tolerated before a run's "
                         "checkpoints are rejected")
    ap.add_argument("--cold", action="store_true",
                    help="train from scratch instead of warm-starting")
    ap.add_argument("--loop", type=int, default=0,
                    help="seconds between cycles; 0 runs one cycle and exits")
    args = ap.parse_args()

    while True:
        print(f"\n=== cycle at {datetime.now(timezone.utc).isoformat()[:19]} ===")
        cycle(args)
        if not args.loop:
            return
        time.sleep(args.loop)


if __name__ == "__main__":
    main()
