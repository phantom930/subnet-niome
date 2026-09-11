#!/usr/bin/env python3
"""seed_width_sweep.py -- train and finetune the seed model at several class WIDTHS, then score
each on the latest N tasks.

The shipped model predicts three 100-seed windows and does not beat chance: walk-forward hits/3 of
0.886 against chance's 1.000, log-loss 2.228 against uniform's ln 9 = 2.197, i.e. it learned to
emit the uniform distribution (seed_model/README.md). The open question this answers is whether
that is a property of the GENERATOR or of the RESOLUTION -- a signal confined to a 20-seed
neighbourhood is invisible to 100-wide classes, because every such class contains five of them.

So each width gets its own model, trained and finetuned from scratch, and its own honest
walk-forward over the last N tasks:

    width  90 75 60 50 30 20 10   ->  classes 10 12 15 18 30 45 90

Each width runs as a separate process under SM_WIDTH, which is what data.py keys N_CLASSES,
FEATURE_DIM and the rebinning off, so the model, its checkpoints and the baselines cannot disagree
about how many classes there are. Width 100 reproduces the shipped nine-class path exactly
(seeds.json's stored fields are reused rather than recomputed) and is included as the control.

Reading the table. Raw hit rates are NOT comparable across widths -- three 100-wide classes cover
300 of 900 seeds while three 10-wide classes cover 30 -- so the columns that matter are:

  lift      model's mean seeds covered / blind chance's, at the same width. 1.00 is chance.
  logloss   against ln(N_CLASSES). Below it means the model found structure; at or above it means
            it learned the uniform distribution, which is the correct answer to noise.
  vs best   model minus the best of the five reference strategies, in mean seeds covered.

A signal that exists at any resolution shows up as lift > 1 and logloss < ln(N) at that width and
nowhere coarser. Uniform seeds predict lift 1.00 and logloss = ln(N) at EVERY width, which is the
null this is measuring against.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ML_PYTHON = ROOT / ".venv-ml" / "bin" / "python"
WIDTHS = [100, 90, 75, 60, 50, 30, 20, 10]


def run(width, n, epochs, refresh):
    tag = "" if width == 100 else f"-w{width}"
    out = ROOT / "seed_model" / f"backtest{tag}.json"
    env = dict(os.environ, SM_WIDTH=str(width),
               OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4")
    # --cold always: update.py retrains pooled.pt through "now", so warm-starting leaks the
    # scored tasks into the weights. Width 100 was the only arm with a checkpoint to load and the
    # only arm to beat chance until this was fixed.
    cmd = [str(ML_PYTHON), "-m", "seed_model.backtest", "--n", str(n),
           "--epochs", str(epochs), "--out", str(out), "--cold"]
    if refresh:
        cmd.append("--refresh")
    print(f"\n=== width {width} ({900 // width} classes) ===", flush=True)
    p = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout[-2000:]); print(p.stderr[-3000:])
        return None
    for line in p.stdout.splitlines():
        if line.strip().startswith(("1_window", "2_windows", "3_windows", "model", "warm", "cold")):
            print("   " + line.strip())
    return json.loads(out.read_text()) if out.exists() else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--widths", default=",".join(str(w) for w in WIDTHS))
    ap.add_argument("--refresh", action="store_true",
                    help="regenerate seeds.json once before the sweep")
    args = ap.parse_args()
    widths = [int(x) for x in args.widths.split(",")]

    if args.refresh:
        print("refreshing seeds.json ...", flush=True)
        env = dict(os.environ, SM_WIDTH="100")
        subprocess.run([str(ML_PYTHON), "-c",
                        "from seed_model.data import refresh_seeds;"
                        "refresh_seeds('2026-08-27T00:00:00','now')"],
                       cwd=ROOT, env=env, check=True)

    results = {}
    for w in widths:
        r = run(w, args.n, args.epochs, False)
        if r:
            results[w] = r
    json.dump(results, open(ROOT / "seed_model" / "width_sweep.json", "w"), indent=1)
    print("\n  wrote seed_model/width_sweep.json")


if __name__ == "__main__":
    main()
