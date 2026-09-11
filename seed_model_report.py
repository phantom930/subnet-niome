#!/usr/bin/env python3
"""seed_model_report.py -- the model's own results at each class width, tested the same way the
baselines were.

`seed_width_report.py` reports the model against same-width chance and against the best reference
strategy. This adds the piece that makes the model comparable to `baseline_cmp.py`'s output: the
uniform-generator null test, run identically -- the model's predictions are held exactly as it made
them and each task's three seeds are redrawn uniformly from 100-999, 20,000 times, scored the way
`score_one` does it (seeds covered, so a predicted duplicate is worth two).

That matters because lift over chance is not on a common scale across widths or methods: at width
10 a coverage of 0.14 against a null of 0.10 reads as 1.40x and is +0.9 sigma, because the null sd
is 0.044. z and p are the comparable columns.

Also reported per width, from the fields backtest.py records:

  d-ll    mean log-loss minus log(N_CLASSES). Below zero means the prediction carried information;
          at or above means the head learned the uniform distribution.
  maxdev  the largest departure of any predicted class probability from 1/N over the scored tasks.
          Small maxdev with d-ll >= 0 is the signature of a model that converged on uniform.

Every backtest read here must have been run with --cold. update.py retrains pooled.pt through
"now", so a warm start loads weights that have already seen the tasks being scored.
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SEEDS_PER_TASK = 3


def load():
    out = {}
    for p in sorted(ROOT.glob("seed_model/backtest*.json")):
        stem = p.stem
        if stem == "backtest":
            w = 100
        elif "-w" in stem:
            w = int(stem.split("-w")[1])
        else:
            continue
        try:
            d = json.loads(p.read_text())
        except Exception:
            continue
        if d.get("model", {}).get("n"):
            out[w] = d
    return dict(sorted(out.items(), reverse=True))


def class_index(name, width):
    return (int(name.split("-")[0]) - 100) // width


def null_test(tasks, width, n_classes, n_iter=20000, seed=0):
    """Same null as baseline_cmp.null_test: predictions fixed, targets resampled uniformly."""
    rng = np.random.default_rng(seed)
    preds = [set(class_index(x, width) for x in t["predicted_windows"]) for t in tasks]
    obs = float(np.mean([t["seeds_covered"] for t in tasks]))
    n = len(tasks)
    draws = rng.integers(100, 1000, size=(n_iter, n, SEEDS_PER_TASK))
    cls = (draws - 100) // width
    null = np.empty(n_iter, dtype=np.float64)
    for k in range(n_iter):
        tot = 0
        for j in range(n):
            pj = preds[j]
            tot += sum(1 for c in cls[k, j] if int(c) in pj)
        null[k] = tot / n
    sd = float(null.std(ddof=1))
    return {"observed": round(obs, 4), "null_mean": round(float(null.mean()), 4),
            "null_sd": round(sd, 4),
            "z": round((obs - float(null.mean())) / sd, 2) if sd else None,
            "p_above": float(((null >= obs).sum() + 1) / (n_iter + 1))}


def main():
    res = load()
    if not res:
        print("no per-width backtests found; run seed_width_sweep.py first")
        return
    n = max(d["model"]["n"] for d in res.values())
    print("=== SeedFormer at each class width, latest %d tasks, cold-start walk-forward ===" % n)
    print("  one model per width, trained and finetuned separately; each task predicted by a model")
    print("  refit on strictly-earlier rounds only\n")
    hdr = ("width", "cls", "n", ">=1", ">=2", "=3", "dup", "win", "seeds",
           "chance", "lift", "null sd", "z", "p", "d-ll", "maxdev")
    print("  " + "".join(f"{h:>9}" for h in hdr))
    zs, dlls = [], []
    for w, d in res.items():
        m, N = d["model"], d["model"]["n"]
        cls = 900 // w
        ch = d["chance"].get("mean_seeds_covered") or 9.0 / cls
        nt = null_test(d["tasks"], w, cls)
        ll = m.get("logloss")
        un = m.get("uniform_logloss") or math.log(cls)
        dll = (ll - un) if ll is not None else None
        zs.append(nt["z"])
        if dll is not None:
            dlls.append(dll)
        print("  " + "".join(f"{c:>9}" for c in (
            str(w), str(cls), str(N),
            "%d/%d" % (m["1_window"], N), "%d/%d" % (m["2_windows"], N),
            "%d/%d" % (m["3_windows"], N),
            "%d/%d" % (m["duplication"]["hit"], m["duplication"]["tasks_with_duplicate"]),
            "%.2f" % m["mean_windows_matched"], "%.3f" % m["mean_seeds_covered"],
            "%.3f" % ch, "%.2fx" % (m["mean_seeds_covered"] / ch if ch else 0),
            "%.3f" % nt["null_sd"], "%+.2f" % nt["z"], "%.3f" % nt["p_above"],
            "%+.4f" % dll if dll is not None else "-",
            "%.4f" % m.get("max_deviation_from_uniform", float("nan")))))
    print("  z across the %d widths: %s" % (len(zs), " ".join("%+.2f" % z for z in zs)))
    print("  mean z %+.2f | widths with p_above < 0.05: %d of %d"
          % (float(np.mean(zs)), sum(1 for z, (w, d) in zip(zs, res.items())
                                     if null_test(d["tasks"], w, 900 // w)["p_above"] < 0.05),
             len(zs)))
    if dlls:
        print("  log-loss minus log(N): %s" % " ".join("%+.4f" % x for x in dlls))
        print("  widths where the model beat uniform on log-loss: %d of %d"
              % (sum(1 for x in dlls if x < 0), len(dlls)))

    cmp_path = ROOT / "seed_model" / "baseline_cmp.json"
    if cmp_path.exists():
        bc = {r["width"]: r for r in json.loads(cmp_path.read_text())}
        if bc and bc[list(bc)[0]]["n"] == n:
            print("\n=== the same z, side by side with the two reference strategies ===")
            print("  %6s %5s %10s %13s %13s" % ("width", "cls", "model z",
                                                "uniform z", "repeat_last z"))
            for (w, d), z in zip(res.items(), zs):
                r = bc.get(w) or {}
                uz = (r.get("uniform_null") or {}).get("z")
                rz = (r.get("repeat_last_null") or {}).get("z")
                print("  %6d %5d %10s %13s %13s" % (
                    w, 900 // w, "%+.2f" % z,
                    "%+.2f" % uz if uz is not None else "-",
                    "%+.2f" % rz if rz is not None else "-"))
        else:
            print("\n  (baseline_cmp.json is at n=%d, not %d -- re-run baseline_cmp.py --n %d "
                  "to line the z columns up)" % (bc[list(bc)[0]]["n"], n, n))


if __name__ == "__main__":
    main()
