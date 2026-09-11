#!/usr/bin/env python3
"""seed_width_report.py -- read the per-width backtests and tabulate them comparably.

Raw hit rates cannot be compared across widths: three 100-wide classes cover 300 of 900 seeds and
three 10-wide classes cover 30, so a coarser model looks better for free. Three columns fix that.

  lift        model mean seeds covered / blind chance's mean, at the SAME width. 1.00 is chance.
              This is the only hit-based number that means anything across the sweep.
  logloss     mean soft cross-entropy of the round's true class multiset under the prediction,
              against log(N_CLASSES) -- what an exactly-uniform head scores. BELOW means structure
              was found; at or above means the model learned the uniform distribution, which is the
              right answer to an unpredictable target.
  vs best     model minus the best reference strategy, in mean seeds covered, same width.

Uniform seeds imply lift 1.00 and logloss = log(N) at every width. A signal confined to a narrow
neighbourhood would appear as lift > 1 and logloss < log(N) at the fine widths only -- which is the
whole reason to run the sweep rather than trust the shipped nine-class result.
"""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load():
    out = {}
    for p in sorted(ROOT.glob("seed_model/backtest*.json")):
        name = p.stem
        w = 100 if name == "backtest" else int(name.split("-w")[1])
        try:
            d = json.loads(p.read_text())
        except Exception:
            continue
        if d.get("model", {}).get("n"):
            out[w] = d
    return dict(sorted(out.items(), reverse=True))


def main():
    res = load()
    if not res:
        print("no per-width backtests found; run seed_width_sweep.py first")
        return
    print("=== seed-window prediction across class widths, latest %d tasks, walk-forward ===\n"
          % max(d["n_tasks"] for d in res.values()))
    hdr = ("width", "cls", "n", ">=1", ">=2", "=3", "seeds", "chance", "lift",
           "logloss", "ln(N)", "d-ll", "best ref", "vs best", "maxdev")
    print("  " + "".join(f"{h:>9}" for h in hdr))
    for w, d in res.items():
        m, ch, bl = d["model"], d["chance"], d["baselines"]
        n = m["n"]
        seeds = m["mean_seeds_covered"]
        cha = ch.get("mean_seeds_covered") or (3.0 * 3.0 / (900 // w))
        best_name, best_val = max(((k, v["mean_seeds_covered"]) for k, v in bl.items()),
                                  key=lambda kv: kv[1])
        ll, un = m.get("logloss"), m.get("uniform_logloss") or math.log(900 // w)
        row = (str(w), str(900 // w), str(n),
               "%d/%d" % (m["1_window"], n), "%d/%d" % (m["2_windows"], n),
               "%d/%d" % (m["3_windows"], n), "%.3f" % seeds, "%.3f" % cha,
               "%.2fx" % (seeds / cha if cha else 0),
               "%.4f" % ll if ll else "-", "%.4f" % un,
               "%+.4f" % (ll - un) if ll else "-",
               best_name[:8], "%+.3f" % (seeds - best_val),
               "%.4f" % m.get("max_deviation_from_uniform", float("nan")))
        print("  " + "".join(f"{c:>9}" for c in row))
    print("\n  lift 1.00 and d-ll >= 0 at every width == the generator is uniform at every")
    print("  resolution tested, and the shipped nine-class result was not a resolution artefact.")
    print("  maxdev is the largest departure of any predicted class probability from 1/N.")


if __name__ == "__main__":
    main()
