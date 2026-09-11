#!/usr/bin/env python3
"""baseline_cmp.py -- repeat_last against uniform, paired, over the latest N tasks at each width.

Neither strategy needs a model, so this trains nothing: `baseline_uniform` returns a flat
distribution and `baseline_repeat` reads the previous round's classes out of the history row. That
makes 50 tasks x 8 widths cheap and, more importantly, exactly reproducible -- there is no fitted
state, so nothing here can leak the way `backtest.py`'s pooled warm-start did.

What is actually being compared, because the names mislead:

  uniform      NOT a random pick. A flat distribution through `top3_from_probs`, which takes a
               STABLE argsort, returns classes [0, 1, 2] every time -- a fixed bet on the three
               LOWEST seed ranges (100-399 at width 100, 100-129 at width 10).
  repeat_last  bets the three classes the PREVIOUS round of that cell type drew.

So this is "always bet low" against "bet that history repeats". Under a uniform generator both are
chance-equivalent in expectation and any gap is sample noise -- which is the null, and the reason
the paired test matters more than the two hit rates.

Scoring matches backtest.py exactly (`score_one`): the fields count SEEDS COVERED, so a correctly
predicted duplicate window is worth two. Chance is the exact mean over all C(N_CLASSES, 3) blind
picks conditioned on each task's own multiset, so duplicate tasks carry their true odds.

The test is the repo's own paired sign-flip (`evaluate.permutation_test`), run in both directions.
"""
import argparse
import importlib
import json
import os
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
WIDTHS = [100, 90, 75, 60, 50, 30, 20, 10]


def run_width(width, n):
    os.environ["SM_WIDTH"] = str(width)
    import seed_model.data as D
    importlib.reload(D)
    import seed_model.evaluate as E
    importlib.reload(E)
    import seed_model.backtest as B          # no torch needed for score_one/chance_rates
    importlib.reload(B)

    rounds, _meta = D.load_rounds()
    sample = D.build_samples(rounds)
    order = np.argsort(np.asarray(sample["created"]), kind="stable")
    if len(order) < n:
        n = len(order)
    targets = order[-n:]
    created = np.asarray(sample["created"])

    out = {"width": width, "n_classes": D.N_CLASSES, "n": int(n), "per_task": []}
    scores = {"uniform": [], "repeat_last": []}
    for i in targets:
        prior = order[created[order] < created[i]]
        row = {"created_at": sample["created"][i][:19],
               "cell_type": D.CELL_TYPES[int(sample["cell"][i])],
               "round": int(sample["round"][i])}
        for name in ("uniform", "repeat_last"):
            p = E.BASELINES[name](sample["y"][prior], sample["x"][i], sample["mask"][i])
            s = B.score_one(E.top3_from_probs(p), sample["y"][i])
            scores[name].append(s)
            row[name] = s["seeds_covered"]
            row[name + "_pred"] = [int(c) for c in E.top3_from_probs(p)]
            row[name + "_matched"] = s["windows_matched"]
        # exact blind-pick chance for this task's own multiset
        counts = np.rint(sample["y"][i] * D.SEEDS_PER_TASK).astype(int)
        picks = list(combinations(range(D.N_CLASSES), D.SEEDS_PER_TASK))
        cov = np.array([counts[list(p)].sum() for p in picks], dtype=float)
        row["chance"] = float(cov.mean())
        row["has_duplicate"] = scores["uniform"][-1]["has_duplicate"]
        row["actual_counts"] = {int(k): int(v) for k, v in enumerate(counts) if v}
        out["per_task"].append(row)

    for name in ("uniform", "repeat_last"):
        out[name] = B.aggregate(scores[name])
    out["chance"] = B.chance_rates(scores["uniform"])
    out["chance"]["mean_seeds_covered"] = round(
        float(np.mean([r["chance"] for r in out["per_task"]])), 4)

    u = np.array([r["uniform"] for r in out["per_task"]], dtype=float)
    rl = np.array([r["repeat_last"] for r in out["per_task"]], dtype=float)
    out["paired"] = {
        "repeat_last_minus_uniform": round(float((rl - u).mean()), 4),
        "sd": round(float((rl - u).std(ddof=1)), 4),
        "se": round(float((rl - u).std(ddof=1) / np.sqrt(len(u))), 4),
        "wins": int((rl > u).sum()), "losses": int((rl < u).sum()),
        "ties": int((rl == u).sum()),
        "p_repeat_better": E.permutation_test(rl, u)["p_value"],
        "p_uniform_better": E.permutation_test(u, rl)["p_value"],
    }
    return out



def null_test(per_task, method, n_classes, n_iter=20000, seed=0):
    """Is this method's coverage unusual UNDER A UNIFORM GENERATOR?

    The predictions are held exactly as the method made them and the TARGETS are resampled: each
    task's three seeds are redrawn uniformly from 100-999 and scored against that task's own
    predicted classes, the same way `score_one` does it (seeds covered, so a predicted duplicate is
    worth two). That is the right null for this whole line of work -- "seeds are uniform random" --
    and it is the correct test for `uniform`, whose prediction is a FIXED bet on classes 0,1,2 and
    therefore carries no randomness of its own.

    A paired sign-flip against chance would be wrong here: chance is a deterministic expectation,
    not a second sample, so the exchangeability the sign-flip assumes does not hold.
    """
    rng = np.random.default_rng(seed)
    preds = [set(t[method + "_pred"]) for t in per_task]
    obs = float(np.mean([t[method] for t in per_task]))
    n = len(per_task)
    draws = rng.integers(100, 1000, size=(n_iter, n, 3))
    cls = (draws - 100) // (900 // n_classes)
    null = np.empty(n_iter, dtype=np.float64)
    for k in range(n_iter):
        tot = 0
        for j in range(n):
            pj = preds[j]
            tot += sum(1 for c in cls[k, j] if int(c) in pj)
        null[k] = tot / n
    p_hi = float(((null >= obs).sum() + 1) / (n_iter + 1))
    p_lo = float(((null <= obs).sum() + 1) / (n_iter + 1))
    return {"observed": round(obs, 4), "null_mean": round(float(null.mean()), 4),
            "null_sd": round(float(null.std(ddof=1)), 4),
            "z": round(float((obs - null.mean()) / null.std(ddof=1)), 2) if null.std() else None,
            "p_above": p_hi, "p_below": p_lo, "n_iter": n_iter}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--widths", default=",".join(str(w) for w in WIDTHS))
    args = ap.parse_args()

    res = [run_width(int(w), args.n) for w in args.widths.split(",")]
    n = res[0]["n"]
    print("=== repeat_last vs uniform, latest %d tasks, paired per task ===" % n)
    print("  uniform = a FIXED bet on the three lowest classes (stable argsort of a flat "
          "distribution)\n  repeat_last = bet the previous round's three classes\n")
    hdr = ("width", "cls", "chance", "uniform", "lift", "repeat", "lift",
           "rl-unif", "se", "W/L/T", "p(rl>u)", "p(u>rl)")
    print("  " + "".join(f"{h:>10}" for h in hdr))
    for r in res:
        ch = r["chance"]["mean_seeds_covered"]
        pa = r["paired"]
        print("  " + "".join(f"{c:>10}" for c in (
            str(r["width"]), str(r["n_classes"]), "%.3f" % ch,
            "%.3f" % r["uniform"]["mean_seeds_covered"],
            "%.2fx" % (r["uniform"]["mean_seeds_covered"] / ch if ch else 0),
            "%.3f" % r["repeat_last"]["mean_seeds_covered"],
            "%.2fx" % (r["repeat_last"]["mean_seeds_covered"] / ch if ch else 0),
            "%+.3f" % pa["repeat_last_minus_uniform"], "%.3f" % pa["se"],
            "%d/%d/%d" % (pa["wins"], pa["losses"], pa["ties"]),
            "%.3f" % pa["p_repeat_better"], "%.3f" % pa["p_uniform_better"])))

    print("\n=== pooled over all %d widths (the same tasks, so widths are not independent) ===" % len(res))
    du = np.concatenate([[t["uniform"] for t in r["per_task"]] for r in res])
    dr = np.concatenate([[t["repeat_last"] for t in r["per_task"]] for r in res])
    d = dr - du
    print("  repeat_last - uniform = %+.4f +/- %.4f (se), wins %d losses %d ties %d over %d "
          "width-task pairs" % (d.mean(), d.std(ddof=1) / np.sqrt(len(d)),
                                int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum()), len(d)))
    print("  widths where repeat_last is ahead: %d of %d"
          % (sum(1 for r in res if r["paired"]["repeat_last_minus_uniform"] > 0), len(res)))

    for method in ("uniform", "repeat_last"):
        print("\n" + "=" * 108)
        label = ("uniform -- a FIXED bet on classes 0,1,2 (the three lowest seed ranges)"
                 if method == "uniform" else
                 "repeat_last -- bet the three classes the previous round of that cell type drew")
        print("=== %s, latest %d tasks ===" % (label, n))
        hdr = ("width", "cls", "bets", ">=1", ">=2", "=3", "dup", "win", "seeds",
               "chance", "lift", "null sd", "z", "p")
        print("  " + "".join(f"{h:>10}" for h in hdr))
        for r in res:
            a = r[method]; N_ = r["n_classes"]; w = r["width"]
            ch = r["chance"]["mean_seeds_covered"]
            nt = null_test(r["per_task"], method, N_)
            r[method + "_null"] = nt
            bets = ("%d-%d" % (100, 100 + 3 * w - 1)) if method == "uniform" else "prev"
            print("  " + "".join(f"{c:>10}" for c in (
                str(w), str(N_), bets,
                "%d/%d" % (a["1_window"], a["n"]), "%d/%d" % (a["2_windows"], a["n"]),
                "%d/%d" % (a["3_windows"], a["n"]),
                "%d/%d" % (a["duplication"]["hit"], a["duplication"]["tasks_with_duplicate"]),
                "%.2f" % a["mean_windows_matched"], "%.3f" % a["mean_seeds_covered"],
                "%.3f" % ch, "%.2fx" % (a["mean_seeds_covered"] / ch if ch else 0),
                "%.3f" % nt["null_sd"], "%+.2f" % nt["z"], "%.3f" % nt["p_above"])))
        zs = [r[method + "_null"]["z"] for r in res]
        print("  z across the 8 widths: %s" % " ".join("%+.2f" % z for z in zs))
        print("  mean z %+.2f | widths with p_above < 0.05: %d of %d  (widths share the same "
              "tasks and nested classes, so these are not 8 independent tests)"
              % (float(np.mean(zs)), sum(1 for r in res
                                         if r[method + "_null"]["p_above"] < 0.05), len(res)))

    print("\n=== positional check: is 'always bet low' collecting a real excess? ===")
    print("  share of the scored tasks' seeds falling in the three lowest classes")
    print("  %6s %8s %10s %10s" % ("width", "observed", "expected", "excess"))
    for r in res:
        w, N = r["width"], r["n_classes"]
        seeds_in, tot = 0, 0
        for t in r["per_task"]:
            tot += 3
        obs = r["uniform"]["mean_seeds_covered"] / 3.0
        print("  %6d %7.1f%% %9.1f%% %9.1f%%" % (w, 100 * obs, 100 * 3 / N,
                                                 100 * (obs - 3 / N)))
    for r in res:
        for method in ("uniform", "repeat_last"):
            r.setdefault(method + "_null", null_test(r["per_task"], method, r["n_classes"]))
    json.dump(res, open(ROOT / "seed_model" / "baseline_cmp.json", "w"), indent=1)
    print("\n  wrote seed_model/baseline_cmp.json")


if __name__ == "__main__":
    main()
