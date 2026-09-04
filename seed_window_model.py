#!/usr/bin/env python3
"""seed_window_model.py — predict which 100-seed windows the next task will draw, and log it live.

Measured edge this is built on: naming the least-frequently-drawn window *for that cell type* hits
one of the round's three seeds 51.8% of the time in-sample and **50.0% on a held-out second half**
(18/36, z = 2.66, p = 0.004), against 29.8% for an i.i.d.-uniform generator. The effect is
frequency-balancing rather than round-robin -- a window repeats before the cycle completes on 22.6%
of draws, and a "longest unseen" rule scores only 42.9% against this rule's 51.8%.

**Not a neural net, deliberately.** There are 71 tasks and 9 classes. A network would fit noise; the
whole signal is carried by one parameter, so this fits exactly that:

    p_w  = softmax(beta * z_w)        z_w = standardised deficit of window w in this cell's history
    P(w drawn among the 3 seeds) = 1 - (1 - p_w)^3

``beta = 0`` reproduces the uniform generator, ``beta > 0`` means under-drawn windows are favoured.
It is fitted by maximum likelihood, walk-forward, on that cell type's own history, and it updates
by recounting -- there is nothing to retrain or fine-tune.

**Shadow mode is the point.** Two accuracy numbers are reported and they are not equivalent:

* *backtest* -- walk-forward over history. No look-ahead (each prediction sees only earlier tasks),
  but the rule was chosen after seeing this data, so it replicates the known result rather than
  testing it.
* *live* -- predictions this script emitted **before** a task existed, resolved once its seeds are
  stamped. This is the only number that can justify reallocating hotkeys.

Break-even for the concentration strategy is ~38% top-1 accuracy; below that, spreading nine
hotkeys as today is better. Do not reallocate on the backtest figure.

Usage:
    python seed_window_model.py            # resolve pending predictions, report, emit new ones
    python seed_window_model.py --no-emit  # report only, do not record a new live prediction
"""
import json
import math
import os
import statistics as st
import sys
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone

TASKS_URL = "https://niome-api.genomes.io/api/v3/tasks"
SINCE = "2026-08-27T00:00:00"
NW = 9                     # windows 100-199 .. 900-999
WARMUP = 15                # tasks before the model is allowed to predict
BETA_GRID = [i / 20 for i in range(-20, 81)]     # -1.0 .. 4.0
LOG = "seed_window_log.json"
CHANCE = 1 - (8 / 9) ** 3  # 0.2977: any one window among three uniform draws


def window(seed):
    return min(NW - 1, max(0, seed // 100 - 1))


def window_label(w):
    return f"{(w + 1) * 100}-{(w + 1) * 100 + 99}"


def load_tasks():
    raw = json.load(urllib.request.urlopen(TASKS_URL, timeout=120))
    raw = raw if isinstance(raw, list) else (raw.get("data") or raw.get("items") or [])
    out = []
    for t in raw:
        c = (t.get("content") or {}).get("contract") or {}
        parts = [x.strip() for x in str(c.get("seed") or "").split(",") if x.strip()]
        at = t.get("created_at") or ""
        if len(parts) == 3 and c.get("cell_type") and at >= SINCE:
            out.append({"id": t["id"], "at": at, "cell": c["cell_type"],
                        "seeds": sorted(int(x) for x in parts)})
    out.sort(key=lambda r: r["at"])
    return out


def counts_for(history, cell):
    c = Counter()
    for t in history:
        if t["cell"] == cell:
            for s in t["seeds"]:
                c[window(s)] += 1
    return [c[w] for w in range(NW)]


def probabilities(counts, beta):
    """Per-draw window probabilities from the standardised deficit."""
    n = sum(counts)
    if n == 0:
        return [1 / NW] * NW
    exp = n / NW
    scale = max(1.0, math.sqrt(exp))
    z = [(exp - c) / scale for c in counts]
    m = max(beta * x for x in z)
    e = [math.exp(beta * x - m) for x in z]
    tot = sum(e)
    return [x / tot for x in e]


def fit_beta(history, cell):
    """Maximum-likelihood beta over that cell type's own walk-forward history."""
    obs = []
    seen = []
    for t in history:
        if t["cell"] != cell:
            continue
        if len(seen) >= 3:                       # need some history to form a deficit
            obs.append((counts_for(seen, cell), [window(s) for s in t["seeds"]]))
        seen.append(t)
    if len(obs) < 5:
        return 0.0
    best, best_ll = 0.0, -1e18
    for beta in BETA_GRID:
        ll = 0.0
        for counts, draws in obs:
            p = probabilities(counts, beta)
            for w in draws:
                ll += math.log(max(p[w], 1e-12))
        if ll > best_ll:
            best, best_ll = beta, ll
    return best


def heuristic(history, cell):
    """The hand-picked rule: least-frequent window for this cell type, ties to the lowest index.

    Kept alongside the fitted model because the two disagree and the disagreement is the whole
    question. It scored 51.8% in-sample and 50.0% held-out, but its edge shrinks to 46.4% when ties
    are broken at random and ties decide 29 of 56 rounds -- while the likelihood fit below puts
    beta at 0, i.e. no balancing at all. Only the live log can separate them.
    """
    counts = counts_for(history, cell)
    return sorted(range(NW), key=lambda w: (counts[w], w))


def predict(history, cell):
    """Ranked windows for the next task of this cell type."""
    counts = counts_for(history, cell)
    beta = fit_beta(history, cell)
    p = probabilities(counts, beta)
    hit = [1 - (1 - x) ** 3 for x in p]           # P(window appears among the three seeds)
    order = sorted(range(NW), key=lambda w: -hit[w])
    return {"beta": round(beta, 3), "counts": counts,
            "ranked": [{"window": w, "label": window_label(w), "p_hit": round(hit[w], 4)}
                       for w in order]}


def backtest(rows):
    """Walk-forward: every prediction sees only strictly earlier tasks."""
    res = defaultdict(lambda: {"n": 0, "top1": 0, "top2": 0, "heur": 0})
    betas = []
    for i, t in enumerate(rows):
        if i < WARMUP:
            continue
        pr = predict(rows[:i], t["cell"])
        if not pr["ranked"]:
            continue
        actual = {window(s) for s in t["seeds"]}
        top = [r["window"] for r in pr["ranked"]]
        hr = heuristic(rows[:i], t["cell"])
        betas.append(pr["beta"])
        for key in ("ALL", t["cell"]):
            d = res[key]
            d["n"] += 1
            d["top1"] += 1 if top[0] in actual else 0
            d["top2"] += 1 if set(top[:2]) & actual else 0
            d["heur"] += 1 if hr[0] in actual else 0
    return res, betas


def wilson(hits, n, z=1.96):
    if not n:
        return (0.0, 0.0)
    p = hits / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    rows = load_tasks()
    log = json.load(open(LOG)) if os.path.exists(LOG) else {"live": []}
    by_id = {t["id"]: t for t in rows}
    now = datetime.now(timezone.utc).isoformat()

    # --- resolve any pending live predictions against tasks that have appeared since ---
    resolved_now = 0
    for entry in log["live"]:
        if entry.get("resolved"):
            continue
        later = [t for t in rows if t["cell"] == entry["cell"] and t["at"] > entry["made_at"]
                 and t["id"] not in entry.get("known_ids", [])]
        if not later:
            continue
        t = later[0]
        actual = sorted({window(s) for s in t["seeds"]})
        top = [r["window"] for r in entry["ranked"]]
        entry["resolved"] = {"task_id": t["id"], "at": t["at"], "seeds": t["seeds"],
                             "actual_windows": actual,
                             "top1_hit": top[0] in actual,
                             "top2_hit": bool(set(top[:2]) & set(actual))}
        resolved_now += 1

    print(f"{len(rows)} three-seed tasks, {rows[0]['at'][:16]} .. {rows[-1]['at'][:16]}")
    if resolved_now:
        print(f"resolved {resolved_now} pending live prediction(s) this run")

    res, betas = backtest(rows)
    a = res["ALL"]
    lo, hi = wilson(a["top1"], a["n"])
    print(f"\n=== BACKTEST (walk-forward; replicates the known result, does NOT test it) ===")
    hl, hh = wilson(a["heur"], a["n"])
    print(f"  {'scope':<12} {'n':>4} {'model':>8} {'95% CI':>16} {'heuristic':>10} "
          f"{'95% CI':>16}   chance {CHANCE:.1%}")
    print(f"  {'ALL':<12} {a['n']:>4} {a['top1']/a['n']:>7.1%} "
          f"{f'[{lo:.1%}, {hi:.1%}]':>16} {a['heur']/a['n']:>9.1%} "
          f"{f'[{hl:.1%}, {hh:.1%}]':>16}")
    for cell in sorted(k for k in res if k != "ALL"):
        d = res[cell]
        print(f"  {cell:<12} {d['n']:>4} {d['top1']/d['n']:>7.1%} {'':>16} "
              f"{d['heur']/d['n']:>9.1%}")
    if betas:
        print(f"  fitted beta: median {st.median(betas):.2f} "
              f"(0 = uniform generator, >0 = under-drawn windows favoured)")

    # --- live scoreboard: the only number that justifies reallocating ---
    done = [e for e in log["live"] if e.get("resolved")]
    print(f"\n=== LIVE (predictions emitted before the task existed) ===")
    if not done:
        print("  no resolved live predictions yet — this is expected until the next task lands.")
        print("  Do not reallocate hotkeys on the backtest figure; break-even is ~38%.")
    else:
        h1 = sum(1 for e in done if e["resolved"]["top1_hit"])
        hh_ = sum(1 for e in done if e["resolved"].get("heuristic_hit"))
        lo, hi = wilson(h1, len(done))
        print(f"  resolved {len(done)}  model {h1}/{len(done)} = {h1/len(done):.1%}  "
              f"95% CI [{lo:.1%}, {hi:.1%}]  |  heuristic {hh_}/{len(done)} = "
              f"{hh_/len(done):.1%}  |  chance {CHANCE:.1%}")
        verdict = ("BELOW break-even (~38%) — keep the current spread" if hi < 0.38
                   else "above chance and break-even" if lo > 0.38
                   else "inconclusive — keep logging")
        print(f"  verdict: {verdict}")
        for e in done[-6:]:
            r = e["resolved"]
            print(f"    {r['at'][:16]}  {e['cell']:<11} predicted "
                  f"{e['ranked'][0]['label']:<9} actual "
                  f"{','.join(window_label(w) for w in r['actual_windows']):<30} "
                  f"{'HIT' if r['top1_hit'] else 'miss'}")

    # --- current prediction for the next task of each cell type ---
    print(f"\n=== PREDICTION for the next task (history through {rows[-1]['at'][:16]}) ===")
    print(f"  {'cell type':<12} {'beta':>5}  top-3 windows by P(drawn)")
    emit = "--no-emit" not in sys.argv
    for cell in sorted({t["cell"] for t in rows}):
        pr = predict(rows, cell)
        top = "  ".join(f"{r['label']} {r['p_hit']:.0%}" for r in pr["ranked"][:3])
        print(f"  {cell:<12} {pr['beta']:>5.2f}  {top}")
        if emit and not any(e["cell"] == cell and not e.get("resolved") for e in log["live"]):
            log["live"].append({"made_at": now, "cell": cell, "beta": pr["beta"],
                                "counts": pr["counts"], "ranked": pr["ranked"][:3],
                                "heuristic": [window_label(w)
                                              for w in heuristic(rows, cell)[:3]],
                                "known_ids": [t["id"] for t in rows if t["cell"] == cell]})
    if emit:
        pend = sum(1 for e in log["live"] if not e.get("resolved"))
        print(f"\n  {pend} live prediction(s) pending; re-run after the next round to resolve them.")
    json.dump(log, open(LOG, "w"), indent=1)


if __name__ == "__main__":
    main()
