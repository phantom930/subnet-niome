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


# instance name -> on-chain hotkey, in miner.sh table order (niome_hotkey == h0). Needed to join
# a logged assignment to what the validator actually scored that hotkey.
INSTANCES = ["niome_hotkey"] + [f"niome_hotkey{i}" for i in range(1, 9)]
SS58 = [
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW", "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU", "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2", "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",
    "5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3", "5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN",
    "5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2",
]
HOTKEY_OF = dict(zip(INSTANCES, SS58))
BAND_HIT = 0.30            # round consistency above this means a seed really landed in the band


def load_log():
    return json.load(open(LOG)) if os.path.exists(LOG) else {"live": []}


def save_log(log):
    tmp = LOG + ".tmp"
    with open(tmp, "w") as handle:
        json.dump(log, handle, indent=1)
    os.replace(tmp, LOG)


def pending_for(log, cell):
    """The unresolved live entry for this cell type, if any — where window_plan records what it
    actually applied, so a resolved round shows predicted *and* applied windows side by side."""
    for entry in log["live"]:
        if entry["cell"] == cell and not entry.get("resolved"):
            return entry
    return None


SCORE_MAX_AGE_H = 48       # past this a still-unscored hotkey was never contacted, stop retrying
_SCORES = None             # the whole score feed, fetched at most once per process


def _all_scores():
    """The backend's score feed, cached for the process. One resolve plus a backfill of several
    stale entries would otherwise download ~38k rows once per entry."""
    global _SCORES
    if _SCORES is None:
        try:
            raw = json.load(urllib.request.urlopen(
                "https://niome-api.genomes.io/api/v3/miners/scores?limit=40000", timeout=180))
            _SCORES = raw if isinstance(raw, list) else (raw.get("data") or raw.get("items") or [])
        except Exception:
            _SCORES = []
    return _SCORES


def fetch_scores(task_id):
    """Our nine hotkeys' scores on one task, plus their rank in that task's full field."""
    raw = _all_scores()
    best = {}
    for x in raw:
        if x["task_id"] != task_id:
            continue
        hk = x["miner_hotkey"]
        if hk not in best or x["final_score"] > best[hk]["final_score"]:
            best[hk] = x
    if not best:
        return {}
    field = sorted(best.values(), key=lambda y: -y["final_score"])
    rank = {y["miner_hotkey"]: i for i, y in enumerate(field, 1)}
    out = {}
    for inst, hk in HOTKEY_OF.items():
        x = best.get(hk)
        if x:
            b = x["breakdown"]
            out[inst] = {"rank": rank[hk], "final": round(x["final_score"], 2),
                         "consistency": round(b["consistency_factor"], 4),
                         "weighted": round(b["total_weighted_score"], 1),
                         "fidelity": round(b["distribution_fidelity_factor"], 4),
                         "band_hit": b["consistency_factor"] >= BAND_HIT}
    return out


WINDOW_USED = "data/inst/{inst}/window_used.json"


def windows_used(task_id):
    """What each instance actually built this task with, from its own per-instance record.

    Preferred over the round plan's assignments: ``window_plan.json`` is rewritten hourly, so a
    round whose build straddled a cron tick was built from a plan that is no longer on disk, and
    reading the plan back would credit a window the hotkey never used. Instances with no record
    (a build that predates this file, or one that never reached all-HDR) are simply absent.
    """
    out = {}
    for inst in INSTANCES:
        try:
            with open(WINDOW_USED.format(inst=inst)) as handle:
                row = json.load(handle).get(task_id)
        except (FileNotFoundError, json.JSONDecodeError, OSError):
            continue
        # A row with no window is still evidence: it means all-HDR was never attempted, so no
        # band existed. Keep it so the reader can say that rather than "unverified".
        if row and (row.get("window") or row.get("all_hdr_built") is False):
            out[inst] = row
    return out


def per_hotkey_block(entry, task, scored):
    """Join one round's seeds to what each hotkey actually built, and what it scored.

    ``window`` comes from the instance's own record where there is one and from the plan only as a
    last resort, marked ``plan_file_unverified`` so an unverifiable row cannot be mistaken for a
    confirmed one. Where the two disagree the plan's value is kept alongside as
    ``planned_window``: that is the straddled-cron case, and hiding it would make the log look
    self-consistent while being wrong.
    """
    applied = (entry.get("applied") or {}).get("assignments") or {}
    roles = (entry.get("applied") or {}).get("roles") or {}
    shipped = windows_used(task["id"])
    per = {}
    for inst in sorted(set(applied) | set(shipped), key=INSTANCES.index):
        rec = shipped.get(inst)
        planned = applied.get(inst)
        # A record with no window of its own (all-HDR never attempted) still needs a window to
        # report containment against, and the plan's assignment is the right one for that: the
        # question "did the plan name the window a seed came from" is separate from "was a band
        # built there". all_hdr_built below carries the second.
        rec_window = (rec or {}).get("window")
        window_ = tuple(rec_window) if rec_window else tuple(planned or (0, 0))
        if not window_ or window_ == (0, 0):
            continue
        per[inst] = {"window": list(window_),
                     "window_source": rec["source"] if rec else "plan_file_unverified",
                     "window_confirmed": bool(rec_window),
                     "role": roles.get(inst),
                     "seeds_in_window": [x for x in task["seeds"]
                                         if window_[0] <= x <= window_[1]],
                     **(scored.get(inst) or {})}
        # Whether all-HDR actually built the band for that window. A decline ships all-cut rows
        # with no band at all, so containment there says something about the plan's accuracy but
        # nothing about whether a spike was even possible.
        if rec and rec.get("all_hdr_built") is not None:
            per[inst]["all_hdr_built"] = rec["all_hdr_built"]
        if planned and rec_window and list(planned) != list(rec_window):
            per[inst]["planned_window"] = list(planned)
    return per


def _fill(per, scored):
    """Merge a score fetch into a per_hotkey block; True if anything landed."""
    got = False
    for inst, row in per.items():
        add = scored.get(inst)
        if add and row.get("consistency") is None:
            row.update(add)
            got = True
    return got


def backfill_scores(log, now):
    """Resolution captures scores at resolve time, but the backend publishes a round's score rows
    minutes to hours after its seeds are stamped — so an entry resolved by the cron tick right
    after its round carries per_hotkey with `consistency: null` throughout, and resolution is
    one-shot. Re-fetch for any resolved entry still missing scores, until the round is old enough
    that a missing hotkey means it was never contacted rather than not yet scored."""
    filled, stale = 0, 0
    for entry in log["live"]:
        res = entry.get("resolved")
        per = (res or {}).get("per_hotkey")
        if not per or all(v.get("consistency") is not None for v in per.values()):
            continue
        at = datetime.fromisoformat(res["at"])          # backend stamps are naive UTC
        if at.tzinfo is None:
            at = at.replace(tzinfo=timezone.utc)
        age = (now - at).total_seconds() / 3600.0
        if age > SCORE_MAX_AGE_H:
            stale += 1
            continue
        if _fill(per, fetch_scores(res["task_id"])):
            res["scores_at"] = now.isoformat()
            res["hotkeys_with_band_hit"] = sorted(
                k for k, v in per.items() if v.get("band_hit"))
            filled += 1
    return filled, stale


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
    log = load_log()
    by_id = {t["id"]: t for t in rows}
    now_dt = datetime.now(timezone.utc)
    now = now_dt.isoformat()

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
        hlab = entry.get("heuristic") or []
        entry["resolved"] = {"task_id": t["id"], "at": t["at"], "seeds": t["seeds"],
                             "actual_windows": actual,
                             "top1_hit": top[0] in actual,
                             "top2_hit": bool(set(top[:2]) & set(actual)),
                             "heuristic_hit": bool(hlab) and hlab[0] in
                             [window_label(w) for w in actual]}
        # What each hotkey was actually running, whether a seed fell in its window, and what the
        # validator scored it. Containment is not a hit: a concentrated hotkey holds ~92% of its
        # 16-seed window as band, a width-100 spread hotkey only ~13% of its 100 -- so band_hit
        # (consistency >= 0.30) is ground truth and seeds_in_window is the part the plan controls.
        scored = fetch_scores(t["id"])
        per = per_hotkey_block(entry, t, scored)
        if per:
            entry["resolved"]["per_hotkey"] = per
            if scored:
                entry["resolved"]["scores_at"] = now
            entry["resolved"]["hotkeys_with_seed_in_window"] = sorted(
                k for k, v in per.items() if v["seeds_in_window"])
            declined = sorted(k for k, v in per.items() if v.get("all_hdr_built") is False)
            if declined:
                entry["resolved"]["hotkeys_without_band"] = declined
            entry["resolved"]["hotkeys_with_band_hit"] = sorted(
                k for k, v in per.items() if v.get("band_hit"))
        resolved_now += 1

    print(f"{len(rows)} three-seed tasks, {rows[0]['at'][:16]} .. {rows[-1]['at'][:16]}")
    if resolved_now:
        print(f"resolved {resolved_now} pending live prediction(s) this run")

    filled, stale = backfill_scores(log, now_dt)
    if filled:
        print(f"backfilled validator scores for {filled} earlier round(s)")
    if stale:
        print(f"{stale} resolved round(s) still unscored past {SCORE_MAX_AGE_H}h — "
              f"those hotkeys were never contacted, not merely unscored")

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
            print(f"    {r['at'][:16]}  {e['cell']:<11} task {r['task_id'][:8]}  predicted "
                  f"{e['ranked'][0]['label']:<9} actual "
                  f"{','.join(window_label(w) for w in r['actual_windows']):<30} "
                  f"{'HIT' if r['top1_hit'] else 'miss'}")
            if r.get("per_hotkey"):
                inw = r.get("hotkeys_with_seed_in_window") or []
                band = r.get("hotkeys_with_band_hit") or []
                n_scored = sum(1 for v in r["per_hotkey"].values()
                               if v.get("consistency") is not None)
                # "none" and "not scored yet" are different facts -- never print the first for the
                # second, or an unpublished round reads as nine hotkeys that all missed the band.
                hit = (f"{', '.join(band)}" if band else
                       "none" if n_scored == len(r["per_hotkey"]) else
                       f"none of {n_scored} scored" if n_scored else "not scored yet")
                print(f"      applied {len(r['per_hotkey'])} hotkeys | seed in window: "
                      f"{', '.join(inw) if inw else 'none'} | scored band hit: {hit}")
                drift = [k for k, v in r["per_hotkey"].items() if v.get("planned_window")]
                # A missing window_source is a round resolved before per-instance recording
                # existed -- unverified for the same reason, so default to that rather than
                # letting an older entry read as confirmed.
                unver = [k for k, v in r["per_hotkey"].items()
                         if not v.get("window_confirmed", False)]
                if drift:
                    print(f"      NOTE {len(drift)} hotkey(s) built a window the plan did not "
                          f"assign: {', '.join(drift)} — the plan was rewritten mid-round")
                nb = r.get("hotkeys_without_band") or []
                if nb:
                    print(f"      NOTE {len(nb)} hotkey(s) shipped without a band (all-HDR "
                          f"declined): {', '.join(nb)} — no spike was possible for them")
                if unver:
                    print(f"      NOTE {len(unver)} hotkey(s) unverified (no per-instance "
                          f"record); their window is the plan's, not a confirmed build")

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
    save_log(log)


if __name__ == "__main__":
    main()
