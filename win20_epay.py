#!/usr/bin/env python3
"""win20_epay.py — price width 20 against the shipped width 100 as expected payout.

Band and score are measured (win20_cells / win20_narrow / win20_fill, 168 builds, 0 declines);
what is not yet known is whether the trade is worth taking. Width 20 buys +1.7 to +1.9 band seeds
on the erythroid types and gives up 2-7 points of total_weighted_score, and ``SCORING_SYSTEM =
"top"`` pays only the top 10 on a fixed curve — so a score change that crosses no rank threshold
pays nothing while the frequency change pays on every round. Score deltas are therefore not
decision-relevant on their own; only E[pay] is.

Method, following the group-80 measurement this repo already trusts:

* **Fields are real and current-regime only.** Each round's field is the other miners' actual
  final scores on that task, taken from /api/v3/miners/scores. Only 3-seed tasks count: the
  single-seed era let the whole field converge on consistency 1.000 and its cutoffs are meaningless
  now (rank-10 medians differ by more than the effect being measured). Our own nine hotkeys are
  removed from the field before our hypothetical scores are inserted, or we would be competing
  against ourselves.
* **Exact enumeration over rounds, not a k=1 point estimate.** A round draws three seeds; each
  lands in at most one hotkey's band because the bands are disjoint. Enumerating the seven possible
  patterns (how many seeds hit, and whether they share a hotkey) gets the k=2 and k=3 tails exactly
  and handles our own hotkeys displacing each other in the ranking when two place on one round.
* **Per contract, not per cell-type mean.** Payout is a step function, so E[pay(mean config)] is
  not mean E[pay]: our weighted spans 193-342 across K562 contracts, and whether a round places
  turns on that spread. Every contract's own measured (band, weighted, fidelity) is crossed with
  every field.
* **All nine hotkeys are ranked, including the misses.** A k=0 hotkey still submits and still
  scores ~29; in a thin field that can place, and dropping it would understate both configs.

consistency at k hits is ``(k + (3 - k) * FLOOR) / 3``, the arithmetic confirmed by the measured
k=2 round: predicted 0.703, measured 0.703.
"""
import json
import math
import os
import random
import statistics as st
import urllib.request
from collections import defaultdict

DIST = [0.3, 0.2, 0.2, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
BURNING_RATE = 0.02
FLOOR = 0.10                      # off-band consistency; measured 0.086-0.106 across these builds
HOTKEYS = 9
SEED_SPACE = 900
ROUNDS = 20_000
CELLS = ["CD34+_HSPC", "HUDEP-2", "K562"]
WIDTHS = [100, 20]
SOURCES = ["win20_cells.json", "win20_narrow.json", "win20_fill.json"]
OURS = {
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW", "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU", "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2", "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",
    "5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3", "5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN",
    "5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2",
}


def configs():
    """Measured (band, weighted, fidelity) per cell type and width, on paired contracts only."""
    recs = []
    for path in SOURCES:
        if os.path.exists(path):
            recs += json.load(open(path))
    seen, uniq = set(), []
    for r in recs:
        key = (r["task_id"], r["width"], r["window"][0])
        if key in seen or "declined" in r:
            continue
        seen.add(key)
        uniq.append(r)
    out = {}
    for cell in CELLS:
        per = {}
        for w in WIDTHS:
            g = [r for r in uniq if r["cell_type"] == cell and r["width"] == w]
            per[w] = {t: [x for x in g if x["task_id"] == t] for t in {x["task_id"] for x in g}}
        common = sorted(set.intersection(*(set(per[w]) for w in WIDTHS)))
        for w in WIDTHS:
            out[(cell, w)] = [
                {"task_id": t,
                 "band": st.mean(r["band_60k"] for r in per[w][t]),
                 "weighted": st.mean(r["total_weighted_score"] for r in per[w][t]),
                 "fid": st.mean(r["distribution_fidelity_factor"] for r in per[w][t])}
                for t in common]
    return out


def fields():
    """Other miners' final scores per current-regime (3-seed) task, keyed by cell type."""
    tl = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks", timeout=120))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    meta = {}
    for t in tl:
        c = (t.get("content") or {}).get("contract") or {}
        seeds = [s for s in str(c.get("seed") or "").split(",") if s.strip()]
        if len(seeds) == 3 and c.get("cell_type"):
            meta[t["id"]] = c["cell_type"]
    sc = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/miners/scores?limit=40000", timeout=180))
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    best = defaultdict(dict)
    for x in sc:
        if x["task_id"] not in meta or x["miner_hotkey"] in OURS:
            continue
        cur = best[x["task_id"]].get(x["miner_hotkey"])
        if cur is None or x["final_score"] > cur:
            best[x["task_id"]][x["miner_hotkey"]] = x["final_score"]
    out = defaultdict(list)
    for tid, hk in best.items():
        if len(hk) >= 10:
            out[meta[tid]].append(sorted(hk.values(), reverse=True))
    return out


def _payout(mine, cut):
    """Payout share for one round: every hotkey ranked against the field and against each other."""
    total = 0.0
    for s in mine:
        rank = 1 + sum(1 for f in cut if f > s) + sum(1 for o in mine if o > s)
        if rank <= 10:
            total += DIST[rank - 1] * (1 - BURNING_RATE)
    return total


def epay(band, weighted, fid, field):
    """Exact expected payout share per round for a 9-hotkey coldkey against one real field.

    Three seeds each land in one of the nine disjoint bands (probability ``p`` each) or in none.
    Given ``m`` seeds landing somewhere, which hotkeys they pick is uniform over the nine, so the
    round collapses to seven patterns whose probabilities are exact.
    """
    p = band / SEED_SPACE
    a = HOTKEYS * p                        # a seed lands in *some* band
    q = 1.0 - a
    s = [weighted * ((k + (3 - k) * FLOOR) / 3) * fid for k in range(4)]
    cut = field[:12]                       # nothing below rank 12 can pay
    z = [s[0]] * HOTKEYS
    def mix(*ks):
        return [s[k] for k in ks] + [s[0]] * (HOTKEYS - len(ks))
    patterns = [
        (q ** 3,                     z),
        (3 * q * q * a,              mix(1)),
        (3 * q * a * a / 9,          mix(2)),
        (3 * q * a * a * 8 / 9,      mix(1, 1)),
        (a ** 3 / 81,                mix(3)),
        (a ** 3 * 24 / 81,           mix(2, 1)),
        (a ** 3 * 56 / 81,           mix(1, 1, 1)),
    ]
    return sum(w * _payout(m, cut) for w, m in patterns)


def main():
    cfg = configs()
    fld = fields()
    print("fields per cell type (current-regime, >=10 miners, our hotkeys removed):")
    for cell in CELLS:
        f = fld.get(cell, [])
        if f:
            print(f"  {cell:<11} {len(f):>3} fields | miners median {st.median(len(x) for x in f):>4.0f}"
                  f" | rank-10 cutoff median {st.median(x[9] for x in f):>6.1f}"
                  f" | rank-1 median {st.median(x[0] for x in f):>6.1f}")
    print(f"\n{'cell':<11} {'w':>4} {'n':>2} {'band':>6} {'weighted':>9} {'fid':>6} "
          f"{'k=1 final':>10} {'E[pay]':>9} {'vs w100':>9}")
    agg = defaultdict(float)
    per_contract = defaultdict(dict)
    for cell in CELLS:
        f = fld.get(cell, [])
        if not f:
            print(f"{cell:<11} no current-regime fields")
            continue
        base = None
        for w in WIDTHS:
            rows = cfg[(cell, w)]
            byc = {c["task_id"]: st.mean(epay(c["band"], c["weighted"], c["fid"], x) for x in f)
                   for c in rows}
            per_contract[cell][w] = byc
            e = st.mean(byc.values())
            agg[w] += e
            k1 = st.mean(c["weighted"] * ((1 + 2 * FLOOR) / 3) * c["fid"] for c in rows)
            rel = "--" if base is None else f"{(e / base - 1) * 100:+8.1f}%"
            if base is None:
                base = e
            print(f"{cell:<11} {w:>4} {len(rows):>2} "
                  f"{st.mean(c['band'] for c in rows):>6.2f} "
                  f"{st.mean(c['weighted'] for c in rows):>9.1f} "
                  f"{st.mean(c['fid'] for c in rows):>6.3f} {k1:>10.1f} {e:>9.5f} {rel:>9}")
    print(f"\n{'AGGREGATE':<11} w100 {agg[100]:>9.5f}   w20 {agg[20]:>9.5f}   "
          f"{(agg[20] / agg[100] - 1) * 100:+.1f}%")

    print("\nper contract, E[pay] w20 vs w100 (a step function: contracts move in jumps):")
    for cell in CELLS:
        if cell not in per_contract:
            continue
        d100, d20 = per_contract[cell][100], per_contract[cell][20]
        deltas = [(d20[t] / d100[t] - 1) * 100 if d100[t] else float("nan") for t in d100]
        wins = sum(1 for x in deltas if x > 0.5)
        losses = sum(1 for x in deltas if x < -0.5)
        print(f"  {cell:<11} {wins} better / {losses} worse / {len(deltas) - wins - losses} flat"
              f"   [{', '.join(f'{x:+.0f}%' for x in sorted(deltas))}]")


if __name__ == "__main__":
    main()
