#!/usr/bin/env python3
"""seed_predict.py — is the next task's seed window predictable from the history of past ones?

The proposal is to predict which 100-seed windows the next round will draw and concentrate 6-8
hotkeys on the most likely one. That is worth doing only if window selection carries exploitable
structure; if the draws are i.i.d. uniform, any model -- deep or otherwise -- will output confident
probabilities that are no better than 1/9, and concentrating on them costs coverage.

So this tests the premise directly, before any model:

1. **Uniformity.** chi-square on window counts, overall and per cell type. A biased generator would
   show here.
2. **Serial dependence.** Does the window set of task t depend on task t-1? Measured as the repeat
   rate against chance, plus a lag-1..5 autocorrelation of the raw seed values.
3. **The balancing hypothesis** -- the specific mechanism proposed: that a window under-selected so
   far becomes more likely next. Measured as the correlation between a window's running deficit and
   whether it is drawn next.
4. **A predictability ceiling.** Every history-based predictor that could be learned -- frequency,
   recency, deficit, per-cell-type frequency, Markov -- scored against the base rate on a proper
   walk-forward split. A neural net can only learn what these carry; if none beats chance, none
   exists to learn.

Base rate: with 3 seeds drawn from 9 windows, any given window appears with probability
1 - (8/9)^3 = 29.9%, and a predictor naming one window is right that often by luck alone.
"""
import json
import math
import statistics as st
import urllib.request
from collections import Counter, defaultdict

LO = "2026-08-27T00:00:00"
NW = 9                                   # windows 100-199 .. 900-999


def win(seed):
    return min(NW - 1, max(0, seed // 100 - 1))


def load():
    tl = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks", timeout=120))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    out = []
    for t in tl:
        c = (t.get("content") or {}).get("contract") or {}
        s = [x.strip() for x in str(c.get("seed") or "").split(",") if x.strip()]
        ca = t.get("created_at") or ""
        if len(s) == 3 and c.get("cell_type") and ca >= LO:
            out.append({"id": t["id"], "at": ca, "cell": c["cell_type"],
                        "seeds": sorted(int(x) for x in s)})
    out.sort(key=lambda r: r["at"])
    return out


def chisq(counts, n_bins):
    n = sum(counts)
    exp = n / n_bins
    x2 = sum((c - exp) ** 2 / exp for c in counts) if exp else 0.0
    # survival of chi-square with df = n_bins-1, via Wilson-Hilferty
    df = n_bins - 1
    z = ((x2 / df) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    p = 0.5 * math.erfc(z / math.sqrt(2))
    return x2, p


def main():
    rows = load()
    print(f"{len(rows)} three-seed tasks from {rows[0]['at'][:16]} to {rows[-1]['at'][:16]}")
    print(f"{len(rows) * 3} seed draws\n")

    print("=== 1. uniformity of window selection ===")
    allc = Counter()
    for r in rows:
        for s in r["seeds"]:
            allc[win(s)] += 1
    counts = [allc[i] for i in range(NW)]
    x2, p = chisq(counts, NW)
    print(f"  overall counts {counts}")
    print(f"  chi-square {x2:.2f}, df {NW-1}, p = {p:.3f}"
          f"   {'-> consistent with uniform' if p > 0.05 else '-> NOT uniform'}")
    for cell in sorted({r["cell"] for r in rows}):
        c = Counter()
        for r in rows:
            if r["cell"] == cell:
                for s in r["seeds"]:
                    c[win(s)] += 1
        cc = [c[i] for i in range(NW)]
        x2c, pc = chisq(cc, NW)
        print(f"  {cell:<11} {cc}  chi2 {x2c:>5.1f}  p = {pc:.3f}")

    print("\n=== 2. serial dependence between consecutive tasks ===")
    rep = same = 0
    for a, b in zip(rows, rows[1:]):
        wa = {win(s) for s in a["seeds"]}
        wb = {win(s) for s in b["seeds"]}
        rep += len(wa & wb)
        same += len(wb)
    exp_overlap = st.mean(len({win(s) for s in a["seeds"]}) for a in rows) * (1 - (8 / 9) ** 3)
    print(f"  mean windows shared with previous task: {rep / (len(rows) - 1):.3f}")
    print(f"  expected under independence:            {exp_overlap:.3f}")
    seq = [s for r in rows for s in r["seeds"]]
    for lag in (1, 2, 3, 4, 5):
        xs, ys = seq[:-lag], seq[lag:]
        mx, my = st.mean(xs), st.mean(ys)
        num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
        den = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
        print(f"  lag-{lag} autocorrelation of raw seeds: {num / den:+.3f}")

    print("\n=== 3. the balancing hypothesis: does a deficit predict selection? ===")
    pairs = []
    for i, r in enumerate(rows):
        if i < 10:
            continue
        hist = Counter()
        for prev in rows[:i]:
            for s in prev["seeds"]:
                hist[win(s)] += 1
        exp = sum(hist.values()) / NW
        drawn = {win(s) for s in r["seeds"]}
        for w in range(NW):
            pairs.append((exp - hist[w], 1.0 if w in drawn else 0.0))
    xs = [a for a, _ in pairs]
    ys = [b for _, b in pairs]
    mx, my = st.mean(xs), st.mean(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    print(f"  corr(running deficit, drawn next) over {len(pairs)} window-rounds: "
          f"r = {num / den:+.4f}")
    hi = [b for a, b in pairs if a > 0]
    lo = [b for a, b in pairs if a <= 0]
    print(f"  P(drawn | under-selected) = {st.mean(hi):.3f}   "
          f"P(drawn | over-selected) = {st.mean(lo):.3f}   base rate {1-(8/9)**3:.3f}")

    print("\n=== 4. walk-forward accuracy of every history-based predictor ===")
    strategies = {}

    def evaluate(name, choose):
        hit = tot = 0
        for i, r in enumerate(rows):
            if i < 15:
                continue
            w = choose(rows[:i], r["cell"])
            if w is None:
                continue
            tot += 1
            hit += 1 if w in {win(s) for s in r["seeds"]} else 0
        strategies[name] = (hit, tot, hit / tot if tot else 0.0)

    def most_frequent(hist, cell):
        c = Counter(win(s) for p in hist for s in p["seeds"])
        return c.most_common(1)[0][0] if c else None

    def least_frequent(hist, cell):
        c = Counter(win(s) for p in hist for s in p["seeds"])
        return min(range(NW), key=lambda w: c[w])

    def most_frequent_cell(hist, cell):
        c = Counter(win(s) for p in hist if p["cell"] == cell for s in p["seeds"])
        return c.most_common(1)[0][0] if c else None

    def least_frequent_cell(hist, cell):
        c = Counter(win(s) for p in hist if p["cell"] == cell for s in p["seeds"])
        return min(range(NW), key=lambda w: c[w]) if c else None

    def recent(hist, cell):
        return win(hist[-1]["seeds"][0]) if hist else None

    def avoid_recent(hist, cell):
        seen = {win(s) for s in hist[-1]["seeds"]} if hist else set()
        return next((w for w in range(NW) if w not in seen), 0)

    def markov(hist, cell):
        trans = Counter()
        for a, b in zip(hist, hist[1:]):
            for wa in {win(s) for s in a["seeds"]}:
                for wb in {win(s) for s in b["seeds"]}:
                    trans[(wa, wb)] += 1
        last = {win(s) for s in hist[-1]["seeds"]} if hist else set()
        best, score = None, -1
        for w in range(NW):
            sc = sum(trans[(x, w)] for x in last)
            if sc > score:
                best, score = w, sc
        return best

    for name, fn in (("most frequent", most_frequent), ("least frequent", least_frequent),
                     ("most frequent (cell)", most_frequent_cell),
                     ("least frequent (cell)", least_frequent_cell),
                     ("repeat last", recent), ("avoid last", avoid_recent),
                     ("markov lag-1", markov)):
        evaluate(name, fn)
    base = 1 - (8 / 9) ** 3
    print(f"  {'predictor':<24} {'hits':>6} {'rounds':>7} {'accuracy':>9} {'vs chance':>10}")
    for name, (h, t, a) in strategies.items():
        se = math.sqrt(base * (1 - base) / t) if t else 1
        z = (a - base) / se if se else 0
        print(f"  {name:<24} {h:>6} {t:>7} {a:>8.1%} {z:>+9.2f} sd")
    print(f"  {'chance (any 1 window)':<24} {'':>6} {'':>7} {base:>8.1%}")


if __name__ == "__main__":
    main()
