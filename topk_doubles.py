#!/usr/bin/env python3
"""topk_doubles.py — top-1..4 accuracy, and how often the prediction lands on a DOUBLE window.

Two questions, walk-forward over every 3-seed task in the range (each prediction sees only
strictly earlier tasks of that cell type):

1. top-k: does any of the k highest-ranked windows contain at least one of the round's seeds?
2. top-k DOUBLE: does any of them contain **two or more** — the case worth having.

Why (2) is the one that matters. A hotkey whose band catches one seed scores consistency 0.40 and
clears ~60% of rank-10 cutoffs; catching two scores 0.70 and clears ~99%. With three seeds over
nine windows a double exists on 30.9% of rounds and, when it does, sits in exactly one window --
so naming it is a 1-in-9 shot per pick, and chance for top-k is 0.3086 * k/9.
"""
import os, sys, json
sys.argv = ["x"]; sys.path.insert(0, "/root/workspace/subnet-niome")
import importlib.util
from collections import Counter, defaultdict

spec = importlib.util.spec_from_file_location("swm", "seed_window_model.py")
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)

SINCE = os.getenv("TK_SINCE", "2026-08-27T00:00:00")
UNTIL = os.getenv("TK_UNTIL", "2026-09-05T23:59:59")
KMAX = 4


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    m = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return max(0.0, c - m), min(1.0, c + m)


def main():
    rows = [t for t in M.load_tasks() if SINCE <= t["at"] <= UNTIL]
    print(f"{len(rows)} three-seed tasks, {rows[0]['at'][:16]} .. {rows[-1]['at'][:16]}")
    print(f"warmup {M.WARMUP} tasks per cell type before the model may predict\n")

    hit = [0] * (KMAX + 1)
    dbl = [0] * (KMAX + 1)
    n = 0
    n_with_double = 0
    per_cell = defaultdict(lambda: {"n": 0, "nd": 0,
                                    "h": [0] * (KMAX + 1), "d": [0] * (KMAX + 1)})
    detail = []
    for i, t in enumerate(rows):
        hist = rows[:i]
        if sum(1 for x in hist if x["cell"] == t["cell"]) < M.WARMUP:
            continue
        pr = M.predict(hist, t["cell"])
        top = [x["window"] for x in pr["ranked"]][:KMAX]
        c = Counter(M.window(s) for s in t["seeds"])
        present = set(c)
        doubles = {w for w, k in c.items() if k >= 2}
        n += 1
        pc = per_cell[t["cell"]]; pc["n"] += 1
        if doubles:
            n_with_double += 1
            pc["nd"] += 1
        for k in range(1, KMAX + 1):
            if set(top[:k]) & present:
                hit[k] += 1; pc["h"][k] += 1
            if set(top[:k]) & doubles:
                dbl[k] += 1; pc["d"][k] += 1
        if doubles:
            detail.append((t["at"][:16], t["cell"], t["id"][:8], t["seeds"],
                           [M.window_label(w) for w in sorted(doubles)],
                           [M.window_label(w) for w in top],
                           min((j + 1 for j, w in enumerate(top) if w in doubles), default=None)))

    print(f"scored {n} predictions ({n_with_double} rounds had a double window "
          f"= {n_with_double/n:.1%}; theory 30.9%)\n")
    print(f"{'k':>2} | {'ANY seed in a top-k window':<38} | {'DOUBLE in a top-k window':<34}")
    print(f"{'':>2} | {'hits':>5} {'rate':>7} {'95% CI':>16} {'chance':>7} | "
          f"{'hits':>5} {'rate':>7} {'chance':>7} {'lift':>6}")
    for k in range(1, KMAX + 1):
        ch = 1 - ((9 - k) / 9) ** 3
        chd = 0.3086 * k / 9
        lo, hi = wilson(hit[k], n)
        print(f"{k:>2} | {hit[k]:>5} {hit[k]/n:>6.1%} [{lo:>5.1%},{hi:>5.1%}] {ch:>7.1%} | "
              f"{dbl[k]:>5} {dbl[k]/n:>6.1%} {chd:>7.1%} {dbl[k]/n/chd:>5.2f}x")

    from math import comb

    def pv(k_, n_, p_):
        return sum(comb(n_, i) * p_ ** i * (1 - p_) ** (n_ - i) for i in range(k_, n_ + 1))

    print("\n" + "=" * 78)
    print("PER CELL TYPE\n")
    for cell, d in sorted(per_cell.items()):
        m = d["n"]
        print(f"  {cell}  (n={m} scored predictions)")
        print(f"    {'k':>2} | {'ANY seed in top-k':<34} | {'DOUBLE in top-k':<30}")
        print(f"    {'':>2} | {'hits':>5} {'rate':>7} {'chance':>7} {'lift':>6} {'p':>6} | "
              f"{'hits':>5} {'rate':>7} {'chance':>7} {'lift':>6} {'p':>6}")
        for k in range(1, KMAX + 1):
            ch = 1 - ((9 - k) / 9) ** 3
            chd = 0.3086 * k / 9
            hk_, dk_ = d["h"][k], d["d"][k]
            print(f"    {k:>2} | {hk_:>5} {hk_/m:>6.1%} {ch:>7.1%} {hk_/m/ch:>5.2f}x "
                  f"{pv(hk_, m, ch):>6.3f} | {dk_:>5} {dk_/m:>6.1%} {chd:>7.1%} "
                  f"{(dk_/m/chd if chd else 0):>5.2f}x {pv(dk_, m, chd):>6.3f}")
        print(f"    doubles occurred on {d['nd']}/{m} = {d['nd']/m:.1%} of this cell's rounds")
        print()

    print(f"\nrounds that HAD a double, and where the model ranked it:")
    print(f"{'when':<17} {'cell':<11} {'task':<9} {'seeds':<18} {'double':<9} {'rank'}")
    for at, cell, tid, seeds, dw, top, rank in detail:
        print(f"{at:<17} {cell:<11} {tid:<9} {str(seeds):<18} {dw[0]:<9} "
              f"{('#' + str(rank)) if rank else 'not in top 4'}")
    ranks = [r for *_, r in detail if r]
    if ranks:
        print(f"\ndouble was inside the top 4 on {len(ranks)}/{len(detail)} of the rounds that had "
              f"one; chance {4/9:.1%}")


if __name__ == "__main__":
    main()
