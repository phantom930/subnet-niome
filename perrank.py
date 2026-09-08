#!/usr/bin/env python3
"""perrank.py — hit rate of EACH ranked window individually, not cumulatively.

top-k cumulative rates are dominated by k itself: naming more windows always hits more, so 78.9%
at k=4 says little. This asks the question that actually tests the ordering: does the window the
model ranks #1 hit more often than the one it ranks #4?

Every individual window has the same baseline -- a specific window appears among three uniform
seeds with probability 1-(8/9)^3 = 29.8%, and is the double window with probability 0.3086/9 =
3.43% -- so the four rates are directly comparable to each other and to one number.
"""
import os, sys, json
sys.argv = ["x"]; sys.path.insert(0, "/root/workspace/subnet-niome")
import importlib.util
from collections import Counter, defaultdict
from math import comb

spec = importlib.util.spec_from_file_location("swm", "seed_window_model.py")
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)

SINCE = os.getenv("PR_SINCE", "2026-08-27T00:00:00")
UNTIL = os.getenv("PR_UNTIL", "2026-09-07T23:59:59")
KMAX = int(os.getenv("PR_KMAX", "4"))
CH_ANY = 1 - (8 / 9) ** 3
CH_DBL = 0.3086 / 9


def pv(k, n, p):
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1)) if n else 1.0


def wilson(k, n, z=1.96):
    if not n:
        return 0.0, 0.0
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    m = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return max(0.0, c - m), min(1.0, c + m)


def main():
    rows = [t for t in M.load_tasks() if SINCE <= t["at"] <= UNTIL]
    hit = [0] * (KMAX + 1)
    dbl = [0] * (KMAX + 1)
    per = defaultdict(lambda: {"n": 0, "h": [0] * (KMAX + 1), "d": [0] * (KMAX + 1)})
    n = 0
    for i, t in enumerate(rows):
        hist = rows[:i]
        if sum(1 for x in hist if x["cell"] == t["cell"]) < M.WARMUP:
            continue
        top = [x["window"] for x in M.predict(hist, t["cell"])["ranked"]][:KMAX]
        c = Counter(M.window(s) for s in t["seeds"])
        doubles = {w for w, k in c.items() if k >= 2}
        n += 1
        p = per[t["cell"]]; p["n"] += 1
        for k in range(1, KMAX + 1):
            w = top[k - 1]
            if w in c:
                hit[k] += 1; p["h"][k] += 1
            if w in doubles:
                dbl[k] += 1; p["d"][k] += 1

    print(f"{len(rows)} three-seed tasks {rows[0]['at'][:10]}..{rows[-1]['at'][:10]}; "
          f"{n} scored predictions (warmup {M.WARMUP}/cell)")
    print(f"baseline for EVERY individual rank: any-seed {CH_ANY:.1%}, double {CH_DBL:.2%}\n")
    print(f"{'rank':>5} | {'that window held a seed':<40} | {'that window was the double':<28}")
    print(f"{'':>5} | {'hits':>5} {'rate':>7} {'95% CI':>16} {'lift':>6} {'p':>6} | "
          f"{'hits':>5} {'rate':>7} {'lift':>6} {'p':>6}")
    for k in range(1, KMAX + 1):
        lo, hi = wilson(hit[k], n)
        print(f"  #{k:<3} | {hit[k]:>5} {hit[k]/n:>6.1%} [{lo:>5.1%},{hi:>5.1%}] "
              f"{hit[k]/n/CH_ANY:>5.2f}x {pv(hit[k], n, CH_ANY):>6.3f} | "
              f"{dbl[k]:>5} {dbl[k]/n:>6.1%} {dbl[k]/n/CH_DBL:>5.2f}x "
              f"{pv(dbl[k], n, CH_DBL):>6.3f}")
    # is rank 1 better than rank 4? the whole question, as a paired sign test
    print(f"\n  monotone in rank?  rates: " + "  ".join(f"#{k} {hit[k]/n:.1%}" for k in
                                                        range(1, KMAX + 1)))
    print(f"\nper cell type — 'held a seed' rate by rank (baseline {CH_ANY:.1%} each):")
    print(f"  {'cell':<12} {'n':>3} " + " ".join(f"{'#'+str(k):>7}" for k in range(1, KMAX + 1))
          + "   | doubles by rank")
    for cell, d in sorted(per.items()):
        m = d["n"]
        print(f"  {cell:<12} {m:>3} " + " ".join(f"{d['h'][k]/m:>7.1%}" for k in range(1, KMAX + 1))
              + "   | " + " ".join(f"{d['d'][k]}" for k in range(1, KMAX + 1)))


if __name__ == "__main__":
    main()
