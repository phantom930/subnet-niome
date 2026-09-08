#!/usr/bin/env python3
"""sd_sweep_all.py — paired comparison of the seed-depend configs across every seed-0 round.

One round decides nothing here: the margins that matter are 1-2 points on a ~300 score, and the
variant tie-break alone moves 0.03-0.14. So each config is compared to ``base`` *within* contract
and the deltas are pooled, the same design the CELL_CONFIG measurements in CLAUDE.md use.

Reports the payout consequence as well as the score, because they are not the same question: a
config that adds a point everywhere is worthless if the point never crosses a rank boundary, and a
config that adds three points on one round can be worth more than the whole rest of the sweep.
"""
import glob
import json
import math
import statistics as st

DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]


def main():
    runs = [json.load(open(p)) for p in sorted(glob.glob("sd_sweep_*.json"))]
    runs = [r for r in runs if "configs" in r]
    if not runs:
        raise SystemExit("no sd_sweep_*.json found")
    names = [n for n in ("base", "hdr", "vps12k", "step1", "greedy1200",
                        "hdr_vps12k", "vps24k")
             if any(n in r["configs"] and "final" in r["configs"][n] for r in runs)]

    print(f"{len(runs)} seed-0 rounds, paired within contract\n")
    print(f"  {'round':<10} {'cell':<11} {'field':>8} " + "".join(f"{n:>11}" for n in names))
    for r in sorted(runs, key=lambda x: x["task"]):
        cells = ""
        for n in names:
            c = r["configs"].get(n, {})
            cells += f"{c['final']:>11.2f}" if "final" in c else f"{'declined':>11}"
        print(f"  {r['task'][:8]:<10} {r['cell']:<11} {r['field_best']:>8.2f} " + cells)

    print(f"\n  {'config':<12} {'mean delta':>11} {'sem':>7} {'t':>7} {'wins':>7} "
          f"{'mean share':>11} {'rank<=4':>8} {'build':>7}")
    base_shares = [r["configs"]["base"]["share"] for r in runs if "share" in r["configs"]["base"]]
    for n in names:
        deltas, shares, ranks, builds = [], [], [], []
        for r in runs:
            c, b = r["configs"].get(n, {}), r["configs"].get("base", {})
            if "final" not in c or "final" not in b:
                continue
            deltas.append(c["final"] - b["final"])
            shares.append(c["share"])
            ranks.append(c["rank"])
            builds.append(c["build_s"])
        if not deltas:
            continue
        sem = st.stdev(deltas) / math.sqrt(len(deltas)) if len(deltas) > 1 else 0.0
        t = deltas and (st.mean(deltas) / sem if sem else float("inf"))
        wins = sum(1 for d in deltas if d > 0)
        print(f"  {n:<12} {st.mean(deltas):>+11.2f} {sem:>7.2f} "
              f"{(t if isinstance(t, float) else 0):>7.2f} {wins:>3}/{len(deltas):<3} "
              f"{st.mean(shares):>11.1%} {sum(1 for x in ranks if x <= 4):>4}/{len(ranks):<3} "
              f"{max(builds):>6.0f}s")

    print(f"\n  base mean share {st.mean(base_shares):.1%} — a config only pays if it moves this.")
    print("  'rank<=4' counts rounds where ONE hotkey places top 4, i.e. where four siblings")
    print("  would still sweep the head of the curve rather than land in its tail.")

    # What the fleet takes, per config: four correlated siblings occupy consecutive ranks.
    print(f"\n  {'config':<12} {'4-sibling share (mean over rounds)':>36}")
    for n in names:
        tot = 0.0
        cnt = 0
        for r in runs:
            c = r["configs"].get(n, {})
            if "rank" not in c:
                continue
            start = c["rank"]
            tot += sum(DIST[i - 1] for i in range(start, start + 4) if i <= len(DIST))
            cnt += 1
        if cnt:
            print(f"  {n:<12} {tot / cnt:>35.1%}")


if __name__ == "__main__":
    main()
