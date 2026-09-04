#!/usr/bin/env python3
"""mingap_detail.py — per-miner consistency vs seed gap, and the window geometry it implies.

mingap.py established the gate: r(min gap, count of cons>=0.60) is -0.85 / -0.96 / -0.90 on
CD34+_HSPC, HUDEP-2 and K562, and the regime is essentially absent above gap 200. HEK293 never
reaches 0.60 at any gap.

This resolves it per miner and infers the window structure. A hotkey reaches k=2 -- consistency
(1.0 + 1.0 + floor)/3 ~ 0.70 -- only if TWO of the round's three seeds lie in its own band, which
requires them inside its window. So every k=2 event pins an interval that hotkey's window must
contain, and the union of those intervals across rounds is a *lower bound* on its window span. A
span past 100 falsifies disjoint 100-seed tiling for that hotkey.

Two coldkey-level signatures separate overlapping windows from disjoint ones:

* two hotkeys of one coldkey spiking on the SAME round -- impossible under disjoint tiling, since
  only one window can contain a given pair;
* the number of simultaneous spikes falling as the gap widens, which estimates the stride: windows
  of width W at stride S cover an interval of length g about (W - g)/S + 1 times.

Which two of the three seeds a k=2 round caught is not observable, so each event contributes all
three candidate pairs and the smallest consistent span is reported.
"""
import json
import statistics as st
import urllib.request
from collections import defaultdict
from itertools import combinations

HIGH = 0.60
K1_LO, K1_HI = 0.33, 0.45
OUT = "mingap_detail.json"


def load():
    tl = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks", timeout=120))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    tasks = {}
    for t in tl:
        c = (t.get("content") or {}).get("contract") or {}
        s = [x.strip() for x in str(c.get("seed") or "").split(",") if x.strip()]
        if len(s) == 3 and c.get("cell_type"):
            seeds = sorted(int(x) for x in s)
            tasks[t["id"]] = {"seeds": seeds, "cell_type": c["cell_type"],
                              "created_at": t.get("created_at", ""),
                              "min_gap": min(seeds[1] - seeds[0], seeds[2] - seeds[1])}
    sc = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/miners/scores?limit=40000", timeout=180))
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    best = defaultdict(dict)
    for x in sc:
        if x["task_id"] not in tasks:
            continue
        d = best[x["miner_hotkey"]]
        if x["task_id"] not in d or x["final_score"] > d[x["task_id"]]["final_score"]:
            d[x["task_id"]] = x
    return tasks, best, json.load(open("hotkey_coldkey.json"))


def corr(xs, ys):
    if len(xs) < 4:
        return None
    mx, my = st.mean(xs), st.mean(ys)
    n = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    d = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    return n / d if d else None


def main():
    tasks, best, ck = load()
    print(f"{len(tasks)} three-seed rounds\n")

    print("=== per-miner: correlation of consistency with the round's min seed gap ===")
    print(f"{'hotkey':<11} {'coldkey':<11} {'n':>3} {'k2 rounds':>10} {'r(gap,cons)':>12} "
          f"{'cons|gap<100':>13} {'cons|gap>=200':>14}")
    rows = []
    for hk, v in best.items():
        if len(v) < 30:
            continue
        g = [tasks[t]["min_gap"] for t in v]
        c = [r["breakdown"]["consistency_factor"] for r in v.values()]
        k2 = [t for t, r in v.items() if r["breakdown"]["consistency_factor"] >= HIGH]
        near = [x for t, x in zip(v, c) if tasks[t]["min_gap"] < 100]
        far = [x for t, x in zip(v, c) if tasks[t]["min_gap"] >= 200]
        rows.append((len(k2), hk, len(v), corr(g, c), st.mean(near) if near else 0,
                     st.mean(far) if far else 0, k2))
    rows.sort(reverse=True)
    for n2, hk, n, r, near, far, _k2 in rows[:12]:
        print(f"{hk[:9]:<11} {ck.get(hk,'?')[:9]:<11} {n:>3} {n2:>10} "
              f"{(f'{r:+.3f}' if r is not None else '-'):>12} {near:>13.3f} {far:>14.3f}")

    print("\n=== window span implied by k=2 events (each pins an interval the window must hold) ===")
    print(f"{'hotkey':<11} {'coldkey':<11} {'k2':>3} {'best-fit interval':>22} {'span':>6}  rounds")
    spans = []
    for n2, hk, n, _r, _a, _b, k2 in rows:
        if n2 < 2:
            continue
        # each k=2 round contributes its three candidate pairs; pick one pair per round so the
        # covering interval is as small as possible (a lower bound on the true window span)
        choices = [[(min(p), max(p)) for p in combinations(tasks[t]["seeds"], 2)] for t in k2]
        bestspan, bestiv = None, None
        def rec(i, lo, hi):
            nonlocal bestspan, bestiv
            if bestspan is not None and hi - lo >= bestspan:
                return
            if i == len(choices):
                bestspan, bestiv = hi - lo, (lo, hi)
                return
            for a, b in choices[i]:
                rec(i + 1, min(lo, a), max(hi, b))
        rec(0, 10**9, -10**9)
        spans.append((bestspan, hk, n2, bestiv, k2))
        print(f"{hk[:9]:<11} {ck.get(hk,'?')[:9]:<11} {n2:>3} "
              f"{f'{bestiv[0]}-{bestiv[1]}':>22} {bestspan:>6}  "
              f"{', '.join(t[:6] for t in k2[:4])}")
    if spans:
        s = [x[0] for x in spans]
        print(f"\n  implied span: median {st.median(s):.0f}, min {min(s)}, max {max(s)}, "
              f"n={len(s)} hotkeys with >=2 k=2 events")
        print(f"  hotkeys whose span exceeds 100 (disjoint 100-tiling impossible): "
              f"{sum(1 for x in s if x > 100)}/{len(s)}")

    print("\n=== same-coldkey simultaneous spikes (impossible under disjoint tiling) ===")
    ev = defaultdict(list)
    for hk, v in best.items():
        for t, r in v.items():
            if r["breakdown"]["consistency_factor"] >= HIGH:
                ev[t].append(hk)
    print(f"{'task':<9} {'cell':<11} {'gap':>5} {'seeds':<20} {'coldkey':<11} {'hotkeys':>8}")
    multi = 0
    for t, hks in sorted(ev.items(), key=lambda kv: tasks[kv[0]]["min_gap"]):
        by = defaultdict(list)
        for hk in hks:
            by[ck.get(hk, "?")].append(hk)
        for c, group in by.items():
            if len(group) >= 2:
                multi += 1
                print(f"{t[:8]:<9} {tasks[t]['cell_type']:<11} {tasks[t]['min_gap']:>5} "
                      f"{str(tasks[t]['seeds']):<20} {c[:9]:<11} {len(group):>8}")
    print(f"  {multi} (round, coldkey) pairs with 2+ simultaneous k=2 spikes")
    json.dump({"spans": [{"hotkey": h, "coldkey": ck.get(h), "k2": n, "interval": iv, "span": s}
                         for s, h, n, iv, _k in spans]}, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
