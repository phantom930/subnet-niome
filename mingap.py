#!/usr/bin/env python3
"""mingap.py — does the gap between a round's three seeds gate the high-consistency regime?

A hotkey can only reach k=2 if two of the round's three seeds fall inside the *same* hotkey's band,
which requires them inside the same window. So the closest pair's gap should gate k=2 directly:
consistency 0.70 is `(1.0 + 1.0 + 0.10)/3`, and the 0.45-0.65 band is a full hit plus a near hit.
With disjoint 100-seed windows a pair only ever shares a window when the gap is small *and* the
pair does not straddle a boundary; with wider, overlapping windows many more hotkeys cover a given
pair, so one coldkey can spike on several hotkeys at once.

Two things this measures that the earlier consistency work did not:

* **min gap against the incidence of consistency >= 0.6**, per cell type, over the latest 3-seed
  rounds -- the direct test of the gate.
* **how many hotkeys of the *same coldkey* spike together**, which is what distinguishes
  overlapping windows from disjoint ones. Under disjoint tiling at most one hotkey can hold both
  seeds, so two or more simultaneous spikes from one coldkey is evidence of overlap, and the count
  estimates the overlap: for windows of width W at stride S, an interval of length g is covered by
  about (W - g)/S + 1 of them.

Coldkeys come from the chain metagraph (hotkey_coldkey.json), not from the scores API, which
exposes only miner_hotkey and miner_uid.
"""
import json
import math
import statistics as st
import sys
import urllib.request
from collections import Counter, defaultdict

LATEST = int(sys.argv[1]) if len(sys.argv) > 1 else 20
HIGH = 0.60
MID = 0.45
OUT = "mingap.json"
OURS = {
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW", "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU", "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2", "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",
    "5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3", "5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN",
    "5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2",
}


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
            gaps = [seeds[1] - seeds[0], seeds[2] - seeds[1], seeds[2] - seeds[0]]
            tasks[t["id"]] = {"seeds": seeds, "min_gap": min(gaps[0], gaps[1]),
                              "cell_type": c["cell_type"], "created_at": t.get("created_at", "")}
    sc = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/miners/scores?limit=40000", timeout=180))
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    best = defaultdict(dict)
    for x in sc:
        if x["task_id"] not in tasks:
            continue
        d = best[x["task_id"]]
        hk = x["miner_hotkey"]
        if hk not in d or x["final_score"] > d[hk]["final_score"]:
            d[hk] = x
    ck = json.load(open("hotkey_coldkey.json"))
    return tasks, best, ck


def corr(xs, ys):
    if len(xs) < 3:
        return 0.0
    mx, my = st.mean(xs), st.mean(ys)
    n = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    d = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    return n / d if d else 0.0


def main():
    tasks, best, ck = load()
    recent = sorted(tasks, key=lambda t: tasks[t]["created_at"], reverse=True)[:LATEST]
    print(f"latest {len(recent)} three-seed rounds "
          f"({tasks[recent[-1]]['created_at'][:16]} .. {tasks[recent[0]]['created_at'][:16]})\n")
    print(f"{'task':<10} {'cell':<11} {'seeds':<20} {'gap':>5} {'miners':>7} {'>=0.60':>7} "
          f"{'>=0.45':>7} {'max':>6} {'ck spiking':>11}")
    rows = []
    for t in recent:
        info = tasks[t]
        rs = best.get(t, {})
        if len(rs) < 30:
            continue
        cons = {hk: r["breakdown"]["consistency_factor"] for hk, r in rs.items()}
        hi = [hk for hk, c in cons.items() if c >= HIGH]
        mid = [hk for hk, c in cons.items() if c >= MID]
        groups = Counter(ck.get(hk, "?") for hk in hi)
        multi = sum(1 for _c, n in groups.items() if n >= 2)
        rows.append({"task_id": t, **info, "miners": len(rs), "n_high": len(hi),
                     "n_mid": len(mid), "max_cons": max(cons.values()),
                     "coldkeys_high": {c: n for c, n in groups.items()},
                     "multi_coldkeys": multi,
                     "max_same_coldkey": max(groups.values()) if groups else 0})
        print(f"{t[:8]:<10} {info['cell_type']:<11} "
              f"{str(info['seeds']):<20} {info['min_gap']:>5} {len(rs):>7} {len(hi):>7} "
              f"{len(mid):>7} {max(cons.values()):>6.3f} "
              f"{(max(groups.values()) if groups else 0):>11}")
    json.dump(rows, open(OUT, "w"), indent=1)

    g = [r["min_gap"] for r in rows]
    print(f"\ncorrelation over these {len(rows)} rounds:")
    print(f"  min_gap vs count(cons>=0.60)  r = {corr(g, [r['n_high'] for r in rows]):+.3f}")
    print(f"  min_gap vs count(cons>=0.45)  r = {corr(g, [r['n_mid'] for r in rows]):+.3f}")
    print(f"  min_gap vs max consistency    r = {corr(g, [r['max_cons'] for r in rows]):+.3f}")

    print(f"\n{'min gap':<12} {'rounds':>7} {'mean >=0.60':>12} {'mean >=0.45':>12} "
          f"{'mean max cons':>14} {'max same coldkey':>17}")
    for lo, hi_ in ((0, 50), (50, 100), (100, 200), (200, 300), (300, 10000)):
        sel = [r for r in rows if lo <= r["min_gap"] < hi_]
        if not sel:
            continue
        label = f"{lo}-{hi_}" if hi_ < 10000 else f"{lo}+"
        print(f"{label:<12} {len(sel):>7} "
              f"{st.mean(r['n_high'] for r in sel):>12.2f} "
              f"{st.mean(r['n_mid'] for r in sel):>12.2f} "
              f"{st.mean(r['max_cons'] for r in sel):>14.3f} "
              f"{max(r['max_same_coldkey'] for r in sel):>17}")

    print("\nper cell type:")
    for cell in sorted({r["cell_type"] for r in rows}):
        sel = [r for r in rows if r["cell_type"] == cell]
        near = [r for r in sel if r["min_gap"] < 150]
        far = [r for r in sel if r["min_gap"] >= 150]
        print(f"  {cell:<11} n={len(sel):>2}  r(gap, >=0.60) = "
              f"{corr([r['min_gap'] for r in sel], [r['n_high'] for r in sel]):+.3f}"
              f"   gap<150: {st.mean(r['n_high'] for r in near) if near else 0:.2f} spikes/round"
              f"   gap>=150: {st.mean(r['n_high'] for r in far) if far else 0:.2f}")


if __name__ == "__main__":
    main()
