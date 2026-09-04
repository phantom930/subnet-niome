#!/usr/bin/env python3
"""mingap_cells.py — the gate per cell type, and the window layout the spikes imply.

mingap_detail.py showed the leaders' windows span 200-500 (median 242, 14 of 15 above 100) and
overlap within a coldkey (11 rounds with 2+ simultaneous k=2 spikes, one with 3). This resolves the
gate per cell type -- HEK293 needs a far smaller gap than the erythroid types -- and estimates the
stride, which is what a tiling would have to copy.

Stride estimate: with windows of width W at stride S, an interval of length g is covered by about
(W - g)/S + 1 of a coldkey's windows. Only hotkeys whose *band* holds both seeds actually spike, so
the observed count is a lower bound on coverage and the fitted S is an upper bound. Reported as a
range rather than a point for that reason.
"""
import json
import statistics as st
import urllib.request
from collections import Counter, defaultdict

HIGH = 0.60


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
                              "min_gap": min(seeds[1] - seeds[0], seeds[2] - seeds[1])}
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
    return tasks, best, json.load(open("hotkey_coldkey.json"))


def main():
    tasks, best, ck = load()
    BUCKETS = [(0, 25), (25, 50), (50, 100), (100, 150), (150, 200), (200, 300), (300, 10000)]

    print("=== k=2 incidence (cons >= 0.60) by min seed gap, per cell type ===")
    print(f"{'gap':<10} " + "".join(f"{c[:10]:>22}" for c in
                                    ["CD34+_HSPC", "HUDEP-2", "K562", "HEK293"]))
    print(f"{'':<10} " + "".join(f"{'rounds  spikes  max':>22}" for _ in range(4)))
    for lo, hi in BUCKETS:
        line = f"{f'{lo}-{hi}' if hi < 10000 else f'{lo}+':<10} "
        for cell in ["CD34+_HSPC", "HUDEP-2", "K562", "HEK293"]:
            ts = [t for t, i in tasks.items()
                  if i["cell_type"] == cell and lo <= i["min_gap"] < hi and len(best.get(t, {})) > 30]
            if not ts:
                line += f"{'-':>22}"
                continue
            sp = [sum(1 for r in best[t].values()
                      if r["breakdown"]["consistency_factor"] >= HIGH) for t in ts]
            mx = max(max(r["breakdown"]["consistency_factor"] for r in best[t].values())
                     for t in ts)
            line += f"{len(ts):>10}{st.mean(sp):>7.1f}{mx:>7.3f}"
        print(line)

    print("\n=== the gap at which each cell type's k=2 regime switches on ===")
    for cell in ["CD34+_HSPC", "HUDEP-2", "K562", "HEK293"]:
        ts = sorted((i["min_gap"], t) for t, i in tasks.items()
                    if i["cell_type"] == cell and len(best.get(t, {})) > 30)
        hit = [(g, sum(1 for r in best[t].values()
                       if r["breakdown"]["consistency_factor"] >= HIGH)) for g, t in ts]
        with_spike = [g for g, n in hit if n > 0]
        without = [g for g, n in hit if n == 0]
        print(f"  {cell:<11} rounds {len(hit):>2} | spikes on {len(with_spike):>2} | "
              f"largest gap WITH a spike {max(with_spike) if with_spike else '-':>4} | "
              f"smallest gap WITHOUT {min(without) if without else '-':>4}")

    print("\n=== coldkey scale vs spike behaviour ===")
    size = Counter(ck.values())
    agg = defaultdict(lambda: {"k2": 0, "multi": 0, "rounds": set()})
    for t, rs in best.items():
        by = defaultdict(int)
        for hk, r in rs.items():
            if r["breakdown"]["consistency_factor"] >= HIGH:
                by[ck.get(hk, "?")] += 1
        for c, n in by.items():
            agg[c]["k2"] += n
            agg[c]["rounds"].add(t)
            if n >= 2:
                agg[c]["multi"] += 1
    print(f"{'coldkey':<11} {'hotkeys':>8} {'k2 events':>10} {'rounds':>7} {'multi-spike':>12}")
    for c, d in sorted(agg.items(), key=lambda kv: -kv[1]["k2"])[:10]:
        print(f"{c[:9]:<11} {size.get(c, 0):>8} {d['k2']:>10} {len(d['rounds']):>7} "
              f"{d['multi']:>12}")

    print("\n=== stride implied by simultaneous spikes:  n ~ (W - gap)/S + 1 ===")
    print(f"{'coldkey':<11} {'hotkeys':>8} {'gap':>5} {'spiking':>8}  implied (W-gap)/S")
    obs = []
    for t, rs in best.items():
        by = defaultdict(int)
        for hk, r in rs.items():
            if r["breakdown"]["consistency_factor"] >= HIGH:
                by[ck.get(hk, "?")] += 1
        for c, n in by.items():
            if n >= 2:
                g = tasks[t]["min_gap"]
                obs.append((c, size.get(c, 0), g, n))
                print(f"{c[:9]:<11} {size.get(c,0):>8} {g:>5} {n:>8}  {n - 1}")
    if obs:
        # W from the measured span lower bounds; S from the counts at each gap
        for W in (250, 300, 400):
            fits = [(W - g) / (n - 1) for _c, _h, g, n in obs if n > 1 and W > g]
            if fits:
                print(f"  assuming W={W}: implied stride S median {st.median(fits):>5.0f} "
                      f"(range {min(fits):.0f}-{max(fits):.0f})  "
                      f"-> {int((900 - W) / st.median(fits)) + 1} hotkeys to tile 100-999")


if __name__ == "__main__":
    main()
