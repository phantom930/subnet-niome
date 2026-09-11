#!/usr/bin/env python3
"""cell_grid_price.py — price every all-HDR mixed config in the grid against that cell's real fields.

`w x fid` alone picks the wrong config: band sets how OFTEN the 1.000 spike lands, and the
band-maximising corner of the grid (small group, narrow window) is the opposite of the
score-maximising one. all-HDR mixed has no elevated floor — every non-band seed sits at the
arithmetic ~0.10 — so a round's consistency is `(k + (3-k)*floor)/3` with k the number of the
round's three seeds inside the band, and expected curve share is an EXACT function of (band, w x fid).
Enumerating the four outcomes k=0..3 gives the expectation with no Monte-Carlo error.

The floor is per cell: 0.101 on the erythroid types, 0.079 on HEK293 (CLAUDE.md, and reproduced in
this session's builds).
"""
import json
import os
import sys
import urllib.request
from collections import defaultdict
from math import comb

DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
SP = os.getenv("CGP_SCORES", "/tmp/claude-0/-root-workspace-subnet-niome/"
                             "d5c07231-737e-45ff-8fd2-aa06c458d312/scratchpad/scores2.json")
OURS = os.getenv("CGP_OURS", "/tmp/claude-0/-root-workspace-subnet-niome/"
                             "d5c07231-737e-45ff-8fd2-aa06c458d312/scratchpad/ours.json")
FLOOR = {"K562": 0.101, "HUDEP-2": 0.101, "CD34+_HSPC": 0.101, "HEK293": 0.079}
FILES = {"K562": "cell_grid_K562.json", "HUDEP-2": "cell_grid_HUDEP_2.json",
         "CD34+_HSPC": "cell_grid_CD34_HSPC.json", "HEK293": "cell_grid_HEK293.json"}
NSEED = 900


def fields_by_cell():
    sc = json.load(open(SP))
    sc = sc if isinstance(sc, list) else (sc.get("items") or sc.get("data"))
    tk = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks?limit=500"))
    tk = tk if isinstance(tk, list) else (tk.get("items") or tk.get("data"))
    meta = {}
    for t in tk:
        tid = t.get("task_id") or t.get("id")
        c = (t.get("content") or {}).get("contract", {})
        s = str(c.get("seed", "") or "")
        meta[tid] = (c.get("cell_type"),
                     len([x for x in s.split(",") if x.strip().isdigit()]))
    ours = set(json.load(open(OURS)).values())
    best = defaultdict(dict)
    for r in sc:
        cell, n = meta.get(r["task_id"], (None, 0))
        if n != 3 or r["miner_hotkey"] in ours:
            continue
        hk = r["miner_hotkey"]
        f = float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    out = defaultdict(list)
    for tid, v in best.items():
        fs = sorted(v.values(), reverse=True)
        if len(fs) >= 10 and fs[0] > 0:
            out[meta[tid][0]].append(fs)
    return out


def share(final, field):
    rank = sum(1 for f in field if f > final) + 1
    return DIST[rank - 1] if rank <= 10 else 0.0


def price(band, wxf, floor, fs):
    p = band / NSEED
    tot = 0.0
    for k in range(4):
        w = comb(3, k) * p ** k * (1 - p) ** (3 - k)
        if w <= 0:
            continue
        cons = (k + (3 - k) * floor) / 3.0
        f = wxf * cons
        tot += w * (sum(share(f, x) for x in fs) / len(fs))
    return tot


def main():
    fbc = fields_by_cell()
    for cell, path in FILES.items():
        fs = fbc.get(cell, [])
        if not fs:
            print(f"{cell}: no fields"); continue
        d = [r for r in json.load(open(path)) if "arms" in r]
        rows = []
        for r in d:
            for light, a in r["arms"].items():
                e = price(r["band"], a["wxfid"], FLOOR[cell], fs)
                rows.append((e, r["width"], r["group"], light, r["band"], a["wxfid"], r["mf"]))
        rows.sort(reverse=True)
        print(f"\n{cell}  ({len(fs)} real 3-seed fields, floor {FLOOR[cell]})")
        print(f"  {'rank':>4}{'E[share]':>10}{'width':>7}{'group':>7}{'light':>7}"
              f"{'band':>6}{'w x fid':>9}{'mf':>5}")
        for i, (e, w, g, l, b, wf, mf) in enumerate(rows[:5], 1):
            print(f"  {i:>4}{e:>10.5f}{w:>7}{g:>7}{l:>7}{b:>6}{wf:>9.1f}{mf:>5}")
        # where the w x fid winner lands
        bw = max(rows, key=lambda x: x[5])
        bi = rows.index(bw) + 1
        print(f"  w x fid winner (w{bw[1]} g{bw[2]} L{bw[3]}, {bw[5]:.1f}) ranks #{bi} "
              f"on payout at {bw[0]:.5f}")
        best = rows[0]
        print(f"  => BEST: width {best[1]}, group {best[2]}, light {best[3]}, mf {best[6]}  "
              f"(band {best[4]}, w x fid {best[5]:.1f})")


if __name__ == "__main__":
    main()
