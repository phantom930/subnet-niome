#!/usr/bin/env python3
"""conj_price_k.py — price the joined-window conjunction's band width against real K562 fields.

The conjunction has three per-seed regimes, unlike every other construction measured here:

    band   (k seeds)        all three stage-4 targets pinned -> consistency 1.000
    clean  (`clean` seeds)  `is_cut` pinned                  -> ~0.24-0.29
    rest                    nothing pinned                   -> ~0.10

so a round's consistency is `(n_band + n_clean*v_clean + n_rest*v_rest) / 3` and the three seeds are
independent uniform draws. That is only 10 distinct (n_band, n_clean, n_rest) outcomes, so the
expectation is computed EXACTLY by enumeration rather than sampled — there is no Monte-Carlo error
in the numbers below, only the measurement error in `v_clean` (3 samples per arm).

Each arm is priced against every current-regime 3-seed K562 field, taking the real rank-10 cutoff
and the real `SCORE_DISTRIBUTION` curve position, one hotkey at a time. A single hotkey avoids the
question of how correlated sibling clean sets would be — they share the cut window, so a fleet
model would need that measured first.

    python conj_price_k.py
"""
import json
import os
import sys
import urllib.request
from collections import defaultdict
from math import factorial

DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
SP = os.getenv("CPK_SCORES", "/tmp/claude-0/-root-workspace-subnet-niome/"
                             "d5c07231-737e-45ff-8fd2-aa06c458d312/scratchpad/scores2.json")
OURS = os.getenv("CPK_OURS", "/tmp/claude-0/-root-workspace-subnet-niome/"
                             "d5c07231-737e-45ff-8fd2-aa06c458d312/scratchpad/ours.json")
CELL = os.getenv("CPK_CELL", "K562")
NSEED = 900


def fields():
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
        cell, nseeds = meta.get(r["task_id"], (None, 0))
        if cell != CELL or nseeds != 3 or r["miner_hotkey"] in ours:
            continue
        hk = r["miner_hotkey"]
        f = float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    out = []
    for tid, v in best.items():
        fs = sorted(v.values(), reverse=True)
        if len(fs) >= 10 and fs[0] > 0:
            out.append(fs)
    return out


def share(final, field):
    rank = sum(1 for f in field if f > final) + 1
    return DIST[rank - 1] if rank <= 10 else 0.0


def price(band, clean, v_clean, v_rest, wxf, fs):
    """Exact expected curve share: enumerate the 10 ways 3 seeds split across the regimes."""
    pb, pc = band / NSEED, clean / NSEED
    pr = 1.0 - pb - pc
    tot, p_place = 0.0, 0.0
    for nb in range(4):
        for nc in range(4 - nb):
            nr = 3 - nb - nc
            w = (factorial(3) / (factorial(nb) * factorial(nc) * factorial(nr))
                 * pb ** nb * pc ** nc * pr ** nr)
            if w <= 0:
                continue
            cons = (nb * 1.0 + nc * v_clean + nr * v_rest) / 3.0
            final = wxf * cons
            s = sum(share(final, f) for f in fs) / len(fs)
            tot += w * s
            p_place += w * (sum(1 for f in fs if share(final, f) > 0) / len(fs))
    return tot, p_place


def main():
    fs = fields()
    cuts = sorted(f[9] for f in fs)
    print(f"{CELL}: {len(fs)} current-regime 3-seed fields, ours excluded")
    print(f"rank-10 cutoff: min {cuts[0]:.1f}  p25 {cuts[len(cuts)//4]:.1f}  "
          f"median {cuts[len(cuts)//2]:.1f}  p75 {cuts[3*len(cuts)//4]:.1f}  max {cuts[-1]:.1f}\n")
    arms = []
    for f in ("conj_joined.json", "conj_joined_k.json", "conj_joined_k12.json"):
        try:
            arms += [r for r in json.load(open(f)) if r.get("light") == 6 and "wxfid" in r]
        except OSError:
            pass
    seen = set()
    rows = []
    for r in arms:
        key = (r["group"], r["k"])
        if key in seen:
            continue
        seen.add(key)
        e, p = price(r["k"], r["clean"], r["clean_cons"], r["dirty_cons"], r["wxfid"], fs)
        rows.append((f"conj g{r['group']} k{r['k']}", r["k"], r["clean"], r["clean_cons"],
                     r["wxfid"], e, p))
    # the shipped-style comparators: a pure band on a ~0.10 floor
    for lbl, band, wxf in (("all-HDR mixed w150 170/80", 12, 233.29),
                           ("all-HDR mixed w30 170/80", 15, 222.74),
                           ("all-HDR shipped w225", 12, 232.0)):
        e, p = price(band, 0, 0.0, 0.101, wxf, fs)
        rows.append((lbl, band, 0, 0.0, wxf, e, p))
    rows.sort(key=lambda x: -x[5])
    print(f"{'arm':<26}{'band':>5}{'clean':>7}{'v_clean':>9}{'w x fid':>9}"
          f"{'E[share]':>10}{'P(place)':>10}")
    for lbl, b, c, vc, w, e, p in rows:
        print(f"{lbl:<26}{b:>5}{c:>7}{vc:>9.3f}{w:>9.2f}{e:>10.5f}{p:>9.2%}")
    best = rows[0]
    base = next(r for r in rows if r[0].startswith("all-HDR mixed w150"))
    print(f"\nbest: {best[0]}  E[share] {best[5]:.5f}  "
          f"({best[5]/base[5]:.2f}x the all-HDR mixed w150 baseline)"
          if base[5] else "")


if __name__ == "__main__":
    main()
