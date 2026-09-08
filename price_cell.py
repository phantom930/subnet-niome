#!/usr/bin/env python3
"""price_cell.py — per cell type, does all-cut or all-HDR earn more? Priced on real rank curves.

E[cons] is NOT the objective. ``SCORING_SYSTEM = "top"`` pays only the top 10 on a fixed curve, so
a construction with a better mean consistency can still earn nothing: what matters is how often its
round score crosses that round's rank-10 cutoff, and how high it lands when it does. all-HDR is
shipped today on exactly that reasoning — it loses ~24/round on the mean and is kept because its
spike can out-RANK a flat score — but CLAUDE.md records that no rank measurement in the repo
supported it. This is that measurement.

Both constructions are three-regime lotteries over a round's three seeds, so the round score is a
small discrete distribution rather than a point:

    all-HDR   each seed lands in the band with p = band/900; cons 1.000 there, floor elsewhere
    all-cut   each seed is cut-clean with p = clean/900;      cons_clean there, cons_dirty elsewhere

    round consistency = mean of the three seeds' values     (stage 12/3-5 average over seeds)
    round final       = total_weighted_score x that x fidelity

Each possible round score is then placed in EVERY real field of that cell type — the actual scored
finals of the other 200+ miners, our own hotkeys excluded — and paid at SCORE_DISTRIBUTION[rank-1].
Averaging over both the lottery and the fields gives expected curve share per round, which is the
quantity the fleet is actually optimising.

Construction parameters come from measurement, never from the module defaults; the sources are named
in PARAMS below.

    python price_cell.py
"""
import glob
import json
import os
import sys
from collections import defaultdict
from itertools import product

DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
OURS_FILE = "sd_task_scores.json"
LISTING = "sd_task_listing.json"

# all-HDR, from CLAUDE.md's measured CELL_CONFIG sweeps and this session's off-band sampling:
#   band      12-16 on the erythroid types, 7 on HEK293 (group 80)
#   spike     consistency exactly 1.000 on a band seed
#   floor     0.1256 measured over 24 uniform off-band seeds (hybrid_floor.py)
#   weighted  the head-to-head legs in CLAUDE.md (230-370 by cell); fidelity from the group sweep
ALL_HDR = {
    "K562":       {"band": 13, "weighted": 258.0, "fid": 0.96, "floor": 0.1256},
    "HUDEP-2":    {"band": 13, "weighted": 240.0, "fid": 0.96, "floor": 0.1256},
    "CD34+_HSPC": {"band": 13, "weighted": 245.0, "fid": 0.96, "floor": 0.1256},
    "HEK293":     {"band": 7,  "weighted": 318.0, "fid": 0.924, "floor": 0.1256},
}


def load_allcut():
    """all-cut parameters measured tonight by conj_test.py / conj300.py at k=0 (pure all-cut)."""
    out = {}
    for path in glob.glob("conj_k0_*.json") + ["conj_test.json", "conj_ctrl_634a512c.json"]:
        if not os.path.exists(path):
            continue
        try:
            d = json.load(open(path))
        except json.JSONDecodeError:
            continue
        arm = (d.get("arms") or {}).get("0") or (d.get("arms") or {}).get(0)
        if not arm or "clean" not in arm:
            continue
        cell = d["cell"]
        # Prefer a whole-window (100-999) measurement; a narrow window is a different construction.
        if d.get("cut_window") and d["cut_window"] != [100, 999]:
            continue
        out[cell] = {"clean": arm["clean"], "cons_clean": arm["clean_cons"],
                     "cons_dirty": arm["dirty_cons"], "weighted": arm["weighted"],
                     "fid": arm["fidelity"], "source": path}
    return out


def fields_by_cell():
    """Real scored finals per three-seed round, excluding our own hotkeys, keyed by cell type."""
    from sd_task import OURS
    rows = json.load(open(OURS_FILE))
    rows = rows if isinstance(rows, list) else (rows.get("data") or rows.get("items") or [])
    listing = json.load(open(LISTING))
    items = listing if isinstance(listing, list) else (listing.get("items") or [])
    cell_of, three = {}, set()
    for t in items:
        c = (t.get("content") or {}).get("contract") or {}
        cell_of[t["id"]] = c.get("cell_type")
        if len([s for s in str(c.get("seed", "")).split(",") if s.strip()]) == 3:
            three.add(t["id"])
    best = defaultdict(dict)                     # task -> hotkey -> final (dedupe validators)
    for r in rows:
        if r.get("task_id") not in three or r.get("miner_hotkey") in OURS:
            continue
        t, h = r["task_id"], r["miner_hotkey"]
        best[t][h] = max(float(r.get("final_score") or 0), best[t].get(h, 0.0))
    out = defaultdict(list)
    for t, per in best.items():
        if len(per) >= 50:
            out[cell_of.get(t)].append(sorted(per.values(), reverse=True))
    return out


def share(score, field):
    """SCORE_DISTRIBUTION share a given round score would earn in one real field."""
    rank = sum(1 for f in field if f > score) + 1
    return DIST[rank - 1] if rank <= len(DIST) else 0.0


def lottery(p_hit, hi, lo):
    """(probability, round consistency) over how many of three seeds land in the good regime."""
    from math import comb
    return [(comb(3, k) * p_hit ** k * (1 - p_hit) ** (3 - k), (k * hi + (3 - k) * lo) / 3.0)
            for k in range(4)]


def main():
    allcut = load_allcut()
    fields = fields_by_cell()
    if not allcut:
        raise SystemExit("no whole-window all-cut measurements found (conj_k0_*.json)")

    print("expected curve share per round, priced on real fields "
          f"({sum(len(v) for v in fields.values())} three-seed rounds)\n")
    print(f"  {'cell':<11} {'rounds':>6} {'construction':<9} {'E[final]':>9} {'P(place)':>9} "
          f"{'E[share]':>9}  {'best case':>10}")
    summary = {}
    for cell in sorted(set(allcut) | set(ALL_HDR)):
        flds = fields.get(cell) or []
        if not flds:
            print(f"  {cell:<11} {'0':>6}  no scored three-seed fields")
            continue
        arms = {}
        if cell in allcut:
            a = allcut[cell]
            arms["all-cut"] = (lottery(a["clean"] / 900.0, a["cons_clean"], a["cons_dirty"]),
                               a["weighted"], a["fid"])
        if cell in ALL_HDR:
            h = ALL_HDR[cell]
            arms["all-HDR"] = (lottery(h["band"] / 900.0, 1.0, h["floor"]),
                               h["weighted"], h["fid"])
        row = {}
        for name, (lot, wtd, fid) in arms.items():
            e_final = sum(p * wtd * c * fid for p, c in lot)
            e_share = 0.0
            p_place = 0.0
            for p, c in lot:
                sc = wtd * c * fid
                s = sum(share(sc, f) for f in flds) / len(flds)
                placed = sum(1 for f in flds if share(sc, f) > 0) / len(flds)
                e_share += p * s
                p_place += p * placed
            best = max(wtd * c * fid for _p, c in lot)
            row[name] = {"e_final": e_final, "e_share": e_share, "p_place": p_place,
                         "best_case": best}
            print(f"  {cell:<11} {len(flds):>6} {name:<9} {e_final:>9.2f} {p_place:>8.1%} "
                  f"{e_share:>9.4f}  {best:>10.1f}")
        if len(row) == 2:
            a, h = row["all-cut"]["e_share"], row["all-HDR"]["e_share"]
            better = "all-cut" if a > h else "all-HDR"
            factor = (max(a, h) / min(a, h)) if min(a, h) > 0 else float("inf")
            print(f"  {'':<11} {'':>6} -> {better} by {factor:.2f}x on expected share")
        summary[cell] = row
        print()

    print("  P(place) is the chance of reaching the top 10 at all; E[share] weights rank as well.")
    print("  all-cut params are measured (conj_*.json); all-HDR params come from CLAUDE.md's")
    print("  CELL_CONFIG sweeps plus this session's off-band floor of 0.1256.")
    json.dump(summary, open("price_cell.json", "w"), indent=1)
    print("\nwrote price_cell.json")


if __name__ == "__main__":
    sys.exit(main())
