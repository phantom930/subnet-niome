#!/usr/bin/env python3
"""conj_price.py — price a conjunction sweep's arms on the real rank curve of their cell type.

E[cons] ranks the arms wrongly, because ``SCORING_SYSTEM = "top"`` is a step function: a band hit
is worth far more than its share of the mean, and a floor round is worth nothing however high the
floor. So each arm is treated as a four-regime lottery over the round's three seeds —

    band          len(band) seeds of 900        consistency 1.000
    in-window clean   clean - len(band)         consistency cons_clean
    in-window dirty   union                     consistency cons_dirty
    off-window        900 - span                consistency cons_off

— every composition of the three seeds is enumerated with its multinomial weight, and the resulting
round final is placed in EVERY real three-seed field of that cell type (our own hotkeys excluded).

    python conj_price.py <arms.json> [more.json ...]
"""
import json
import sys
from collections import defaultdict
from itertools import product
from math import factorial

DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]


def field_of_task(task_id):
    """The one real field that played THIS contract — the matched comparison.

    Pooling fields across contracts measures field softness, not the construction: the K562
    rank-10 cutoff swings 45.9-135.1, and a score measured on a high-weighted contract placed
    against another contract's softer field scores a placement it never earned.
    """
    from sd_task import OURS
    rows = json.load(open("sd_task_scores.json"))
    rows = rows if isinstance(rows, list) else (rows.get("data") or rows.get("items") or [])
    best = {}
    for r in rows:
        if r.get("task_id") == task_id and r.get("miner_hotkey") not in OURS:
            h = r["miner_hotkey"]
            best[h] = max(float(r.get("final_score") or 0), best.get(h, 0.0))
    return sorted(best.values(), reverse=True)


def fields_for(cell):
    from sd_task import OURS
    rows = json.load(open("sd_task_scores.json"))
    rows = rows if isinstance(rows, list) else (rows.get("data") or rows.get("items") or [])
    listing = json.load(open("sd_task_listing.json"))
    items = listing if isinstance(listing, list) else (listing.get("items") or [])
    keep = set()
    for t in items:
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") == cell and len(
                [s for s in str(c.get("seed", "")).split(",") if s.strip()]) == 3:
            keep.add(t["id"])
    best = defaultdict(dict)
    for r in rows:
        if r.get("task_id") in keep and r.get("miner_hotkey") not in OURS:
            t, h = r["task_id"], r["miner_hotkey"]
            best[t][h] = max(float(r.get("final_score") or 0), best[t].get(h, 0.0))
    return [sorted(p.values(), reverse=True) for p in best.values() if len(p) >= 50]


def share(sc, field):
    rank = sum(1 for f in field if f > sc) + 1
    return DIST[rank - 1] if rank <= len(DIST) else 0.0


def lottery(counts, values):
    """(probability, round consistency) over every composition of 3 seeds across the regimes.

    Generic in the number of regimes: a whole-window build has no off-window regime and a k=0 arm
    has no band, so callers drop empty ones. Enumerating a fixed four parts against a shorter
    ``counts`` would let ``zip`` silently discard a part, and the probabilities would not sum to 1.
    """
    tot = float(sum(counts))
    ps = [n / tot for n in counts]
    out = []
    for ks in product(range(4), repeat=len(counts)):
        if sum(ks) != 3:
            continue
        p = factorial(3)
        for k, pi in zip(ks, ps):
            p *= pi ** k / factorial(k)
        if p > 0:
            out.append((p, sum(k * v for k, v in zip(ks, values)) / 3.0))
    assert abs(sum(p for p, _c in out) - 1.0) < 1e-9, "regime lottery must be a distribution"
    return out


def arm_rows(path):
    d = json.load(open(path))
    span_lo, span_hi = (d.get("cut_window") or [100, 999])
    span = span_hi - span_lo + 1
    off = 900 - span
    for k, a in sorted((d.get("arms") or {}).items(), key=lambda kv: int(kv[0])):
        if "clean" not in a or a.get("weighted") is None:
            continue
        nb = len(a.get("band") or [])
        counts = (nb, a["clean"] - nb, a["union"], off)
        values = (a.get("band_cons") or 0.0, a["clean_cons"], a["dirty_cons"],
                  a.get("offwindow_cons") or 0.0)
        counts, values = zip(*[(n, v) for n, v in zip(counts, values) if n > 0])
        yield (f"{path.split('.')[0]} k={k} (g{d.get('group')}, w{span})", d["cell"],
               lottery(counts, values), a["weighted"], a["fidelity"], d["task"])


def main(paths):
    own = "--own" in paths
    paths = [p for p in paths if not p.startswith("--")]
    cache = {}
    print(f"  {'construction':<44} {'E[final]':>9} {'P(place)':>9} {'E[share]':>9} "
          f"{'best case':>10}")
    best = (None, -1.0)
    for path in paths:
        for name, cell, lot, wtd, fid, tid in arm_rows(path):
            flds = ([field_of_task(tid)] if own
                    else cache.setdefault(cell, fields_for(cell)))
            flds = [f for f in flds if len(f) >= 50]
            if not flds:
                print(f"  {name:<44}  no scored three-seed {cell} fields")
                continue
            e_final = sum(p * wtd * c * fid for p, c in lot)
            e_share = sum(p * sum(share(wtd * c * fid, f) for f in flds) / len(flds)
                          for p, c in lot)
            p_place = sum(p * sum(1 for f in flds if share(wtd * c * fid, f) > 0) / len(flds)
                          for p, c in lot)
            top = max(wtd * c * fid for _p, c in lot)
            print(f"  {name:<44} {e_final:>9.2f} {p_place:>8.1%} {e_share:>9.4f} {top:>10.1f}")
            if e_share > best[1]:
                best = (name, e_share)
    if best[0]:
        print(f"\n  best on expected share: {best[0]}  ({best[1]:.4f})")
    if own:
        print("  each arm priced against the ONE field that played its own contract.")
    else:
        print(f"  priced on pooled fields: "
              + ", ".join(f"{c} {len(f)}" for c, f in sorted(cache.items()))
              + " — cross-contract, so it measures field softness too; prefer --own.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["conj_test.json"]))
