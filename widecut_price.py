#!/usr/bin/env python3
"""widecut_price.py — is the 900-seed CUT worth keeping on the erythroid cells?

`conj_wall.py` settled the feasibility half on 2026-09-17: the wide cut lowers HEK293's band wall
9 -> 8 (so HEK293 was excluded from `joined_window.WIDE_CUT_CELLS`) and leaves CD34+/K562/HUDEP-2 at
12, where the live `band_k` of 11 sits one step below either way. That says the wide cut costs those
three cells nothing on AVAILABILITY. It does not say the wide cut EARNS anything, and nothing has
measured that: every k-arm valuation in CLAUDE.md's band-depth table was priced at the narrow cut.

So this is the paired comparison the file flags as open, at the live arm only:

    narrow  cut min-union over the 300-seed band space, band carved by `sub_window` (pre-09-17)
    wide    cut min-union over the full 900, band from explicit `band_candidates` -- SHIPS today

Everything else is the cell's live `CELL_CONFIG`: k=11, group 100, width 150, light 6, cell-aware.
k is deliberately NOT swept here. The wall sweep showed it does not move on these cells, so the live
value is not in doubt, and holding it fixed makes this one question instead of two.

**Pricing method is `hud_resolve.py`'s, not `band_hit.py`'s, and the reason is not cosmetic.**
`band_hit.py` derives its regimes from the config's seed space -- it samples "dirty" as
`SPACE - joined` and reads the cut-clean count out of `meta["clean"]`, which counts within the CUT
space. With a 900-seed cut there is no "outside the cut space" and `meta["clean"]` is a count over
900 rather than over 300, so those two regimes stop denoting the same thing across the two arms and
the comparison would be biased by construction. `hud_resolve.evaluate` instead enumerates both sets
from the SHIPPED ROWS over all 900 -- band = every seed where every row satisfies `hdr`, cut-clean =
the same on `cut` -- and scores each regime at its measured consistency. That is identical work for
both arms regardless of what seed space built them, which is the only way the pairing is honest.

Also recorded per contract: the band seeds each arm actually chose. CLAUDE.md currently asserts the
two arms "produce the identical band" -- true of the 150-seed CANDIDATE set, but `choose_band` picks
from the guides that survive the cut min-union, and that pool differs between arms. Whether the
chosen bands coincide is therefore a measurement, not a given, and it is printed as `band overlap`.

    WP_CELLS=K562,CD34+_HSPC WP_N=4 python widecut_price.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "wcprice")

import dataclasses
import json
import logging
import random
import time
from collections import defaultdict
from math import factorial

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from niome_subnet.genomics import mt19937 as MT          # noqa: E402
from conj_stageb import DIST, OURS, API                  # noqa: E402
from sd_task import score, fetch                         # noqa: E402

CELLS = [c for c in os.getenv("WP_CELLS", "K562,CD34+_HSPC").split(",") if c]
N_CONTRACTS = int(os.getenv("WP_N", "4"))
NS = int(os.getenv("WP_NS", "8"))
HK = os.getenv("WP_HK", "niome_hotkey")
BUDGET = float(os.getenv("WP_BUDGET", "1800"))
SPACE = list(range(100, 1000))
OUT = os.getenv("WP_JSON", "widecut_price.json")


def records_of(rows):
    return [{"guide": r["guideRNA"], "start": int(r["target_alignment_start"]),
             "length": int(r["target_alignment_end"]) - int(r["target_alignment_start"]),
             "strand": r["strand"], "mutation": r["mutation"],
             "cas_system": r["cas_system"]} for r in rows]


def compliant(rows, contract, cell_types, ctx, rule):
    """Every seed of 100-999 on which EVERY row satisfies `rule`."""
    ok = CJ.hdr_compliance(records_of(rows), contract, cell_types, ctx, SPACE, rule)
    MT.free_gpu_memory()
    if not ok:
        return []
    keep = set(SPACE)
    for i in sorted(ok):
        keep &= ok[i]
        if not keep:
            break
    return sorted(keep)


def price(nb, nc, v_band, v_clean, v_rest, wxf, field):
    """E[own-field curve share] over a 3-seed draw at the measured per-regime values."""
    pb, pc = nb / 900.0, max(0, nc) / 900.0
    pr = max(0.0, 1.0 - pb - pc)
    tot = 0.0
    for j in range(4):
        for i in range(4 - j):
            r = 3 - j - i
            w = (factorial(3) / (factorial(j) * factorial(i) * factorial(r))
                 * pb ** j * pc ** i * pr ** r)
            if w <= 0:
                continue
            final = wxf * (j * v_band + i * v_clean + r * v_rest) / 3.0
            rank = 1 + sum(1 for x in field if x > final)
            if rank <= 10:
                tot += w * DIST[rank - 1]
    return tot


def own_fields(cell):
    """rank-ordered final_scores per task for this cell, our own hotkeys excluded by ss58."""
    sc = fetch(f"{API}/miners/scores?limit=40000", cache="sd_task_scores.json")
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    tk = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    tk = tk if isinstance(tk, list) else (tk.get("items") or tk.get("data") or [])
    meta = {}
    for t in tk:
        c = (t.get("content") or {}).get("contract", {})
        s = str(c.get("seed", "") or "")
        meta[t.get("task_id") or t.get("id")] = (
            c.get("cell_type"), len([x for x in s.split(",") if x.strip().isdigit()]))
    best = defaultdict(dict)
    for r in sc:
        c, n = meta.get(r.get("task_id"), (None, 0))
        if c != cell or n != 3 or r.get("miner_hotkey") in OURS:
            continue
        hk, f = r["miner_hotkey"], float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    return {tid: sorted(v.values(), reverse=True)
            for tid, v in best.items() if len(v) >= 10 and max(v.values()) > 0}


def evaluate(rows, contract, reference, cell_types, ctx, field, label, t0):
    """hud_resolve.py's method verbatim: both sets enumerated over 900 from the shipped rows."""
    band = compliant(rows, contract, cell_types, ctx, "hdr")
    cut = compliant(rows, contract, cell_types, ctx, "cut")
    if not band:
        print(f"    {label}: no hdr-clean seed over 900 — unpriceable", flush=True)
        return None
    cutonly = [s for s in cut if s not in set(band)]
    sb = score(rows, contract, reference, cell_types, seed=band[0])
    wxf = sb["weighted"] * sb["fidelity"]
    rng = random.Random(0)
    vc = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(cutonly, min(NS, len(cutonly)))] if cutonly else []
    rest = [s for s in SPACE if s not in set(cut)]
    vr = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(rest, min(NS, len(rest)))]
    v_band = sb["consistency"]
    v_clean = float(np.mean(vc)) if vc else float(np.mean(vr))
    v_rest = float(np.mean(vr))
    e = price(len(band), len(cutonly), v_band, v_clean, v_rest, wxf, field)
    print(f"    {label:6s} band {len(band):3d} cut-clean {len(cut):3d} (+{len(cutonly):3d}) "
          f"wxf {wxf:6.1f} | v_band {v_band:.4f} v_clean {v_clean:.4f} v_rest {v_rest:.4f} "
          f"| E[own] {e:.6f}  ({time.monotonic()-t0:.0f}s)", flush=True)
    return {"label": label, "band": len(band), "cut": len(cut), "cutonly": len(cutonly),
            "wxf": wxf, "v_band": v_band, "v_clean": v_clean, "v_rest": v_rest, "share": e}


def arm_cfg(base, band_space, arm):
    """(cfg, shape) for one arm at the cell's LIVE band_k — k is not swept here."""
    if arm == "wide":
        cands = CJ.sub_window(band_space, base.band_width, JW.band_offset_frac(HK) or 0.0)
        return (dataclasses.replace(base, seed_list=tuple(SPACE),
                                    band_candidates=tuple(cands)),
                f"cut 900 / band cands {len(cands)}")
    return (dataclasses.replace(base, seed_list=tuple(band_space), band_candidates=()),
            f"cut {len(band_space)} / sub_window {base.band_width}")


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])

    recs, agg = [], defaultdict(lambda: defaultdict(list))
    for cell in CELLS:
        base = CJ.config_for(cell)
        if base is None:
            print(f"=== {cell}: no conjunction config, skipped ===", flush=True)
            continue
        ftasks = own_fields(cell)
        band_space = sorted({s for a, b in JW.fallback_for(cell, HK) for s in range(a, b + 1)})
        cands = []
        for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
            c = (t.get("content") or {}).get("contract") or {}
            s = str(c.get("seed", "") or "")
            if c.get("cell_type") != cell:
                continue
            if len([x for x in s.split(",") if x.strip().isdigit()]) != 3:
                continue
            if (t.get("task_id") or t["id"]) not in ftasks:
                continue
            cands.append(t)
            if len(cands) >= N_CONTRACTS:
                break
        print(f"\n=== {cell}  live k={base.band_k} group={base.group_size} "
              f"width={base.band_width} light={base.light_cell_rows} ca={base.band_cell_aware}"
              f" | band space {len(band_space)} | {len(cands)} contracts with own fields ===",
              flush=True)

        for t in cands:
            tid = (t.get("task_id") or t["id"])
            contract = dict(t["content"]["contract"])
            reference = t["content"]["hbb_reference"]
            field = ftasks[tid]
            print(f"  {tid[:8]}  own cut10 {field[9]:.1f}  top {field[0]:.1f}", flush=True)
            ctx = G.build_context(contract, reference, cell_types)
            bands = {}
            for arm in ("narrow", "wide"):
                cfg, shape = arm_cfg(base, band_space, arm)
                t0 = time.monotonic()
                rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                                 budget_s=BUDGET)
                if not rows:
                    print(f"    {arm:6s} DECLINED — {meta.get('reason')} "
                          f"({time.monotonic()-t0:.0f}s)", flush=True)
                    agg[cell][arm].append(0.0)
                    recs.append({"cell": cell, "task": tid[:8], "arm": arm, "built": False,
                                 "reason": meta.get("reason"), "share": 0.0})
                    continue
                bands[arm] = sorted(int(x) for x in meta["band_seeds"])
                r = evaluate(rows, contract, reference, cell_types, ctx, field, arm, t0)
                agg[cell][arm].append(r["share"] if r else 0.0)
                recs.append(dict(r or {}, cell=cell, task=tid[:8], arm=arm, built=True,
                                 shape=shape, pinned=len(bands[arm]),
                                 cas9=meta.get("cas9_pool"), bank=meta.get("bank"),
                                 clean_in_cut=meta.get("clean"),
                                 band_seeds=bands[arm], cut10=field[9]))
            if len(bands) == 2:
                a, b = set(bands["narrow"]), set(bands["wide"])
                print(f"    band overlap {len(a & b)} of {len(a)}/{len(b)} "
                      f"{'IDENTICAL' if a == b else 'DIFFERENT'}", flush=True)
                recs.append({"cell": cell, "task": tid[:8], "arm": "_overlap",
                             "overlap": len(a & b), "n_narrow": len(a), "n_wide": len(b),
                             "identical": a == b})

    with open(OUT, "w") as fh:
        json.dump(recs, fh, indent=1)

    print("\n=== E[own-field share], paired within contract ===")
    print(f"{'cell':12s} {'n':>3s} {'narrow':>10s} {'wide':>10s} {'wide/narrow':>12s} "
          f"{'wide wins':>10s}")
    for cell in CELLS:
        nar, wid = agg[cell]["narrow"], agg[cell]["wide"]
        if not nar or not wid:
            continue
        n = min(len(nar), len(wid))
        mn, mw = float(np.mean(nar[:n])), float(np.mean(wid[:n]))
        wins = sum(1 for i in range(n) if wid[i] > nar[i])
        ties = sum(1 for i in range(n) if wid[i] == nar[i])
        ratio = (mw / mn) if mn > 0 else float("nan")
        print(f"{cell:12s} {n:3d} {mn:10.6f} {mw:10.6f} {ratio:11.2f}x "
              f"{wins:4d}/{n} ({ties} tie)")
    print("\nPaired per contract (wide - narrow):")
    for cell in CELLS:
        nar, wid = agg[cell]["narrow"], agg[cell]["wide"]
        n = min(len(nar), len(wid))
        if n:
            print(f"  {cell:12s} " + "  ".join(f"{wid[i]-nar[i]:+.6f}" for i in range(n)))


if __name__ == "__main__":
    main()
