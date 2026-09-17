#!/usr/bin/env python3
"""conj_wall.py — where the Cas12a band-formation wall sits under the 2026-09-17 wide-cut layout.

CLAUDE.md brackets the wall at 9 (HEK293) / 12 (erythroid), but every one of those measurements was
taken with the cut min-union over the same ~300-seed space the band came from. Since 2026-09-17 a
`joined_window.BAND_HK` hotkey min-unions the cut over the FULL 900 while the band is a 150-wide
slice of the plan's 300, and `CELL_CONFIG`'s `band_k` (8 / 11) was set one step below the OLD wall
by hand rather than re-measured. The wall is a property of the Cas12a pool that reaches
`choose_band`, and the wider cut changes which guides survive to be in it, so it has to be re-read.

Method, and why it is cheap. `choose_band` is a greedy PREFIX that breaks as soon as no remaining
candidate keeps `group_size` guides alive, so the depth it reaches does not depend on the `k` it was
asked for — ask for more than the wall and it returns the wall. `build_submission` then declines at
"band reached N of K" BEFORE `scan_cas9` runs, so one probe build per contract costs a bank scan
plus `hdr_compliance` and skips the Cas9 scan and the assembly entirely. `meta["band"]` is the wall.

Two arms per contract, and the second one is the point:

  wide    the live path — cut over `conjunction_cut_seeds` (900), band from explicit
          `band_candidates` (150-wide slice of the 300-seed band space at this hotkey's offset)
  narrow  the control — cut over the 300-seed band space, band carved by `sub_window` as before.
          This reproduces the configuration CLAUDE.md measured at 9/12. If it does not come back
          9/12 the probe is wrong, not the wall, and nothing in the wide column can be trusted.

    CW_N=4 CW_CELLS=HEK293,HUDEP-2 python conj_wall.py
    CW_ARMS=wide python conj_wall.py            # skip the control once it has been established

No network: reads sd_task_listing.json / test_hek/cell_types.json, same as conj_widecut_check.py.
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "wall")

import dataclasses
import json
import logging
import time

logging.basicConfig(level=logging.ERROR)

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402

N_PER_CELL = int(os.getenv("CW_N", "4"))
CELLS = [c for c in os.getenv("CW_CELLS", "HEK293,HUDEP-2,K562,CD34+_HSPC").split(",") if c]
ARMS = [a for a in os.getenv("CW_ARMS", "wide,narrow").split(",") if a]
# Deeper than any wall this file has ever recorded (26 on `not_mhnhej`, 12 on `hdr`), so the greedy
# always stops on its own rather than on `k`. Costs nothing: the loop breaks at the wall either way.
PROBE_K = int(os.getenv("CW_PROBE_K", "40"))
HK = os.getenv("CW_HK", "niome_hotkey")
BUDGET = float(os.getenv("CW_BUDGET", "1800"))


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def tasks_for(items, cell, n):
    """The n most recent 3-seed contracts for this cell, newest first."""
    out = []
    for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        s = str(c.get("seed", "") or "")
        if c.get("cell_type") != cell:
            continue
        if len([x for x in s.split(",") if x.strip().isdigit()]) != 3:
            continue
        out.append(t)
        if len(out) >= n:
            break
    return out


def band_space_for(cell):
    return sorted({s for a, b in JW.fallback_for(cell, HK) for s in range(a, b + 1)})


def arm_cfg(base, cell, arm):
    """(cfg, label) for one arm, at PROBE_K so the band step reports the wall and then declines."""
    space = band_space_for(cell)
    if arm == "wide":
        # The full 900 explicitly, NOT via `conjunction_cut_seeds` -- that now applies the
        # per-cell `WIDE_CUT_CELLS` gate this arm exists to measure against, so routing
        # through it would silently turn the wide arm into a second narrow one on any
        # excluded cell (HEK293 since 2026-09-17).
        cut = list(range(100, 1000))
        cands = CJ.sub_window(space, base.band_width, JW.band_offset_frac(HK) or 0.0)
        cfg = dataclasses.replace(base, seed_list=tuple(cut), band_candidates=tuple(cands),
                                  band_k=PROBE_K)
        return cfg, f"cut {len(cut)} / band cands {len(cands)}"
    # narrow: the pre-2026-09-17 shape — cut and band share the 300-seed space, `sub_window` carves.
    cfg = dataclasses.replace(base, seed_list=tuple(space), band_candidates=(),
                              band_k=PROBE_K)
    return cfg, f"cut {len(space)} / band cands {base.band_width} (sub_window)"


def main():
    cell_types = load_json("test_hek/cell_types.json")
    G.load_sequence()
    items = load_json("sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])

    rows = []
    for cell in CELLS:
        base = CJ.config_for(cell)
        if base is None:
            print(f"=== {cell}: no conjunction config, skipped ===", flush=True)
            continue
        cands = tasks_for(items, cell, N_PER_CELL)
        print(f"\n=== {cell}  live band_k={base.band_k} group={base.group_size} "
              f"width={base.band_width}  {len(cands)} contracts ===", flush=True)
        for t in cands:
            tid = (t.get("task_id") or t["id"])[:8]
            contract = dict(t["content"]["contract"])
            reference = t["content"]["hbb_reference"]
            for arm in ARMS:
                cfg, shape = arm_cfg(base, cell, arm)
                t0 = time.monotonic()
                _rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                                  budget_s=BUDGET)
                dt = time.monotonic() - t0
                wall = meta.get("band")
                bank = meta.get("bank")
                p_row = meta.get("probe_p_row")
                # Self-consistency: the wall should sit where bank * P(comply)**k falls under the
                # group, so back out the per-seed decay the observed wall implies and print it. A
                # value far from the ~0.57 erythroid / ~0.37 HEK293 per-row rate means the wall is
                # being set by something other than pool decay.
                implied = None
                if wall and bank and wall > 0 and bank > 0:
                    implied = (base.group_size / bank) ** (1.0 / wall)
                rows.append({"cell": cell, "task": tid, "arm": arm, "wall": wall,
                             "bank": bank, "p_row": p_row, "implied": implied,
                             "reason": meta.get("reason"), "s": round(dt)})
                w = "n/a" if wall is None else str(wall)
                print(f"  {tid} {arm:6s} {shape:34s} wall {w:>3s}  bank {bank}  "
                      f"p_row {p_row if p_row is None else round(p_row, 3)}  "
                      f"implied {implied if implied is None else round(implied, 3)}  ({dt:.0f}s)",
                      flush=True)
                if wall is None:
                    print(f"       no band step reached — {meta.get('reason')}", flush=True)

    with open("conj_wall.json", "w") as fh:
        json.dump(rows, fh, indent=1)

    print("\n=== summary: wall per cell per arm (CLAUDE.md records 9 HEK293 / 12 erythroid, "
          "measured at the NARROW cut) ===")
    print(f"{'cell':12s} {'arm':7s} {'walls seen':22s} {'bank (mean)':>12s} {'live band_k':>12s}")
    for cell in CELLS:
        base = CJ.config_for(cell)
        if base is None:
            continue
        for arm in ARMS:
            sel = [r for r in rows if r["cell"] == cell and r["arm"] == arm and r["wall"]]
            if not sel:
                continue
            walls = sorted({r["wall"] for r in sel})
            banks = [r["bank"] for r in sel if r["bank"]]
            shown = ", ".join(f"{w}x{sum(1 for r in sel if r['wall'] == w)}" for w in walls)
            print(f"{cell:12s} {arm:7s} {shown:22s} "
                  f"{(sum(banks)/len(banks) if banks else 0):12.0f} {base.band_k:12d}")


if __name__ == "__main__":
    main()
