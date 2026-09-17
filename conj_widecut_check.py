#!/usr/bin/env python3
"""conj_widecut_check.py — one-off sanity/timing check for the 2026-09-17 wide-cut fleet layout.

Builds the conjunction exactly the way `neurons/miner.py`'s new wide-cut path now does for a
`joined_window.BAND_HK` hotkey: cut `seed_list` = the full 900 (`joined_window.conjunction_cut_seeds`),
band drawn via explicit `band_candidates` from a 150-wide slice of the (still 300-seed) band space at
this hotkey's usual offset. Reports whether it builds and how long the cold Cas12a bank scan takes,
since `Miner.CONJUNCTION_MIN_BUDGET_S` (380/430/450s) was calibrated at the OLD 300-seed cut window
and this is a materially bigger scan.

No network: reads sd_task_listing.json / test_hek/cell_types.json directly, same as conj_nomh.py.

    python conj_widecut_check.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "wctest")

import json
import logging
import time

logging.basicConfig(level=logging.ERROR)

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def cell_types_table():
    return load_json("test_hek/cell_types.json")


def pick(items, cell):
    for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        s = str(c.get("seed", "") or "")
        if c.get("cell_type") == cell and len([x for x in s.split(",") if x.strip().isdigit()]) == 3:
            return t
    return None


def main():
    cell_types = cell_types_table()
    G.load_sequence()
    items = load_json("sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])

    band_space = sorted({s for a, b in JW.fallback_for("HEK293", "niome_hotkey")
                         for s in range(a, b + 1)})

    for cell, hk in (("CD34+_HSPC", "niome_hotkey2"), ("HUDEP-2", "niome_hotkey3")):
        t = pick(items, cell)
        if not t:
            print(f"{cell}: no cached 3-seed task found, skipping")
            continue
        tid = t.get("task_id") or t["id"]
        contract = dict(t["content"]["contract"])
        reference = t["content"]["hbb_reference"]
        cj_cfg = CJ.config_for(cell)
        conj_offset = JW.band_offset_frac(hk)
        cell_band_space = sorted({s for a, b in JW.fallback_for(cell, "niome_hotkey")
                                  for s in range(a, b + 1)})
        band_candidates = CJ.sub_window(cell_band_space, cj_cfg.band_width, conj_offset)
        wide_cut = JW.conjunction_cut_seeds(hk, cell)
        print(f"=== {cell} {tid[:8]} hotkey={hk} band_k={cj_cfg.band_k} "
              f"group={cj_cfg.group_size} cut_seeds={len(wide_cut)} "
              f"band_candidates={len(band_candidates)} ===", flush=True)
        t0 = time.monotonic()
        rows, meta = CJ.build_for_cell(contract, reference, cell_types, budget_s=1200.0,
                                       seed_list=wide_cut, band_candidates=band_candidates)
        dt = time.monotonic() - t0
        if not rows:
            print(f"  DECLINED — {meta.get('reason')}  ({dt:.0f}s)", flush=True)
            continue
        print(f"  built: band {meta['band']} clean {meta['clean']}/{len(wide_cut)} "
              f"pool {meta.get('pool')} cas9_pool {meta.get('cas9_pool')} "
              f"rows {meta['rows']} cells {meta['cells']}/8  "
              f"elapsed(meta) {meta['elapsed_s']}s  wall {dt:.0f}s", flush=True)


if __name__ == "__main__":
    main()
