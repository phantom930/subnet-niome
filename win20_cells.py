#!/usr/bin/env python3
"""win20_cells.py — does the width-20 band gain hold across cell types and contracts?

Everything measured so far is one K562 contract: band 15 of 20 against 13 of 100 at the shipped
load cap, flat across widths 20-50, with the wide-width decline traced to bank starvation rather
than width. This asks whether that transfers, over several contracts for each of the four cell
types the backend issues.

The comparison is shipped-against-shipped. Width 100 runs exactly what a hotkey runs today — the
cell type's own ``CELL_CONFIG`` window and ``main_max_fail``, ``load_bank`` at its 60,000 default.
Width 20 is that same config with the window narrowed to 20 seeds at the shipped window's ``lo``
and ``lo + 40``, two windows per contract so window-to-window noise is averaged rather than assumed
away. Both are also solved at limit 150,000, because the cap binds asymmetrically: it cost a band
seed at width 100 and none at width 20, and on one width-30 window a *larger* pool made the band
*smaller* (FastGreedy is a heuristic min-union, so it is not monotone in pool size).

``main_max_fail`` cannot be shared across cell types. HEK293 sits at P(HDR) ~ 0.37 against the
erythroid types' ~0.57, so it fails ~12.6 seeds of 20 where they fail ~8.6, and the erythroid
``mf6`` would starve its bank. The nominal values below are z-matched to ``mf6``; any bank landing
under the 60,000 load cap is rebuilt at ``mf + 2``, because that starvation is exactly what
corrupted the first width sweep (36k records at width 50 against 150k at width 20, read as a width
effect until the screen was loosened).

Resumable: re-running skips any (contract, width, window) already in the output.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import collections
import dataclasses
import hashlib
import json
import sys
import time
import urllib.request
from collections import Counter

import numpy as np

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import fastgreedy as FG
from niome_subnet.genomics import mt19937 as MT
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

TASKS_URL = "https://niome-api.genomes.io/api/v3/tasks"
CONTRACTS_PER_CELL = 5
LIMITS = [60_000, 150_000]
NARROW_OFFSETS = [0, 40]          # 20-seed windows at lo and lo+40 of the shipped window
NOMINAL_MF20 = {"HEK293": 10}     # z-matched to the erythroid mf6; default below
DEFAULT_MF20 = 6
MIN_BANK = 60_000                 # under the load cap => the screen, not the width, is the limit
OUT = "win20_cells.json"


def pick_contracts():
    """Distinct contracts per cell type, spread across history rather than clustered."""
    tl = json.load(urllib.request.urlopen(TASKS_URL, timeout=60))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    by = collections.defaultdict(dict)
    for t in sorted(tl, key=lambda x: x.get("created_at") or ""):
        content = t.get("content") or {}
        c = content.get("contract") or {}
        cell = c.get("cell_type")
        if not cell or not content.get("hbb_reference"):
            continue
        key = hashlib.sha256(json.dumps(
            {k: c.get(k) for k in ("active_mutations", "mutation_weights",
                                   "mutation_regions", "rules")},
            sort_keys=True).encode()).hexdigest()[:12]
        by[cell].setdefault(key, t)
    out = []
    for cell, d in sorted(by.items()):
        tasks = list(d.values())
        step = max(1, len(tasks) // CONTRACTS_PER_CELL)
        out.extend(tasks[::step][:CONTRACTS_PER_CELL])
    return out


def bank_path_for(contract, cell_types, cfg):
    return os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")


def ensure_bank(contract, reference, cell_types, ctx, sites, cfg, adapt):
    """Bank for this config, loosening the screen while it lands under the load cap."""
    for bump in (0, 2, 4) if adapt else (0,):
        use = dataclasses.replace(cfg, main_max_fail=cfg.main_max_fail + bump)
        path = bank_path_for(contract, cell_types, use)
        if not os.path.exists(path):
            bank = AH.build_bank(contract, reference, cell_types, ctx, sites, use, None)
            MT.free_gpu_memory()
            if not bank:
                continue
            save_bank(path, bank)
        n = int(np.load(path, allow_pickle=False)["fails"].shape[0])
        if n >= MIN_BANK or bump == (4 if adapt else 0):
            return path, n, use
    return None, 0, cfg


def one_build(contract, reference, cell_types, ctx, sites, cfg, adapt, probe, seeds):
    started = time.monotonic()
    path, disk, cfg = ensure_bank(contract, reference, cell_types, ctx, sites, cfg, adapt)
    rec = {"max_fail": cfg.main_max_fail, "disk": disk,
           "window": [cfg.start_seed, cfg.end_seed],
           "width": cfg.end_seed - cfg.start_seed + 1}
    if path is None:
        return {**rec, "declined": "bank scan produced nothing"}
    bands = {}
    for limit in LIMITS:
        records = load_bank(path, limit=limit)
        if len(records) < cfg.group_size:
            continue
        sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                            per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg))
        index, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
        bad = set()
        for i in index:
            bad.update(int(x) for x in records[i]["fails"])
        bands[limit] = (records, index,
                        sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad))
    if not bands:
        return {**rec, "declined": f"bank {disk} short of group {cfg.group_size}"}
    rec.update({f"band_{k // 1000}k": len(v[2]) for k, v in bands.items()})
    best_limit = max(bands, key=lambda L: (len(bands[L][2]), -L))
    records, index, band = bands[best_limit]
    rec.update(best_limit=best_limit, band=band, clean=len(band),
               clean_fraction=round(len(band) / rec["width"], 4),
               k=len([s for s in seeds if s in set(band)]))
    if not band:
        return {**rec, "declined": "failures cover the whole window"}

    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    want = n_rows - cfg.group_size
    cas9 = AH.scan_cas9(np.array(band, dtype=np.int64), contract, cell_types, ctx, sites, cfg,
                        want, None)
    MT.free_gpu_memory()
    rec["cas9_pool"] = len(cas9)
    cells4 = len({(r["mutation"], r["strand"]) for r in cas9})
    if len(cas9) < want or cells4 < 4:
        return {**rec, "declined": f"Cas9 pool {len(cas9)} over {cells4} cells short of {want}",
                "build_s": round(time.monotonic() - started, 1)}
    rows = AH.assemble([records[i] for i in index], cas9, contract, ctx, cfg, n_rows)
    build_s = time.monotonic() - started
    cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    run_stage3(seed=probe)
    run_stage4(seed=probe)
    res = run_stage5()
    rec.update(rows=len(rows), cells=cells, build_s=round(build_s, 1),
               total_weighted_score=res["total_weighted_score"],
               distribution_fidelity_factor=res["distribution_fidelity_factor"])
    if len(rows) < n_rows or cells < 8:
        rec["declined"] = f"assembled {len(rows)} rows over {cells} cells"
    return rec


def main():
    tasks = pick_contracts()
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    done, out = set(), []
    if os.path.exists(OUT):
        out = json.load(open(OUT))
        done = {(r["task_id"], r["width"], r["window"][0]) for r in out}
    counts = collections.Counter(
        (t["content"]["contract"]["cell_type"]) for t in tasks)
    print(f"{len(tasks)} contracts: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    print(f"{len(done)} builds already done\n")
    print(f"{'cell':<11} {'task':<9} {'w':>4} {'mf':>3} {'window':>9} {'disk':>7} {'b@60k':>6} "
          f"{'b@150k':>7} {'band':>5} {'frac':>6} {'cas9':>6} {'rows':>5} {'weighted':>9} "
          f"{'fid':>6} {'build':>7}")

    started = time.monotonic()
    total = len(tasks) * (1 + len(NARROW_OFFSETS))
    step = 0
    for task in tasks:
        contract = task["content"]["contract"]
        reference = task["content"]["hbb_reference"]
        cell = contract["cell_type"]
        seeds = [int(x) for x in str(contract.get("seed") or "0").split(",") if x.strip()]
        probe = next((s for s in seeds if s >= 100), 500)
        base = AH.config_for(cell)
        if base is None:
            continue
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        lo = base.start_seed
        plan = [(dataclasses.replace(base), False)]
        mf20 = NOMINAL_MF20.get(cell, DEFAULT_MF20)
        for off in NARROW_OFFSETS:
            plan.append((dataclasses.replace(base, hdr_range=(lo + off, lo + off + 19),
                                             main_max_fail=mf20), True))
        for cfg, adapt in plan:
            step += 1
            key = (task["id"], cfg.end_seed - cfg.start_seed + 1, cfg.start_seed)
            if key in done:
                continue
            rec = one_build(contract, reference, cell_types, ctx, sites, cfg, adapt, probe, seeds)
            rec.update(task_id=task["id"], cell_type=cell)
            out.append(rec)
            json.dump(out, open(OUT, "w"), indent=1)
            eta = (time.monotonic() - started) / step * (total - step) / 60
            tail = (f"{rec.get('cas9_pool', 0):>6} {rec.get('rows', 0):>5} "
                    f"{rec.get('total_weighted_score', 0):>9.1f} "
                    f"{rec.get('distribution_fidelity_factor', 0):>6.3f} "
                    f"{rec.get('build_s', 0):>6.1f}s")
            if "declined" in rec:
                tail = f"DECLINED: {rec['declined']}"
            print(f"{cell:<11} {task['id'][:8]:<9} {rec['width']:>4} {rec['max_fail']:>3} "
                  f"{rec['window'][0]}-{rec['window'][1]:<4} {rec['disk']:>7} "
                  f"{rec.get('band_60k', 0):>6} {rec.get('band_150k', 0):>7} "
                  f"{rec.get('clean', 0):>5} {rec.get('clean_fraction', 0) * 100:>5.1f}% {tail}"
                  f"   eta {eta:.0f}m", flush=True)

    print("\nband by cell type and width (mean over contracts, declines excluded):")
    for cell in sorted({r["cell_type"] for r in out}):
        for width in (100, 20):
            got = [r for r in out if r["cell_type"] == cell and r["width"] == width
                   and "declined" not in r]
            all_n = [r for r in out if r["cell_type"] == cell and r["width"] == width]
            if not got:
                print(f"  {cell:<11} w{width:<4} no builds ({len(all_n)} attempted)")
                continue
            print(f"  {cell:<11} w{width:<4} band {sum(r['clean'] for r in got) / len(got):>5.2f} "
                  f"weighted {sum(r['total_weighted_score'] for r in got) / len(got):>6.1f} "
                  f"fid {sum(r['distribution_fidelity_factor'] for r in got) / len(got):>5.3f} "
                  f"built {len(got)}/{len(all_n)}")


if __name__ == "__main__":
    main()
