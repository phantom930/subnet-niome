#!/usr/bin/env python3
"""band_rebuild.py — reconstruct each hotkey's band on tasks we have no archived submission for.

Fills the gap band_scan.py cannot: the fleet only began archiving 2026-09-01, so most three-seed
tasks in the window have no stored submission. This rebuilds the band under the CURRENT config
(group_size 80, light_cell_rows 6) and asks how often it would have caught a drawn seed.

Two economies, both load-bearing:

* Only the three windows that own a drawn seed are scanned. The other six hotkeys have no
  opportunity on that round and cannot change the hit count.
* The band is the complement of the Cas12a min-union's failed-seed set. ``scan_cas9`` fills rows
  onto that band but cannot alter it, so the Cas9 scan and assemble are skipped.

**Rejected approach, do not retry:** one 900-seed bank per contract, windows derived by filtering
each guide's fail list. It is 3x faster and wrong — ``main_max_fail`` scaled by rate is far
stricter over a wide window (45/100 is -1 SD, 405/900 is -3 SD), so the bank collapses from 300k
guides to ~1k and every band comes out at 4 instead of 12-13.

Caveat recorded in the output: skipping scan_cas9 means these are the bands the construction WOULD
produce, not proof the build would have completed. A contract whose Cas9 pool came up short would
have declined and submitted nothing. That has not happened on any contract built at group 80.
"""
import json, os, sys, time, dataclasses, urllib.request

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging; logging.basicConfig(level=logging.ERROR)
import genExp as G
from niome_subnet.genomics import all_hdr as AH, mt19937 as MT, fastgreedy as FG

LO, HI = "2026-08-27T00:00:00", "2026-09-03T23:59:59"
WINDOWS = [(lo, lo + 99) for lo in range(100, 1000, 100)]
LABEL = {lo: f"h{i}" for i, (lo, _hi) in enumerate(WINDOWS)}
OUT = "band_rebuilt.json"


def band_for(contract, reference, cell_types, window):
    """Clean band for one window: bank -> min-union -> complement of the failed-seed union."""
    cell = contract.get("cell_type")
    base = AH.config_for(cell)
    if base is None:
        return None, "no all-HDR config for this cell type"
    cfg = dataclasses.replace(base, hdr_range=window)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
    MT.free_gpu_memory()
    if not bank:
        return None, "bank scan produced nothing"
    if len(bank) < cfg.group_size:
        return None, f"bank {len(bank)} short of group {cfg.group_size}"
    sel = FG.FastGreedy(bank, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg))
    idx, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
    bad = set()
    for i in idx:
        bad.update(int(x) for x in bank[i]["fails"])
    clean = sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad)
    return clean, None


def main():
    tl = json.load(urllib.request.urlopen("https://niome-api.genomes.io/api/v3/tasks", timeout=60))
    tl = tl if isinstance(tl, list) else (tl.get("data") or tl.get("items") or [])
    tasks = {}
    for t in tl:
        ca = t.get("created_at") or ""
        if not (LO <= ca <= HI):
            continue
        c = (t.get("content") or {}).get("contract") or {}
        seeds = [x.strip() for x in str(c.get("seed") or "").split(",") if x.strip()]
        if len(seeds) == 3:
            tasks[t["id"]] = (t, [int(x) for x in seeds])

    done = set()
    if os.path.exists(OUT):
        prev = json.load(open(OUT))
        done = {r["task_id"] for r in prev["tasks"]}
        out = prev
    else:
        out = {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "note": ("bands rebuilt under the current config (group 80, light_cell_rows 6); "
                        "only the 3 windows owning a drawn seed are scanned; scan_cas9 skipped "
                        "so a build that would have DECLINED is not detected"),
               "tasks": []}

    scanned = json.load(open("band_efficiency.json"))["tasks"] if os.path.exists("band_efficiency.json") else []
    archived = {t["task_id"] for t in scanned}
    todo = [k for k in tasks if k not in archived and k not in done]
    todo.sort(key=lambda k: tasks[k][0].get("created_at") or "")
    print(f"{len(tasks)} three-seed tasks | {len(archived)} already measured from archives | "
          f"{len(todo)} to rebuild", flush=True)

    cell_types = G.fetch_cell_types(); G.load_sequence()
    t_start = time.time()
    for n, tid in enumerate(todo, 1):
        t, seeds = tasks[tid]
        contract = dict(t["content"]["contract"]); contract["seed"] = 500
        reference = t["content"]["hbb_reference"]
        rec = {"task_id": tid, "created_at": t["created_at"],
               "cell_type": contract.get("cell_type"), "seeds": seeds,
               "source": "rebuilt", "hotkeys": {}, "hits": 0, "opportunities": 0}
        # Group the seeds by window FIRST: two of a round's three seeds can land in the same
        # hotkey's window, and writing rec["hotkeys"][tag] per seed would overwrite the earlier
        # record (losing a hit from the detail while still counting it) and rescan the same
        # window twice for an identical band.
        by_window = {}
        for seed in seeds:
            by_window.setdefault((seed // 100) * 100, []).append(seed)
        for lo, owned in sorted(by_window.items()):
            window = (lo, lo + 99)
            tag = LABEL[lo]
            t0 = time.time()
            band, err = band_for(contract, reference, cell_types, window)
            if band is None:
                rec["hotkeys"][tag] = {"window": list(window), "error": err,
                                       "owned_seeds": owned, "hits": []}
                print(f"  [{n}/{len(todo)}] {tid[:8]} {tag} seeds {owned}: DECLINED ({err})",
                      flush=True)
                rec["opportunities"] += len(owned)
                continue
            hits = [s for s in owned if s in band]
            rec["hotkeys"][tag] = {
                "window": list(window), "band": band, "band_size": len(band),
                "owned_seeds": owned, "hits": hits,
                "nearest": {str(s): min(band, key=lambda b: abs(b - s)) for s in owned}
                           if band else {},
                "scan_s": round(time.time() - t0, 1)}
            rec["hits"] += len(hits)
            rec["opportunities"] += len(owned)
            print(f"  [{n}/{len(todo)}] {tid[:8]} {rec['cell_type']:11} {tag} seeds {owned} "
                  f"band {len(band):>2} -> {len(hits)} hit  ({time.time()-t0:.0f}s)", flush=True)
        out["tasks"].append(rec)
        json.dump(out, open(OUT, "w"), indent=2)
        el = time.time() - t_start
        print(f"    ({n}/{len(todo)} tasks, {el/60:.0f}m elapsed, "
              f"~{el/n*(len(todo)-n)/60:.0f}m left)", flush=True)
    h = sum(r["hits"] for r in out["tasks"]); o = sum(r["opportunities"] for r in out["tasks"])
    print(f"\nrebuilt: {h}/{o} hits over {len(out['tasks'])} tasks")


if __name__ == "__main__":
    main()
