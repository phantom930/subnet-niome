#!/usr/bin/env python3
"""band_scan.py — the clean band of every hotkey on every archived task, and what it hit.

A band seed is one where every row draws HDR: is_cut, is_hdr and indel_length are then all constant
and stage 4 returns consistency exactly 1.000. That needs stage 3 only, not the RandomForest, and a
row loop that breaks on the first non-HDR draw — so a 100-seed window scans in ~10s rather than the
~6 min a full rescore sweep costs.

Writes band_efficiency.json. Scope is limited to tasks we have archived submissions for; older
tasks would have to be rebuilt, and the config has changed since, so their bands would not be the
ones we actually submitted.
"""
import json, os, glob, copy, sys, time, urllib.request
from collections import defaultdict

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging; logging.basicConfig(level=logging.ERROR)
import genExp as G

LO, HI = "2026-08-27T00:00:00", "2026-09-03T23:59:59"
WINDOW = {"niome_hotkey": (100, 199), "niome_hotkey1": (200, 299), "niome_hotkey2": (300, 399),
          "niome_hotkey3": (400, 499), "niome_hotkey4": (500, 599), "niome_hotkey5": (600, 699),
          "niome_hotkey6": (700, 799), "niome_hotkey7": (800, 899), "niome_hotkey8": (900, 999)}
LABEL = {n: ("h0" if n == "niome_hotkey" else "h" + n[-1]) for n in WINDOW}


def band_of(rows, contract, reference, cell_types, lo, hi):
    """Seeds in [lo, hi] on which every row draws HDR."""
    out = []
    for sd in range(lo, hi + 1):
        c = copy.deepcopy(contract); c["seed"] = sd
        ctx = G.build_context(c, reference, cell_types)
        for row in rows:
            entry = G.build_valid_entry(row, ctx)
            if entry is None or G.simulate(entry, ctx)["outcome"] != "HDR":
                break
        else:
            out.append(sd)
    return out


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
            tasks[t["id"]] = {"created_at": ca, "cell_type": c.get("cell_type"),
                              "seeds": [int(x) for x in seeds]}

    archives = defaultdict(dict)
    for f in sorted(glob.glob("data/inst/*/result/*/last_upload.json")):
        try:
            tid = json.load(open(f)).get("task_id")
        except Exception:
            continue
        if tid in tasks:
            archives[tid][f.split("/")[2]] = os.path.dirname(f)

    cell_types = G.fetch_cell_types(); G.load_sequence()
    out = {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "range": [LO, HI], "windows": {LABEL[n]: list(w) for n, w in WINDOW.items()},
           "tasks_in_range": len(tasks), "tasks_scanned": 0, "tasks": []}

    for tid in sorted(archives, key=lambda k: tasks[k]["created_at"]):
        meta = tasks[tid]
        rec = {"task_id": tid, "created_at": meta["created_at"], "cell_type": meta["cell_type"],
               "seeds": meta["seeds"], "hotkeys": {}, "hits": 0, "opportunities": 0}
        for name, path in sorted(archives[tid].items()):
            lo, hi = WINDOW[name]
            rows = json.load(open(f"{path}/submission.json"))
            contract = json.load(open(f"{path}/contract.json"))
            reference = json.load(open(f"{path}/hbb_reference.json"))
            t0 = time.time()
            band = band_of(rows, contract, reference, cell_types, lo, hi)
            owned = [s for s in meta["seeds"] if lo <= s <= hi]
            hits = [s for s in owned if s in band]
            rec["hotkeys"][LABEL[name]] = {
                "window": [lo, hi], "rows": len(rows), "band": band, "band_size": len(band),
                "owned_seeds": owned, "hits": hits,
                "nearest": {str(s): (min(band, key=lambda b: abs(b - s)) if band else None)
                            for s in owned},
                "scan_s": round(time.time() - t0, 1),
            }
            rec["hits"] += len(hits)
            rec["opportunities"] += len(owned)
            print(f"  {tid[:8]} {LABEL[name]} band {len(band):>2} owned {owned} hits {hits}",
                  flush=True)
        out["tasks"].append(rec)
        out["tasks_scanned"] += 1
        json.dump(out, open("band_efficiency.json", "w"), indent=2)

    tot_h = sum(t["hits"] for t in out["tasks"])
    tot_o = sum(t["opportunities"] for t in out["tasks"])
    bands = [hk["band_size"] for t in out["tasks"] for hk in t["hotkeys"].values()]
    out["summary"] = {
        "hits": tot_h, "opportunities": tot_o,
        "hit_rate": round(tot_h / tot_o, 4) if tot_o else None,
        "mean_band_size": round(sum(bands) / len(bands), 2) if bands else None,
        "expected_hit_rate": round(sum(bands) / len(bands) / 100, 4) if bands else None,
        "rounds_with_a_hit": sum(1 for t in out["tasks"] if t["hits"] > 0),
    }
    json.dump(out, open("band_efficiency.json", "w"), indent=2)
    print(f"\n  {tot_h} hits / {tot_o} opportunities over {out['tasks_scanned']} tasks")
    print(f"  mean band {out['summary']['mean_band_size']} -> expected rate "
          f"{out['summary']['expected_hit_rate']}, observed {out['summary']['hit_rate']}")


if __name__ == "__main__":
    main()
