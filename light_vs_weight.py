#!/usr/bin/env python3
"""light_vs_weight.py — does the best `light_cell_rows` track the contract's mutation weight spread?

`light_cell_rows` replaces `assemble`'s smooth mutation apportionment with a hard quota: every
light (mutation, Cas9, strand) cell gets exactly N rows and the heaviest mutation takes the rest.
CLAUDE.md records the effect as SPREAD-DEPENDENT — paired over 105 configs, light=6 measured
-0.0008 E[pay] at weight spread <= 2.0 but +0.0026 at 2.0-2.6 — which implies the shipped constant
6 is wrong for some contracts. K562's 28 current-regime tasks span spread 1.50-2.54, straddling
that boundary, so the relation can be measured directly rather than inferred.

One all-HDR mixed build per task (joined-150 band window, group 80), then `assemble` re-run at each
`light_cell_rows` value. Everything upstream of the apportionment — bank, min-union, band, Cas9 scan
— is shared across the light arms, because `light_cell_rows` only governs how the Cas9 rows already
found are distributed. That is what makes 28 tasks affordable.

`total_weighted_score` and `distribution_fidelity_score` are design-only and so seed-independent
(verified in fidelity_window.py), so one scoring pass per arm at an arbitrary seed gives the whole
comparison; the band is identical across light values by construction.

    python light_vs_weight.py
    LW_TASKS=7306c626,f2d36f4d LW_LIGHT=6,12 python light_vs_weight.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "lvw")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402
import urllib.request                # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank   # noqa: E402
from sd_task import score                                   # noqa: E402

CLASSES = [tuple(int(y) for y in x.split("-"))
           for x in os.getenv("LW_JOINED", "100-149,400-449,700-749").split(",")]
GROUP = int(os.getenv("LW_GROUP", "80"))
LIGHTS = [None if x.lower() in ("none", "") else int(x)
          for x in os.getenv("LW_LIGHT", "none,6,12,25").split(",")]
ONLY = [t for t in os.getenv("LW_TASKS", "").split(",") if t]
SEED = int(os.getenv("LW_SEED", "500"))
OUT = os.getenv("LW_OUT", "light_vs_weight.json")
MAXLOAD = float(os.getenv("LW_MAXLOAD", "4"))


def k562_tasks():
    tk = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/tasks?limit=500"))
    tk = tk if isinstance(tk, list) else (tk.get("items") or tk.get("data"))
    out = []
    for t in tk:
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != "K562":
            continue
        s = str(c.get("seed", "") or "")
        if len([x for x in s.split(",") if x.strip().isdigit()]) != 3:
            continue
        w = c.get("mutation_weights") or {}
        if len(w) != 2:
            continue
        tid = (t.get("task_id") or t["id"])
        if ONLY and not any(tid.startswith(p) for p in ONLY):
            continue
        out.append((tid, t, max(w.values()) / min(w.values())))
    out.sort(key=lambda r: r[2])
    return out


def wait_for_idle():
    while True:
        with open("/proc/loadavg") as fh:
            if float(fh.read().split()[0]) < MAXLOAD:
                return
        time.sleep(20)


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    tasks = k562_tasks()
    print(f"K562, all-HDR mixed, joined-{len(joined)} band window, group {GROUP}")
    print(f"{len(tasks)} tasks, light_cell_rows {LIGHTS}\n")
    print(f"{'task':<10}{'spread':>7}{'band':>6}{'cas9':>7}"
          + "".join(f"{('L'+str(l)):>18}" for l in LIGHTS) + f"{'best':>7}{'s':>6}")
    out = []
    for tid, task, spread in tasks:
        wait_for_idle()
        t0 = time.monotonic()
        content = task["content"]
        contract, reference = content["contract"], content["hbb_reference"]
        try:
            ctx = G.build_context(contract, reference, cell_types)
            sites = G.enumerate_sites(ctx, 3000, (20, 23))
            n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
            cfg = AH.config_for("K562")
            cfg = _dc.replace(cfg, hdr_range=(min(joined), max(joined)),
                              seed_list=tuple(joined), group_size=GROUP,
                              main_max_fail=AH._scaled_max_fail("K562", 45, len(joined)),
                              variants=min(cfg.variants, AH.WIDE_WINDOW_VARIANTS))
            path = os.path.join(AH.HDR_BANK_DIR,
                                f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
            if not os.path.exists(path):
                bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg)
                MT.free_gpu_memory()
                if not bank:
                    print(f"{tid[:8]:<10}{spread:>7.2f}   bank empty"); continue
                save_bank(path, bank)
            records = load_bank(path)
            if len(records) < GROUP:
                print(f"{tid[:8]:<10}{spread:>7.2f}   bank {len(records)} < group"); continue
            sel = FG.FastGreedy(records, window_lo=min(joined), window_hi=max(joined),
                                per_cell_min=cfg.per_cell_min,
                                seeds=np.asarray(joined, dtype=np.int64))
            idx, _u = sel.best(GROUP, restarts=cfg.restarts)
            group = [records[i] for i in idx]
            bad = set()
            for rec in group:
                bad.update(int(x) for x in rec["fails"])
            clean = np.array(sorted(set(joined) - bad), dtype=np.int64)
            if clean.size == 0:
                print(f"{tid[:8]:<10}{spread:>7.2f}   band empty"); continue
            cas9 = AH.scan_cas9(clean, contract, cell_types, ctx, sites, cfg, n_rows - GROUP)
            if len(cas9) < n_rows - GROUP:
                print(f"{tid[:8]:<10}{spread:>7.2f}{clean.size:>6}{len(cas9):>7}"
                      f"   Cas9 short of {n_rows-GROUP}")
                out.append({"task": tid[:8], "spread": spread, "band": int(clean.size),
                            "cas9": len(cas9), "reason": "cas9 short"})
                json.dump(out, open(OUT, "w"), indent=1); continue
            rec = {"task": tid[:8], "spread": spread, "heavy": max(
                       contract["mutation_weights"].values()),
                   "light_w": min(contract["mutation_weights"].values()),
                   "band": int(clean.size), "cas9": len(cas9), "arms": {}}
            line = f"{tid[:8]:<10}{spread:>7.2f}{clean.size:>6}{len(cas9):>7}"
            for light in LIGHTS:
                rows = AC.assemble(group, cas9, contract, ctx,
                                   _dc.replace(cfg, light_cell_rows=light), n_rows)
                s = score(rows, contract, reference, cell_types, seed=SEED)
                wf = s["weighted"] * s["fidelity"]
                heavy_m = max(contract["mutation_weights"],
                              key=lambda m: contract["mutation_weights"][m])
                rec["arms"][str(light)] = {
                    "weighted": s["weighted"], "fidelity": s["fidelity"], "wxfid": wf,
                    "rows": len(rows),
                    "heavy_rows": sum(1 for r in rows if r["mutation"] == heavy_m)}
                line += f"{wf:>10.1f}/{s['fidelity']:.3f}"
            best = max(rec["arms"], key=lambda k: rec["arms"][k]["wxfid"])
            rec["best_light"] = best
            rec["build_s"] = round(time.monotonic() - t0, 1)
            print(line + f"{best:>7}{rec['build_s']:>6.0f}")
            out.append(rec)
            json.dump(out, open(OUT, "w"), indent=1)
            MT.free_gpu_memory()
        except Exception as exc:
            print(f"{tid[:8]:<10}{spread:>7.2f}   ERROR {type(exc).__name__}: {exc}")
            out.append({"task": tid[:8], "spread": spread,
                        "error": f"{type(exc).__name__}: {exc}"})
            json.dump(out, open(OUT, "w"), indent=1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
