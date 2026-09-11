#!/usr/bin/env python3
"""allhdr_mf.py — all-HDR's clean band against `main_max_fail`, at width 225 on K562.

`main_max_fail` is the screen that decides which Cas12a guides enter the bank: a guide qualifies if
it reaches HDR on all but `mf` of the window's seeds. It is calibrated at 45 of 100 and scaled to
wider windows by `_scaled_max_fail` -- linearly for K562, which puts width 225 at **101**, the value
the live fleet ships and the only one ever measured there (band 12, `joined300.json`).

The screen is not what produces the band; the MIN-UNION is. So the sweep has two opposing forces
and an optimum somewhere between:

  * looser mf -> more candidates -> more freedom to find guides whose failures COINCIDE, which is
    what a small union needs. CLAUDE.md's note at `main_max_fail` argues this is why 45/100 is
    deliberately loose.
  * looser mf -> the extra candidates each carry more failures, and `load_bank` truncates to the
    best 60,000 by fail count -- so past some point the additions are discarded and the band
    saturates. Measured on HEK293 at span 300: band FLAT at 8 for mf 164/184/192/204 while build
    time ran 57/123/190/255s.

And a wider band is not automatically better, which is why every arm here runs the FULL build
rather than just the min-union: the Cas9 half must reach HDR on every clean seed, a requirement
that decays as P(HDR)**band, so a band bought with a loose screen can leave the Cas9 pool unable to
fill and the build declines. An arm that declines is worth zero, not its band.

Reported per arm: raw bank (pre-truncation), loaded bank, group union, clean band in-window, the
true band over 100-999, Cas9 pool, and whether it assembled 250 rows over 8/8 cells.

    python allhdr_mf.py
    AM_MF=45,101,170 AM_WINDOW=200-424 python allhdr_mf.py
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "ahmf")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import random                        # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank   # noqa: E402
from joined300 import fetch_tasks, true_band                # noqa: E402
from sd_task import score                                   # noqa: E402

CELL = os.getenv("AM_CELL", "K562")
MFS = [int(x) for x in os.getenv("AM_MF", "45,70,85,101,120,140,170,200").split(",")]
WINDOW = tuple(int(x) for x in os.getenv("AM_WINDOW", "100-324").split("-"))
TASK = os.getenv("AM_TASK", "")           # blank -> newest stamped task for this cell
BUDGET = float(os.getenv("AM_BUDGET", "2400"))
OUT = os.getenv("AM_OUT", "allhdr_mf.json")


def pick(items):
    for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL:
            continue
        if TASK:
            if (t.get("task_id") or t.get("id", "")).startswith(TASK):
                return t
            continue
        if str(c.get("seed", "0")).strip() not in ("0", ""):
            return t
    raise SystemExit(f"no task found for {CELL} {TASK}")


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    rng = random.Random(20260911)
    lo, hi = WINDOW
    span = hi - lo + 1
    base = AH.config_for(CELL)
    linear = AH._scaled_max_fail(CELL, base.main_max_fail, span)
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  "
          f"{task.get('created_at','')[:16]}  window {lo}-{hi} (span {span})  "
          f"group {base.group_size}  shipped mf {linear}\n")
    print(f"{'mf':>5}{'raw bank':>10}{'loaded':>8}{'union':>7}{'band':>6}{'of 900':>8}"
          f"{'cas9 pool':>11}{'rows':>6}{'cells':>7}{'weighted':>10}{'fid':>8}{'s':>6}")
    out = []
    for mf in MFS:
        t0 = time.monotonic()
        cfg = _dc.replace(base, hdr_range=(lo, hi), main_max_fail=mf,
                          variants=min(base.variants, AH.WIDE_WINDOW_VARIANTS))
        rec = {"cell": CELL, "task": (task.get("task_id") or task["id"])[:8],
               "window": [lo, hi], "span": span, "mf": mf, "shipped_mf": linear,
               "group_size": base.group_size}
        # Build the bank explicitly so the RAW candidate count is visible: load_bank truncates to
        # the best 60,000 by fail count, which is the mechanism behind any saturation here.
        path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
        if not os.path.exists(path):
            ctx = G.build_context(contract, reference, cell_types)
            sites = G.enumerate_sites(ctx, 3000, (20, 23))
            bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg)
            MT.free_gpu_memory()
            rec["raw_bank"] = len(bank)
            if not bank:
                rec["declined"] = f"no Cas12a guide qualifies at mf {mf}"
                print(f"{mf:>5}{0:>10}   —  no qualifying guide")
                out.append(rec); json.dump(out, open(OUT, "w"), indent=1); continue
            save_bank(path, bank)
        rec["loaded_bank"] = len(load_bank(path))

        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg,
                                         budget_s=BUDGET)
        rec["build_s"] = round(time.monotonic() - t0, 1)
        rec["meta"] = {k: v for k, v in meta.items() if k != "bank_path"}
        rec.update(union=meta.get("union"), clean_in_window=meta.get("clean"),
                   cas9_pool=meta.get("cas9_pool"), cells=meta.get("cells"))
        if not rows:
            rec["declined"] = meta.get("reason", "declined")
            print(f"{mf:>5}{rec.get('raw_bank','-'):>10}{rec['loaded_bank']:>8}"
                  f"{str(meta.get('union','-')):>7}{str(meta.get('clean','-')):>6}"
                  f"{'':>8}{str(meta.get('cas9_pool','-')):>11}   DECLINED — {rec['declined']}")
            MT.free_gpu_memory()
            out.append(rec); json.dump(out, open(OUT, "w"), indent=1); continue

        band, _n = true_band(rows, contract, reference, cell_types)
        hit = score(rows, contract, reference, cell_types,
                    seed=band[len(band) // 2] if band else lo)
        rec.update(rows=len(rows), band=len(band),
                   band_in_window=len([s for s in band if lo <= s <= hi]),
                   leak=len(band) - len([s for s in band if lo <= s <= hi]),
                   band_seeds=band, weighted=hit["weighted"], fidelity=hit["fidelity"],
                   wxfid=hit["weighted"] * hit["fidelity"], spike_cons=hit["consistency"])
        MT.free_gpu_memory()
        print(f"{mf:>5}{rec.get('raw_bank','cached'):>10}{rec['loaded_bank']:>8}"
              f"{rec['union']:>7}{rec['clean_in_window']:>6}{rec['band']:>8}"
              f"{rec['cas9_pool']:>11}{rec['rows']:>6}{rec['cells']:>5}/8"
              f"{rec['weighted']:>10.1f}{rec['fidelity']:>8.4f}{rec['build_s']:>6.0f}")
        out.append(rec)
        json.dump(out, open(OUT, "w"), indent=1)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
