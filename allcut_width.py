#!/usr/bin/env python3
"""allcut_width.py — all-cut's CLEAN SEED COUNT against window width, per cell type.

all-cut min-unions its Cas12a group over a seed window; the group's failed-seed union is what the
submission cannot save, and the complement is the clean set. Narrowing the window confines those
failures to a smaller space, so the clean *fraction* rises and the clean *set* shrinks -- and the
set is what pays, because a round draws its three seeds from the whole 100-999.

What is already measured is partial and on old contracts: K562 at width 225 (clean 225 of 225, four
tiles, `clean225.json`) and 300 (280-281 of 300, `clean300.json`), HEK293 at 225 (102 of 225,
`acut225.json`), and the whole-window figures CLAUDE.md records (K562 549-570, HEK293 137). HUDEP-2
and CD34+_HSPC have never been measured at any narrow width, and nothing has been measured at 150.
This fills the table on one current contract per cell type, so the four cells are comparable.

Three things are reported per arm, and the difference between them is the point:

  * **clean in-window** -- what the min-union leaves inside its own window, the number the earlier
    scripts report.
  * **clean over 100-999** -- seeds where every one of the 250 ASSEMBLED rows cuts, scanned over
    the whole space. Outside its window the Cas12a group is unconstrained, so this is the in-window
    count plus whatever leaks; it is the number a round actually draws against.
  * **E[final]** -- `weighted x fidelity x E[cons]` at the measured clean/dirty per-seed
    consistencies, which is what decides whether a wider clean fraction was worth the smaller set.

`cas12a_max_fail` is a screen over the window and must scale with the span or a narrow window
demands the impossible (100 of 900 is 11%; held fixed at 225 it would demand 56% clean per guide).
Linear scaling is used, matching `clean225.py`/`acut225.py`. Note this scales the SAFE way here:
narrowing makes a fixed rate a shallower tail cut, so the bank grows rather than starves -- the
opposite of the span-900 failure `joined300_mf.py` found for all-HDR.

    python allcut_width.py                      # four cell types, widths 900/300/225/150
    AW_CELLS=K562 AW_WIDTHS=225 AW_TILES=4 python allcut_width.py
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "acw")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import random                        # noqa: E402
import statistics as st              # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.validation import stage3         # noqa: E402
from niome_subnet.genomics.validation.stage12 import run_stage12   # noqa: E402
from niome_subnet.utils import settings                     # noqa: E402
from joined300 import fetch_tasks, pick_tasks               # noqa: E402
from sd_task import score                                   # noqa: E402

CELLS = [c for c in os.getenv("AW_CELLS", "K562,HUDEP-2,CD34+_HSPC,HEK293").split(",") if c]
WIDTHS = [int(w) for w in os.getenv("AW_WIDTHS", "900,300,225,150").split(",") if w]
TILES = int(os.getenv("AW_TILES", "1"))          # how many tiles of each width to build
N_CLEAN = int(os.getenv("AW_N_CLEAN", "2"))      # clean seeds sampled for the clean-leg value
N_DIRTY = int(os.getenv("AW_N_DIRTY", "2"))
BUDGET = float(os.getenv("AW_BUDGET", "3000"))
OUT = os.getenv("AW_OUT", "allcut_width.json")
LO, HI = 100, 999


def windows_for(width):
    """The first `TILES` tiles of `width` covering 100-999 from the bottom."""
    n = max(1, (HI - LO + 1) // width)
    return [(LO + i * width, min(HI, LO + (i + 1) * width - 1)) for i in range(min(TILES, n))]


def true_clean(rows, contract, reference, cell_types):
    """Seeds in 100-999 where EVERY assembled row cuts, and the valid-row count.

    Bails on the first row that does not cut, so the 900-seed sweep is seconds rather than minutes.
    """
    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
    clean = [s for s in range(LO, HI + 1)
             if all(stage3.simulate(e, s)["outcome"] != "no_cut" for e in valid)]
    return clean, len(valid)


def measure(cell, task, window, contract, reference, cell_types, rng):
    lo, hi = window
    span = hi - lo + 1
    base = AC.config_for(cell) or AC.AllCutConfig()
    mf = base.cas12a_max_fail if span == 900 else max(1, round(base.cas12a_max_fail * span / 900))
    cfg = _dc.replace(base, start_seed=lo, end_seed=hi, cas12a_max_fail=mf)
    rec = {"cell": cell, "task": (task.get("task_id") or task["id"])[:8], "width": span,
           "window": [lo, hi], "max_fail": mf, "group_size": cfg.group_size}
    t0 = time.monotonic()
    rows, meta = AC.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=BUDGET)
    rec["build_s"] = round(time.monotonic() - t0, 1)
    rec["meta"] = {k: v for k, v in meta.items() if k != "bank_path"}
    if not rows:
        rec["declined"] = meta.get("reason", "declined")
        print(f"  w{span:<4} {lo}-{hi:<4} DECLINED after {rec['build_s']:.0f}s — {rec['declined']}")
        MT.free_gpu_memory()
        return rec

    clean, n_valid = true_clean(rows, contract, reference, cell_types)
    inside = [s for s in clean if lo <= s <= hi]
    rec.update(rows=len(rows), valid=n_valid, bank=meta.get("bank"),
               group_union=meta.get("union"), clean_in_window=meta.get("clean"),
               clean_fraction=meta.get("clean_fraction"), clean_total=len(clean),
               clean_inside=len(inside), leak=len(clean) - len(inside),
               cas9_pool=meta.get("cas9_pool"), cells=meta.get("cells"))

    dirty_pool = [s for s in range(LO, HI + 1) if s not in set(clean)]
    cl = [score(rows, contract, reference, cell_types, seed=s)
          for s in rng.sample(clean, min(N_CLEAN, len(clean)))] if clean else []
    dr = [score(rows, contract, reference, cell_types, seed=s)
          for s in rng.sample(dirty_pool, min(N_DIRTY, len(dirty_pool)))]
    base_score = (cl or dr)[0]
    p = len(clean) / (HI - LO + 1)
    cons_clean = st.mean(x["consistency"] for x in cl) if cl else 0.0
    cons_dirty = st.mean(x["consistency"] for x in dr) if dr else 0.0
    e_cons = p * cons_clean + (1 - p) * cons_dirty
    rec.update(weighted=base_score["weighted"], fidelity=base_score["fidelity"],
               wxfid=base_score["weighted"] * base_score["fidelity"],
               cons_clean=cons_clean, cons_dirty=cons_dirty, e_cons=e_cons,
               e_final=base_score["weighted"] * base_score["fidelity"] * e_cons)
    MT.free_gpu_memory()
    print(f"  w{span:<4} {lo}-{hi:<4} mf {mf:<4} bank {rec['bank']:<6} "
          f"clean {rec['clean_in_window']:>3}/{span:<4} ({rec['clean_fraction']:.1%})  "
          f"of 900: {rec['clean_total']:<3} (leak {rec['leak']})  w {rec['weighted']:6.1f} "
          f"fid {rec['fidelity']:.4f}  cons {rec['cons_clean']:.3f}/{rec['cons_dirty']:.3f}  "
          f"E[final] {rec['e_final']:6.2f}  cells {rec['cells']}/8  {rec['build_s']:.0f}s")
    return rec


def main():
    tasks = pick_tasks(fetch_tasks())
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    rng = random.Random(20260911)
    results = []
    for cell in CELLS:
        task = tasks.get(cell)
        if task is None:
            print(f"{cell}: no stamped task in the feed, skipped")
            continue
        content = task["content"]
        contract, reference = content["contract"], content["hbb_reference"]
        print(f"\n{cell}  task {(task.get('task_id') or task['id'])[:8]}  "
              f"{task.get('created_at', '')[:16]}  seeds {contract.get('seed')}")
        for width in WIDTHS:
            for window in windows_for(width):
                try:
                    results.append(measure(cell, task, window, contract, reference,
                                           cell_types, rng))
                except Exception as exc:
                    print(f"  w{width:<4} ERROR {type(exc).__name__}: {exc}")
                    results.append({"cell": cell, "width": width,
                                    "error": f"{type(exc).__name__}: {exc}"})
                json.dump(results, open(OUT, "w"), indent=1)
    json.dump(results, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
