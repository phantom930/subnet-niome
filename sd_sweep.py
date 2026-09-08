#!/usr/bin/env python3
"""sd_sweep.py — hunt the ~1.8 weighted points seed-depend now gives up to the field.

On a66f01fa (2026-09-08) the field's best reached weighted 281.5 / fid 0.9342 against our 279.7 /
0.9334, which moved the four siblings from ranks 1-4 (85% of the curve) to 6-9 (9%). The gap is
tiny and it is on BOTH terms, so it is not the mutation split — ``allocate`` already searches that,
refining to granularity 2 around the optimum — it is the candidate pool the split draws from.

Configs tested, one knob each, paired within contract at a fixed variant so the tie-break noise
(0.03-0.14 points) does not confound the comparison:

  base        what ships today
  hdr         the other rule. Chosen as "mh" on a single measurement (331.70 vs 330.40 on
              8f02f1a4, before the k-mer fix); the two admit different guide populations, so the
              weighted available differs per contract and one round is not evidence.
  vps12k      3x the guide variants per site — more chances at a high-GC, near-mutation guide that
              still satisfies the rule. The most direct route to a better pool.
  step1       finest split search. Cheap to test, and the note that the k-mer price moves the peak
              by "a couple of rows" says granularity 2 may be leaving something.
  greedy1200  3x the candidates the greedy considers per cell.

Build time is reported because it is a constraint, not a curiosity: the prefetch path allows 900s
but SEED_DEPEND_MIN_BUDGET_S is 190s, so a config that cannot finish in ~190s is unavailable on any
round whose prefetch failed.

    python sd_sweep.py <task_id> [config ...]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

_TASK = _ARGV[1] if len(_ARGV) > 1 else "a66f01fa-9795-4044-a243-15e33cdf9f81"
os.environ["NIOME_INSTANCE"] = f"sdsweep_{_TASK[:8]}"   # before settings resolves DATA_DIR

import dataclasses  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import time  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G  # noqa: E402
from niome_subnet.genomics import seed_depend as SD  # noqa: E402
from sd_task import DIST, field, score, task_content  # noqa: E402

TASK = _TASK
VARIANT = int(os.getenv("SD_VARIANT", "1"))
BUDGET = float(os.getenv("SD_BUDGET", "1800"))          # research budget, not the miner's 900
OUT = f"sd_sweep_{TASK[:8]}.json"

CONFIGS = {
    "base": {},
    "hdr": {"rule": "hdr"},
    "vps12k": {"variants_per_site": 12000},
    "step1": {"alloc_step": 1},
    "greedy1200": {"greedy_window": 1200},
    # The two winners together. hdr (+0.53) and vps12k (+1.47) act on different stages -- which
    # guides qualify, and how many variants of each site are tried -- so they may be additive.
    # On a66f01fa vps12k alone lands 0.02 below the field's best, and that round is the whole
    # question: rank 5 pays 5% where rank 1 starts a four-sibling sweep worth 85%.
    "hdr_vps12k": {"rule": "hdr", "variants_per_site": 12000},
    # Does the pool keep paying, or is 12k already the plateau? Build time is the constraint:
    # 12k costs ~300s against SEED_DEPEND_MIN_BUDGET_S of 190s, so this only ever runs on the
    # prefetch path.
    "vps24k": {"variants_per_site": 24000},
}


def rows_for(name, overrides, contract, reference, cell_types):
    # The default config at variant 1 is exactly what sd_task.py already built and cached.
    cache = (f"sd_task_rows_{TASK[:8]}.json" if name == "base" and VARIANT == 1
             else f"sd_sweep_rows_{TASK[:8]}_{name}_v{VARIANT}.json")
    if os.path.exists(cache):
        got = json.load(open(cache))
        return got["rows"], got["meta"], 0.0
    cfg = dataclasses.replace(SD.SeedDependConfig(), variant=VARIANT, **overrides)
    t0 = time.monotonic()
    rows, meta = SD.build(contract, reference, cell_types, seed=0, cfg=cfg, budget_s=BUDGET)
    secs = time.monotonic() - t0
    if not rows:
        return None, meta, secs
    json.dump({"rows": rows, "meta": meta}, open(cache, "w"))
    return rows, meta, secs


def main():
    names = _ARGV[2:] or list(CONFIGS)
    task, contract, reference = task_content(TASK)
    seeds = str(contract.get("seed", "")).split(",")
    if not (len(seeds) == 1 and seeds[0].strip() in ("0", "")):
        raise SystemExit(f"{TASK} is stamped ({contract.get('seed')}); not a seed-depend round")

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    fld = field(TASK)                       # excludes our own hotkeys
    best_other = fld[0][0] if fld else 0.0

    print(f"task {TASK[:8]}  {contract.get('cell_type')}  "
          f"opened {task.get('created_at','?')[:16]}  variant {VARIANT}")
    print(f"  field best (excluding ours) {best_other:.2f}, "
          f"{sum(1 for _f, c in [(r[0], r[1]) for r in fld] if c >= 0.999)} at cons 1.000\n")
    print(f"  {'config':<12} {'weighted':>9} {'fid':>7} {'final':>8} {'vs base':>8} "
          f"{'rank':>5} {'share':>6} {'build':>7}")

    out, base_final = {}, None
    for name in names:
        rows, meta, secs = rows_for(name, CONFIGS[name], contract, reference, cell_types)
        if not rows:
            print(f"  {name:<12} declined: {meta.get('reason', 'unknown')}  ({secs:.0f}s)")
            out[name] = {"declined": meta.get("reason", "unknown")}
            continue
        s = score(rows, contract, reference, cell_types, seed=0)
        if name == "base":
            base_final = s["final"]
        rank = sum(1 for r in fld if r[0] > s["final"]) + 1
        share = DIST[rank - 1] if rank <= len(DIST) else 0.0
        delta = "" if base_final is None else f"{s['final'] - base_final:+8.2f}"
        out[name] = {**{k: s[k] for k in ("weighted", "consistency", "fidelity", "final")},
                     "rank": rank, "share": share, "build_s": round(secs, 1),
                     "heavy": meta.get("heavy")}
        print(f"  {name:<12} {s['weighted']:>9.1f} {s['fidelity']:>7.4f} {s['final']:>8.2f} "
              f"{delta:>8} {rank:>5} {share:>6.1%} {secs:>6.0f}s")

    # Merge, don't overwrite: a run naming a subset of configs would otherwise erase the rest and
    # they would read as "declined" in the aggregate — a measurement that never happened looking
    # exactly like one that failed.
    try:
        prior = json.load(open(OUT)).get("configs", {})
    except (FileNotFoundError, json.JSONDecodeError):
        prior = {}
    json.dump({"task": TASK, "cell": contract.get("cell_type"),
               "field_best": best_other, "variant": VARIANT, "configs": {**prior, **out}},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
