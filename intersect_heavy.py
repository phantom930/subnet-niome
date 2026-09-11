#!/usr/bin/env python3
"""intersect_heavy.py — mutation-weight skew applied to the two-group all-HDR intersection.

`intersect_score.py` measured the three 250-row splits at weighted 188-195 against a field median
of 249.6 on the same round, and the whole deficit was `total_weighted_score` rather than fidelity
(0.9833, the best number in that field). The obvious suspect is the mutation mix: `assemble`'s
apportionment is INERT in these builds -- with `need = 0` the rows are exactly the two min-union
groups -- and CLAUDE.md records that the min-union is "blind to mutation_weight and lands near
50/50", while the shipped build skews its Cas9 half to the heavy mutation via `light_cell_rows = 6`
(the 79/79/6/6 split).

So the skew has to move into the group SELECTION, through `FastGreedy`'s caps: bound each light
(mutation, cas, strand) cell and the remainder lands on the heavy mutation. Floors still win, so a
cap can never empty a stage-5 cell.

The cap is swept rather than fixed at the shipped 6, because 6 was tuned for a ~170-row Cas9 half
and these groups are 80-170 rows of a DIFFERENT construction. Two things pull against each other:

  * skewing raises mean `mutation_weight` and with it `total_weighted_score`
  * constraining the min-union costs it freedom, so the union grows and the BAND shrinks -- and
    CLAUDE.md found group caps "strictly harmful everywhere" in the shipped context, where a
    separate Cas9 half did the apportionment instead

Reports band, weighted, fidelity and w x fid for every (split, cap), against the uncapped baseline.

    python intersect_heavy.py
    IH_LIGHT=none,6 IH_SPLITS=125/125 python intersect_heavy.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "isectH")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.all_cut import bank_key, load_bank   # noqa: E402
from niome_subnet.genomics.validation import stage3          # noqa: E402
from niome_subnet.genomics.validation.stage12 import run_stage12   # noqa: E402
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity  # noqa: E402
from niome_subnet.utils import settings                      # noqa: E402
from allhdr_mf import pick                                   # noqa: E402
from joined300 import fetch_tasks                            # noqa: E402
from sd_task import score                                    # noqa: E402

WINDOW = tuple(int(x) for x in os.getenv("IH_WINDOW", "100-129").split("-"))
MF9, MF12 = int(os.getenv("IH_MF9", "10")), int(os.getenv("IH_MF12", "12"))
SPLITS = [tuple(int(x) for x in s.split("/"))
          for s in os.getenv("IH_SPLITS", "125/125,80/170,170/80").split(",")]
LIGHTS = [None if x.lower() in ("none", "") else int(x)
          for x in os.getenv("IH_LIGHT", "none,6,12,25").split(",")]
OUT = os.getenv("IH_OUT", "intersect_heavy.json")


def caps_for(contract, ctx, cas, light):
    """Cap every light (mutation, cas, strand) cell at `light`; the heavy mutation takes the rest."""
    if light is None:
        return None
    weights = contract.get("mutation_weights", {})
    heavy = max(ctx.mutations, key=lambda m: weights.get(m, 1.0))
    return {(m, cas, s): light for m in ctx.mutations if m != heavy for s in ("+", "-")}


def group_from(path, size, cfg, cas, contract, ctx, light):
    records = load_bank(path, limit=300_000)
    for rec in records:              # load_bank hardcodes Cas12a; see intersect_score.py
        rec["cas_system"] = cas
    sel = FG.FastGreedy(records, window_lo=WINDOW[0], window_hi=WINDOW[1],
                        per_cell_min=cfg.per_cell_min,
                        caps=caps_for(contract, ctx, cas, light))
    idx, _u = sel.best(size, restarts=cfg.restarts)
    group = [records[i] for i in idx]
    bad = set()
    for rec in group:
        bad.update(int(x) for x in rec["fails"])
    return group, sorted(set(range(WINDOW[0], WINDOW[1] + 1)) - bad)


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context(contract, reference, cell_types)
    n_rows = contract["rules"].get("max_experiments") or 250
    base = AH.config_for("K562")
    cfg9 = _dc.replace(base, hdr_range=WINDOW, main_max_fail=MF9)
    cfg12 = _dc.replace(base, hdr_range=WINDOW, main_max_fail=MF12)
    p9 = os.path.join(AH.HDR_BANK_DIR, f"cas9-{bank_key(contract, cell_types, cfg9)}.npz")
    p12 = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg12)}.npz")
    weights = contract.get("mutation_weights", {})
    heavy = max(ctx.mutations, key=lambda m: weights.get(m, 1.0))
    print(f"K562  task {(task.get('task_id') or task['id'])[:8]}  window {WINDOW[0]}-{WINDOW[1]}"
          f"  mutation_weights {weights}  heavy = {heavy}\n")
    print(f"{'split':>9}{'light cap':>11}{'heavy rows':>12}{'band':>6}{'weighted':>10}"
          f"{'fidelity':>10}{'w x fid':>10}{'mut ratio':>11}")
    out = []
    for g9, g12 in SPLITS:
        for light in LIGHTS:
            grp9, b9 = group_from(p9, g9, cfg9, "Cas9", contract, ctx, light)
            grp12, b12 = group_from(p12, g12, cfg12, "Cas12a", contract, ctx, light)
            rows = AC.assemble(grp9 + grp12, [], contract, ctx, cfg12, n_rows)
            os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
            json.dump(contract, open(settings.CONTRACT_PATH, "w"))
            json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
            json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
            run_stage12(cell_types)
            valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
            band = [s for s in range(100, 1000)
                    if all(stage3.simulate(e, s)["outcome"] == "HDR" for e in valid)]
            off = [s for s in range(100, 1000) if s not in set(band)]
            dirty = score(rows, contract, reference, cell_types, seed=off[0])
            mix = Counter(r["mutation"] for r in rows)
            n_heavy = mix.get(heavy, 0)
            # `sd_task.score` computes fidelity in-process; stage 5 is never run here, so
            # DISTRIBUTION_FIDELITY_PATH does not exist. Recompute the terms the same way.
            d = compute_distribution_fidelity(
                valid, json.load(open(settings.STAGE3_DATASET)), contract, k=12)
            w, fid = dirty["weighted"], dirty["fidelity"]
            print(f"{f'{g9}/{g12}':>9}{str(light):>11}{n_heavy:>6}/{len(rows):<5}"
                  f"{len(band):>6}{w:>10.1f}{fid:>10.4f}{w*fid:>10.2f}"
                  f"{d.get('mutation_coverage_entropy_ratio', float('nan')):>11.3f}")
            out.append({"split": f"{g9}/{g12}", "light_cap": light, "rows": len(rows),
                        "heavy_rows": n_heavy, "cells": len(Counter(
                            (r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
                        "band_cas9": len(b9), "band_cas12a": len(b12), "band": len(band),
                        "band_seeds": band, "weighted": w, "fidelity": fid, "wxfid": w * fid,
                        "mutation_ratio": d.get("mutation_coverage_entropy_ratio"),
                        "cas_ratio": d.get("cas_system_coverage_entropy_ratio"),
                        "joint_ratio": d.get("joint_coverage_entropy_ratio"),
                        "cons_off_band": dirty["consistency"]})
            MT.free_gpu_memory()
            json.dump(out, open(OUT, "w"), indent=1)
        print()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
