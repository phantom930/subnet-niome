#!/usr/bin/env python3
"""mixed_split.py — the SHIPPED all-HDR construction at width 30, at three Cas9/Cas12a splits.

`intersect_score.py` / `intersect_heavy.py` built the two-group version: min-union Cas9 and Cas12a
independently and ship both, so the band is the INTERSECTION of two independent bands (10 seeds,
matching the independence prediction). This runs the construction all_hdr actually ships, at the
same width and the same three splits, so the two are comparable:

    two-group   Cas9 min-union  AND  Cas12a min-union     -> band = intersection
    mixed       Cas12a min-union, Cas9 required to reach HDR on THAT band -> band preserved

The mixed form should be strictly wider by construction -- the Cas9 half is a conditional fill, not
a second independent draw -- at the cost of a Cas9 pool that decays as P(HDR) per band seed, which
is the constraint CLAUDE.md's "the band cannot be widened" entry is about. Whether that fill
survives a 14-15 seed band at width 30 is the open question here; at the shipped width 225 the band
is 12 and the pool 1363.

The split is just `group_size`, since the Cas9 half is `n_rows - group_size`:

    125/125 (Cas9/Cas12a) -> group_size 125      170/80 -> group_size 80      80/170 -> group_size 170

`light_cell_rows` is LIVE in this path (unlike the two-group builds, where `need = 0` made
`assemble`'s apportionment inert), so both the shipped 6 and None are run.

    python mixed_split.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "mixsplit")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import time                          # noqa: E402
from collections import Counter      # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_hdr as AH             # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics.validation import stage3         # noqa: E402
from niome_subnet.genomics.validation.stage12 import run_stage12   # noqa: E402
from niome_subnet.genomics.validation.stage5 import compute_distribution_fidelity  # noqa: E402
from niome_subnet.utils import settings                     # noqa: E402
from allhdr_mf import pick                                  # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402
from sd_task import score                                   # noqa: E402

WINDOW = tuple(int(x) for x in os.getenv("MS_WINDOW", "100-129").split("-"))
MF = int(os.getenv("MS_MF", "12"))
# (label, group_size) -- group_size is the Cas12a half; Cas9 takes n_rows - group_size.
# MS_GROUPS lists Cas12a group sizes; the label is derived as <Cas9>/<Cas12a> at 250 rows.
SPLITS = [("%d/%d" % (250 - g, g), g)
          for g in (int(x) for x in os.getenv("MS_GROUPS", "125,80,170").split(","))]
LIGHTS = [None if x.lower() in ("none", "") else int(x)
          for x in os.getenv("MS_LIGHT", "6,none").split(",")]
BUDGET = float(os.getenv("MS_BUDGET", "1200"))
OUT = os.getenv("MS_OUT", "mixed_split.json")
TERMS = ["mutation_coverage_entropy_ratio", "cas_system_coverage_entropy_ratio",
         "strand_coverage_entropy_ratio", "joint_coverage_entropy_ratio",
         "kmer_diversity_entropy_ratio", "distinct_guide_ratio"]


def main():
    task = pick(fetch_tasks())
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    n_rows = contract["rules"].get("max_experiments") or 250
    base = AH.config_for("K562")
    weights = contract.get("mutation_weights", {})
    print(f"K562  task {(task.get('task_id') or task['id'])[:8]}  window {WINDOW[0]}-{WINDOW[1]}"
          f"  mf {MF}  rows {n_rows}  weights {weights}\n")
    print(f"{'split':>9}{'grp(C12a)':>11}{'light':>7}{'band':>6}{'cas9 pool':>11}{'mix':>12}"
          f"{'heavy':>7}{'weighted':>10}{'fid':>9}{'w x fid':>10}{'s':>6}")
    out = []
    for label, g12 in SPLITS:
        for light in LIGHTS:
            t0 = time.monotonic()
            cfg = _dc.replace(base, hdr_range=WINDOW, main_max_fail=MF, group_size=g12,
                              light_cell_rows=light)
            rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=BUDGET)
            el = time.monotonic() - t0
            rec = {"split": label, "group_cas12a": g12, "light_cell_rows": light, "mf": MF,
                   "meta": {k: v for k, v in meta.items() if k != "bank_path"},
                   "build_s": round(el, 1)}
            if not rows:
                rec["declined"] = meta.get("reason")
                print(f"{label:>9}{g12:>11}{str(light):>7}   DECLINED — {rec['declined']}"
                      f"  [{el:.0f}s]")
                out.append(rec); json.dump(out, open(OUT, "w"), indent=1)
                MT.free_gpu_memory(); continue
            os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
            json.dump(contract, open(settings.CONTRACT_PATH, "w"))
            json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
            json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
            run_stage12(cell_types)
            valid = json.load(open(settings.VALID_EXPERIMENTS_PATH))
            band = [s for s in range(100, 1000)
                    if all(stage3.simulate(e, s)["outcome"] == "HDR" for e in valid)]
            off = [s for s in range(100, 1000) if s not in set(band)]
            hit = score(rows, contract, reference, cell_types, seed=band[0]) if band else None
            dirty = score(rows, contract, reference, cell_types, seed=off[0])
            d = compute_distribution_fidelity(
                valid, json.load(open(settings.STAGE3_DATASET)), contract, k=12)
            mix = Counter(r["cas_system"] for r in rows)
            heavy = max(contract.get("mutation_weights", {}), key=lambda m: weights.get(m, 1.0)) \
                if weights else None
            n_heavy = sum(1 for r in rows if r["mutation"] == heavy)
            w, fid = dirty["weighted"], dirty["fidelity"]
            rec.update(rows=len(rows), valid=len(valid), band=len(band), band_seeds=band,
                       cas9_pool=meta.get("cas9_pool"), cells=meta.get("cells"),
                       cas_mix=dict(mix), heavy_rows=n_heavy, weighted=w, fidelity=fid,
                       wxfid=w * fid, terms={t: d.get(t) for t in TERMS},
                       cons_on_band=hit["consistency"] if hit else None,
                       final_on_band=hit["final"] if hit else None,
                       cons_off_band=dirty["consistency"])
            mixs = "%d/%d" % (mix.get("Cas9", 0), mix.get("Cas12a", 0))
            print(f"{label:>9}{g12:>11}{str(light):>7}{len(band):>6}"
                  f"{meta.get('cas9_pool', 0):>11}{mixs:>12}{n_heavy:>7}"
                  f"{w:>10.1f}{fid:>9.4f}{w*fid:>10.2f}{el:>6.0f}")
            out.append(rec)
            json.dump(out, open(OUT, "w"), indent=1)
            MT.free_gpu_memory()
        print()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
