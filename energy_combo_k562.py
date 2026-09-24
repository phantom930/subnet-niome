#!/usr/bin/env python3
"""energy_combo_k562.py — does {consistency, gc, dist, energy, mh} beat energy alone for is_hdr?

**Prediction, from reading `stage3.repair_mode(cas, energy, mh, rng)` directly**: no combination of
these 5 features can beat `energy` alone (conditional on `not_mhnhej`), because the function only
takes `cas`, `energy`, `mh` and a fresh draw as arguments -- `gc`/`distance`/`dist_score`/
`consistency` reach the outcome ONLY through `energy` (a deterministic transform of them), and
conditioning on `not_mhnhej` (missing the `mh_nhej` draw segment) makes the HDR-vs-BLUNT split
`hdr/(hdr+blunt)`, in which `mh_nhej`'s width -- and therefore `mh` -- cancels out exactly. So
`energy` is a sufficient statistic for all four of the others, for this specific question, by
construction rather than by measurement.

This is the population-scale check of that prediction, not a 240-row submission-scale one (that
was `energy_hdr_k562.py`, which already found a small submission-scale signal was invisible or
harmful once selected into an actual build). Reusing the same not_mhnhej candidate pool
enumeration, at proper sample size (thousands of candidates, not 240), fits:

    A  logistic regression on energy alone
    B  logistic regression on all 5 named features (consistency, gc, distance, energy, mh)
    C  the CLOSED-FORM formula hdr/(hdr+blunt) directly, as the ceiling neither model can beat

and reports cross-validated AUC/log-loss for A vs B (expect indistinguishable, per the derivation
above) against C's theoretical ceiling (expect A and B to approach it, since there IS a real
relationship, just a narrow one on K562 -- see `energy_hdr_k562.py`'s finding that the achievable
energy range on this cell only spans a ~3-4 point P(HDR) shift).

Memory: reuses `seed_depend.enumerate_candidates` at `variants_per_site=300` (see
`energy_hdr_k562.py`'s docstring for why the shipped 12000 OOM-killed this box's shared fleet).

    python energy_combo_k562.py [task_id]

RESULT (2026-09-20, K562, task 0441422d, 193,135 not_mhnhej-compliant candidates, all 8 cells):
the prediction holds exactly. AUC: A (energy alone) 0.5226 +/- 0.0015, B (all 5 features)
0.5229 +/- 0.0016 -- a gap of +0.0003, statistically indistinguishable, against a closed-form
theoretical ceiling C of only 0.5271 (barely above chance, since K562's energy sits mean 0.9149
with 57% of the pool clamped at 1.0 -- the same saturation `energy_hdr_k562.py` found on this
cell). Model B's own coefficients confirm the mechanism rather than just the aggregate score:
distance fits to EXACTLY 0.0000, gc/consistency get small wrong-signed coefficients
(-0.0523/-0.0333), mh fits to +0.0009 (indistinguishable from its algebraic zero), and energy
absorbs essentially the whole signal (+1.0181, matching model A's own +0.9015). So there is no
sample size or feature combination that beats energy alone here -- the ceiling was already energy
alone, confirmed at population scale (800x the 240-row submission scale) rather than assumed from
the formula. See CLAUDE.md's falsified-table entry for the mathematical derivation and the
companion HEK293/K562 energy-only results this extends.
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "energycombo")

import logging                                          # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                        # noqa: E402
from sklearn.linear_model import LogisticRegression       # noqa: E402
from sklearn.model_selection import cross_val_score, StratifiedKFold  # noqa: E402
from sklearn.metrics import roc_auc_score                 # noqa: E402

import genExp as G                                        # noqa: E402
from niome_subnet.genomics import seed_depend as SD       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds  # noqa: E402
from sd_task import fetch, task_content                   # noqa: E402
from energy_hdr_k562 import candidates_by_cell, pick_task  # noqa: E402

CELL = os.getenv("EH_CELL", "K562")


def theoretical_p_hdr(cas, energy):
    hdr_base = 0.32 if cas == "Cas9" else 0.24
    hdr = hdr_base + 0.35 * energy
    return hdr / (hdr + 0.35)


def main():
    task_id = _ARGV[1] if len(_ARGV) > 1 else pick_task()
    task, contract, reference = task_content(task_id)
    seeds = _parse_seeds(contract["seed"])
    seed = seeds[0]
    print(f"task {task_id[:8]}  {CELL}  real seeds {seeds}  testing at seed {seed}\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    cfg = SD.SeedDependConfig(rule="not_mhnhej", variants_per_site=300)
    ctx, by_cell = candidates_by_cell(contract, reference, cell_types, cfg, seed)

    rows = []
    for cell, recs in by_cell.items():
        for r in recs:
            rec = r["record"]
            feat = rec["features"]
            rows.append({
                "consistency": feat["consistency"],
                "gc": feat["gc"],
                "distance": feat["distance"],
                "energy": rec["energy"],
                "mh": int(rec["mh"]),
                "cas": rec["cas"],
                "is_hdr": 1 if rec["outcome"] == "HDR" else 0,
            })
    n = len(rows)
    print(f"pool size: {n} not_mhnhej-compliant candidates across {len(by_cell)}/8 cells\n",
          flush=True)

    energy = np.array([r["energy"] for r in rows])
    is_hdr = np.array([r["is_hdr"] for r in rows])
    X5 = np.array([[r["consistency"], r["gc"], r["distance"], r["energy"], r["mh"]]
                   for r in rows])
    theo = np.array([theoretical_p_hdr(r["cas"], r["energy"]) for r in rows])

    print(f"observed P(HDR | not_mhnhej) overall: {is_hdr.mean():.4f}")
    print(f"theoretical P(HDR | not_mhnhej) (formula, mean over pool): {theo.mean():.4f}")
    print(f"energy: mean {energy.mean():.4f} std {energy.std():.4f} "
          f"min {energy.min():.4f} max {energy.max():.4f}\n", flush=True)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # A: energy alone
    Xa = energy.reshape(-1, 1)
    model_a = LogisticRegression(max_iter=1000)
    auc_a = cross_val_score(model_a, Xa, is_hdr, cv=cv, scoring="roc_auc")
    model_a.fit(Xa, is_hdr)

    # B: all 5 named features
    model_b = LogisticRegression(max_iter=1000)
    auc_b = cross_val_score(model_b, X5, is_hdr, cv=cv, scoring="roc_auc")
    model_b.fit(X5, is_hdr)

    # C: the closed-form theoretical probability itself, scored as a "prediction"
    auc_c = roc_auc_score(is_hdr, theo)

    print("=== cross-validated AUC (5-fold, is_hdr as target) ===")
    print(f"  A energy alone            : {auc_a.mean():.4f} +/- {auc_a.std():.4f}  {auc_a}")
    print(f"  B all 5 (cons,gc,dist,e,mh): {auc_b.mean():.4f} +/- {auc_b.std():.4f}  {auc_b}")
    print(f"  C closed-form formula (ceiling, no fitting): {auc_c:.4f}")
    print(f"\n  B - A = {auc_b.mean() - auc_a.mean():+.4f} "
          f"(prediction: ~0, since mh/gc/distance/consistency add nothing beyond energy)")

    print("\n=== model B's fitted coefficients (standardized would be better, raw shown) ===")
    for name, coef in zip(["consistency", "gc", "distance", "energy", "mh"], model_b.coef_[0]):
        print(f"  {name:<12} {coef:+.4f}")
    print(f"  (model A's energy coefficient alone: {model_a.coef_[0][0]:+.4f})")


if __name__ == "__main__":
    main()
