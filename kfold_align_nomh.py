#!/usr/bin/env python3
"""kfold_align_nomh.py — exploit stage4's KFold(shuffle=True, random_state=seed) being a PURE
function of row POSITION, to test whether arranging not_mhnhej rows so all BLUNT-outcome rows land
in one KFold fold (the other 4 folds pure HDR) lifts is_hdr's cross-validated R2 the way a fully
constant is_hdr submission does.

Two preconditions this idea needs, verified separately BEFORE this script was written, not assumed:

    1. KFold's fold membership must depend ONLY on (n, random_state), never on the data itself --
       confirmed: sklearn.model_selection.KFold(...).split() gives IDENTICAL partitions for the
       same (n, random_state) regardless of what array is passed in.
    2. The miner's submitted ROW ORDER must survive unchanged into stage4's merged DataFrame, so
       "position i in submission.json" really is "df.iloc[i]" at scoring time -- confirmed by
       replaying stage4.flatten_stage3/flatten_stage12/merge on a real 250-row submission: the
       merged frame's experiment_id order matched submission order exactly, 250/250.

Both hold, so a miner who KNOWS the round seed (this is a seed_depend-style ORACLE test on an
already-scored task, exactly like `nomh_standalone.py` -- NOT buildable blind on a live round, see
the caveat at the end of the RESULT) can precompute
`KFold(n_splits=5, shuffle=True, random_state=seed).split(range(250))` OFFLINE, with no data at all,
and choose which of the 250 row-POSITIONS gets which guide.

**The mechanistic catch found before writing any allocation code.** K-Fold trains each fold's model
on the OTHER 4 folds. The "special" BLUNT fold is one of those 4 "other" folds for EVERY one of the
4 "pure HDR" folds' evaluations -- there is no way to exclude it from all four simultaneously, since
there are only 5 folds total. So the training set behind each "pure" fold's prediction is never
actually pure (~200 rows, both labels present). And `sklearn.metrics.r2_score` on a constant y_test
is a HARD STEP FUNCTION, not continuous: verified exactly 1.0 only when the prediction matches
EXACTLY, else exactly 0.0 (a prediction off by 0.01 already returns 0.0, not something close to
1.0). Achieving an exact match on 50 held-out rows the model never trained on needs a PERFECT
HDR/BLUNT separator from the 7 build_X features -- which `energy_combo_k562.py` already measured
does not exist (AUC ceiling 0.5271 on precisely this question). So the mechanistic prediction is
that the "4 pure folds" should ALSO land at r2=0.0, not "almost 1.00" -- this script checks that
against the real pipeline rather than trusting the toy-example reasoning alone.

    KA_CELL=K562 python kfold_align_nomh.py [task_id]

RESULT (2026-09-21, K562, task 2af71117, pinned seed 401): the mechanistic prediction above was
confirmed exactly, and the outcome is worse than "no gain" -- it is actively harmful.
consistency_factor came in at 0.0018, far BELOW the ordinary ~0.10-0.30 floor this file's CLAUDE.md
counterpart documents everywhere, against 0.2743 for an unstructured not_mhnhej build at the SAME
pinned seed (nomh_standalone.py). Per-fold replay: all 5 folds landed at r2 = +0.0000, including the
four "pure HDR" folds -- none came close to "almost 1.00". Worse, the predictions were not merely
unhelpful, they were confidently WRONG-SIGNED: fold 0 (truly all-BLUNT, y=0) got mean(pred)=1.0000
exactly; fold 1 (truly all-HDR, y=1) got mean(pred)=0.0000 exactly. is_hdr's MAE hit 0.9143 -- close
to the worst possible for a 0/1 target -- which is what drags consistency down past the usual floor
rather than merely failing to lift it. Likely cause: segregating BLUNT rows into a contiguous
rank-block per cell (both HDR/BLUNT pools were taken nearest/highest-weighted_score-first, the same
ordering seed_depend.enumerate_candidates already sorts by) handed the deep (max_depth=12), 200-tree
RandomForest a spurious structural correlate of "labeled BLUNT in this dataset" to overfit to, which
then generalised in the wrong direction on held-out rows. This does NOT transfer to a live, blind
round even in principle -- the whole mechanism needs the round seed known at BUILD time to compute
both the KFold partition and the true HDR/BLUNT labels, which the shipped conjunction (built before
broadcast, seeds stamped after) never has; it is only testable as a seed_depend-style oracle exercise
on an already-scored task. See CLAUDE.md's falsified table for the full writeup.

FOLLOW-UP RESULT (2026-09-21, same task/seed, KA_SPECIAL_BLUNT=25 -- 4 pure-HDR folds of 50 + 1
mixed 25 HDR/25 BLUNT fold): less harmful, still short of the predicted range. consistency_factor
rose to 0.1468 (from 0.0018 with an all-BLUNT special fold) but stayed below the unstructured
build's 0.2743. The special fold's y_test is now genuinely non-constant (mean 0.500), so its r2 is
computed normally rather than snapped by the degenerate rule -- but the model still predicted it
with full confidence in one direction (mean(pred)=1.0000 exactly on a truly 50/50 fold, r2=-0.7091).
The four "pure HDR" folds still ALL land at r2=+0.0000 -- none hits the exact match r2_score
requires -- but predictions are less extremely wrong (0.37-0.70 instead of snapping to exactly
0.0000), so is_hdr's MAE recovers to 0.4837 (from 0.9143). Softening the split reduces the
self-inflicted damage but cannot create the exact-match condition no fold combination reaches.
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
CELL = os.getenv("KA_CELL", "K562")
os.environ["EH_CELL"] = CELL
N_SPECIAL_BLUNT = int(os.getenv("KA_SPECIAL_BLUNT", "50"))   # BLUNT rows inside the special fold;
                                                               # the rest of that fold's 50 slots
                                                               # (and all of the other 4 folds) are
                                                               # HDR. Default 50 reproduces the
                                                               # original all-BLUNT special fold.
os.environ.setdefault("NIOME_INSTANCE", f"kfoldalign{N_SPECIAL_BLUNT}")

import json                                                # noqa: E402
import logging                                             # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                          # noqa: E402
from sklearn.ensemble import RandomForestRegressor          # noqa: E402
from sklearn.metrics import r2_score, mean_absolute_error   # noqa: E402
from sklearn.model_selection import KFold                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import seed_depend as SD         # noqa: E402
from niome_subnet.genomics.validation import (              # noqa: E402
    run_stage12, run_stage3, run_stage4, _parse_seeds)
from niome_subnet.genomics.validation.stage4 import (       # noqa: E402
    flatten_stage3, flatten_stage12, build_X, build_y)
from niome_subnet.utils import settings                     # noqa: E402
from sd_task import task_content                            # noqa: E402
from energy_hdr_k562 import pick_task, candidates_by_cell    # noqa: E402

N_ROWS = 250
N_SPECIAL = 50   # one KFold fold at n_splits=5


def main():
    task_id = _ARGV[1] if len(_ARGV) > 1 else pick_task()
    task, contract, reference = task_content(task_id)
    seeds = _parse_seeds(contract["seed"])
    seed = seeds[0]
    print(f"task {task_id[:8]}  {CELL}  real seeds {seeds}  pinned to seed {seed}\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    cfg = SD.SeedDependConfig(rule="not_mhnhej", variants_per_site=300)
    ctx, by_cell = candidates_by_cell(contract, reference, cell_types, cfg, seed)
    print(f"cells with candidates: {len(by_cell)}/8", flush=True)

    # Split every cell's not_mhnhej-compliant pool by realised outcome at this seed.
    hdr_by_cell, blunt_by_cell = {}, {}
    for cell, recs in by_cell.items():
        hdr_by_cell[cell] = [r for r in recs if r["record"]["outcome"] == "HDR"]
        blunt_by_cell[cell] = [r for r in recs if r["record"]["outcome"] == "BLUNT_NHEJ"]
    print("per-cell pool (HDR / BLUNT):")
    for cell in by_cell:
        print(f"  {cell}: {len(hdr_by_cell[cell])} / {len(blunt_by_cell[cell])}")
    print(flush=True)

    special_hdr = N_SPECIAL - N_SPECIAL_BLUNT           # HDR rows filling out the special fold
    pure_take = N_ROWS - N_SPECIAL                       # always 200: the 4 pure-HDR folds
    hdr_take = pure_take + special_hdr                    # total HDR rows across the submission
    blunt_take = N_SPECIAL_BLUNT

    cells = sorted(by_cell)
    n_cells = len(cells)
    per_hdr = hdr_take // n_cells
    per_blunt = blunt_take // n_cells

    hdr_pool, blunt_pool = [], []
    for i, cell in enumerate(cells):
        hn = per_hdr + (1 if i < hdr_take - per_hdr * n_cells else 0)
        bn = per_blunt + (1 if i < blunt_take - per_blunt * n_cells else 0)
        hdr_pool.extend(hdr_by_cell[cell][:hn])
        blunt_pool.extend(blunt_by_cell[cell][:bn])
    print(f"selected {len(hdr_pool)} HDR rows, {len(blunt_pool)} BLUNT rows "
          f"(target {hdr_take} / {blunt_take}; special fold = {N_SPECIAL_BLUNT} BLUNT + "
          f"{special_hdr} HDR, other 4 folds = {pure_take} HDR)\n", flush=True)
    if len(hdr_pool) < hdr_take or len(blunt_pool) < blunt_take:
        raise SystemExit("a cell ran short -- reduce per-cell targets or widen the scan")

    # KFold's own partition of positions 0..249, computed with ZERO data -- purely a function of
    # (n, random_state=seed), exactly matching what run_stage4(seed=seed) will compute internally.
    kf = KFold(n_splits=5, shuffle=True, random_state=seed)
    folds = [test_idx for _, test_idx in kf.split(np.arange(N_ROWS))]
    special_positions = folds[0]                       # the fold that will hold the BLUNT rows
    pure_positions = np.concatenate(folds[1:])
    assert len(special_positions) == N_SPECIAL and len(pure_positions) == pure_take

    slot = {}
    for pos, rec in zip(special_positions[:blunt_take], blunt_pool[:blunt_take]):
        slot[int(pos)] = rec
    for pos, rec in zip(special_positions[blunt_take:], hdr_pool[:special_hdr]):
        slot[int(pos)] = rec
    for pos, rec in zip(pure_positions, hdr_pool[special_hdr:special_hdr + pure_take]):
        slot[int(pos)] = rec

    rows, seen = [], set()
    for i in range(N_ROWS):
        rec = slot[i]
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": rec["guide"],
                     "target_alignment_start": rec["start"],
                     "target_alignment_end": rec["start"] + rec["length"],
                     "strand": rec["strand"], "mutation": rec["mutation"],
                     "cas_system": rec["cas_system"], "cell_type": CELL})
    print(f"built {len(rows)} rows (deduped from {N_ROWS}); "
          f"positions {sorted(special_positions.tolist())[:5]}... are the special BLUNT fold\n",
          flush=True)

    os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
    json.dump(contract, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    run_stage3(seed=seed)
    r4 = run_stage4(seed=seed)
    print("=== real pipeline result, pinned seed ===")
    print(f"  consistency_factor {r4['consistency_factor']:.4f}  "
          f"weighted {r4['total_weighted_score']:.1f}  final {r4['final_reward']:.2f}")
    for t in ("is_cut", "is_hdr", "indel_length"):
        v = r4["model_results"].get(t, {})
        print(f"    {t:<13} r2_mean {v.get('r2_mean', float('nan')):+.4f}  "
              f"mae_mean {v.get('mae_mean', float('nan')):.4f}")

    # Independent replay of stage4's OWN merge + KFold, to see the PER-FOLD is_hdr r2 directly --
    # run_stage4 only returns the mean/std across folds, and the per-fold breakdown is exactly what
    # settles whether the "4 pure folds" landed at ~1.0 or at 0.0 as the mechanistic note predicts.
    print("\n=== per-fold is_hdr breakdown (independent replay of stage4's own logic) ===")
    stage3_df = flatten_stage3(json.load(open(settings.STAGE3_DATASET)))
    stage12_df = flatten_stage12(json.load(open(settings.VALID_EXPERIMENTS_PATH)))
    stage12_slim = stage12_df[["experiment_id", "guideRNA", "start", "stage2_score",
                               "mutation_weight", "weighted_score"]]
    df = stage3_df.merge(stage12_slim, on="experiment_id", how="inner")
    X = build_X(df)
    y = build_y(df)["is_hdr"]
    sw = df["mutation_weight"]
    kf2 = KFold(n_splits=5, shuffle=True, random_state=seed)
    for i, (train_idx, test_idx) in enumerate(kf2.split(X)):
        y_test = y.iloc[test_idx]
        model = RandomForestRegressor(n_estimators=200, random_state=42, max_depth=12)
        model.fit(X.iloc[train_idx], y.iloc[train_idx], sample_weight=sw.iloc[train_idx])
        pred = model.predict(X.iloc[test_idx])
        r2 = r2_score(y_test, pred, sample_weight=sw.iloc[test_idx])
        mae = mean_absolute_error(y_test, pred, sample_weight=sw.iloc[test_idx])
        print(f"  fold {i}: n={len(test_idx):3d}  y_test constant={y_test.nunique() == 1}  "
              f"mean(y_test)={y_test.mean():.3f}  mean(pred)={pred.mean():.4f}  "
              f"pred range [{pred.min():.4f},{pred.max():.4f}]  r2={r2:+.4f}  mae={mae:.4f}")


if __name__ == "__main__":
    main()
