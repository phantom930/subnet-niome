#!/usr/bin/env python3
"""energy_hdr_k562.py — can `energy` distinguish HDR from BLUNT_NHEJ inside not_mhnhej?

Cell type is `EH_CELL` (default K562, env-overridable) -- the filename is K562-specific because
that is where the question started, but the mechanism and the script are not cell-specific.

**RESULT on K562 (2026-09-20, one contract, real stamped seed): no — the dependence is real but the
effect size is unusably small on this cell.** The natural (nearest-first) not_mhnhej candidate pool
saturates `energy` to exactly 1.0000 with ZERO variance across all 240 rows (matching CLAUDE.md's
"energy is CLAMPED at 1.0 on every erythroid row"), so `Pearson(energy, is_hdr) = 0` by mathematical
necessity. Deliberately forcing variance (half nearest/half farthest, energy 0.7259-1.0000) did not
help: `Pearson` came out **-0.0362** (wrong sign) and is_hdr's R2 got WORSE (-0.1528 vs the
zero-variance baseline's -0.1092), because the achievable range only implies a ~3.5-4.2 point
theoretical P(HDR) shift against ~5.3 points of sampling noise at n=80/tercile. Not a sample-size
problem — more seeds sharpen the estimate toward that same small shift, not reveal a bigger one.
See CLAUDE.md's falsified table for the full writeup; do not re-run without first widening
`max_distance` well past 600 to reach genuinely lower-energy sites (untested, and it costs
`total_weighted_score` — the moderate widening already tried here dropped it 216 -> 112).

**RESULT on HEK293 (2026-09-20, one contract, real stamped seed): still no, for a DIFFERENT
reason.** HEK293's low accessibility (0.35) means energy never saturates -- the whole pool spans
0.32-0.64 with zero candidates at the 1.0 clamp -- so real variance exists even unforced, and the
`spread` arm finally shows the theoretically-predicted POSITIVE correlation (`Pearson` **+0.1229**,
tercile means trending the right way). It still fails: is_hdr's R2 gets WORSE (-0.1004 -> -0.5070),
`indel_length`'s collapses alongside it, and `total_weighted_score` more than halves (270.9 ->
161.8) -- reaching for energy variance drags in farther/lower-GC candidates that vary on several
OTHER features too, and that heterogeneity gives the cross-validated RandomForest more to overfit
on than the energy signal gives it to generalise from. `default` (no energy engineering) wins on
`consistency_factor` on BOTH cells (K562 0.3058, HEK293 0.2896). Same verdict, different failure
mode -- do not re-run this direction on either cell without a plan for the overfitting problem, not
just the variance problem.

Mechanistic derivation, straight from `stage3.repair_mode`:

    hdr_w   = hdr_base + 0.35*energy      (hdr_base 0.32 Cas9 / 0.24 Cas12a)
    blunt_w = 0.35                        (constant, no energy term)
    P(HDR | not_mhnhej) = hdr_w / (hdr_w + blunt_w)

which is a genuine, deterministic function of energy: 0.41-0.48 at energy=0 up to 0.63-0.66 at
energy=1 (Cas12a/Cas9). That is real, not spurious — and `energy` is literally one of stage4's
seven `build_X` input features (`gc, distance, gc_score, dist_score, consistency, energy, mh`), so
in principle the RandomForest could exploit it to lift is_hdr's R2 above the near-zero/negative
floor CLAUDE.md documents everywhere else for that target.

The catch, also mechanistic: a feature that does not VARY across a submission's own 250 rows
cannot explain variance in a target, however real the underlying dependence is. Erythroid cells'
high accessibility already saturates `energy` toward its 1.0 clamp for most rows (CLAUDE.md,
"energy is CLAMPED at 1.0 on every erythroid row") -- confirmed on K562 (accessibility 0.77) below.
**HEK293 is the interesting counter-case**: accessibility 0.35 means `energy = 0.35 * (1.8*gc +
0.6*exp(-dist/1500) + offset)` rarely reaches the clamp at all -- a close, high-GC guide lands
around 0.5-0.6, not 1.0 -- so the natural (unforced) candidate pool may already carry real energy
variance where K562's cannot. That is the reason to run this file again with `EH_CELL=HEK293`
rather than assume the K562 verdict transfers.

So the question this asks is not "does energy affect
is_hdr" (it provably does) but "does a not_mhnhej-compliant submission carry enough ENERGY VARIANCE,
correlated with the REALISED outcome at a real seed, for stage4's cross-validated RandomForest to
turn that into a measurably better is_hdr R2" — and it tests three row-selection strategies to
separate the two:

    default       nearest-first per cell (what a real build does) — some natural energy spread
    spread        half nearest (high energy) + half farthest-but-still-compliant (low energy)
                  per cell — deliberately maximises energy variance
    uniform_high  all nearest per cell — deliberately ZERO energy variance, as a negative control

All three are built via `seed_depend`'s not_mhnhej rule (added alongside this script) at the
contract's REAL, already-stamped seed — oracle knowledge, exactly like the earlier seed-depend
pricing work, used here purely to get ground truth for the mechanism question, not as a claim that
this is buildable on a live round.

    EH_CELL=HEK293 python energy_hdr_k562.py [task_id]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "energyhdr")

import json                                             # noqa: E402
import logging                                          # noqa: E402
from collections import defaultdict                     # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA    # noqa: E402
from niome_subnet.genomics import seed_depend as SD      # noqa: E402
from niome_subnet.genomics.validation import (           # noqa: E402
    run_stage12, run_stage3, run_stage4, _parse_seeds)
from niome_subnet.utils import settings                  # noqa: E402
from sd_task import fetch, task_content                  # noqa: E402

API = "https://niome-api.genomes.io/api/v3"
CELL = os.getenv("EH_CELL", "K562")
N_ROWS_PER_CELL = 30  # x 8 cells = 240, close to the real 250 and evenly apportioned on purpose --
                      # this is a mechanism test, not a weighted_score optimisation.


def pick_task():
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") == CELL and len(_parse_seeds(c.get("seed"))) == 3:
            return t.get("task_id") or t["id"]
    raise SystemExit("no stamped 3-seed K562 task found")


def candidates_by_cell(contract, reference, cell_types, cfg, seed):
    """Every not_mhnhej-compliant (site, mutation, guide) at the real seed, keyed by stage-5 cell."""
    pinned = dict(contract)
    pinned["seed"] = seed
    ctx = G.build_context(pinned, reference, cell_types)
    sites = G.enumerate_sites(ctx, cfg.flank, cfg.lengths)
    by_cell = SD.enumerate_candidates(ctx, sites, cfg, None)
    return ctx, by_cell


def make_rows(by_cell, n_per_cell, mode):
    """Select n_per_cell rows per cell under one of the three strategies, from records already
    carrying `record["energy"]`/`record["record"]["outcome"]` at the real seed."""
    rows, kept = [], []
    for cell, recs in by_cell.items():
        if len(recs) < n_per_cell:
            return None, None  # this cell can't fill its quota; caller reports the shortfall
        ordered = sorted(recs, key=lambda r: r["record"]["energy"])  # ascending energy
        if mode == "uniform_high":
            chosen = ordered[-n_per_cell:]
        elif mode == "spread":
            half = n_per_cell // 2
            chosen = ordered[:half] + ordered[-(n_per_cell - half):]
        else:  # "default": nearest-first, i.e. as scanned (already appended in distance order)
            chosen = recs[:n_per_cell]
        kept.extend(chosen)
    seen = set()
    for i, rec in enumerate(kept):
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": rec["guide"],
                     "target_alignment_start": rec["start"],
                     "target_alignment_end": rec["start"] + rec["length"],
                     "strand": rec["strand"], "mutation": rec["mutation"],
                     "cas_system": rec["cas_system"], "cell_type": CELL})
    return rows, kept


def pearson(xs, ys):
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    if xs.std() == 0 or ys.std() == 0:
        return 0.0
    return float(np.corrcoef(xs, ys)[0, 1])


def main():
    task_id = _ARGV[1] if len(_ARGV) > 1 else pick_task()
    task, contract, reference = task_content(task_id)
    seeds = _parse_seeds(contract["seed"])
    seed = seeds[0]
    print(f"task {task_id[:8]}  {CELL}  real seeds {seeds}  testing at seed {seed}\n", flush=True)

    cell_types = G.fetch_cell_types()
    G.load_sequence()
    # variants_per_site is the shipped 12000 by default. not_mhnhej admits ~70-85% of draws (far
    # more than hdr's ~40-55%), and `enumerate_candidates` keeps a full dict (entry + kmers +
    # record) per SURVIVOR, not per guide scanned -- at 12000 that is enough surviving candidates
    # across ~50-200 in-range sites to reach tens of GB of RSS. This box's live 10-hotkey fleet
    # already holds ~26 GB of a 49 GB total, so a mechanism check does not get to spend the
    # `hdr_saturate7.py`-scale budget; 300 is ample for the ~30/cell this needs.
    cfg = SD.SeedDependConfig(rule="not_mhnhej", variants_per_site=300)

    ctx, by_cell = candidates_by_cell(contract, reference, cell_types, cfg, seed)
    print(f"cells with candidates: {len(by_cell)}/8, sizes "
          f"{sorted(len(v) for v in by_cell.values())} "
          f"(total {sum(len(v) for v in by_cell.values())})", flush=True)
    all_energy = [r["record"]["energy"] for recs in by_cell.values() for r in recs]
    print(f"whole pool energy: mean {np.mean(all_energy):.4f} std {np.std(all_energy):.4f} "
          f"min {min(all_energy):.4f} max {max(all_energy):.4f} "
          f"(<1.0 count: {sum(1 for e in all_energy if e < 0.9999)}/{len(all_energy)})\n",
          flush=True)

    for mode in ("default", "spread", "uniform_high"):
        rows, kept = make_rows(by_cell, N_ROWS_PER_CELL, mode)
        if rows is None:
            print(f"=== {mode}: a cell has fewer than {N_ROWS_PER_CELL} candidates, skipped ===\n",
                  flush=True)
            continue

        energies = [r["record"]["energy"] for r in kept]
        is_hdr = [1.0 if r["record"]["outcome"] == "HDR" else 0.0 for r in kept]
        corr = pearson(energies, is_hdr)
        print(f"=== {mode}: {len(rows)} rows ({len(kept)} pre-dedup) ===")
        print(f"  energy: mean {np.mean(energies):.4f} std {np.std(energies):.4f} "
              f"min {min(energies):.4f} max {max(energies):.4f}")
        print(f"  raw P(HDR) overall: {np.mean(is_hdr):.4f}")
        # bin by energy tercile to show the theoretical relationship against the observed one
        order = np.argsort(energies)
        n = len(order)
        for name, idx in (("low third", order[:n // 3]), ("mid third", order[n // 3:2 * n // 3]),
                          ("high third", order[2 * n // 3:])):
            e = [energies[i] for i in idx]
            h = [is_hdr[i] for i in idx]
            print(f"    {name}: n={len(idx):3d} mean energy {np.mean(e):.4f} "
                  f"observed P(HDR) {np.mean(h):.4f}")
        print(f"  Pearson(energy, is_hdr) = {corr:+.4f}", flush=True)

        os.makedirs(os.path.dirname(settings.CONTRACT_PATH), exist_ok=True)
        json.dump(contract, open(settings.CONTRACT_PATH, "w"))
        json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
        json.dump(rows, open(settings.MINER_SUBMISSION_PATH, "w"))
        run_stage12(cell_types)
        run_stage3(seed=seed)
        r4 = run_stage4(seed=seed)
        print(f"  stage4: consistency_factor {r4['consistency_factor']:.4f}  "
              f"weighted {r4['total_weighted_score']:.1f}")
        mr = r4.get("model_results") or {}
        for t in ("is_cut", "is_hdr", "indel_length"):
            v = mr.get(t, {})
            print(f"    {t:<13} r2_mean {v.get('r2_mean', float('nan')):+.4f}  "
                  f"mae_mean {v.get('mae_mean', float('nan')):.4f}")
        print(flush=True)


if __name__ == "__main__":
    main()
