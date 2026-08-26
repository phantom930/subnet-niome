#!/usr/bin/env python3
"""assemble.py — build a submission from the seed-agnostic banks and test it across the seed lottery.

Composes 250 rows from the strict-900 Cas9 bank (guides that cut under every seed) and a min-union
Cas12a group (guides chosen so their combined no_cut seeds are few), then scores the result across
seeds 100..999 **split into clean vs bad** — clean = a seed no Cas12a group member no_cuts under, so
the whole submission cuts and is_cut is constant. The hypothesis to test: on clean seeds is_cut's
stage-4 r2 recovers to 1.0 (lifting consistency well above the ~0.10 the vanilla build gets), while
on bad seeds a few Cas12a rows no_cut and it collapses as usual.

    python assemble.py --task test/task.json --cell-types test/cell_types.json \\
        --cas9-bank test/guides_cas9_k562_gc/guides.jsonl \\
        --cas12a-group test/guides_cas12a_k562/group_mf22.json --seeds 120
"""

from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import statistics as st
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import niome_subnet.utils.settings as settings  # noqa: E402
import genExp as G  # noqa: E402

STATE: dict = {}


def to_experiment(rec: dict, cell_type: str, exp_id: str) -> dict:
    """A bank record carries a design; stage 1 needs the full experiment row around it."""
    start = rec["target_alignment_start"]
    return {
        "experiment_id": exp_id,
        "guideRNA": rec["guide"],
        "target_alignment_start": start,
        "target_alignment_end": start + rec["length"],
        "strand": rec["strand"],
        "mutation": rec["mutation"],
        "cas_system": rec["cas_system"],
        "cell_type": cell_type,
    }


def cas9_quota(contract, mutations, n_cas9, weight_skew):
    """Rows per (mutation, +/-) Cas9 cell: mutation-weighted like the miner, even across strands."""
    w = contract.get("mutation_weights", {})
    share = {m: max(w.get(m, 1.0), 1e-9) ** weight_skew for m in mutations}
    tot = sum(share.values())
    per_mut = {m: round(n_cas9 * share[m] / tot) for m in mutations}
    quota = {}
    for m in mutations:
        half = per_mut[m] // 2
        quota[(m, "+")] = half
        quota[(m, "-")] = per_mut[m] - half
    return quota


def assemble(cas9_strict, cas12a_group, contract, cas_mix, max_experiments=250):
    """Cas12a rows = the whole min-union group; Cas9 rows fill the rest from the strict bank.

    Ranked by weighted_score within each cell so the strongest rows go in, and capped at the bank's
    availability. Dedup on (cas, start, strand, guide) mirrors stage 1.
    """
    mutations = sorted({r["mutation"] for r in cas9_strict + cas12a_group})
    cell_type = contract.get("cell_type")
    n_cas12a = len(cas12a_group)
    n_cas9 = max_experiments - n_cas12a

    by_cell = {}
    for r in cas9_strict:
        by_cell.setdefault((r["mutation"], r["strand"]), []).append(r)
    for rows in by_cell.values():
        rows.sort(key=lambda r: -r["weighted_score"])

    quota = cas9_quota(contract, mutations, n_cas9, weight_skew=1.25)
    chosen = list(cas12a_group)
    for cell, want in quota.items():
        chosen.extend(by_cell.get(cell, [])[:want])

    # Backfill any shortfall (a cell short of its quota) from the strongest unused Cas9 rows.
    if len(chosen) < max_experiments:
        used = {(r["cas_system"], r["target_alignment_start"], r["strand"], r["guide"])
                for r in chosen}
        spare = sorted((r for r in cas9_strict
                        if (r["cas_system"], r["target_alignment_start"], r["strand"], r["guide"])
                        not in used), key=lambda r: -r["weighted_score"])
        chosen.extend(spare[:max_experiments - len(chosen)])

    rows, seen = [], set()
    for i, rec in enumerate(chosen[:max_experiments]):
        key = (rec["cas_system"], rec["target_alignment_start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(to_experiment(rec, cell_type, f"exp-{i:05d}"))
    return rows


def worker_init(payload):
    logging.getLogger().setLevel(logging.ERROR)
    sys.stdout = open(os.devnull, "w")
    task = json.load(open(payload["task"]))
    STATE["contract"] = task["content"]["contract"]
    STATE["reference"] = task["content"]["hbb_reference"]
    STATE["cell_types"] = payload["cell_types"]
    STATE["subs"] = payload["subs"]        # {"assembled": rows, "vanilla": rows}


def evaluate(job):
    label, seed = job
    sc = copy.deepcopy(STATE["contract"])
    sc["seed"] = seed
    ctx = G.build_context(sc, STATE["reference"], STATE["cell_types"])
    report, valid, results = G.rescore_under(STATE["subs"][label], ctx)
    s4 = report  # score_rows output
    outcomes = Counter(r["outcome"] for r in results)
    # per-target r2 needs the in-memory stage 4; recompute for the is_cut recovery story
    detail = G.stage4_in_memory(valid, results, fold_seed=seed).get("per_target", {})
    return {
        "label": label, "seed": seed,
        "final_score": s4.get("final_score", 0.0),
        "consistency_factor": s4.get("consistency_factor", 0.0),
        "fidelity": s4.get("distribution_fidelity_factor", 0.0),
        "no_cut": outcomes.get("no_cut", 0),
        "r2_is_cut": detail.get("is_cut", {}).get("r2"),
        "r2_is_hdr": detail.get("is_hdr", {}).get("r2"),
        "r2_indel": detail.get("indel_length", {}).get("r2"),
    }


def stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return {}
    return {"min": min(vals), "mean": st.mean(vals), "median": st.median(vals), "max": max(vals)}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--task", default="test/task.json")
    p.add_argument("--cell-types", default="test/cell_types.json")
    p.add_argument("--cas9-bank", default="test/guides_cas9_k562_gc/guides.jsonl")
    p.add_argument("--cas12a-group", default="test/guides_cas12a_k562/group_mf22.json")
    p.add_argument("--cas-mix", default="70/30")
    p.add_argument("--seeds", type=int, default=120, help="total seeds sampled, half clean/half bad")
    p.add_argument("--jobs", type=int, default=12)
    p.add_argument("--out-dir", default="test/assembled")
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.ERROR, format="%(levelname)s %(message)s")
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    task = json.load(open(args.task))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell_types = json.load(open(args.cell_types))

    cas9 = [json.loads(l) for l in open(args.cas9_bank)]
    strict = [r for r in cas9 if r["cut_seeds"] == r["total_seeds"]]
    grp_doc = json.load(open(args.cas12a_group))
    group = grp_doc["group"]
    bad_seeds = set(grp_doc["failed_seed_union"])
    lo, hi = grp_doc["seed_range"]

    print("=" * 92)
    print(f"  task {task.get('id')}  cell_type {contract.get('cell_type')}")
    print(f"  strict-900 Cas9 available: {len(strict)}   Cas12a group: {len(group)} "
          f"(clean {grp_doc['clean_seeds']}/{hi - lo + 1})")
    print("=" * 92)

    assembled = assemble(strict, group, contract, args.cas_mix)
    # Vanilla baseline: the miner's own build at seed 0, the fragile reference point.
    from neurons.miner import Miner
    from dataclasses import replace
    bc = copy.deepcopy(contract); bc["seed"] = 0
    ctx0 = G.build_context(bc, reference, cell_types)
    vcfg = Miner.gen_config_for(bc)
    vsites = G.enumerate_sites(ctx0, vcfg.flank, tuple(sorted(set(vcfg.lengths))))
    vcfg = replace(vcfg, weight_skew=G.choose_weight_skew(ctx0, vsites, vcfg))
    vrows, vvalid, vres = G._generate(ctx0, vsites, vcfg)
    vordered = G.order_rows(vrows, vvalid)
    vanilla = vordered

    cas = Counter((r["cas_system"], r["strand"]) for r in assembled)
    print(f"  assembled {len(assembled)} rows: cas/strand {dict(cas)}")
    with open(out / "submission.json", "w") as h:
        json.dump(assembled, h, indent=2)

    # Sample equal numbers of clean and bad seeds so both regimes are estimated well.
    all_seeds = list(range(lo, hi + 1))
    clean = [s for s in all_seeds if s not in bad_seeds]
    bad = [s for s in all_seeds if s in bad_seeds]
    k = args.seeds // 2
    step_c = max(1, len(clean) // k); step_b = max(1, len(bad) // k)
    sample_clean = clean[::step_c][:k]
    sample_bad = bad[::step_b][:k]

    payload = {"task": str(Path(args.task).resolve()), "cell_types": cell_types,
               "subs": {"assembled": assembled, "vanilla": vanilla}}
    jobs = [(lab, s) for lab in ("assembled", "vanilla")
            for s in sample_clean + sample_bad]

    recs = []
    with Pool(args.jobs, initializer=worker_init, initargs=(payload,)) as pool:
        for r in pool.imap_unordered(evaluate, jobs, chunksize=4):
            r["regime"] = "clean" if r["seed"] not in bad_seeds else "bad"
            recs.append(r)

    print(f"\n  scored {len(sample_clean)} clean + {len(sample_bad)} bad seeds per submission\n")
    hdr = f"  {'submission':<12}{'regime':<7}{'final.mean':>11}{'cons.mean':>11}{'r2_cut':>9}{'r2_hdr':>9}{'r2_indel':>10}{'nocut':>7}"
    print(hdr)
    for lab in ("assembled", "vanilla"):
        for reg in ("clean", "bad"):
            sub = [r for r in recs if r["label"] == lab and r["regime"] == reg]
            if not sub:
                continue
            print(f"  {lab:<12}{reg:<7}"
                  f"{stats([r['final_score'] for r in sub])['mean']:>11.2f}"
                  f"{stats([r['consistency_factor'] for r in sub])['mean']:>11.4f}"
                  f"{stats([r['r2_is_cut'] for r in sub]).get('mean', 0):>9.3f}"
                  f"{stats([r['r2_is_hdr'] for r in sub]).get('mean', 0):>9.3f}"
                  f"{stats([r['r2_indel'] for r in sub]).get('mean', 0):>10.3f}"
                  f"{stats([r['no_cut'] for r in sub])['mean']:>7.2f}")

    # The headline: assembled mean weighted by the true clean fraction vs vanilla overall.
    p_clean = grp_doc["clean_seeds"] / (hi - lo + 1)
    a_clean = stats([r["final_score"] for r in recs
                     if r["label"] == "assembled" and r["regime"] == "clean"])["mean"]
    a_bad = stats([r["final_score"] for r in recs
                   if r["label"] == "assembled" and r["regime"] == "bad"])["mean"]
    a_blend = p_clean * a_clean + (1 - p_clean) * a_bad
    v_mean = stats([r["final_score"] for r in recs if r["label"] == "vanilla"])["mean"]
    print(f"\n  in-range expectation (clean fraction {p_clean:.3f}):")
    print(f"    assembled  {a_blend:8.2f}   (clean {a_clean:.2f} x {p_clean:.3f} + "
          f"bad {a_bad:.2f} x {1 - p_clean:.3f})")
    print(f"    vanilla    {v_mean:8.2f}")
    print(f"    lift       {a_blend - v_mean:+8.2f}  ({a_blend / max(v_mean, 1e-9):.2f}x)")

    with open(out / "records.jsonl", "w") as h:
        for r in sorted(recs, key=lambda r: (r["label"], r["regime"], r["seed"])):
            h.write(json.dumps(r) + "\n")
    print(f"\n  artifacts in {out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
