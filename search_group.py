#!/usr/bin/env python3
"""search_group.py — pick the group whose combined failed-seed set is smallest.

Each Cas12a candidate no_cuts under some seeds (at most --max-fail of them). For a *group* of guides,
the "failed-group-seed count" is the union of those seeds — each seed counted once, however many
members fail it. The 900 - union seeds where no member fails are exactly the seeds under which the
whole group cuts cleanly, so minimising the union maximises the seeds for which the entire group is
seed-agnostic. This selects --group-size of them minimising that union.

Min-union of k sets chosen from n is NP-hard, but n=150 with tiny sets is easy to get near-optimal:
removal-greedy (drop the member covering the most seeds nothing else covers) then swap local search,
from several randomised orders, keeping the best.

    python search_group.py --bank test/guides_cas12a_k562/guides.jsonl --max-fail 16 --group-size 75
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path


def cut_draw(seed, mutation, cas, guide, start, strand):
    key = f"{seed}|{mutation}|{cas}|{guide}|{start}|{strand}"
    rng = random.Random(int.from_bytes(hashlib.sha256(key.encode()).digest()[-4:], "big"))
    rng.random()
    return rng.random()


def failed_seeds(cand, seeds):
    """The exact seeds this candidate no_cuts under — recomputed, not stored in the bank."""
    cp = cand["cut_p"]
    return frozenset(s for s in seeds
                     if cut_draw(s, cand["mutation"], cand["cas_system"], cand["guide"],
                                 cand["target_alignment_start"], cand["strand"]) > cp)


def union_size(indices, fails):
    u = set()
    for i in indices:
        u |= fails[i]
    return len(u)


def removal_greedy(fails, keep, order, cell=None, floor=None):
    """Start with all, drop the member removing the most union each step until `keep` remain.

    Removing a member shrinks the union only by the seeds it *uniquely* covers, so a per-seed
    coverage count makes each step O(sets). `order` breaks ties, giving different runs to compare.
    With `cell`/`floor`, a member is only removable while its cell is above its floor, so the
    per-cell minimums the submission needs for stage-5 coverage are never violated.
    """
    alive = set(range(len(fails)))
    cover = Counter()
    for i in alive:
        cover.update(fails[i])
    cell_count = Counter(cell[i] for i in alive) if cell else None
    while len(alive) > keep:
        best, best_gain = None, -1
        for i in order:
            if i not in alive:
                continue
            if cell and cell_count[cell[i]] <= floor[cell[i]]:
                continue      # dropping i would starve its cell below the floor
            gain = sum(1 for s in fails[i] if cover[s] == 1)
            if gain > best_gain:
                best, best_gain = i, gain
        if best is None:      # every remaining cell is at its floor; cannot shrink further
            break
        alive.discard(best)
        if cell:
            cell_count[cell[best]] -= 1
        for s in fails[best]:
            cover[s] -= 1
    return alive


def swap_improve(chosen, fails, pool, cell=None, floor=None):
    """Local search: swap an in-member for an out-member whenever it shrinks the union.

    A swap is allowed only if it respects the per-cell floors: dropping i either leaves its cell
    above the floor, or is offset by adding a j from the same cell.
    """
    chosen = set(chosen)
    cell_count = Counter(cell[i] for i in chosen) if cell else None
    improved = True
    while improved:
        improved = False
        base = union_size(chosen, fails)
        for i in list(chosen):
            for j in pool - chosen:
                if cell and cell[i] != cell[j] and cell_count[cell[i]] <= floor[cell[i]]:
                    continue
                trial = (chosen - {i}) | {j}
                if union_size(trial, fails) < base:
                    chosen = trial
                    base = union_size(chosen, fails)
                    if cell:
                        cell_count[cell[i]] -= 1
                        cell_count[cell[j]] += 1
                    improved = True
                    break
            if improved:
                break
    return chosen


def forward_greedy(fails, keep, cell, floor, rng, jitter):
    """Add the member that grows the union least, until `keep`, respecting per-cell floors.

    O(keep * pool) per build via a seed->coverage-count map, so it scales to the thousands-strong
    pools the removal method chokes on. Floor feasibility: a member from a cell whose floor is
    already met is addable only while there are free slots beyond the still-unmet floors, so the
    minimums are always reachable.
    """
    n = len(fails)
    chosen: set[int] = set()
    cover = Counter()
    union = 0
    cell_count = Counter()
    floors = floor or {}
    total_need = sum(floors.values())
    while len(chosen) < keep:
        slots_left = keep - len(chosen)
        unmet = sum(max(0, floors.get(c, 0) - cell_count[c]) for c in floors) if floors else 0
        free_slot = slots_left > unmet
        best_j, best_cost = None, 10 ** 9
        for j in range(n):
            if j in chosen:
                continue
            c = cell[j] if cell else None
            if floors and not (cell_count[c] < floors.get(c, 0) or free_slot):
                continue
            cost = 0
            for sd in fails[j]:
                if cover[sd] == 0:
                    cost += 1
            if cost < best_cost or (jitter and cost == best_cost and rng.random() < 0.5):
                best_cost, best_j = cost, j
        chosen.add(best_j)
        if cell:
            cell_count[cell[best_j]] += 1
        for sd in fails[best_j]:
            if cover[sd] == 0:
                union += 1
            cover[sd] += 1
    return chosen, cover, union


def swap_cover(chosen, cover, union, fails, cell, floor):
    """Swap-improve using the coverage map, so each trial is O(set) not O(all sets)."""
    chosen = set(chosen)
    cover = Counter(cover)
    cell_count = Counter(cell[i] for i in chosen) if cell else None
    floors = floor or {}
    improved = True
    while improved:
        improved = False
        for i in list(chosen):
            ci = cell[i] if cell else None
            saved = sum(1 for sd in fails[i] if cover[sd] == 1)   # seeds i alone covers
            for j in range(len(fails)):
                if j in chosen:
                    continue
                cj = cell[j] if cell else None
                if cell and ci != cj and cell_count[ci] <= floors.get(ci, 0):
                    continue
                # new seeds j adds after i is gone: cover minus i's contribution hits 0
                added = sum(1 for sd in fails[j]
                            if cover[sd] - (1 if sd in fails[i] else 0) == 0)
                if added - saved < 0:
                    for sd in fails[i]:
                        cover[sd] -= 1
                    for sd in fails[j]:
                        cover[sd] += 1
                    union += added - saved
                    chosen.discard(i)
                    chosen.add(j)
                    if cell:
                        cell_count[ci] -= 1
                        cell_count[cj] += 1
                    improved = True
                    break
            if improved:
                break
    return chosen, union


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bank", default="test/guides_cas12a_k562/guides.jsonl")
    p.add_argument("--max-fail", type=int, default=16,
                   help="only candidates that no_cut under at most this many seeds")
    p.add_argument("--group-size", type=int, default=75)
    p.add_argument("--start-seed", type=int, default=100)
    p.add_argument("--end-seed", type=int, default=999)
    p.add_argument("--restarts", type=int, default=40)
    p.add_argument("--per-cell-min", type=int, default=0,
                   help="minimum members from each (mutation, cas, strand) cell; capped at what "
                        "the pool has, so a thin cell keeps all of its candidates")
    p.add_argument("--out", default="test/guides_cas12a_k562/group.json")
    return p.parse_args()


def main():
    args = parse_args()
    seeds = list(range(args.start_seed, args.end_seed + 1))
    n_seeds = len(seeds)

    bank = [json.loads(l) for l in open(args.bank)]
    pool = [c for c in bank if c["cut_seeds"] >= n_seeds - args.max_fail]
    print(f"pool: {len(pool)} candidates with <= {args.max_fail} no_cut, choosing "
          f"{args.group_size}; seeds {args.start_seed}..{args.end_seed}")
    if len(pool) < args.group_size:
        print(f"  ! only {len(pool)} candidates; need {args.group_size}")
        return 1

    fails = [failed_seeds(c, seeds) for c in pool]
    all_idx = set(range(len(pool)))
    print(f"  union of all {len(pool)}: {union_size(all_idx, fails)} seeds; "
          f"mean fails/candidate {sum(len(f) for f in fails) / len(fails):.1f}")

    cell = [(c["mutation"], c["cas_system"], c["strand"]) for c in pool]
    available = Counter(cell)
    floor = {k: min(args.per_cell_min, available[k]) for k in available}
    if args.per_cell_min:
        capped = {k: v for k, v in floor.items() if v < args.per_cell_min}
        print(f"  per-cell floor {args.per_cell_min}; available {dict(available)}")
        if capped:
            print(f"  ! floor capped by the pool for {capped} (fewer than {args.per_cell_min})")
        if sum(floor.values()) > args.group_size:
            print(f"  ! floors sum to {sum(floor.values())} > group size {args.group_size}")
            return 1
    use_floor = args.per_cell_min > 0

    # Forward-greedy (cover-count, O(keep*pool)) scales to the large pools; removal-greedy is only
    # affordable when the pool is small. Run what fits, refine each with a cover-count swap, keep
    # the smallest union across randomised restarts.
    cf = cell if use_floor else None
    fl = floor if use_floor else None
    best, best_u = None, 10 ** 9
    rng = random.Random(0xC0FFEE)
    small = len(pool) <= 400
    for r in range(args.restarts):
        chosen, cover, union = forward_greedy(fails, args.group_size, cf, fl, rng, jitter=(r > 0))
        chosen, union = swap_cover(chosen, cover, union, fails, cf, fl)
        if union < best_u:
            best, best_u = set(chosen), union
        if small:
            order = list(range(len(pool)))
            if r:
                rng.shuffle(order)
            rc = removal_greedy(fails, args.group_size, order, cf, fl)
            rc = swap_improve(rc, fails, all_idx, cf, fl)
            u = union_size(rc, fails)
            if u < best_u:
                best, best_u = set(rc), u
    chosen = best

    # Baseline: the individually-best group-size guides (fewest no_cut each), for comparison.
    top = set(sorted(range(len(pool)), key=lambda i: len(fails[i]))[:args.group_size])
    clean = n_seeds - best_u
    print(f"\n  min-union group: union {best_u} -> {clean} clean seeds "
          f"({100 * clean / n_seeds:.1f}% of the range fully cut by the whole group)")
    print(f"  baseline (individually fewest-fail 75): union {union_size(top, fails)} "
          f"-> {n_seeds - union_size(top, fails)} clean")

    # What actually hits the score under a bad seed: how many rows no_cut at once.
    per_seed = Counter()
    for i in chosen:
        per_seed.update(fails[i])
    hist = Counter(per_seed.values())
    worst = max(per_seed.values()) if per_seed else 0
    exp_nocut = sum(per_seed.values()) / n_seeds
    print(f"  under a random in-range seed: E[no_cut rows] {exp_nocut:.2f}, "
          f"worst-case simultaneous {worst}")
    print(f"  bad-seed fail multiplicity: {dict(sorted(hist.items()))}")

    group = [pool[i] for i in sorted(chosen)]
    cells = Counter((c["mutation"][:20], c["strand"]) for c in group)
    print(f"  per (mutation, strand) in the chosen 75: {dict(cells)}")

    doc = {
        "bank": args.bank, "group_size": args.group_size, "max_fail": args.max_fail,
        "seed_range": [args.start_seed, args.end_seed],
        "union_failed_seeds": best_u, "clean_seeds": clean,
        "baseline_top_union": union_size(top, fails),
        "expected_no_cut_per_seed": exp_nocut, "worst_case_simultaneous_no_cut": worst,
        "per_cell": {f"{m}|{s}": n for (m, s), n in cells.items()},
        "failed_seed_union": sorted(set().union(*[fails[i] for i in chosen])),
        "group": group,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as h:
        json.dump(doc, h, indent=2)
    print(f"\n  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
