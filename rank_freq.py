#!/usr/bin/env python3
"""rank_freq.py — backtest the rank-frequency seed-window predictor, per cell type.

The scheme, exactly as specified:

  1. seed_window_cumulative_counts[r][w] — how often window w (of nine width-100 windows) was
     drawn, over rounds < r.
  2. seed_value_rank[r][0:3] — the rank of each of round r's three seed windows, on the counts
     from rounds < r.
  3. Row 0 is 0/0/0: no prior data to rank.
  4. The nine windows get nine DISTINCT ranks. Ties on count break by most-recently-updated first,
     then by lower window index.
  5. rank_cumulative_counts[r][k] — how often rank k (0..9, so ten columns) appeared in
     seed_value_rank over rounds < r.
  6. The three highest-frequency ranks. Ties break by most-recently-updated position, then by
     lower rank number.
  7. Those three ranks map to the three windows currently holding them -> a joined 300-seed window.
  8. Eleven hotkeys tile it at stride 10, width 200 (offsets 0,10,...,100 in the joined space).

Everything a round's prediction uses comes from strictly earlier rounds, so this is a real
out-of-sample backtest and not a fit. The number that matters is how many of a round's three seeds
land in the predicted 300-seed window against the 1.00 a uniform draw gives (300/900 x 3).

    python rank_freq.py [--from 2026-08-27T00:00:00]
"""
import json
import os
import sys
from collections import defaultdict

FROM = "2026-08-27T00:00:00"
for i, a in enumerate(sys.argv):
    if a == "--from" and i + 1 < len(sys.argv):
        FROM = sys.argv[i + 1]

N_WIN = 9
STRIDE, HK_WIDTH, N_HK = 10, 200, 11


def win_of(seed):
    """Nine width-100 windows over 100-999: window 0 is 100-199, window 8 is 900-999."""
    return (int(seed) - 100) // 100


def rank_windows(counts, last_upd):
    """Nine distinct ranks 1..9. Highest count is rank 1; ties -> latest updated, then lower index."""
    order = sorted(range(N_WIN), key=lambda w: (-counts[w], -last_upd[w], w))
    rank_of = {w: i + 1 for i, w in enumerate(order)}
    return rank_of, order


def top3_ranks(rank_counts, rank_last_upd):
    """The three most frequent ranks; ties -> latest updated position, then lower rank number."""
    cand = [k for k in range(1, N_WIN + 1)]
    return sorted(cand, key=lambda k: (-rank_counts[k], -rank_last_upd[k], k))[:3]


def joined_windows(top3, rank_of):
    """Map ranks back to the windows holding them, in seed order."""
    holder = {r: w for w, r in rank_of.items()}
    return sorted(holder[r] for r in top3 if r in holder)


def hotkey_windows(joined):
    """Eleven width-200 slices at stride 10 over the concatenated 300 seeds of `joined`."""
    space = [w * 100 + 100 + i for w in joined for i in range(100)]   # actual seed values
    out = []
    for h in range(N_HK):
        off = h * STRIDE
        out.append(set(space[off:off + HK_WIDTH]))
    return out


def main():
    listing = json.load(open("sd_task_listing.json"))
    items = listing if isinstance(listing, list) else (listing.get("items") or [])
    rounds = defaultdict(list)
    for t in items:
        c = (t.get("content") or {}).get("contract") or {}
        seeds = [int(x) for x in str(c.get("seed", "")).split(",") if x.strip().isdigit()]
        at = t.get("created_at", "")
        if len(seeds) == 3 and at >= FROM and all(100 <= s <= 999 for s in seeds):
            rounds[c.get("cell_type")].append((at, t["id"], seeds))
    print(f"rank-frequency backtest, three-seed rounds from {FROM}\n")
    grand = {"n": 0, "hits": 0, "any": 0, "hk": 0, "hk_n": 0}
    per_cell = {}
    for cell in sorted(rounds):
        rs = sorted(rounds[cell])
        counts = [0] * N_WIN
        last_upd = [-1] * N_WIN
        rank_counts = defaultdict(int)
        rank_last_upd = defaultdict(lambda: -1)
        svr = []
        hits_seq, any_seq, hk_seq = [], [], []
        print(f"  {cell}  ({len(rs)} rounds)")
        print(f"    {'#':>3} {'task':<9} {'seeds':<18} {'top3 ranks':<12} {'joined windows':<22} "
              f"{'hit':>3} {'hk cover':>9}")
        for r, (at, tid, seeds) in enumerate(rs):
            rank_of, _order = rank_windows(counts, last_upd)
            if r == 0:
                pred_top3, joined = [], []
            else:
                pred_top3 = top3_ranks(rank_counts, rank_last_upd)
                joined = joined_windows(pred_top3, rank_of)
            jset = {w * 100 + 100 + i for w in joined for i in range(100)}
            hit = sum(1 for s in seeds if s in jset)
            # how many of the 11 hotkeys hold at least one of the round's seeds
            cover = 0
            if joined:
                cover = sum(1 for hw in hotkey_windows(joined)
                            if any(s in hw for s in seeds))
            if r > 0:
                hits_seq.append(hit); any_seq.append(1 if hit else 0); hk_seq.append(cover)
            print(f"    {r:>3} {tid[:8]:<9} {str(seeds):<18} "
                  f"{str(pred_top3) if pred_top3 else '-':<12} "
                  f"{','.join(f'{w*100+100}-{w*100+199}' for w in joined) if joined else '-':<22} "
                  f"{hit:>3} {cover:>6}/11")
            # ---- now fold round r in ----
            ranks_this = [0, 0, 0] if r == 0 else [rank_of[win_of(s)] for s in seeds]
            svr.append(ranks_this)
            for k in ranks_this:
                rank_counts[k] += 1
                rank_last_upd[k] = r
            for s in seeds:
                w = win_of(s)
                counts[w] += 1
                last_upd[w] = r
        n = len(hits_seq)
        if n:
            mh = sum(hits_seq) / n
            pa = sum(any_seq) / n
            mc = sum(hk_seq) / n
            per_cell[cell] = (n, mh, pa, mc)
            grand["n"] += n; grand["hits"] += sum(hits_seq)
            grand["any"] += sum(any_seq); grand["hk"] += sum(hk_seq)
            print(f"    -> {n} predicted rounds: mean seeds in joined window {mh:.3f} "
                  f"(chance 1.000)  P(>=1) {pa:.1%} (chance 70.4%)  mean hotkeys covering "
                  f"{mc:.2f}/11\n")
    if grand["n"]:
        n = grand["n"]
        print(f"  ALL CELLS: {n} predicted rounds")
        print(f"    mean seeds in the joined 300-seed window : {grand['hits']/n:.3f}  "
              f"(uniform chance 1.000)")
        print(f"    P(at least one seed in it)               : {grand['any']/n:.1%}  "
              f"(uniform chance 70.4%)")
        print(f"    mean hotkeys (of 11) holding a seed      : {grand['hk']/n:.2f}")
        lift = (grand['hits']/n)
        print(f"\n    lift over chance: {lift:.3f}x — a value near 1.00 means the predictor is")
        print(f"    doing nothing beyond naming any 300 seeds.")
    json.dump({"from": FROM, "per_cell": {k: {"rounds": v[0], "mean_hits": v[1],
               "p_any": v[2], "mean_hk": v[3]} for k, v in per_cell.items()}},
              open("rank_freq.json", "w"), indent=1)
    print("\nwrote rank_freq.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
