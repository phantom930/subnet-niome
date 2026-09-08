#!/usr/bin/env python3
"""sd_fleet.py — what the 4-hotkey seed-depend fleet would have scored on the seed-0 rounds.

Reads the per-round results sd_variants.py wrote and puts them side by side with each round's
actual field: how far our four siblings sit above (or below) the best miner there, which ranks they
take once they compete with each other, and what share of SCORE_DISTRIBUTION that is against a
single hotkey.

Only rounds where at least one miner reached consistency 1.000 are used. A round whose listing says
``seed: 0`` but where nobody hit 1.000 is ambiguous — CLAUDE.md records that the listing reports 0
for rounds that were in fact stamped — and scoring a seed-0 build against it would invent a result.

    python sd_fleet.py [task_id ...]
"""
import glob
import json
import sys

DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
UNSTAMPED_RATE = None          # filled from the listing if it is on disk


def load(task_ids):
    out = []
    paths = ([f"sd_variants_{t[:8]}.json" for t in task_ids] if task_ids
             else sorted(glob.glob("sd_variants_*.json")))
    for p in paths:
        try:
            out.append(json.load(open(p)))
        except FileNotFoundError:
            print(f"  (missing {p} — run sd_variants.py for that task first)", file=sys.stderr)
    return out


def field_finals(task_id):
    """That round's scored finals, from the cached score feed, excluding our own hotkeys.

    From 2026-09-07 18:46 our own submissions are in this feed — h0 placed rank 1 on the 19:06 and
    21:29 rounds — so a replay compared against an unfiltered field would be racing itself.
    """
    from sd_task import OURS
    rows = json.load(open("sd_task_scores.json"))
    rows = rows if isinstance(rows, list) else (rows.get("data") or rows.get("items") or [])
    out = []
    for r in rows:
        if r.get("task_id") != task_id or r.get("miner_hotkey") in OURS:
            continue
        b = r.get("breakdown") or {}
        out.append((float(r.get("final_score") or 0),
                    float(b.get("consistency_factor") or 0), r.get("miner_hotkey")))
    # Deduplicate: the feed carries one row per validator per miner.
    best = {}
    for row in out:
        if row[2] not in best or row[0] > best[row[2]][0]:
            best[row[2]] = row
    return sorted(((f, c) for f, c, _h in best.values()), reverse=True)


def opened(task_id):
    listing = json.load(open("sd_task_listing.json"))
    listing = listing if isinstance(listing, list) else (listing.get("items") or [])
    for t in listing:
        if t["id"] == task_id:
            return (t.get("created_at") or "")[:16]
    return "?"


def main():
    results = load(sys.argv[1:])
    if not results:
        raise SystemExit("no sd_variants_*.json found")

    print(f"{'round opened':<17} {'cell':<11} {'tied':>4} {'field best':>10} "
          f"{'our best':>9} {'margin':>8} {'ranks':>9} {'4hk':>6} {'1hk':>6}")
    tot4 = tot1 = 0.0
    rows = []
    for res in sorted(results, key=lambda r: opened(r["task"]), reverse=True):
        tid = res["task"]
        finals = sorted((v["final"] for v in res["variants"].values()), reverse=True)
        fld = field_finals(tid)
        if not fld:
            continue
        best = fld[0][0]
        tied = sum(1 for _f, c in fld if c >= 0.999)
        combined = sorted([(f, True) for f in finals] + [(f, False) for f, _c in fld],
                          key=lambda x: (-x[0], not x[1]))
        ours_ranks = [i for i, (_f, mine) in enumerate(combined, 1) if mine]
        share4 = sum(DIST[i - 1] for i in ours_ranks if i <= len(DIST))
        share1 = DIST[ours_ranks[0] - 1] if ours_ranks[0] <= len(DIST) else 0.0
        tot4 += share4
        tot1 += share1
        rows.append((tid, finals, best, tied, ours_ranks, share4, share1))
        rng = f"{ours_ranks[0]}-{ours_ranks[-1]}"
        print(f"{opened(tid):<17} {res['cell']:<11} {tied:>4} {best:>10.2f} "
              f"{finals[0]:>9.2f} {finals[0] - best:>+8.2f} {rng:>9} "
              f"{share4:>6.1%} {share1:>6.1%}")

    n = len(rows)
    print(f"\n{'mean over ' + str(n) + ' rounds':<17} {'':<11} {'':>4} {'':>10} {'':>9} {'':>8} "
          f"{'':>9} {tot4 / n:>6.1%} {tot1 / n:>6.1%}")
    print(f"\n  four hotkeys are worth {tot4 / tot1:.2f}x one on these rounds")
    for tid, finals, best, tied, ranks, s4, s1 in rows:
        spread = max(finals) - min(finals)
        print(f"\n  {opened(tid)}  {tid[:8]}  {tied} miners at cons 1.000")
        print(f"    ours: " + ", ".join(f"{f:.2f}" for f in finals)
              + f"   (spread {spread:.2f})")
        print(f"    field top 5: " + ", ".join(f"{f:.2f}" for f, _c in field_finals(tid)[:5]))
    json.dump({"rounds": [{"task": t, "ours": f, "field_best": b, "n_tied": c,
                           "ranks": r, "share_4hk": s4, "share_1hk": s1}
                          for t, f, b, c, r, s4, s1 in rows],
               "mean_share_4hk": tot4 / n, "mean_share_1hk": tot1 / n},
              open("sd_fleet.json", "w"), indent=1)
    print("\nwrote sd_fleet.json")


if __name__ == "__main__":
    main()
