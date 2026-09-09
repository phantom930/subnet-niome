"""Count how often each miner landed in ranks 1..10 over a window of tasks.

Ranking mirrors weight_utils.process_scores_top: within one task, only positive
final_score rows are eligible, sorted descending, and the first TOP_MINER_COUNT
of them get ranks 1..10 (ties broken by miner_uid so the output is stable).

Each task is scored by exactly one validator, so a task's score list is already
a single ranking - no cross-validator aggregation is needed.

Writes fetchRank.json in the repo root, and - when --uid is given - a
per-task detail file (special_score.json) for those uids: task_id, cell_type,
consistency_factor, distribution fidelity and rank for every task they were
scored on in the window.
"""

import argparse
import json
import time
import sys
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from niome_subnet.utils.settings import BASE_URL, MINER_SCORE_URL, TOP_MINER_COUNT

TASKS_URL = f"{BASE_URL}/api/v3/tasks"
OUT_PATH = Path(__file__).resolve().parents[1] / "fetchRank.json"
SPECIAL_OUT_PATH = Path(__file__).resolve().parents[1] / "special_score.json"


def fetch(url, timeout=180, attempts=5, base_delay=3.0):
    """GET with backoff - the API intermittently refuses connections, 504s, or
    truncates a response mid-read, and a whole multi-page crawl should not die
    with it."""
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    delay = base_delay
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except Exception as exc:
            if attempt == attempts:
                raise
            print(f"  fetch failed ({type(exc).__name__}: {exc}); "
                  f"retry {attempt}/{attempts - 1} in {delay:.0f}s", flush=True)
            time.sleep(delay)
            delay = min(delay * 2, 30.0)


def fetch_items(url):
    """Collect `items` across pages; both endpoints expose page/size + a total."""
    sep = "&" if "?" in url else "?"
    items, page = [], 1
    while True:
        payload = fetch(f"{url}{sep}page={page}&size=200")
        batch = payload["items"] if isinstance(payload, dict) else payload
        items.extend(batch)
        meta = payload.get("pagination") or payload if isinstance(payload, dict) else {}
        pages = meta.get("pages") or meta.get("total_pages") or 1
        if page >= pages or not batch:
            break
        page += 1
    return items


def parse_ts(value):
    if value == "now":
        return datetime.now(timezone.utc).replace(tzinfo=None)
    return datetime.fromisoformat(value.replace("Z", "+00:00")).replace(tzinfo=None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default=None,
                    help="ISO date, inclusive. Default: first task in the API.")
    ap.add_argument("--end", default="2026-08-14",
                    help="ISO date, inclusive (whole day), or 'now'. "
                         "Default: 2026-08-14.")
    ap.add_argument("--out", default=str(OUT_PATH))
    ap.add_argument("--uid", action="append", default=[],
                    help="miner uid to dump per-task detail for; repeatable "
                         "or comma-separated (e.g. --uid 39,42)")
    ap.add_argument("--special-out", default=str(SPECIAL_OUT_PATH))
    args = ap.parse_args()

    special_uids = [int(u) for chunk in args.uid for u in chunk.split(",") if u.strip()]

    tasks = fetch_items(TASKS_URL)
    scores = fetch_items(MINER_SCORE_URL)
    print(f"fetched {len(tasks)} tasks, {len(scores)} score rows")

    task_created = {t["id"]: parse_ts(t["created_at"]) for t in tasks}
    task_cell_type = {
        t["id"]: t.get("content", {}).get("contract", {}).get("cell_type")
        for t in tasks
    }

    start = parse_ts(args.start) if args.start else min(task_created.values())
    end = parse_ts(args.end)
    if end.hour == end.minute == end.second == 0:
        end = end.replace(hour=23, minute=59, second=59, microsecond=999999)

    by_task = defaultdict(list)
    for row in scores:
        by_task[row["task_id"]].append(row)

    # rank_counts[uid][rank] -> times, plus hotkeys / scores for context
    rank_counts = defaultdict(lambda: defaultdict(int))
    hotkey_rank_counts = defaultdict(lambda: defaultdict(int))
    hotkeys = defaultdict(set)
    best_scores = defaultdict(float)
    special_records = defaultdict(list)
    ranked_tasks, skipped_no_task_date, orphan_tasks = 0, 0, 0

    for task_id, rows in by_task.items():
        created = task_created.get(task_id)
        if created is None:
            # score rows whose task is no longer served by /tasks
            created = min(parse_ts(r["created_at"]) for r in rows)
            orphan_tasks += 1
        if not (start <= created <= end):
            continue
        eligible = [r for r in rows if r["final_score"] > 0]
        if not eligible:
            skipped_no_task_date += 1
            continue
        eligible.sort(key=lambda r: (-r["final_score"], r["miner_uid"]))
        ranked_tasks += 1
        for rank, row in enumerate(eligible[:TOP_MINER_COUNT], start=1):
            uid = row["miner_uid"]
            rank_counts[uid][rank] += 1
            hotkey_rank_counts[(uid, row["miner_hotkey"])][rank] += 1
            hotkeys[uid].add(row["miner_hotkey"])
            best_scores[uid] = max(best_scores[uid], row["final_score"])

        # per-task detail for the requested uids - ranks past 10 are kept so a
        # near-miss task still shows where the miner actually landed
        for uid in special_uids:
            for position, row in enumerate(eligible, start=1):
                if row["miner_uid"] != uid:
                    continue
                breakdown = row.get("breakdown") or {}
                special_records[uid].append({
                    "task_id": task_id,
                    "created_at": created.isoformat(),
                    "cell_type": task_cell_type.get(task_id),
                    "rank": position,
                    "in_top10": position <= TOP_MINER_COUNT,
                    "consistency_factor": breakdown.get("consistency_factor"),
                    "consistency_score": breakdown.get("consistency_score"),
                    "distribution_fidelity_factor":
                        breakdown.get("distribution_fidelity_factor"),
                    "distribution_fidelity_score":
                        breakdown.get("distribution_fidelity_score"),
                    "final_score": row["final_score"],
                    "miner_hotkey": row["miner_hotkey"],
                    "validator_uid": row["validator_uid"],
                    "ranked_miners_in_task": len(eligible),
                })

    miners = []
    for uid, counts in rank_counts.items():
        per_rank = {str(r): counts.get(r, 0) for r in range(1, TOP_MINER_COUNT + 1)}
        total = sum(counts.values())
        weighted = sum(r * n for r, n in counts.items())
        # a uid can be re-registered under a new hotkey, so keep the split too
        by_hotkey = {}
        for hk in sorted(hotkeys[uid]):
            hk_counts = hotkey_rank_counts[(uid, hk)]
            by_hotkey[hk] = {
                "top10_count": sum(hk_counts.values()),
                "rank_counts": {str(r): hk_counts.get(r, 0)
                                for r in range(1, TOP_MINER_COUNT + 1)},
            }
        miners.append({
            "miner_uid": uid,
            "miner_hotkeys": sorted(hotkeys[uid]),
            "top10_count": total,
            "rank_counts": per_rank,
            "best_rank": min(counts),
            "mean_rank": round(weighted / total, 4),
            "best_final_score": round(best_scores[uid], 6),
            "by_hotkey": by_hotkey,
        })

    # best performers first: most rank-1s, then rank-2s, ... then fewer appearances
    miners.sort(key=lambda m: tuple(-m["rank_counts"][str(r)]
                                    for r in range(1, TOP_MINER_COUNT + 1)))

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sources": {"scores": MINER_SCORE_URL, "tasks": TASKS_URL},
        "window": {"start": start.isoformat(), "end": end.isoformat()},
        "ranking_rule": (
            "per task: positive final_score only, sorted desc "
            f"(ties by miner_uid), first {TOP_MINER_COUNT} get ranks 1..10"
        ),
        "totals": {
            "tasks_fetched": len(tasks),
            "score_rows_fetched": len(scores),
            "tasks_with_scores": len(by_task),
            "tasks_ranked_in_window": ranked_tasks,
            "tasks_in_window_without_positive_scores": skipped_no_task_date,
            "tasks_dated_from_scores_only": orphan_tasks,
            "miners_ranked_at_least_once": len(miners),
        },
        "miners": miners,
    }

    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
    print(f"ranked {ranked_tasks} tasks in {start.date()}..{end.date()}, "
          f"{len(miners)} miners -> {args.out}")

    if not special_uids:
        return

    special_miners = []
    for uid in special_uids:
        records = sorted(special_records[uid], key=lambda r: r["created_at"])
        counts = defaultdict(int)
        for r in records:
            if r["in_top10"]:
                counts[r["rank"]] += 1
        factors = [r["consistency_factor"] for r in records
                   if r["consistency_factor"] is not None]
        fidelities = [r["distribution_fidelity_factor"] for r in records
                      if r["distribution_fidelity_factor"] is not None]
        special_miners.append({
            "miner_uid": uid,
            "miner_hotkeys": sorted({r["miner_hotkey"] for r in records}),
            "tasks_scored": len(records),
            "top10_count": sum(counts.values()),
            "rank_counts": {str(r): counts.get(r, 0)
                            for r in range(1, TOP_MINER_COUNT + 1)},
            "best_rank": min((r["rank"] for r in records), default=None),
            "mean_consistency_factor": (round(sum(factors) / len(factors), 6)
                                        if factors else None),
            "mean_distribution_fidelity_factor":
                (round(sum(fidelities) / len(fidelities), 6) if fidelities else None),
            "cell_types": sorted({r["cell_type"] for r in records
                                  if r["cell_type"]}),
            "records": records,
        })

    special = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sources": {"scores": MINER_SCORE_URL, "tasks": TASKS_URL},
        "window": {"start": start.isoformat(), "end": end.isoformat()},
        "ranking_rule": out["ranking_rule"] + "; rank here is the full position "
                        "among positive scores, so values above 10 are possible",
        "uids": special_uids,
        "miners": special_miners,
    }
    Path(args.special_out).write_text(json.dumps(special, indent=2) + "\n")
    print("special detail: " + ", ".join(
        f"uid {m['miner_uid']} -> {m['tasks_scored']} tasks "
        f"({m['top10_count']} in top10)" for m in special_miners)
        + f" -> {args.special_out}")


if __name__ == "__main__":
    main()
