"""Replay auto_rank on a frozen task history and save latest-10/latest-20 reports.

Both tables use the current 10/20/30 average-rank rule; 10 and 20 specify the
number of tasks displayed. Production files are read only.
"""
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import sys

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))
import strategy_rank as SR
import seed_model.data as D


def fmt(values):
    return "/".join(map(str, values)) if values is not None else "—"


def summary(rows):
    ready = [r for r in rows if r["matched_seed_count"] is not None]
    return {
        "tasks": len(rows), "scored_tasks": len(ready),
        "matched_seed_count": sum(r["matched_seed_count"] for r in ready),
        "matched_distinct_window_count": sum(r["matched_distinct_window_count"] for r in ready),
        "actual_distinct_windows": sum(len(set(r["actual_windows"])) for r in ready),
        "tasks_with_match": sum(r["matched_seed_count"] > 0 for r in ready),
        "method_counts": dict(Counter(r["picked_method"] for r in rows)),
    }


def markdown_table(rows):
    lines = [
        "| No | Created (UTC) | Cell type | Hist | Picked method | Predicted windows | Actual windows | Matched seeds |",
        "|---:|---|---|---:|---|---|---|---:|",
    ]
    for i, r in enumerate(rows, 1):
        lines.append(f"| {i} | {r['created_at'][:16].replace('T', ' ')} | {r['cell_type']} | "
                     f"{r['history_count']} | {r['picked_method']} | {fmt(r['predicted_windows'])} | "
                     f"{fmt(r['actual_windows'])} | {r['matched_seed_count']} |")
    return lines


def write_csv(path, rows):
    fields = ["No", "Created", "Cell type", "Hist", "Picked method", "Predicted windows",
              "Actual windows", "Matched seeds", "Matched distinct windows", "Average rank", "Task ID"]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for i, r in enumerate(rows, 1):
            writer.writerow([i, r["created_at"], r["cell_type"], r["history_count"], r["picked_method"],
                             fmt(r["predicted_windows"]), fmt(r["actual_windows"]),
                             r["matched_seed_count"], r["matched_distinct_window_count"],
                             r["picked_average_rank"], r["task_id"]])


def main():
    assert D.WIDTH == 100
    snapshot = json.loads((OUT / "task_history.json").read_text())
    rows = snapshot["rows"]
    # Historical model predictions must come from the per-task walk record,
    # never today's next_prediction.json.
    SR.WALK_RECORD = OUT / "walk_record.json"
    model_records = SR._walk_record()
    configs = {"average_rank_10_20_30": SR.SIZES}
    assert SR.SIZES == (10, 20, 30)
    doc = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": snapshot["source"], "fetched_at": snapshot["fetched_at"],
        "since": snapshot["since"], "history_task_count": len(rows),
        "latest_task_created_at": rows[-1]["at"], "window_width": D.WIDTH,
        "protocol": "Historical replay: target excluded from both strategy ranking and prediction; Hist counts prior stamped tasks of that cell since the configured cutoff.",
        "matching": "Matched seeds counts each actual seed covered by any predicted 100-wide window; distinct matches also saved separately.",
        "code_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest()
                        for f in ("strategy_rank.py", "seed_model/evaluate.py", "seed_model/data.py")},
        "configs": {},
    }
    for label, sizes in configs.items():
        results = []
        for target in rows[-20:]:
            before = [r for r in rows if r["at"] < target["at"]]
            cell = target["cell"]
            tbl = SR.rank_table(before, cell, sizes)
            winner = tbl["winner"]
            assert winner is not None, (target, tbl)
            if winner == "model":
                pred = model_records.get(target["at"][:19])
                note = "Historical model prediction from the saved walk-forward record."
            else:
                pred, note = SR.next_windows(before, cell, winner)
            actual = sorted((seed // 100) * 100 for seed in target["seeds"])
            predicted = sorted(100 + 100*int(c) for c in pred) if pred is not None else None
            covered = sum(w in predicted for w in actual) if predicted is not None else None
            unique = len(set(actual) & set(predicted)) if predicted is not None else None
            if predicted is not None:
                assert len(set(predicted)) == 3
                assert covered == SR._seeds_covered(pred, SR._drawn(target["seeds"]) / 3)
            for k, table in tbl["sizes"].items():
                if table:
                    assert k <= tbl["available"]
            results.append({
                "task_id": target["id"], "created_at": target["at"], "cell_type": cell,
                "history_count": sum(r["cell"] == cell for r in before),
                "history_through": before[-1]["at"], "picked_method": winner,
                "picked_average_rank": tbl["avg_rank"][winner],
                "predicted_windows": predicted, "actual_windows": actual,
                "actual_seeds": target["seeds"], "matched_seed_count": covered,
                "matched_distinct_window_count": unique, "prediction_note": note,
                "ranking": tbl,
            })
        results.reverse()  # newest first
        doc["configs"][label] = {
            "ranking_lookbacks": list(sizes), "latest_20": results,
            "summary_latest_10": summary(results[:10]), "summary_latest_20": summary(results),
        }
        for n in (10, 20):
            write_csv(OUT / f"{label}_latest_{n}.csv", results[:n])
        print(label, "latest10", summary(results[:10]), "latest20", summary(results), flush=True)
    (OUT / "analysis.json").write_text(json.dumps(doc, indent=2) + "\n")
    current = doc["configs"]["average_rank_10_20_30"]
    lines = ["# Seed-window method selections", "",
             f"Latest stamped task: **{doc['latest_task_created_at']} UTC**. "
             f"Snapshot fetched {doc['fetched_at']}. [Task feed]({doc['source']}).", "",
             "These are historical replays of the current **auto_rank** selector. For each task, "
             "strategies are ranked using only earlier tasks, separately for that cell, over lookbacks "
             "of **10, 20, and 30 tasks**; the lowest average rank wins. "
             "The six candidates are model, uniform, marginal, hot_hand, cold_hand, and repeat_last. "
             "Ties use the larger lookback's mean, then method name. "
             "Historical model predictions come from the saved walk record. "
             "These tables do not assert that every historical miner used the replayed selection: "
             "live plan timing or configuration can differ.", "",
             "**Hist** is the number of earlier stamped tasks of that cell since "
             f"{doc['since'][:10]}, excluding the target. "
             "**200** denotes the window 200–299; windows are sorted numerically. "
             "**Matched seeds** follows the selector's metric: actual seed draws covered by the predicted "
             "windows, from 0 to 3. Repeated actual windows count once per seed (100/100/100 can score 3). "
             "Distinct-window matches are included in the CSV and JSON outputs.", ""]
    for n in (10, 20):
        s = current[f"summary_latest_{n}"]
        lines += [f"## Latest {n} tasks", "",
                  f"**{s['matched_seed_count']}/{3*s['scored_tasks']} seeds matched**, "
                  f"{s['matched_seed_count']/s['scored_tasks']:.3f} per task; "
                  f"{s['tasks_with_match']}/{s['scored_tasks']} tasks had a match. "
                  f"Distinct windows matched: {s['matched_distinct_window_count']}/{s['actual_distinct_windows']}.", ""]
        lines += markdown_table(current["latest_20"][:n]) + [""]
    lines += ["## Artifacts", "",
              "- [Latest 10 CSV](average_rank_10_20_30_latest_10.csv)",
              "- [Latest 20 CSV](average_rank_10_20_30_latest_20.csv)",
              "- [Full results, task IDs, per-strategy ranks and means](analysis.json)",
              "- [Frozen task history](task_history.json) and [model walk record](walk_record.json)", "",
              "Reproduce from the repository root: "
              "`.venv/bin/python reports/seed_window_methods_20260918/generate.py`.", ""]
    (OUT / "report.md").write_text("\n".join(lines))
    print("\n" + "\n".join(markdown_table(current["latest_20"])))


if __name__ == "__main__":
    main()
