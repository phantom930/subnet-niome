"""Compare production auto_rank with the five-baseline specification, read only.

Run from the repository root with .venv/bin/python followed by this file's path.
Inputs are frozen beside this script. MIN_HISTORY=10 follows auto_rank_alt.py.
"""
import csv
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))
import auto_rank_alt as ALT
import strategy_rank as SR

MIN_HISTORY = 10
ALT.MIN_HISTORY = MIN_HISTORY


def spec_table(before, cell):
    scores, _ = SR.per_task_scores(before, cell, 30)
    tables = {}
    for horizon in ALT.HORIZONS:
        take = min(horizon, len(scores))
        if take < MIN_HISTORY:
            continue
        means = {s: sum(r[s] for r in scores[-take:]) / take for s in ALT.CANDIDATES}
        ranks = ALT._frac_rank(means, ALT.CANDIDATES)
        assert sum(ranks.values()) == 15
        tables[horizon] = {"tasks_used": take, "means": means, "ranks": ranks}
    average = {s: sum(t["ranks"][s] for t in tables.values()) / len(tables)
               for s in ALT.CANDIDATES} if tables else {}
    winner = min(ALT.CANDIDATES, key=lambda s: (
        average[s], -tables[max(tables)]["means"][s], s)) if tables else "uniform"
    fallback = not tables
    partial = any(t["tasks_used"] < k for k, t in tables.items())
    assert (winner, fallback, partial, len(scores)) == ALT.alt_select(before, cell)
    return {"winner": winner, "pool": list(ALT.CANDIDATES), "available": len(scores),
            "sizes": tables, "avg_rank": average, "fallback": fallback, "partial": partial}


def predict(before, target, table, model_record):
    method = table["winner"]
    if method == "model":
        pred = model_record[target["at"][:19]]
    elif table.get("fallback"):
        pred = [0, 1, 2]  # Stable top-three of equal probabilities.
    else:
        pred, _ = SR.next_windows(before, target["cell"], method)
    assert pred is not None and len(set(pred)) == 3
    windows = sorted(100 + 100 * int(c) for c in pred)
    actual = [s // 100 * 100 for s in target["seeds"]]
    matched = sum(w in windows for w in actual)
    assert matched == SR._seeds_covered(pred, SR._drawn(target["seeds"]) / 3)
    return {"method": method, "windows": windows, "matched_seeds": matched, "ranking": table}


def fmt(windows):
    return "/".join(map(str, windows))


def summary(rows):
    return {"tasks": len(rows), "total_seeds": 3 * len(rows),
            "current_matched": sum(r["current"]["matched_seeds"] for r in rows),
            "spec_matched": sum(r["spec"]["matched_seeds"] for r in rows),
            "changed_methods": sum(r["current"]["method"] != r["spec"]["method"] for r in rows),
            "spec_fallbacks": sum(r["spec"]["ranking"]["fallback"] for r in rows),
            "spec_partial_horizons": sum(r["spec"]["ranking"]["partial"] for r in rows),
            "current_methods": dict(Counter(r["current"]["method"] for r in rows)),
            "spec_methods": dict(Counter(r["spec"]["method"] for r in rows))}


def detail_table(rows):
    lines = ["| No | Created (UTC) | Cell | Hist | Current method | Current windows | Proposed method | Proposed windows | Actual windows | Current matches | Proposed matches |",
             "|---:|---|---|---:|---|---|---|---|---|---:|---:|"]
    for i, r in enumerate(rows, 1):
        a, b = r["current"], r["spec"]
        lines.append(f"| {i} | {r['created_at'][:16].replace('T', ' ')} | {r['cell']} | "
                     f"{r['history_count']} | {a['method']} | {fmt(a['windows'])} | {b['method']} | "
                     f"{fmt(b['windows'])} | {fmt(r['actual_windows'])} | {a['matched_seeds']} | {b['matched_seeds']} |")
    return lines


def main():
    snapshot = json.loads((OUT / "task_history.json").read_text())
    rows = snapshot["rows"]
    SR.WALK_RECORD = OUT / "walk_record.json"
    model_record = SR._walk_record()
    results = []
    for target in reversed(rows[-20:]):
        before = [r for r in rows if r["at"] < target["at"]]
        cell = target["cell"]
        samecell = [r for r in before if r["cell"] == cell]
        current = SR.rank_table(before, cell)
        proposed = spec_table(before, cell)
        strict = spec_table(samecell, cell)
        assert current["winner"] is not None  # No production fallback in this comparison period.
        assert len(current["pool"]) == 6
        assert current["used_sizes"] == [10, 20, 30]
        results.append({
            "task_id": target["id"], "created_at": target["at"], "cell": cell,
            "history_count": len(samecell), "scored_history_count": max(0, len(samecell) - 1),
            "actual_seeds": target["seeds"],
            "actual_windows": sorted(s // 100 * 100 for s in target["seeds"]),
            "current": predict(before, target, current, model_record),
            "spec": predict(before, target, proposed, model_record),
            "strict_cell_baselines": predict(samecell, target, strict, model_record),
        })
    # Check fallback/short-history behavior over the complete frozen record.
    warmup = Counter()
    for target in rows:
        before = [r for r in rows if r["at"] < target["at"]]
        a = SR.rank_table(before, target["cell"])
        b = spec_table(before, target["cell"])
        warmup["tasks"] += 1
        warmup["current_no_qualifying_horizon"] += a["winner"] is None
        warmup["spec_uniform_fallbacks"] += b["fallback"]
        warmup["spec_partial_horizon_tasks"] += b["partial"]
    strict_changes = sum((r["spec"]["method"], r["spec"]["windows"]) !=
                         (r["strict_cell_baselines"]["method"], r["strict_cell_baselines"]["windows"])
                         for r in results)
    doc = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": snapshot["source"], "fetched_at": snapshot["fetched_at"],
        "latest_task": rows[-1]["at"], "since": snapshot["since"],
        "history_tasks": len(rows), "min_history": MIN_HISTORY,
        "candidate_assumption": "Five existing baselines: uniform, marginal, hot_hand, cold_hand, repeat_last.",
        "history_rule": "created_at strictly before target; excludes unstamped tasks; one initial task per cell is not a scored prediction.",
        "interpretation": "Selector-only comparison keeps baseline definitions; also checks all baseline inputs restricted to the target cell.",
        "code_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in
                        ("strategy_rank.py", "auto_rank_alt.py", "seed_model/evaluate.py", "seed_model/data.py", "window_plan.py")},
        "summary_latest_10": summary(results[:10]), "summary_latest_20": summary(results),
        "warmup_counts_entire_history": dict(warmup),
        "strict_cell_baselines_changed_picks_or_windows_latest_20": strict_changes,
        "latest_20": results,
    }
    (OUT / "analysis.json").write_text(json.dumps(doc, indent=2) + "\n")
    for n in (10, 20):
        with (OUT / f"latest_{n}.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["No", "Created UTC", "Cell", "Hist", "Current method", "Current windows",
                             "Proposed method", "Proposed windows", "Actual windows", "Current matched seeds",
                             "Proposed matched seeds", "Task ID"])
            for i, r in enumerate(results[:n], 1):
                a, b = r["current"], r["spec"]
                writer.writerow([i, r["created_at"], r["cell"], r["history_count"], a["method"], fmt(a["windows"]),
                                 b["method"], fmt(b["windows"]), fmt(r["actual_windows"]),
                                 a["matched_seeds"], b["matched_seeds"], r["task_id"]])
    lines = ["# Current auto_rank versus the proposed selector", "",
             f"Latest stamped task: **{doc['latest_task']} UTC**. Snapshot fetched {doc['fetched_at']}. "
             f"[Task feed]({doc['source']}).", "",
             "Both use the target cell's earlier scored tasks, mean seeds covered, lookbacks 10/20/30, "
             "fractional ranks for tied means, and equal-weight average ranks. Both break an average-rank "
             "tie by the longest used horizon's mean, then by name. Neither selector fits parameters.", "",
             "| Rule | Current implementation | Proposed specification |",
             "|---|---|---|",
             "| Candidates | Five baselines plus model, if model records cover the scoring history | Five baselines only |",
             "| Incomplete horizon | Drop unless all k scored tasks exist | Use min(k, available), provided at least MIN_HISTORY |",
             "| No qualifying horizon | No winner; window_plan uses rank_freq | Uniform; count fallback |",
             "| Fitted model candidate | Includes saved model walk-forward predictions | No model candidate |", "",
             "Assumptions: the five candidates are **uniform, marginal, hot_hand, cold_hand, repeat_last**. "
             "**MIN_HISTORY=10**, following the existing auto_rank_alt.py comparator; your text leaves its value unspecified. "
             "For example, 15 available scored tasks gives current horizons [10] and proposed horizons [10, 15, 15]. "
             "Thus the two truncated horizons receive separate votes. At 25 available, the horizons are [10, 20] versus [10, 20, 25].", "",
             "The main comparison preserves existing baseline predictions. The selector scores only cell C, "
             "but marginal's underlying probabilities use pooled earlier labels across cells. "
             f"A separate replay restricting baseline inputs to C changes **{strict_changes}/20** proposed picks or predicted window sets here. "
             "Full ranks for that variant are also saved in analysis.json.", "",
             "**Hist** counts earlier stamped tasks of the cell since " + snapshot["since"][:10] + ". "
             "The first task per cell initializes history and has no scored prediction; available scored history is Hist − 1. "
             "The ranking functions retain at most the newest 30 scored tasks. "
             "A window labelled 200 means 200–299. Matched seeds counts every covered draw; duplicate actual windows can contribute multiple hits. "
             "Uniform uses the code's stable top-three tie rule, predicting 100/200/300.", "",
             "| Report period | Current seeds covered | Proposed seeds covered | Changed picks | Partial horizons | Proposed fallbacks |",
             "|---|---:|---:|---:|---:|---:|"]
    for n in (10, 20):
        s = doc[f"summary_latest_{n}"]
        lines.append(f"| Latest {n} | {s['current_matched']}/{s['total_seeds']} | {s['spec_matched']}/{s['total_seeds']} | "
                     f"{s['changed_methods']}/{n} | {s['spec_partial_horizons']} | {s['spec_fallbacks']} |")
    lines += ["", "All latest-20 rankings include model and have all three complete horizons. "
              "Therefore every pick difference in this period is caused by candidate-pool changes.", "",
              "## Why the one changed pick differs", ""]
    changed = [r for r in results if r["current"]["method"] != r["spec"]["method"]]
    for r in changed:
        lines += [f"Task **{r['task_id']}**, {r['created_at']} UTC, **{r['cell']}**.", "",
                  "| Candidate | Current ranks 10/20/30 | Current average | Proposed ranks 10/20/30 | Proposed average | 30-task mean |",
                  "|---|---|---:|---|---:|---:|"]
        a, b = r["current"]["ranking"], r["spec"]["ranking"]
        for name in a["pool"]:
            ar = "/".join(f"{a['sizes'][k][name]['rank']:g}" for k in SR.SIZES)
            br = "/".join(f"{b['sizes'][k]['ranks'][name]:g}" for k in SR.SIZES) if name in b["pool"] else "excluded"
            ba = f"{b['avg_rank'][name]:.6f}" if name in b["pool"] else "—"
            lines.append(f"| {name} | {ar} | {a['avg_rank'][name]:.6f} | {br} | {ba} | {a['sizes'][30][name]['mean_seeds']:.6f} |")
        lines += ["", "Current cold_hand and repeat_last tie at average rank 2.166667; "
                  "repeat_last wins on its higher 30-task mean (1.2 versus 1.066667). "
                  "Removing model changes cold_hand's average to 1.833333 and repeat_last's to 2.0, so cold_hand wins. "
                  "Actual windows are **300/500/600**. Current predicts **200/400/600**, covering **1** seed; "
                  "proposed predicts **300/600/900**, covering **2** seeds.", ""]
    for n in (10, 20):
        lines += [f"## Latest {n} tasks", ""] + detail_table(results[:n]) + [""]
    lines += ["## Early-history behavior", "",
              f"Across all {len(rows)} stamped tasks since the cutoff, current has no qualifying horizon for "
              f"**{warmup['current_no_qualifying_horizon']}** tasks. The proposed rule has "
              f"**{warmup['spec_uniform_fallbacks']}** uniform fallbacks and uses a partial horizon on "
              f"**{warmup['spec_partial_horizon_tasks']}** tasks. These are initialization counts over the whole record, "
              "not the latest-20 totals. Current no-winner cases route to rank_freq; this report does not score those production fallbacks.", "",
              "These are historical replays ordered by task creation time, using only strictly earlier targets for each choice. "
              "They do not establish when historical seed stamps became available or which plan each live miner actually received. "
              "Model predictions are read from the frozen walk record; this comparison does not re-audit model training. "
              "The one-seed gain on 20 tasks does not establish a lasting advantage.", "",
              "## Reproduce and inspect", "",
              "Run `.venv/bin/python reports/rank_method_comparison_20260918/generate.py` from the repository root. "
              "The report generator cross-checks proposed picks against auto_rank_alt.alt_select and seed counts against the production metric.", "",
              "- [Latest 10 CSV](latest_10.csv)", "- [Latest 20 CSV](latest_20.csv)",
              "- [Full ranks, means, task IDs and strict-cell replay](analysis.json)",
              "- [Frozen task history](task_history.json)", "- [Frozen model predictions](walk_record.json)", ""]
    (OUT / "report.md").write_text("\n".join(lines))
    print(json.dumps({k: v for k, v in doc.items() if k.startswith(('summary', 'warmup', 'strict_cell'))}, indent=2))


if __name__ == "__main__":
    main()
