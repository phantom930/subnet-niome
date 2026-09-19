# Current auto_rank versus the proposed selector

Latest stamped task: **2026-09-18T03:03:14.557809 UTC**. Snapshot fetched 2026-09-18T04:49:10.033652+00:00. [Task feed](https://niome-api.genomes.io/api/v3/tasks).

Both use the target cell's earlier scored tasks, mean seeds covered, lookbacks 10/20/30, fractional ranks for tied means, and equal-weight average ranks. Both break an average-rank tie by the longest used horizon's mean, then by name. Neither selector fits parameters.

| Rule | Current implementation | Proposed specification |
|---|---|---|
| Candidates | Five baselines plus model, if model records cover the scoring history | Five baselines only |
| Incomplete horizon | Drop unless all k scored tasks exist | Use min(k, available), provided at least MIN_HISTORY |
| No qualifying horizon | No winner; window_plan uses rank_freq | Uniform; count fallback |
| Fitted model candidate | Includes saved model walk-forward predictions | No model candidate |

Assumptions: the five candidates are **uniform, marginal, hot_hand, cold_hand, repeat_last**. **MIN_HISTORY=10**, following the existing auto_rank_alt.py comparator; your text leaves its value unspecified. For example, 15 available scored tasks gives current horizons [10] and proposed horizons [10, 15, 15]. Thus the two truncated horizons receive separate votes. At 25 available, the horizons are [10, 20] versus [10, 20, 25].

The main comparison preserves existing baseline predictions. The selector scores only cell C, but marginal's underlying probabilities use pooled earlier labels across cells. A separate replay restricting baseline inputs to C changes **0/20** proposed picks or predicted window sets here. Full ranks for that variant are also saved in analysis.json.

**Hist** counts earlier stamped tasks of the cell since 2026-08-27. The first task per cell initializes history and has no scored prediction; available scored history is Hist − 1. The ranking functions retain at most the newest 30 scored tasks. A window labelled 200 means 200–299. Matched seeds counts every covered draw; duplicate actual windows can contribute multiple hits. Uniform uses the code's stable top-three tie rule, predicting 100/200/300.

| Report period | Current seeds covered | Proposed seeds covered | Changed picks | Partial horizons | Proposed fallbacks |
|---|---:|---:|---:|---:|---:|
| Latest 10 | 14/30 | 14/30 | 0/10 | 0 | 0 |
| Latest 20 | 26/60 | 27/60 | 1/20 | 0 | 0 |

All latest-20 rankings include model and have all three complete horizons. Therefore every pick difference in this period is caused by candidate-pool changes.

## Why the one changed pick differs

Task **3fd0e6e3-5fee-4a8c-9547-edbe5059546c**, 2026-09-17T00:26:18.340389 UTC, **HEK293**.

| Candidate | Current ranks 10/20/30 | Current average | Proposed ranks 10/20/30 | Proposed average | 30-task mean |
|---|---|---:|---|---:|---:|
| model | 6/2.5/2.5 | 3.666667 | excluded | — | 1.133333 |
| uniform | 2/4.5/6 | 4.166667 | 2/3.5/5 | 3.500000 | 1.000000 |
| marginal | 4/4.5/4.5 | 4.333333 | 4/3.5/3.5 | 3.666667 | 1.066667 |
| hot_hand | 5/6/2.5 | 4.500000 | 5/5/2 | 4.000000 | 1.133333 |
| cold_hand | 1/1/4.5 | 2.166667 | 1/1/3.5 | 1.833333 | 1.066667 |
| repeat_last | 3/2.5/1 | 2.166667 | 3/2/1 | 2.000000 | 1.200000 |

Current cold_hand and repeat_last tie at average rank 2.166667; repeat_last wins on its higher 30-task mean (1.2 versus 1.066667). Removing model changes cold_hand's average to 1.833333 and repeat_last's to 2.0, so cold_hand wins. Actual windows are **300/500/600**. Current predicts **200/400/600**, covering **1** seed; proposed predicts **300/600/900**, covering **2** seeds.

## Latest 10 tasks

| No | Created (UTC) | Cell | Hist | Current method | Current windows | Proposed method | Proposed windows | Actual windows | Current matches | Proposed matches |
|---:|---|---|---:|---|---|---|---|---|---:|---:|
| 1 | 2026-09-18 03:03 | CD34+_HSPC | 51 | repeat_last | 100/200/900 | repeat_last | 100/200/900 | 300/700/800 | 0 | 0 |
| 2 | 2026-09-18 00:35 | HEK293 | 52 | cold_hand | 200/300/900 | cold_hand | 200/300/900 | 300/700/800 | 1 | 1 |
| 3 | 2026-09-17 22:08 | CD34+_HSPC | 50 | repeat_last | 300/400/800 | repeat_last | 300/400/800 | 100/900/900 | 0 | 0 |
| 4 | 2026-09-17 19:44 | K562 | 46 | cold_hand | 200/800/900 | cold_hand | 200/800/900 | 200/300/700 | 1 | 1 |
| 5 | 2026-09-17 17:19 | HEK293 | 51 | cold_hand | 200/300/900 | cold_hand | 200/300/900 | 100/800/900 | 1 | 1 |
| 6 | 2026-09-17 14:54 | HEK293 | 50 | cold_hand | 200/300/900 | cold_hand | 200/300/900 | 200/300/800 | 2 | 2 |
| 7 | 2026-09-17 12:29 | K562 | 45 | cold_hand | 100/800/900 | cold_hand | 100/800/900 | 100/100/100 | 3 | 3 |
| 8 | 2026-09-17 10:04 | K562 | 44 | cold_hand | 600/800/900 | cold_hand | 600/800/900 | 500/600/900 | 2 | 2 |
| 9 | 2026-09-17 07:40 | K562 | 43 | cold_hand | 600/800/900 | cold_hand | 600/800/900 | 400/800/900 | 2 | 2 |
| 10 | 2026-09-17 05:15 | HUDEP-2 | 49 | uniform | 100/200/300 | uniform | 100/200/300 | 100/300/400 | 2 | 2 |

## Latest 20 tasks

| No | Created (UTC) | Cell | Hist | Current method | Current windows | Proposed method | Proposed windows | Actual windows | Current matches | Proposed matches |
|---:|---|---|---:|---|---|---|---|---|---:|---:|
| 1 | 2026-09-18 03:03 | CD34+_HSPC | 51 | repeat_last | 100/200/900 | repeat_last | 100/200/900 | 300/700/800 | 0 | 0 |
| 2 | 2026-09-18 00:35 | HEK293 | 52 | cold_hand | 200/300/900 | cold_hand | 200/300/900 | 300/700/800 | 1 | 1 |
| 3 | 2026-09-17 22:08 | CD34+_HSPC | 50 | repeat_last | 300/400/800 | repeat_last | 300/400/800 | 100/900/900 | 0 | 0 |
| 4 | 2026-09-17 19:44 | K562 | 46 | cold_hand | 200/800/900 | cold_hand | 200/800/900 | 200/300/700 | 1 | 1 |
| 5 | 2026-09-17 17:19 | HEK293 | 51 | cold_hand | 200/300/900 | cold_hand | 200/300/900 | 100/800/900 | 1 | 1 |
| 6 | 2026-09-17 14:54 | HEK293 | 50 | cold_hand | 200/300/900 | cold_hand | 200/300/900 | 200/300/800 | 2 | 2 |
| 7 | 2026-09-17 12:29 | K562 | 45 | cold_hand | 100/800/900 | cold_hand | 100/800/900 | 100/100/100 | 3 | 3 |
| 8 | 2026-09-17 10:04 | K562 | 44 | cold_hand | 600/800/900 | cold_hand | 600/800/900 | 500/600/900 | 2 | 2 |
| 9 | 2026-09-17 07:40 | K562 | 43 | cold_hand | 600/800/900 | cold_hand | 600/800/900 | 400/800/900 | 2 | 2 |
| 10 | 2026-09-17 05:15 | HUDEP-2 | 49 | uniform | 100/200/300 | uniform | 100/200/300 | 100/300/400 | 2 | 2 |
| 11 | 2026-09-17 02:51 | HEK293 | 49 | cold_hand | 200/300/900 | cold_hand | 200/300/900 | 500/600/700 | 0 | 0 |
| 12 | 2026-09-17 00:26 | HEK293 | 48 | repeat_last | 200/400/600 | cold_hand | 300/600/900 | 300/500/600 | 1 | 2 |
| 13 | 2026-09-16 22:01 | HUDEP-2 | 48 | uniform | 100/200/300 | uniform | 100/200/300 | 200/600/800 | 1 | 1 |
| 14 | 2026-09-16 19:37 | HUDEP-2 | 47 | uniform | 100/200/300 | uniform | 100/200/300 | 200/700/900 | 1 | 1 |
| 15 | 2026-09-16 17:12 | K562 | 42 | cold_hand | 600/800/900 | cold_hand | 600/800/900 | 400/500/700 | 0 | 0 |
| 16 | 2026-09-16 14:47 | HUDEP-2 | 46 | uniform | 100/200/300 | uniform | 100/200/300 | 100/200/400 | 2 | 2 |
| 17 | 2026-09-16 12:22 | CD34+_HSPC | 49 | repeat_last | 100/200/500 | repeat_last | 100/200/500 | 300/400/800 | 0 | 0 |
| 18 | 2026-09-16 09:58 | CD34+_HSPC | 48 | repeat_last | 100/200/700 | repeat_last | 100/200/700 | 100/100/500 | 2 | 2 |
| 19 | 2026-09-16 07:33 | K562 | 41 | cold_hand | 200/800/900 | cold_hand | 200/800/900 | 200/800/800 | 3 | 3 |
| 20 | 2026-09-16 05:08 | HUDEP-2 | 45 | uniform | 100/200/300 | uniform | 100/200/300 | 300/300/700 | 2 | 2 |

## Early-history behavior

Across all 202 stamped tasks since the cutoff, current has no qualifying horizon for **44** tasks. The proposed rule has **44** uniform fallbacks and uses a partial horizon on **80** tasks. These are initialization counts over the whole record, not the latest-20 totals. Current no-winner cases route to rank_freq; this report does not score those production fallbacks.

These are historical replays ordered by task creation time, using only strictly earlier targets for each choice. They do not establish when historical seed stamps became available or which plan each live miner actually received. Model predictions are read from the frozen walk record; this comparison does not re-audit model training. The one-seed gain on 20 tasks does not establish a lasting advantage.

## Reproduce and inspect

Run `.venv/bin/python reports/rank_method_comparison_20260918/generate.py` from the repository root. The report generator cross-checks proposed picks against auto_rank_alt.alt_select and seed counts against the production metric.

- [Latest 10 CSV](latest_10.csv)
- [Latest 20 CSV](latest_20.csv)
- [Full ranks, means, task IDs and strict-cell replay](analysis.json)
- [Frozen task history](task_history.json)
- [Frozen model predictions](walk_record.json)
