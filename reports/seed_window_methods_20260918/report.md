# Seed-window method selections

Latest stamped task: **2026-09-18T03:03:14.557809 UTC**. Snapshot fetched 2026-09-18T04:40:30.787299+00:00. [Task feed](https://niome-api.genomes.io/api/v3/tasks).

These are historical replays of the current **auto_rank** selector. For each task, strategies are ranked using only earlier tasks, separately for that cell, over lookbacks of **10, 20, and 30 tasks**; the lowest average rank wins. The six candidates are model, uniform, marginal, hot_hand, cold_hand, and repeat_last. Ties use the larger lookback's mean, then method name. Historical model predictions come from the saved walk record. These tables do not assert that every historical miner used the replayed selection: live plan timing or configuration can differ.

**Hist** is the number of earlier stamped tasks of that cell since 2026-08-27, excluding the target. **200** denotes the window 200–299; windows are sorted numerically. **Matched seeds** follows the selector's metric: actual seed draws covered by the predicted windows, from 0 to 3. Repeated actual windows count once per seed (100/100/100 can score 3). Distinct-window matches are included in the CSV and JSON outputs.

## Latest 10 tasks

**14/30 seeds matched**, 1.400 per task; 8/10 tasks had a match. Distinct windows matched: 12/27.

| No | Created (UTC) | Cell type | Hist | Picked method | Predicted windows | Actual windows | Matched seeds |
|---:|---|---|---:|---|---|---|---:|
| 1 | 2026-09-18 03:03 | CD34+_HSPC | 51 | repeat_last | 100/200/900 | 300/700/800 | 0 |
| 2 | 2026-09-18 00:35 | HEK293 | 52 | cold_hand | 200/300/900 | 300/700/800 | 1 |
| 3 | 2026-09-17 22:08 | CD34+_HSPC | 50 | repeat_last | 300/400/800 | 100/900/900 | 0 |
| 4 | 2026-09-17 19:44 | K562 | 46 | cold_hand | 200/800/900 | 200/300/700 | 1 |
| 5 | 2026-09-17 17:19 | HEK293 | 51 | cold_hand | 200/300/900 | 100/800/900 | 1 |
| 6 | 2026-09-17 14:54 | HEK293 | 50 | cold_hand | 200/300/900 | 200/300/800 | 2 |
| 7 | 2026-09-17 12:29 | K562 | 45 | cold_hand | 100/800/900 | 100/100/100 | 3 |
| 8 | 2026-09-17 10:04 | K562 | 44 | cold_hand | 600/800/900 | 500/600/900 | 2 |
| 9 | 2026-09-17 07:40 | K562 | 43 | cold_hand | 600/800/900 | 400/800/900 | 2 |
| 10 | 2026-09-17 05:15 | HUDEP-2 | 49 | uniform | 100/200/300 | 100/300/400 | 2 |

## Latest 20 tasks

**26/60 seeds matched**, 1.300 per task; 15/20 tasks had a match. Distinct windows matched: 21/54.

| No | Created (UTC) | Cell type | Hist | Picked method | Predicted windows | Actual windows | Matched seeds |
|---:|---|---|---:|---|---|---|---:|
| 1 | 2026-09-18 03:03 | CD34+_HSPC | 51 | repeat_last | 100/200/900 | 300/700/800 | 0 |
| 2 | 2026-09-18 00:35 | HEK293 | 52 | cold_hand | 200/300/900 | 300/700/800 | 1 |
| 3 | 2026-09-17 22:08 | CD34+_HSPC | 50 | repeat_last | 300/400/800 | 100/900/900 | 0 |
| 4 | 2026-09-17 19:44 | K562 | 46 | cold_hand | 200/800/900 | 200/300/700 | 1 |
| 5 | 2026-09-17 17:19 | HEK293 | 51 | cold_hand | 200/300/900 | 100/800/900 | 1 |
| 6 | 2026-09-17 14:54 | HEK293 | 50 | cold_hand | 200/300/900 | 200/300/800 | 2 |
| 7 | 2026-09-17 12:29 | K562 | 45 | cold_hand | 100/800/900 | 100/100/100 | 3 |
| 8 | 2026-09-17 10:04 | K562 | 44 | cold_hand | 600/800/900 | 500/600/900 | 2 |
| 9 | 2026-09-17 07:40 | K562 | 43 | cold_hand | 600/800/900 | 400/800/900 | 2 |
| 10 | 2026-09-17 05:15 | HUDEP-2 | 49 | uniform | 100/200/300 | 100/300/400 | 2 |
| 11 | 2026-09-17 02:51 | HEK293 | 49 | cold_hand | 200/300/900 | 500/600/700 | 0 |
| 12 | 2026-09-17 00:26 | HEK293 | 48 | repeat_last | 200/400/600 | 300/500/600 | 1 |
| 13 | 2026-09-16 22:01 | HUDEP-2 | 48 | uniform | 100/200/300 | 200/600/800 | 1 |
| 14 | 2026-09-16 19:37 | HUDEP-2 | 47 | uniform | 100/200/300 | 200/700/900 | 1 |
| 15 | 2026-09-16 17:12 | K562 | 42 | cold_hand | 600/800/900 | 400/500/700 | 0 |
| 16 | 2026-09-16 14:47 | HUDEP-2 | 46 | uniform | 100/200/300 | 100/200/400 | 2 |
| 17 | 2026-09-16 12:22 | CD34+_HSPC | 49 | repeat_last | 100/200/500 | 300/400/800 | 0 |
| 18 | 2026-09-16 09:58 | CD34+_HSPC | 48 | repeat_last | 100/200/700 | 100/100/500 | 2 |
| 19 | 2026-09-16 07:33 | K562 | 41 | cold_hand | 200/800/900 | 200/800/800 | 3 |
| 20 | 2026-09-16 05:08 | HUDEP-2 | 45 | uniform | 100/200/300 | 300/300/700 | 2 |

## Artifacts

- [Latest 10 CSV](average_rank_10_20_30_latest_10.csv)
- [Latest 20 CSV](average_rank_10_20_30_latest_20.csv)
- [Full results, task IDs, per-strategy ranks and means](analysis.json)
- [Frozen task history](task_history.json) and [model walk record](walk_record.json)

Reproduce from the repository root: `.venv/bin/python reports/seed_window_methods_20260918/generate.py`.
