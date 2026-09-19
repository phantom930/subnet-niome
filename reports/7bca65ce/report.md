# Submission audit: 7bca65ce-693d-4ad9-bf60-6297014bed83

HEK293 · task opened 2026-09-17T17:19:11.656625 UTC · scoring seeds **117, 842, 989**.

**Best of h0–h5: h0, 24.337843, rank 176/248. Top miner: UID 83, 157.829909.** h0 reached 15.42% of the top score, a gap of 133.492066 points. All six submissions were uploaded successfully, contained 250 valid rows, and received zero reward weight.

Scores are from the [task-specific public API](https://niome-api.genomes.io/api/v3/miners/scores?task_id=7bca65ce-693d-4ad9-bf60-6297014bed83&limit=40000); 248 records, 248 distinct miners, one validator. Score timestamp: 2026-09-17T19:43:54.537196 UTC. Audit generated: 2026-09-18T03:58:34.976957+00:00.

## Official score comparison

| Miner | UID this round | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 83 | 1 | 157.829909 | 292.737406 | 0.561675 | 0.959900 | 0.294000 |
| h0 | 248 | 176 | 24.337843 | 306.420997 | 0.086152 | 0.921927 | 0.000000 |
| h1 | 41 | 197 | 23.418563 | 308.007712 | 0.082579 | 0.920725 | 0.000000 |
| h2 | 163 | 227 | 20.022675 | 311.645380 | 0.070241 | 0.914686 | 0.000000 |
| h3 | 118 | 210 | 22.247743 | 307.017767 | 0.078707 | 0.920681 | 0.000000 |
| h4 | 190 | 222 | 20.594702 | 309.883338 | 0.072664 | 0.914614 | 0.000000 |
| h5 | 234 | 215 | 21.571634 | 307.121364 | 0.076726 | 0.915442 | 0.000000 |

Top miner hotkey: `5H8x24ZVrSWQenChH96fRc52CVvw2YP5zC1vpvKLLVkCa5Q1`. The rank-10 cutoff was **119.453361**.

Final score is the mean of the per-seed scores. For these submissions the weighted score and fidelity do not change across seeds, so final = weighted score × mean consistency × fidelity.

## Joined window, clean sets, and bands

The original 17:19 builds used joined band space **200–399 ∪ 900–999** (300 seeds), with a separate **100–999 cut search** (900 seeds). Band candidate slices were 150 seeds wide, at offsets 0/50/100/150/200/250 for h0–h5. All six shipped the conjunction build with 80 Cas12a and 170 Cas9 rows, pinning eight HDR seeds.

**Clean** means all 250 submitted rows cut at that seed. **Band** means all 250 return HDR; the HDR band is a subset of the clean set. Both sets below were recovered by replaying the exact uploaded rows on every seed 100–999. The full-row clean sets equal the Cas12a group clean sets and their sizes match the original build logs.

| Hotkey | Candidate slice within joined space | Clean / 900 | HDR band seeds | Hits among 117, 842, 989 |
|---|---|---:|---|---|
| h0 | 200–349 | 19 | 207, 221, 263, 273, 279, 289, 299, 336 | 0 clean; 0 band |
| h1 | 250–399 | 16 | 274, 293, 321, 335, 349, 353, 357, 361 | 0 clean; 0 band |
| h2 | 300–399, 900–949 | 18 | 336, 343, 357, 360, 901, 926, 929, 949 | 0 clean; 0 band |
| h3 | 350–399, 900–999 | 18 | 357, 382, 925, 928, 929, 938, 956, 974 | 0 clean; 0 band |
| h4 | 200–249, 900–999 | 17 | 207, 238, 240, 243, 903, 922, 935, 993 | 0 clean; 0 band |
| h5 | 200–299, 950–999 | 19 | 207, 208, 220, 234, 263, 954, 986, 987 | 0 clean; 0 band |

Candidate slices are reconstructed from the original logged joined space and the six-hotkey offset layout; all recovered bands lie inside their corresponding slices.

The fleet covers **99 distinct clean seeds** and **41 distinct HDR band seeds** (48 band memberships before overlap). **None of 117, 842, or 989 lies in either union.** Seeds 117 and 842 are outside the joined band space. Seed 989 is inside the candidate slices of h3/h4/h5 but was not selected into any band and is not clean for any hotkey.

Exact clean sets (including band seeds):

- **h0 (19):** 159, 207, 221, 253, 259, 263, 273, 279, 289, 299, 336, 430, 473, 512, 637, 766, 854, 869, 894.
- **h1 (16):** 139, 168, 196, 274, 293, 297, 321, 335, 349, 353, 357, 361, 386, 587, 789, 906.
- **h2 (18):** 209, 336, 343, 357, 360, 404, 445, 611, 639, 666, 886, 901, 926, 929, 941, 949, 976, 994.
- **h3 (18):** 232, 322, 331, 357, 382, 446, 519, 656, 678, 753, 863, 925, 928, 929, 938, 943, 956, 974.
- **h4 (17):** 190, 197, 207, 238, 240, 243, 479, 484, 658, 663, 682, 695, 903, 922, 935, 992, 993.
- **h5 (19):** 181, 207, 208, 220, 223, 234, 263, 374, 492, 583, 704, 747, 749, 782, 903, 936, 954, 986, 987.

## What happened at the three real seeds

Each cell below shows **per-seed final score / no-cut rows**. Every hotkey had at least one no-cut row on every drawn seed, losing the all-cut consistency benefit.

| Hotkey | Seed 117 | Seed 842 | Seed 989 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 22.4826 / 19 | 27.3916 / 10 | 23.1393 / 11 | 24.337843 |
| h1 | 25.6767 / 14 | 25.2531 / 13 | 19.3259 / 22 | 23.418563 |
| h2 | 22.8654 / 14 | 17.7354 / 19 | 19.4672 / 22 | 20.022675 |
| h3 | 21.1196 / 11 | 22.5595 / 17 | 23.0641 / 14 | 22.247743 |
| h4 | 18.7829 / 18 | 21.6042 / 11 | 21.3970 / 19 | 20.594702 |
| h5 | 21.6297 / 15 | 22.1442 / 14 | 20.9410 / 20 | 21.571634 |

## Why the scores trail the leader

All six weighted scores (306.42–311.65) exceed the top miner's 292.74. Fidelity is somewhat lower (0.9146–0.9219 versus 0.9599), but the dominant difference is consistency: **0.0702–0.0862 versus 0.561675**. Missing both the clean sets and HDR bands leaves all six on the low-consistency outcomes for this round. This establishes what happened on this task; it does not estimate long-run strategy performance. The top miner's bands and clean sets cannot be inferred from its aggregate score, and were not reconstructed.

## Submission history caveat

The uploaded rows came from the **17:19–17:23 UTC builds** and were submitted between 17:23:59 and 17:55:56 UTC. All six used prepared submissions with 283–286 seconds of upload TTL remaining. Later, at **18:54 UTC**, each hotkey rebuilt this same task after HEK293's cut search was narrowed to 300 seeds. Those later builds overwrite the task's `window_used.json` band entry and have different bands and clean counts. They were **not the uploaded submissions**: each hotkey has one recorded upload for this task, before the rebuild. This report therefore uses the uploaded archives and original logs, not the later overwritten band entries. h0's UID for this round was **248**, as shown by its upload target and the official score record.

## Verification and artifacts

- All six archived contracts match the published task, allowing for the subsequently stamped seed; their reference files also match after that same seed update.
- Each archived stage-1/2 valid-row list matches the submitted designs, with accessibility 0.35.
- All six archived three-seed final scores match the official API exactly.
- Stage-3 replay matches every row's archived outcome and indel length at all three real seeds.
- Scanned all 900 seeds for each hotkey: 1,350,000 row simulations; no new build or upload was performed.
- [Machine-readable analysis and full seed sets](analysis.json), [official score snapshot](scores.json), [published task snapshot](task.json), and [sanitized build/upload evidence](build_evidence.json).
- Reproduce with `.venv/bin/python reports/7bca65ce/audit.py` then `.venv/bin/python reports/7bca65ce/render_report.py` from the repository root.

Uploaded archives:

- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — [data/inst/niome_hotkey/result/2026-09-17T17:25:08](../../data/inst/niome_hotkey/result/2026-09-17T17:25:08/submission.json)
- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-17T17:34:01](../../data/inst/niome_hotkey1/result/2026-09-17T17:34:01/submission.json)
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-17T17:29:12](../../data/inst/niome_hotkey2/result/2026-09-17T17:29:12/submission.json)
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-17T17:29:10](../../data/inst/niome_hotkey3/result/2026-09-17T17:29:10/submission.json)
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — [data/inst/niome_hotkey4/result/2026-09-17T17:56:27](../../data/inst/niome_hotkey4/result/2026-09-17T17:56:27/submission.json)
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — [data/inst/niome_hotkey5/result/2026-09-17T17:24:31](../../data/inst/niome_hotkey5/result/2026-09-17T17:24:31/submission.json)
