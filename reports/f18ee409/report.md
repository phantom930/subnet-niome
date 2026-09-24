# Submission audit: f18ee409-bdb0-4e80-9818-e25f9f456f46

**HEK293**, created **2026-09-23T15:43:17.518522 UTC**. Scoring seeds: **214, 249, 784**.

**Best of h0–h9: h4, 17.296429, rank 31/248. Leader: UID 64, 118.530507.** h4 reached 14.59% of the leader, a gap of 101.234078. Top-10 cutoff: 79.718066.

All ten hotkeys uploaded **250 valid rows: 80 Cas12a + 170 Cas9**. **h0–h6 used prepared conjunctions; h7–h9 used prepared all-HDR submissions.** No ordinary fallback occurred. Every uploaded submission has **eight actual HDR-band seeds**. All ten locally reproduced final scores match the official scores within 1e-10.

Fleet scoring-seed clean hits: **214**. Fleet HDR-band hits: **none**.

**Only h0 was clean at seed 214**: all 250 rows cut, with 111 HDR rows, for a per-seed score of **21.597457**. No uploaded clean set contained 249 or 784. None of the three scoring seeds hit an uploaded all-HDR band.

[Official score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=f18ee409-bdb0-4e80-9818-e25f9f456f46&limit=40000), fetched 2026-09-23T18:08:55.223885+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-23T17:55:32.548885 UTC. Audit generated 2026-09-23T18:13:19.303439+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 64 | 1 | 118.530507 | 215.528826 | 0.570253 | 0.964400 | 250 | 0.294000 |
| h0 | 226 | 34 | 17.213616 | 203.443335 | 0.091002 | 0.929777 | 250 | 0.000000 |
| h1 | 189 | 108 | 15.913998 | 205.159995 | 0.084416 | 0.918890 | 250 | 0.000000 |
| h2 | 45 | 173 | 14.098758 | 205.528707 | 0.074328 | 0.922908 | 250 | 0.000000 |
| h3 | 136 | 170 | 14.233496 | 205.586921 | 0.074900 | 0.924341 | 250 | 0.000000 |
| h4 | 228 | 31 | 17.296429 | 207.526121 | 0.090600 | 0.919927 | 250 | 0.000000 |
| h5 | 49 | 179 | 13.488907 | 205.964789 | 0.070893 | 0.923808 | 250 | 0.000000 |
| h6 | 227 | 182 | 13.250904 | 202.571313 | 0.070347 | 0.929868 | 250 | 0.000000 |
| h7 | 116 | 142 | 15.183356 | 198.940856 | 0.081416 | 0.937418 | 250 | 0.000000 |
| h8 | 7 | 177 | 13.690927 | 195.798587 | 0.074330 | 0.940717 | 250 | 0.000000 |
| h9 | 143 | 168 | 14.376825 | 196.557337 | 0.078067 | 0.936930 | 250 | 0.000000 |

Ranks count strictly higher scores, so ties share rank. Identities were matched by full public hotkeys using this task's UIDs. All ten requested hotkeys have zero reward weight in this snapshot. Final scores and breakdown factors are means across three seeds; the product of mean factors need not reproduce the mean final score.

## Joined windows and construction

The task used a fixed layout independent of the prediction plan. **h0–h6** used seven **15-seed candidate windows at stride 15 over 200–299**, every hotkey on **loop 1**, with a common **200-seed cut search: 200–299 ∪ 400–499**. h6 wraps to **200–204 ∪ 290–299**, sharing five candidate seeds with h0.

**h7–h9** used the all-HDR builder directly, over **34-seed candidate windows at stride 34 over 400–499**. h9 wraps to **400–401 ∪ 468–499**, sharing two candidate seeds with h7. This builder selects Cas12a rows to share HDR seeds inside that window, then selects Cas9 rows that return HDR on those same seeds. It does not use the conjunction's separate 200-seed cut search.

**The current layout was moved after these builds.** Its present constants target 100–199 and 500–599. This report uses the task's historical build logs, archived rows and recovered bands, which establish the 200–299 and 400–499 windows actually uploaded.

| Hotkey | Uploaded construction | Joined band candidates | Candidate count | Conjunction cut search | Logged clean statistic | Full clean /900 |
|---|---|---|---:|---|---|---:|
| h0 | conjunction | 200–214 | 15 | 200–299 ∪ 400–499 | 22/200 cut-clean | 22 |
| h1 | conjunction | 215–229 | 15 | 200–299 ∪ 400–499 | 21/200 cut-clean | 21 |
| h2 | conjunction | 230–244 | 15 | 200–299 ∪ 400–499 | 21/200 cut-clean | 21 |
| h3 | conjunction | 245–259 | 15 | 200–299 ∪ 400–499 | 20/200 cut-clean | 20 |
| h4 | conjunction | 260–274 | 15 | 200–299 ∪ 400–499 | 19/200 cut-clean | 19 |
| h5 | conjunction | 275–289 | 15 | 200–299 ∪ 400–499 | 20/200 cut-clean | 20 |
| h6 | conjunction | 200–204 ∪ 290–299 | 15 | 200–299 ∪ 400–499 | 22/200 cut-clean | 22 |
| h7 | all-HDR | 400–433 | 34 | not applicable | 8/34 HDR-clean | 8 |
| h8 | all-HDR | 434–467 | 34 | not applicable | 8/34 HDR-clean | 8 |
| h9 | all-HDR | 400–401 ∪ 468–499 | 34 | not applicable | 8/34 HDR-clean | 8 |

The logged clean counts have different definitions: conjunction logs count the Cas12a group's **no-cut-free** seeds inside the cut window; all-HDR logs count its **all-HDR** seeds inside the candidate window. The report's full clean set always requires **every one of the 250 uploaded rows to cut**, and is recovered over all 900 seeds.

## Actual uploaded bands and hits

An **HDR-band seed** makes all 250 uploaded rows return HDR. Each actual band below was recovered by replay and is a subset of that hotkey's clean set and candidate window. A clean seed need not be all-HDR.

**h7–h9's `window_used.json` band strings describe the 34-seed candidate windows, not 34 actual HDR seeds.** Their build logs report `clean 8/34`; replay identifies the exact eight HDR seeds. h0–h6's JSON band lists already contain their actual eight HDR seeds and match replay. The bounding `window` fields on h6 and h9 omit the gaps in their joined sets.

| Hotkey | Full clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits |
|---|---:|---|---|---|
| h0 | 22 | 201, 202, 203, 205, 206, 207, 211, 213 | 214 | none |
| h1 | 21 | 216, 219, 221, 222, 223, 224, 226, 227 | none | none |
| h2 | 21 | 231, 232, 235, 236, 237, 241, 243, 244 | none | none |
| h3 | 20 | 248, 250, 252, 254, 255, 257, 258, 259 | none | none |
| h4 | 19 | 261, 262, 266, 268, 269, 270, 271, 273 | none | none |
| h5 | 20 | 276, 277, 280, 282, 283, 285, 286, 288 | none | none |
| h6 | 22 | 202, 203, 204, 292, 293, 295, 296, 298 | none | none |
| h7 | 8 | 404, 405, 414, 423, 426, 429, 430, 433 | none | none |
| h8 | 8 | 434, 438, 440, 441, 445, 463, 464, 465 | none | none |
| h9 | 8 | 473, 474, 480, 481, 483, 487, 491, 498 | none | none |

Fleet union: **121 distinct clean seeds** and **78 distinct HDR-band seeds** from 80 band memberships. Overlapping HDR seeds: 202 (h0/h6); 203 (h0/h6).

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.

| Hotkey | Seed 214 | Seed 249 | Seed 784 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 21.597457 / 111 / 0 | 14.921023 / 103 / 18 | 15.122367 / 112 / 15 | 17.213616 |
| h1 | 17.299701 / 104 / 12 | 14.355290 / 117 / 15 | 16.087001 / 113 / 18 | 15.913998 |
| h2 | 14.092849 / 104 / 18 | 17.445274 / 115 / 13 | 10.758152 / 96 / 22 | 14.098758 |
| h3 | 12.553075 / 100 / 19 | 16.273848 / 128 / 13 | 13.873564 / 99 / 19 | 14.233496 |
| h4 | 14.780370 / 103 / 21 | 18.882396 / 100 / 15 | 18.226522 / 99 / 13 | 17.296429 |
| h5 | 10.113151 / 108 / 28 | 15.236897 / 102 / 16 | 15.116673 / 99 / 21 | 13.488907 |
| h6 | 14.603598 / 105 / 20 | 13.807963 / 103 / 13 | 11.341150 / 107 / 18 | 13.250904 |
| h7 | 15.065654 / 108 / 16 | 15.311683 / 124 / 16 | 15.172731 / 121 / 20 | 15.183356 |
| h8 | 14.601201 / 109 / 15 | 10.959704 / 113 / 19 | 15.511876 / 119 / 15 | 13.690927 |
| h9 | 13.177656 / 96 / 21 | 15.029540 / 95 / 28 | 14.923280 / 94 / 21 | 14.376825 |

## Comparison with the leader

h4's weighted score was **207.526121**, versus the leader's **215.528826**. Mean consistency was **0.090600 versus 0.570253**, and fidelity was **0.919927 versus 0.964400**. The largest factor gap is consistency.

Seeds 214 and 249 were within the conjunction fleet's candidate class, specifically h0 and h3's windows, but candidate membership is not a guarantee of clean or HDR coverage. Seed 784 was outside both construction groups' candidate windows and the conjunction cut search. The all-HDR group targeted class 400–499; this task drew no seed in that class.

Leader hotkey: `5E2GMrert3uDw1VhhLJqyHkayvew2nnPvHPnKEi1qLGBx6av`. The leader's submitted rows, exact clean sets, HDR bands and per-seed outcomes were unavailable in the inspected data. Its aggregate score does not identify those sets. One task does not establish which construction performs better across tasks.

## Exact uploaded clean sets

Each list covers all seeds 100–999 and includes the corresponding uploaded HDR band.

### h0: 22 clean seeds

201, 202, 203, 205, 206, 207, 211, 213, 214, 228, 261, 292, 400, 405, 406, 439, 441, 442, 454, 459, 460, 478.

### h1: 21 clean seeds

216, 219, 221, 222, 223, 224, 226, 227, 231, 277, 283, 290, 409, 413, 427, 430, 453, 479, 492, 494, 498.

### h2: 21 clean seeds

201, 211, 222, 227, 231, 232, 235, 236, 237, 238, 241, 243, 244, 250, 264, 285, 292, 435, 438, 445, 452.

### h3: 20 clean seeds

205, 222, 229, 243, 248, 250, 252, 254, 255, 257, 258, 259, 272, 281, 298, 421, 426, 429, 435, 464.

### h4: 19 clean seeds

216, 242, 253, 261, 262, 266, 268, 269, 270, 271, 273, 277, 401, 403, 414, 424, 439, 468, 494.

### h5: 20 clean seeds

205, 209, 232, 233, 239, 245, 250, 276, 277, 280, 282, 283, 285, 286, 288, 402, 405, 424, 473, 488.

### h6: 22 clean seeds

202, 203, 204, 210, 235, 245, 280, 282, 292, 293, 295, 296, 298, 407, 409, 410, 415, 442, 445, 452, 488, 493.

### h7: 8 clean seeds

404, 405, 414, 423, 426, 429, 430, 433.

### h8: 8 clean seeds

434, 438, 440, 441, 445, 463, 464, 465.

### h9: 8 clean seeds

473, 474, 480, 481, 483, 487, 491, 498.

## Verification and provenance

- All ten successful upload records identify this task. Contracts and references match the published task after substituting its scoring seed stamp.
- All 2,500 rows passed validation with HEK293 accessibility 0.35. All ten local three-seed final scores match the official scores within 1e-10.
- All validation artifacts were generated under this report directory; miner archives were read without modification.
- Replayed **900 × 250 × 10 = 2,250,000 row/seed combinations**. At all three scoring seeds, every row's outcome and indel length matches local stage-3 detail.
- Conjunction cut-clean counts match all seven build logs. All-HDR HDR-clean counts match the three 8/34 logs and the actual eight-seed bands. Every upload used its completed prepared build.
- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [historical layout evidence](layout_evidence.json).
- Reproduce with `.venv/bin/python reports/f18ee409/validate.py`, then `.venv/bin/python reports/f18ee409/audit.py`, then `.venv/bin/python reports/f18ee409/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.

| Hotkey | Upload UTC | Prepared build completed UTC | Public hotkey | Submission archive |
|---|---|---|---|---|
| h0 | 2026-09-23 16:35:07,891 | 2026-09-23 15:45:31,248 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | [data/inst/niome_hotkey/result/2026-09-23T16:35:36](../../data/inst/niome_hotkey/result/2026-09-23T16:35:36/submission.json) |
| h1 | 2026-09-23 16:16:32,704 | 2026-09-23 15:45:32,353 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | [data/inst/niome_hotkey1/result/2026-09-23T16:17:01](../../data/inst/niome_hotkey1/result/2026-09-23T16:17:01/submission.json) |
| h2 | 2026-09-23 16:44:21,791 | 2026-09-23 15:45:31,407 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | [data/inst/niome_hotkey2/result/2026-09-23T16:44:50](../../data/inst/niome_hotkey2/result/2026-09-23T16:44:50/submission.json) |
| h3 | 2026-09-23 16:09:12,193 | 2026-09-23 15:45:32,092 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | [data/inst/niome_hotkey3/result/2026-09-23T16:09:41](../../data/inst/niome_hotkey3/result/2026-09-23T16:09:41/submission.json) |
| h4 | 2026-09-23 16:18:54,878 | 2026-09-23 15:45:35,621 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | [data/inst/niome_hotkey4/result/2026-09-23T16:19:24](../../data/inst/niome_hotkey4/result/2026-09-23T16:19:24/submission.json) |
| h5 | 2026-09-23 16:11:24,316 | 2026-09-23 15:45:32,525 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | [data/inst/niome_hotkey5/result/2026-09-23T16:11:53](../../data/inst/niome_hotkey5/result/2026-09-23T16:11:53/submission.json) |
| h6 | 2026-09-23 16:07:16,119 | 2026-09-23 15:45:31,549 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | [data/inst/niome_hotkey6/result/2026-09-23T16:07:45](../../data/inst/niome_hotkey6/result/2026-09-23T16:07:45/submission.json) |
| h7 | 2026-09-23 15:51:05,551 | 2026-09-23 15:44:22,380 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | [data/inst/niome_hotkey7/result/2026-09-23T15:51:35](../../data/inst/niome_hotkey7/result/2026-09-23T15:51:35/submission.json) |
| h8 | 2026-09-23 15:55:09,917 | 2026-09-23 15:44:22,609 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | [data/inst/niome_hotkey8/result/2026-09-23T15:55:39](../../data/inst/niome_hotkey8/result/2026-09-23T15:55:39/submission.json) |
| h9 | 2026-09-23 16:07:41,882 | 2026-09-23 15:44:22,077 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | [data/inst/niome_hotkey9/result/2026-09-23T16:08:11](../../data/inst/niome_hotkey9/result/2026-09-23T16:08:11/submission.json) |
