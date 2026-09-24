# Submission audit: 35e41a29-5658-494b-9baa-7c83a18c25bf

**HEK293**, created **2026-09-24T10:54:37.147759 UTC**. Scoring seeds: **872, 548, 939**.

**Best of h0–h9: h3, 20.139513, rank 37/248. Leader: UID 22, 116.362948.** h3 reached 17.31% of the leader, a gap of 96.223435. Top-10 cutoff: 93.933806.

All ten hotkeys uploaded **250 valid rows**, comprising **80 Cas12a and 170 Cas9** rows each. **h0–h8 used prepared conjunctions; h9 used a prepared all-HDR submission.** All ten locally reproduced final scores match the official scores within 1e-10.

[Official score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=35e41a29-5658-494b-9baa-7c83a18c25bf&limit=40000), fetched 2026-09-24T13:13:18.012127+00:00. 248 records, 248 miners, 1 validator. Leader score timestamp: 2026-09-24T13:05:21.266908 UTC. Audit generated 2026-09-24T13:30:46.848095+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 22 | 1 | 116.362948 | 202.134784 | 0.592497 | 0.971600 | 250 | 0.29400003 |
| h0 | 226 | 74 | 19.145704 | 249.965507 | 0.083333 | 0.919122 | 250 | 0.00000000 |
| h1 | 189 | 89 | 18.822711 | 251.041319 | 0.081569 | 0.919204 | 250 | 0.00000000 |
| h2 | 45 | 79 | 19.068843 | 248.769222 | 0.083139 | 0.921985 | 250 | 0.00000000 |
| h3 | 136 | 37 | 20.139513 | 248.444691 | 0.087654 | 0.924799 | 250 | 0.00000000 |
| h4 | 228 | 119 | 18.315747 | 249.110594 | 0.079698 | 0.922543 | 250 | 0.00000000 |
| h5 | 49 | 165 | 16.457426 | 249.351088 | 0.071655 | 0.921095 | 250 | 0.00000000 |
| h6 | 227 | 161 | 16.697850 | 252.930799 | 0.073082 | 0.903338 | 250 | 0.00000000 |
| h7 | 116 | 80 | 19.062518 | 249.594845 | 0.083234 | 0.917581 | 250 | 0.00000000 |
| h8 | 7 | 81 | 19.017503 | 251.818408 | 0.082987 | 0.910033 | 250 | 0.00000000 |
| h9 | 143 | 114 | 18.363996 | 254.752417 | 0.082129 | 0.877711 | 250 | 0.00000000 |

Ranks count strictly higher scores; ties share rank. Identities were matched by full public hotkey using this task's UIDs. All ten hotkeys have zero reported task reward weight. These weights come from the task score feed. Final scores and factors are means across three seeds; multiplying mean factors need not reproduce the mean final score.

## Joined windows and construction

The fixed layout used **nine 12-seed candidate windows at stride 12 over 500–599** for h0–h8, with the common **200-seed cut search 100–199 ∪ 500–599**. h8 wraps to **500–507 ∪ 596–599**, sharing eight candidate seeds with h0. h9 used the all-HDR builder directly across **100–199**, without a separate conjunction cut search.

The cut search was reconstructed by testing all 36 pairs of 100-seed classes: only 100–199 ∪ 500–599 matches every hotkey’s logged Cas12a clean count. The live layout has since changed and was not used to infer these historical windows.

Build logs record loop 1 for h0–h6 and loop None for h7–h8. The conjunction function defaults to loop 1 and only overrides it for a non-None input; h7–h8 therefore also use the first loop. The recorded source evidence and replay checks are saved with the report.

| Hotkey | Construction | Joined band candidates | Conjunction cut search | Logged clean statistic | Full clean /900 |
|---|---|---|---|---|---:|
| h0 | conjunction | 500–511 | 100–199 ∪ 500–599 | 22/200 cut-clean | 22 |
| h1 | conjunction | 512–523 | 100–199 ∪ 500–599 | 22/200 cut-clean | 22 |
| h2 | conjunction | 524–535 | 100–199 ∪ 500–599 | 23/200 cut-clean | 23 |
| h3 | conjunction | 536–547 | 100–199 ∪ 500–599 | 20/200 cut-clean | 20 |
| h4 | conjunction | 548–559 | 100–199 ∪ 500–599 | 21/200 cut-clean | 21 |
| h5 | conjunction | 560–571 | 100–199 ∪ 500–599 | 22/200 cut-clean | 22 |
| h6 | conjunction | 572–583 | 100–199 ∪ 500–599 | 19/200 cut-clean | 19 |
| h7 | conjunction | 584–595 | 100–199 ∪ 500–599 | 22/200 cut-clean | 22 |
| h8 | conjunction | 500–507 ∪ 596–599 | 100–199 ∪ 500–599 | 19/200 cut-clean | 19 |
| h9 | all-HDR | 100–199 | not applicable | 7/100 HDR-clean | 7 |

Conjunction logs count the Cas12a group's no-cut-free seeds inside its cut search. The all-HDR log counts its shared HDR seeds inside its candidate window. A **full clean seed requires all 250 uploaded rows to cut**, replayed over seeds 100–999. It may include seeds outside the construction window. Cas12a-only clean seeds are excluded when any Cas9 row fails to cut.

## Actual uploaded HDR bands and hits

An **HDR-band seed** makes all 250 uploaded rows return HDR. Each recovered band is a subset of its clean set and candidate window. A clean hit may contain mixed HDR and NHEJ outcomes.

**h9's recorded band string 100–199 describes its candidate window. Replay recovers its seven actual HDR seeds below.** h0–h8's recorded JSON band lists match replay. h8's bounding window field 500–599 omits the gap in its joined 12-seed candidate window.

| Hotkey | Clean /900 | HDR-band size | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits |
|---|---:|---:|---|---|---|
| h0 | 22 | 8 | 500, 501, 504, 506, 507, 508, 509, 511 | none | none |
| h1 | 22 | 8 | 512, 515, 516, 518, 519, 521, 522, 523 | none | none |
| h2 | 23 | 8 | 524, 525, 526, 527, 532, 533, 534, 535 | none | none |
| h3 | 20 | 8 | 538, 540, 541, 542, 543, 544, 545, 546 | none | none |
| h4 | 21 | 8 | 549, 551, 552, 555, 556, 557, 558, 559 | none | none |
| h5 | 22 | 8 | 560, 561, 563, 564, 566, 567, 569, 570 | none | none |
| h6 | 19 | 8 | 572, 573, 577, 578, 580, 581, 582, 583 | none | none |
| h7 | 22 | 8 | 584, 585, 586, 589, 590, 592, 593, 595 | none | none |
| h8 | 19 | 8 | 500, 501, 502, 503, 504, 507, 597, 599 | none | none |
| h9 | 7 | 7 | 103, 105, 107, 119, 126, 154, 187 | none | none |

Fleet union: **140 distinct clean seeds** and **75 distinct HDR-band seeds**, from 79 band memberships. Repeated HDR seeds: 500 (h0/h8); 501 (h0/h8); 504 (h0/h8); 507 (h0/h8).

| Scoring seed | Hotkeys clean at this seed | Hotkeys all-HDR at this seed |
|---:|---|---|
| 872 | none | none |
| 548 | none | none |
| 939 | none | none |

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.

| Hotkey | Seed 872 | Seed 548 | Seed 939 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 20.265282 / 96 / 18 | 19.586835 / 105 / 14 | 17.584996 / 114 / 12 | 19.145704 |
| h1 | 18.485411 / 109 / 14 | 21.234478 / 116 / 12 | 16.748243 / 96 / 12 | 18.822711 |
| h2 | 17.602145 / 101 / 22 | 18.289830 / 107 / 15 | 21.314554 / 116 / 14 | 19.068843 |
| h3 | 17.648857 / 101 / 11 | 20.179355 / 123 / 15 | 22.590326 / 107 / 16 | 20.139513 |
| h4 | 20.790696 / 115 / 13 | 17.553247 / 102 / 13 | 16.603299 / 100 / 13 | 18.315747 |
| h5 | 16.498036 / 111 / 20 | 17.633003 / 96 / 16 | 15.241240 / 95 / 24 | 16.457426 |
| h6 | 20.348519 / 115 / 21 | 14.501454 / 108 / 20 | 15.243576 / 97 / 20 | 16.697850 |
| h7 | 17.254992 / 104 / 14 | 21.110588 / 110 / 12 | 18.821974 / 99 / 15 | 19.062518 |
| h8 | 18.774970 / 102 / 15 | 17.129250 / 103 / 22 | 21.148287 / 117 / 9 | 19.017503 |
| h9 | 18.401934 / 111 / 18 | 17.304725 / 108 / 21 | 19.385330 / 103 / 18 | 18.363996 |

## Comparison with the leader

h3's weighted score was **248.444691**, versus the leader's **202.134784**. Mean consistency was **0.087654 versus 0.592497**, and fidelity was **0.924799 versus 0.971600**.

Seed 548 was in h4's candidate window and the conjunctions' cut search, but was absent from h4's uploaded HDR band. Seeds 872 and 939 were outside both construction groups' candidate windows and the conjunction cut search. Actual coverage is determined by the replay above.

Leader hotkey: `5GLMBXc73iueF8vPSooSmujNVjs3HPuuAw7BDTEx9458j5eX`. The inspected data contains the leader's aggregate score, but does not provide its uploaded rows, clean sets, HDR bands or per-seed outcomes.

## Exact uploaded clean sets

Each list covers all seeds 100–999 and includes its uploaded HDR band.

### h0: 22 clean seeds

103, 104, 108, 112, 135, 140, 165, 166, 180, 500, 501, 504, 506, 507, 508, 509, 511, 542, 554, 575, 588, 593.

### h1: 22 clean seeds

112, 123, 126, 132, 134, 154, 174, 510, 512, 515, 516, 518, 519, 521, 522, 523, 532, 534, 535, 550, 575, 597.

### h2: 23 clean seeds

100, 128, 130, 162, 171, 181, 197, 524, 525, 526, 527, 532, 533, 534, 535, 558, 560, 573, 586, 590, 592, 594, 598.

### h3: 20 clean seeds

102, 117, 158, 172, 185, 193, 507, 529, 536, 538, 540, 541, 542, 543, 544, 545, 546, 554, 567, 572.

### h4: 21 clean seeds

109, 110, 115, 134, 167, 173, 179, 500, 518, 537, 547, 549, 551, 552, 555, 556, 557, 558, 559, 561, 597.

### h5: 22 clean seeds

125, 137, 145, 165, 175, 182, 196, 520, 525, 544, 550, 552, 560, 561, 563, 564, 566, 567, 569, 570, 575, 587.

### h6: 19 clean seeds

129, 133, 149, 153, 175, 183, 198, 507, 520, 521, 547, 572, 573, 577, 578, 580, 581, 582, 583.

### h7: 22 clean seeds

133, 137, 153, 155, 176, 181, 187, 515, 517, 544, 556, 570, 584, 585, 586, 588, 589, 590, 592, 593, 595, 596.

### h8: 19 clean seeds

106, 132, 169, 500, 501, 502, 503, 504, 507, 509, 530, 547, 549, 555, 562, 574, 577, 597, 599.

### h9: 7 clean seeds

103, 105, 107, 119, 126, 154, 187.

## Verification and provenance

- All ten successful upload records identify this task. Contracts and references match the published task after substituting its scoring seed stamp.
- All 2,500 rows pass validation at HEK293 accessibility 0.35. Every locally reproduced three-seed final score matches its official score within 1e-10.
- Replayed **900 × 250 × 10 = 2,250,000 row/seed combinations**. At each scoring seed, every row's outcome and indel length matches local stage-3 detail.
- Conjunction cut-clean counts match the build logs. The all-HDR HDR-clean count matches the recovered band. All uploads used completed prepared builds.
- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).
- Reproduce from saved snapshots with `.venv/bin/python reports/35e41a29/validate.py`, then `.venv/bin/python reports/35e41a29/audit.py`, then `.venv/bin/python reports/35e41a29/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.

| Hotkey | Upload UTC | Prepared build completed UTC | Public hotkey | Submission archive |
|---|---|---|---|---|
| h0 | 2026-09-24 11:30:39,862 | 2026-09-24 10:56:32,241 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | [data/inst/niome_hotkey/result/2026-09-24T11:31:09](../../data/inst/niome_hotkey/result/2026-09-24T11:31:09/submission.json) |
| h1 | 2026-09-24 11:44:38,393 | 2026-09-24 10:56:30,202 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | [data/inst/niome_hotkey1/result/2026-09-24T11:45:08](../../data/inst/niome_hotkey1/result/2026-09-24T11:45:08/submission.json) |
| h2 | 2026-09-24 11:29:24,279 | 2026-09-24 10:56:29,661 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | [data/inst/niome_hotkey2/result/2026-09-24T11:30:08](../../data/inst/niome_hotkey2/result/2026-09-24T11:30:08/submission.json) |
| h3 | 2026-09-24 11:44:20,487 | 2026-09-24 10:56:29,952 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | [data/inst/niome_hotkey3/result/2026-09-24T11:44:50](../../data/inst/niome_hotkey3/result/2026-09-24T11:44:50/submission.json) |
| h4 | 2026-09-24 11:07:18,158 | 2026-09-24 10:56:28,125 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | [data/inst/niome_hotkey4/result/2026-09-24T11:07:47](../../data/inst/niome_hotkey4/result/2026-09-24T11:07:47/submission.json) |
| h5 | 2026-09-24 10:59:35,646 | 2026-09-24 10:56:27,678 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | [data/inst/niome_hotkey5/result/2026-09-24T11:00:04](../../data/inst/niome_hotkey5/result/2026-09-24T11:00:04/submission.json) |
| h6 | 2026-09-24 11:45:16,424 | 2026-09-24 10:56:30,369 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | [data/inst/niome_hotkey6/result/2026-09-24T11:45:45](../../data/inst/niome_hotkey6/result/2026-09-24T11:45:45/submission.json) |
| h7 | 2026-09-24 11:17:41,544 | 2026-09-24 10:56:30,008 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | [data/inst/niome_hotkey7/result/2026-09-24T11:18:11](../../data/inst/niome_hotkey7/result/2026-09-24T11:18:11/submission.json) |
| h8 | 2026-09-24 11:24:53,576 | 2026-09-24 10:56:30,247 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | [data/inst/niome_hotkey8/result/2026-09-24T11:25:23](../../data/inst/niome_hotkey8/result/2026-09-24T11:25:23/submission.json) |
| h9 | 2026-09-24 10:56:01,141 | 2026-09-24 10:54:55,881 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | [data/inst/niome_hotkey9/result/2026-09-24T10:56:30](../../data/inst/niome_hotkey9/result/2026-09-24T10:56:30/submission.json) |
