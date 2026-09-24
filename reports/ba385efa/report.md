# Submission audit: ba385efa-d25a-48ae-82dc-73473543bd01

**CD34+_HSPC**, created **2026-09-21T08:31:16.562618 UTC**. Scoring seeds: **263, 486, 269**.

**Best of h0–h9: h8, 87.923075, rank 13/248. Leader: UID 62, 149.014247.** h8 reached 59.00% of the leader, a gap of 61.091172. Top-10 cutoff: 116.884749.

All ten hotkeys uploaded **250 valid rows** each. Nine used prepared conjunctions with **100 Cas12a + 150 Cas9**, each with an **11-seed HDR band**. **h7 uploaded an ordinary fallback with 76 Cas12a + 174 Cas9.** Its prepared conjunction finished after upload, so the late 11-seed band in its window record is not submitted coverage.

[Official score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=ba385efa-d25a-48ae-82dc-73473543bd01&limit=40000), fetched 2026-09-21T10:55:26.143691+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-21T10:51:54.526607 UTC. Audit generated 2026-09-21T10:59:20.492871+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 62 | 1 | 149.014247 | 230.168414 | 0.668747 | 0.968100 | 250 | 0.294000 |
| h0 | 248 | 219 | 25.439591 | 247.524763 | 0.113619 | 0.904566 | 250 | 0.000000 |
| h1 | 189 | 179 | 35.921393 | 245.716267 | 0.160373 | 0.911568 | 250 | 0.000000 |
| h2 | 45 | 225 | 25.003760 | 245.681212 | 0.111534 | 0.912482 | 250 | 0.000000 |
| h3 | 136 | 220 | 25.419753 | 245.661275 | 0.112052 | 0.923457 | 250 | 0.000000 |
| h4 | 190 | 180 | 28.457368 | 244.413985 | 0.126190 | 0.922663 | 250 | 0.000000 |
| h5 | 234 | 161 | 41.587316 | 245.936704 | 0.186061 | 0.908831 | 250 | 0.000000 |
| h6 | 175 | 191 | 27.284422 | 246.095681 | 0.121150 | 0.915138 | 250 | 0.000000 |
| h7 | 65 | 238 | 21.607719 | 229.811577 | 0.101243 | 0.928695 | 250 | 0.000000 |
| h8 | 211 | 13 | 87.923075 | 247.680904 | 0.393374 | 0.902411 | 250 | 0.000000 |
| h9 | 245 | 234 | 23.374141 | 244.223966 | 0.104431 | 0.916468 | 250 | 0.000000 |

Ranks count strictly higher scores; ties share rank. Identities were matched using full public hotkeys. All ten have zero reward weight in this snapshot. Final scores and breakdown factors are means across three seeds; multiplying mean factors need not reproduce the mean final score.

## Joined seed windows and cut search

The round plan generated **08:17:19 UTC** predicted width-100 windows **100 / 200 / 300**, joined as **100–399**. The actual conjunction layout instead uses fixed windows **100 / 200 / 500**, joined as **100–299 ∪ 500–599**: 300 seeds shared by all loops. Both band candidates and narrow cut search use this fixed shared set. The plan does not steer these band candidates. The plan and shared set both contain scoring seeds **263 and 269**, while **486** is outside. The actual width-100 window labels are **200 / 400 / 200**.

h0–h9 were assigned loops 1–10, each excluding earlier loops' band seeds. The uploaded conjunctions are loops **1–7 and 9–10**. h7's assigned loop 8 completed too late and its uploaded fallback has no conjunction candidate or cut-search window.

| Hotkey | Uploaded construction | Uploaded loop | Band candidates and cut search | Clean in cut search | Full clean set /900 |
|---|---|---:|---|---:|---:|
| h0 | conjunction | 1 | 100–299 ∪ 500–599 | 108/300 | 113 |
| h1 | conjunction | 2 | 100–299 ∪ 500–599 | 111/300 | 115 |
| h2 | conjunction | 3 | 100–299 ∪ 500–599 | 101/300 | 104 |
| h3 | conjunction | 4 | 100–299 ∪ 500–599 | 107/300 | 110 |
| h4 | conjunction | 5 | 100–299 ∪ 500–599 | 104/300 | 105 |
| h5 | conjunction | 6 | 100–299 ∪ 500–599 | 107/300 | 111 |
| h6 | conjunction | 7 | 100–299 ∪ 500–599 | 105/300 | 109 |
| h7 | ordinary fallback | — | not applicable | — | 5 |
| h8 | conjunction | 9 | 100–299 ∪ 500–599 | 107/300 | 110 |
| h9 | conjunction | 10 | 100–299 ∪ 500–599 | 101/300 | 104 |

The full clean counts exceed the logged narrow-window counts because replay also finds some clean seeds outside the 300 construction seeds. Both counts require every uploaded row to cut; Cas12a-only clean seeds outside that window are excluded where any Cas9 row fails to cut.

## Uploaded HDR bands and scoring-seed hits

**Clean:** all 250 rows cut. **HDR band:** all 250 rows return HDR. Both sets were recovered from actual archived submissions by replay over all seeds 100–999. A clean hit can have mixed HDR/NHEJ outcomes. Candidate membership alone guarantees neither clean nor HDR.

| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits |
|---|---:|---|---|---|
| h0 | 113 | 105, 118, 197, 251, 253, 291, 533, 546, 576, 588, 593 | 269 | none |
| h1 | 115 | 130, 134, 139, 147, 175, 245, 293, 537, 569, 596, 598 | 269 | none |
| h2 | 104 | 123, 164, 199, 224, 250, 507, 525, 559, 564, 568, 584 | 269 | none |
| h3 | 110 | 122, 133, 149, 156, 166, 202, 249, 262, 500, 582, 591 | 263 | none |
| h4 | 105 | 113, 124, 170, 174, 178, 222, 229, 283, 295, 513, 527 | 269 | none |
| h5 | 111 | 182, 194, 198, 264, 274, 280, 519, 524, 560, 586, 592 | 263, 269 | none |
| h6 | 109 | 119, 129, 144, 235, 248, 254, 281, 284, 290, 532, 599 | 263 | none |
| h7 | 5 | none | none | none |
| h8 | 110 | 154, 186, 221, 240, 259, 261, 269, 536, 543, 558, 565 | 269 | 269 |
| h9 | 104 | 131, 153, 181, 210, 298, 504, 505, 521, 549, 578, 580 | none | none |

The uploaded fleet union contains **328 distinct clean seeds** and **99 distinct HDR-band seeds**. Of the clean union, **295 are inside the shared 300-seed window** and **33 are outside**. All nine uploaded bands are mutually disjoint. Fleet clean hits: 263, 269; HDR hits: 269.

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.

| Hotkey | Seed 263 | Seed 486 | Seed 269 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 26.315124 / 138 / 9 | 19.763018 / 127 / 6 | 30.240630 / 141 / 0 | 25.439591 |
| h1 | 21.973702 / 141 / 4 | 24.680032 / 134 / 3 | 61.110445 / 139 / 0 | 35.921393 |
| h2 | 23.476940 / 143 / 3 | 22.677583 / 100 / 5 | 28.856757 / 135 / 0 | 25.003760 |
| h3 | 28.835224 / 134 / 0 | 22.192190 / 127 / 6 | 25.231845 / 136 / 1 | 25.419753 |
| h4 | 22.706437 / 137 / 4 | 25.403445 / 125 / 2 | 37.262223 / 129 / 0 | 28.457368 |
| h5 | 51.396805 / 124 / 0 | 20.937062 / 125 / 10 | 52.428080 / 129 / 0 | 41.587316 |
| h6 | 39.045792 / 147 / 0 | 20.283516 / 135 / 9 | 22.523959 / 128 / 3 | 27.284422 |
| h7 | 22.691785 / 132 / 7 | 21.040383 / 124 / 8 | 21.090988 / 139 / 6 | 21.607719 |
| h8 | 23.066512 / 115 / 7 | 17.192637 / 122 / 3 | 223.510074 / 250 / 0 | 87.923075 |
| h9 | 21.605702 / 130 / 3 | 25.110091 / 135 / 3 | 23.406628 / 130 / 4 | 23.374141 |

## Comparison with the leader

h8's weighted score was **247.680904**, versus the leader's **230.168414**. Mean consistency was **0.393374 versus 0.668747**, and fidelity was **0.902411 versus 0.968100**. The largest factor gap is consistency.

h8 returned HDR on all 250 rows at **269**, scoring **223.510074** at that seed. Its final score averages this result with the other two seed scores shown above.

Leader hotkey: `5E4SsYpA8TdqeYEJFue7EfH8m9mCq8vraFcRmF5QiAFLWgRC`. The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were unavailable in the inspected data; its aggregate score does not identify those sets.

## h7 fallback and late prepared band

h7 received the task at **08:31:31 UTC** while prefetch was still building. After waiting, at **08:34:35** the runtime reported the prepared round unusable (still building at 187 seconds). It skipped hedges because the prepared build still held the GPU and generated an ordinary submission, which uploaded at **08:34:40.849**. The prepared conjunction finished only at **08:37:06.650**, with **99/300 clean seeds and an 11-seed band**. Those are preparation statistics, not statistics of h7's upload.

Late, unsubmitted h7 band: **115, 141, 158, 207, 214, 233, 297, 512, 526, 534, 574**. The actual fallback has **5 clean seeds over 100–999** and **0 all-HDR band seeds**. Its window record combines `source: all_hdr_not_attempted` with `all_hdr_built: true` and the late band. This audit uses the upload archive and timing to distinguish them.

## Exact uploaded clean sets

Each list covers all 100–999 and includes the corresponding uploaded HDR band. The full-row definition excludes seeds where only the Cas12a group cuts cleanly.

### h0: 113 clean seeds

100, 101, 103, 105, 107, 111, 112, 113, 115, 118, 122, 123, 128, 129, 133, 134, 140, 141, 151, 152, 153, 154, 169, 173, 174, 178, 181, 182, 185, 186, 188, 189, 192, 197, 198, 202, 207, 210, 212, 213, 220, 221, 225, 229, 232, 234, 238, 241, 242, 243, 247, 248, 249, 251, 253, 255, 260, 264, 267, 269, 270, 273, 274, 276, 277, 279, 285, 286, 291, 394, 410, 506, 507, 512, 513, 519, 521, 522, 523, 529, 530, 533, 535, 537, 540, 542, 546, 549, 551, 552, 554, 555, 558, 561, 563, 565, 568, 571, 573, 574, 576, 581, 583, 586, 587, 588, 592, 593, 594, 597, 795, 944, 995.

Cas12a-only clean seeds excluded from the full-row set: 338, 362, 707, 780, 832, 986.

### h1: 115 clean seeds

100, 101, 105, 109, 110, 111, 114, 116, 117, 118, 122, 124, 126, 130, 132, 134, 135, 137, 139, 141, 145, 146, 147, 149, 150, 153, 154, 155, 157, 158, 159, 163, 165, 169, 170, 171, 172, 173, 175, 179, 184, 190, 191, 192, 194, 199, 201, 207, 214, 216, 219, 232, 236, 238, 241, 242, 243, 245, 247, 250, 253, 256, 258, 261, 265, 267, 268, 269, 274, 277, 278, 279, 283, 286, 293, 294, 314, 502, 505, 507, 513, 515, 516, 517, 519, 521, 523, 524, 529, 537, 542, 547, 549, 551, 556, 566, 569, 570, 572, 574, 575, 577, 583, 584, 587, 588, 589, 590, 593, 596, 598, 599, 746, 826, 848.

Cas12a-only clean seeds excluded from the full-row set: 455, 608, 647, 784, 896.

### h2: 104 clean seeds

100, 102, 105, 106, 107, 108, 111, 114, 115, 120, 123, 125, 126, 127, 129, 134, 135, 139, 143, 149, 153, 157, 159, 161, 164, 165, 168, 170, 190, 195, 197, 199, 202, 207, 219, 222, 224, 226, 234, 235, 237, 239, 242, 244, 245, 250, 256, 261, 269, 274, 284, 287, 289, 290, 291, 293, 295, 298, 313, 447, 506, 507, 508, 509, 510, 512, 513, 514, 516, 518, 519, 524, 525, 526, 528, 537, 538, 539, 543, 545, 548, 553, 554, 555, 556, 557, 559, 561, 562, 564, 565, 568, 569, 570, 573, 577, 578, 579, 584, 587, 588, 591, 597, 780.

Cas12a-only clean seeds excluded from the full-row set: 301, 353, 354, 458, 717, 776, 779, 879, 955, 999.

### h3: 110 clean seeds

100, 104, 105, 106, 107, 109, 113, 114, 117, 118, 122, 127, 129, 131, 133, 139, 141, 142, 145, 146, 149, 151, 155, 156, 157, 158, 166, 167, 169, 174, 176, 178, 180, 181, 183, 190, 195, 199, 201, 202, 203, 204, 206, 207, 212, 214, 216, 223, 226, 229, 231, 237, 242, 243, 247, 248, 249, 250, 252, 257, 260, 261, 262, 263, 266, 271, 272, 277, 278, 279, 280, 281, 282, 284, 286, 290, 292, 294, 295, 304, 378, 500, 507, 510, 519, 520, 522, 523, 527, 529, 530, 538, 541, 552, 554, 559, 564, 568, 572, 577, 578, 579, 580, 582, 583, 584, 591, 593, 596, 613.

Cas12a-only clean seeds excluded from the full-row set: 310, 344, 404, 419, 440, 442, 446, 821, 834, 855, 882, 932.

### h4: 105 clean seeds

101, 103, 106, 107, 110, 113, 114, 116, 119, 124, 126, 138, 139, 140, 151, 152, 153, 159, 161, 163, 168, 169, 170, 171, 174, 175, 178, 181, 187, 188, 190, 191, 192, 197, 199, 201, 206, 210, 216, 222, 227, 229, 231, 234, 235, 241, 242, 243, 245, 249, 262, 264, 266, 267, 268, 269, 273, 278, 283, 284, 286, 288, 291, 293, 295, 298, 299, 501, 505, 513, 514, 516, 518, 521, 525, 526, 527, 528, 531, 532, 538, 543, 547, 553, 556, 557, 559, 560, 565, 566, 570, 571, 572, 574, 577, 582, 583, 588, 592, 593, 594, 595, 596, 597, 972.

Cas12a-only clean seeds excluded from the full-row set: 347, 409, 453, 487, 675, 722, 732, 815, 890, 971, 986.

### h5: 111 clean seeds

100, 102, 103, 108, 113, 118, 131, 132, 134, 136, 140, 141, 143, 148, 152, 154, 158, 159, 161, 162, 165, 174, 180, 182, 184, 185, 188, 189, 190, 194, 195, 196, 198, 200, 203, 205, 211, 215, 219, 222, 226, 227, 231, 242, 244, 247, 248, 250, 257, 260, 261, 262, 263, 264, 265, 267, 268, 269, 271, 272, 274, 276, 277, 280, 281, 288, 291, 297, 327, 330, 411, 500, 502, 505, 509, 510, 511, 514, 515, 517, 518, 519, 521, 524, 525, 526, 527, 534, 536, 540, 541, 547, 549, 550, 554, 555, 560, 562, 567, 569, 574, 576, 578, 581, 582, 586, 587, 592, 593, 595, 881.

Cas12a-only clean seeds excluded from the full-row set: 346, 463, 469, 643, 714, 724, 789, 816, 861, 940, 973, 991, 998.

### h6: 109 clean seeds

100, 105, 106, 109, 114, 115, 117, 119, 120, 121, 124, 125, 129, 131, 133, 137, 143, 144, 147, 150, 151, 153, 154, 155, 158, 160, 162, 166, 168, 169, 175, 181, 185, 186, 194, 195, 198, 201, 204, 213, 215, 220, 222, 225, 228, 229, 230, 235, 237, 241, 242, 243, 247, 248, 254, 255, 259, 261, 263, 265, 270, 274, 276, 278, 279, 281, 282, 283, 284, 288, 289, 290, 292, 297, 298, 474, 482, 500, 501, 503, 504, 506, 507, 518, 519, 520, 527, 530, 531, 532, 537, 543, 548, 555, 556, 557, 560, 562, 573, 580, 582, 586, 587, 589, 593, 594, 599, 656, 919.

Cas12a-only clean seeds excluded from the full-row set: 353, 376, 400, 494, 647, 863, 977, 998.

### h7: 5 clean seeds

237, 586, 689, 719, 841.

Cas12a-only clean seeds excluded from the full-row set: 138, 180, 224, 302, 305, 318, 324, 329, 373, 439, 455, 458, 482, 500, 506, 511, 532, 559, 622, 667, 730, 742, 783, 837, 846, 917, 929, 932, 970.

### h8: 110 clean seeds

100, 103, 108, 110, 126, 128, 129, 132, 136, 138, 139, 140, 148, 150, 152, 154, 155, 160, 163, 165, 169, 170, 172, 175, 176, 177, 183, 184, 186, 190, 192, 193, 196, 198, 201, 202, 204, 206, 212, 213, 215, 217, 219, 220, 221, 222, 227, 231, 233, 234, 235, 238, 240, 243, 246, 248, 255, 257, 258, 259, 261, 264, 269, 271, 272, 274, 279, 291, 294, 298, 299, 340, 431, 501, 502, 504, 505, 508, 515, 522, 524, 526, 527, 529, 531, 532, 534, 536, 537, 543, 546, 548, 550, 553, 555, 556, 557, 558, 565, 567, 568, 574, 581, 582, 583, 587, 590, 593, 597, 854.

Cas12a-only clean seeds excluded from the full-row set: 323, 441, 673, 684, 724, 801, 809, 826, 888.

### h9: 104 clean seeds

102, 104, 105, 112, 117, 122, 127, 131, 136, 138, 139, 140, 143, 144, 146, 148, 150, 151, 153, 155, 167, 168, 169, 171, 172, 175, 177, 178, 180, 181, 183, 185, 191, 196, 198, 199, 201, 210, 211, 212, 213, 220, 222, 225, 226, 239, 244, 245, 249, 251, 262, 266, 267, 268, 272, 273, 274, 275, 279, 284, 289, 294, 295, 298, 403, 500, 502, 504, 505, 506, 510, 513, 517, 519, 520, 521, 522, 524, 525, 528, 530, 537, 538, 542, 544, 549, 558, 560, 565, 566, 570, 571, 572, 576, 577, 578, 580, 581, 586, 588, 589, 598, 713, 742.

Cas12a-only clean seeds excluded from the full-row set: 327, 448, 915, 977.

## Verification and provenance

- All ten successful upload archives identify this task. Contracts and references match the published task after substituting its seed stamp.
- All 2,500 rows passed validation with CD34+_HSPC accessibility 0.87. All ten local three-seed final scores match official scores within 1e-10.
- Generated all ten validations under this report directory without modifying miner archives.
- Replayed **900 × 250 × 10 = 2,250,000 row/seed combinations**. Every outcome and indel length at the three scoring seeds matches local stage-3 detail.
- For nine conjunctions, Cas12a and full-row clean counts inside the cut-search window agree with build logs; recovered uploaded bands match the window records. h7's late prepared band is excluded.
- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).
- Reproduce with `.venv/bin/python reports/ba385efa/validate.py`, then `.venv/bin/python reports/ba385efa/audit.py`, then `.venv/bin/python reports/ba385efa/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.

| Hotkey | Upload UTC | Prepared build completed UTC | Prepared build uploaded | Public hotkey | Submission archive |
|---|---|---|---|---|---|
| h0 | 2026-09-21 08:54:22,970 | 2026-09-21 08:36:12,130 | True | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | [data/inst/niome_hotkey/result/2026-09-21T08:54:53](../../data/inst/niome_hotkey/result/2026-09-21T08:54:53/submission.json) |
| h1 | 2026-09-21 08:45:16,239 | 2026-09-21 08:36:21,905 | True | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | [data/inst/niome_hotkey1/result/2026-09-21T08:45:46](../../data/inst/niome_hotkey1/result/2026-09-21T08:45:46/submission.json) |
| h2 | 2026-09-21 08:43:12,386 | 2026-09-21 08:36:29,039 | True | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | [data/inst/niome_hotkey2/result/2026-09-21T08:43:42](../../data/inst/niome_hotkey2/result/2026-09-21T08:43:42/submission.json) |
| h3 | 2026-09-21 08:55:45,858 | 2026-09-21 08:36:44,639 | True | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | [data/inst/niome_hotkey3/result/2026-09-21T08:56:15](../../data/inst/niome_hotkey3/result/2026-09-21T08:56:15/submission.json) |
| h4 | 2026-09-21 08:45:16,451 | 2026-09-21 08:36:51,092 | True | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | [data/inst/niome_hotkey4/result/2026-09-21T08:45:46](../../data/inst/niome_hotkey4/result/2026-09-21T08:45:46/submission.json) |
| h5 | 2026-09-21 09:00:21,552 | 2026-09-21 08:36:57,833 | True | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | [data/inst/niome_hotkey5/result/2026-09-21T09:00:51](../../data/inst/niome_hotkey5/result/2026-09-21T09:00:51/submission.json) |
| h6 | 2026-09-21 09:17:28,824 | 2026-09-21 08:37:04,102 | True | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | [data/inst/niome_hotkey6/result/2026-09-21T09:17:58](../../data/inst/niome_hotkey6/result/2026-09-21T09:17:58/submission.json) |
| h7 | 2026-09-21 08:34:40,849 | 2026-09-21 08:37:06,650 | False | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | [data/inst/niome_hotkey7/result/2026-09-21T08:35:33](../../data/inst/niome_hotkey7/result/2026-09-21T08:35:33/submission.json) |
| h8 | 2026-09-21 09:18:58,886 | 2026-09-21 08:37:12,805 | True | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | [data/inst/niome_hotkey8/result/2026-09-21T09:19:28](../../data/inst/niome_hotkey8/result/2026-09-21T09:19:28/submission.json) |
| h9 | 2026-09-21 09:07:13,336 | 2026-09-21 08:37:13,411 | True | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | [data/inst/niome_hotkey9/result/2026-09-21T09:07:42](../../data/inst/niome_hotkey9/result/2026-09-21T09:07:42/submission.json) |
