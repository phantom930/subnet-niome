# Submission audit: e6caf784-f2af-48cc-938a-1377c441a1dd

**HUDEP-2**, created **2026-09-24T03:43:20.026079 UTC**. Scoring seeds: **481, 208, 809**.

**Best of h0–h9: h0, 115.215325, rank 9/248. Leader: UID 239, 153.196375.** h0 reached 75.21% of the leader, a gap of 37.981050. Top-10 cutoff: 108.499136.

All ten hotkeys uploaded **250 valid rows**. **h0–h6 used prepared conjunctions; h7–h9 used prepared all-HDR submissions.** No ordinary fallback occurred. All ten locally reproduced final scores match the official scores within 1e-10.

[Official score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=e6caf784-f2af-48cc-938a-1377c441a1dd&limit=40000), fetched 2026-09-24T09:32:28.020557+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-24T05:25:13.471431 UTC. Audit generated 2026-09-24T09:32:31.496970+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 239 | 1 | 153.196375 | 310.681567 | 0.506885 | 0.972800 | 250 | 0.29400000 |
| h0 | 226 | 9 | 115.215325 | 304.660064 | 0.403248 | 0.937826 | 250 | 0.01470000 |
| h1 | 189 | 20 | 31.493406 | 306.780090 | 0.110769 | 0.926779 | 250 | 0.00000000 |
| h2 | 45 | 13 | 47.707291 | 304.811067 | 0.167390 | 0.935026 | 250 | 0.00000000 |
| h3 | 136 | 29 | 28.509311 | 303.308210 | 0.100829 | 0.932221 | 250 | 0.00000000 |
| h4 | 228 | 38 | 26.664934 | 297.435487 | 0.095569 | 0.938058 | 250 | 0.00000000 |
| h5 | 49 | 36 | 26.896929 | 305.016567 | 0.095334 | 0.924977 | 250 | 0.00000000 |
| h6 | 227 | 16 | 33.235374 | 287.536086 | 0.121074 | 0.954676 | 250 | 0.00000000 |
| h7 | 116 | 40 | 26.462338 | 292.762529 | 0.099790 | 0.905786 | 250 | 0.00000000 |
| h8 | 7 | 54 | 23.881695 | 295.980660 | 0.088904 | 0.907570 | 250 | 0.00000000 |
| h9 | 143 | 10 | 108.499136 | 297.008617 | 0.408638 | 0.893961 | 250 | 0.00980000 |

Ranks count strictly higher scores; ties share rank. Identities were matched by full public hotkeys using this task's UIDs. Nonzero reported task weights: **h0: 0.014700 (1.47%)**; **h9: 0.009800 (0.98%)**. These are task weights from the score feed, not a calculation of on-chain payout. Final scores and factors are means across three seeds; multiplying mean factors need not reproduce the mean final score.

## Joined windows and construction

The task used a fixed layout independent of the prediction plan. **h0–h6** used seven **15-seed band-candidate windows at stride 15 over 200–299**, all on **loop 1**, with a common **200-seed cut search: 200–299 ∪ 400–499**. h6 wraps to **200–204 ∪ 290–299**, sharing five candidate seeds with h0.

**h7–h9** used the all-HDR builder directly over **34-seed candidate windows at stride 34 over 400–499**. h9 wraps to **400–401 ∪ 468–499**, sharing two candidate seeds with h7. This builder chooses Cas12a rows with a shared HDR set inside its window, then chooses Cas9 rows that return HDR on the same set. It does not use the conjunction's separate 200-seed cut search. The windows below come from this task's build logs and archived rows.

| Hotkey | Construction | Cas12a / Cas9 | Joined band candidates | Conjunction cut search | Logged clean statistic | Full clean /900 |
|---|---|---:|---|---|---|---:|
| h0 | conjunction | 100 / 150 | 200–214 | 200–299 ∪ 400–499 | 80/200 cut-clean | 81 |
| h1 | conjunction | 100 / 150 | 215–229 | 200–299 ∪ 400–499 | 74/200 cut-clean | 74 |
| h2 | conjunction | 100 / 150 | 230–244 | 200–299 ∪ 400–499 | 74/200 cut-clean | 77 |
| h3 | conjunction | 100 / 150 | 245–259 | 200–299 ∪ 400–499 | 80/200 cut-clean | 83 |
| h4 | conjunction | 100 / 150 | 260–274 | 200–299 ∪ 400–499 | 77/200 cut-clean | 79 |
| h5 | conjunction | 100 / 150 | 275–289 | 200–299 ∪ 400–499 | 77/200 cut-clean | 81 |
| h6 | conjunction | 100 / 150 | 200–204 ∪ 290–299 | 200–299 ∪ 400–499 | 70/200 cut-clean | 73 |
| h7 | all-HDR | 80 / 170 | 400–433 | not applicable | 14/34 HDR-clean | 19 |
| h8 | all-HDR | 80 / 170 | 434–467 | not applicable | 14/34 HDR-clean | 25 |
| h9 | all-HDR | 80 / 170 | 400–401 ∪ 468–499 | not applicable | 14/34 HDR-clean | 23 |

Conjunction logs count the Cas12a group's no-cut-free seeds inside its cut search. All-HDR logs count its all-HDR seeds inside the candidate window. **The full clean set always requires all 250 uploaded rows to cut**, replayed over all seeds 100–999. It may include seeds outside the construction window; Cas12a-only clean seeds are excluded where any Cas9 row fails to cut.

## Actual uploaded HDR bands and hits

An **HDR-band seed** makes all 250 uploaded rows return HDR. Each band below was recovered by replay and is a subset of its clean set and candidate window. A clean hit can contain mixed HDR and NHEJ outcomes.

**h7–h9's recorded band strings describe their 34-seed search windows, not their actual HDR seeds.** Replay recovers the exact bands below. h0–h6's recorded JSON band lists already identify the actual HDR seeds and match replay. The bounding `window` fields on h6/h9 omit the gaps in their joined windows.

| Hotkey | Clean /900 | HDR-band size | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits |
|---|---:|---:|---|---|---|
| h0 | 81 | 11 | 200, 202, 203, 204, 207, 208, 209, 210, 211, 213, 214 | 208 | 208 |
| h1 | 74 | 11 | 215, 219, 220, 221, 222, 223, 224, 225, 227, 228, 229 | 208 | none |
| h2 | 77 | 11 | 230, 231, 233, 234, 236, 237, 238, 240, 241, 242, 244 | 208 | none |
| h3 | 83 | 11 | 246, 248, 250, 251, 252, 253, 255, 256, 257, 258, 259 | 481 | none |
| h4 | 79 | 11 | 260, 261, 262, 263, 264, 266, 267, 269, 270, 271, 272 | none | none |
| h5 | 81 | 11 | 275, 276, 277, 279, 282, 283, 284, 285, 287, 288, 289 | none | none |
| h6 | 73 | 11 | 200, 201, 202, 203, 291, 293, 294, 295, 296, 298, 299 | 208, 481 | none |
| h7 | 19 | 14 | 400, 404, 405, 407, 408, 412, 413, 415, 416, 417, 420, 423, 425, 428 | none | none |
| h8 | 25 | 14 | 435, 440, 441, 442, 443, 444, 445, 446, 447, 448, 454, 456, 461, 462 | none | none |
| h9 | 23 | 14 | 469, 471, 474, 475, 477, 478, 481, 484, 485, 487, 488, 489, 490, 493 | 481 | 481 |

Fleet union: **235 distinct clean seeds** and **116 distinct HDR-band seeds**, from 119 band memberships. Repeated HDR seeds: 200 (h0/h6); 202 (h0/h6); 203 (h0/h6).

| Scoring seed | Hotkeys clean at this seed | Hotkeys all-HDR at this seed |
|---:|---|---|
| 481 | h3, h6, h9 | h9 |
| 208 | h0, h1, h2, h6 | h0 |
| 809 | none | none |

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.

| Hotkey | Seed 481 | Seed 208 | Seed 809 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 28.950690 / 132 / 5 | 285.718071 / 250 / 0 | 30.977214 / 132 / 6 | 115.215325 |
| h1 | 32.756702 / 122 / 7 | 33.026712 / 136 / 0 | 28.696803 / 139 / 8 | 31.493406 |
| h2 | 30.607092 / 146 / 5 | 83.083331 / 138 / 0 | 29.431450 / 132 / 6 | 47.707291 |
| h3 | 30.985369 / 123 / 0 | 26.936711 / 142 / 2 | 27.605853 / 128 / 7 | 28.509311 |
| h4 | 28.643332 / 132 / 4 | 29.583047 / 137 / 7 | 21.768423 / 119 / 11 | 26.664934 |
| h5 | 26.108427 / 124 / 6 | 32.304965 / 128 / 9 | 22.277394 / 130 / 4 | 26.896929 |
| h6 | 34.702409 / 145 / 0 | 37.788727 / 129 / 0 | 27.214985 / 127 / 6 | 33.235374 |
| h7 | 28.527213 / 123 / 8 | 25.568846 / 135 / 4 | 25.290956 / 127 / 8 | 26.462338 |
| h8 | 22.794398 / 130 / 3 | 20.983053 / 133 / 8 | 27.867633 / 139 / 6 | 23.881695 |
| h9 | 265.514136 / 250 / 0 | 32.523142 / 132 / 2 | 27.460129 / 129 / 5 | 108.499136 |

## Comparison with the leader

**h0 hit HDR seed 208 on all 250 rows**, scoring **285.718071** at that seed. Its three-seed mean was **115.215325**, rank **9**.

**h9 hit HDR seed 481 on all 250 rows**, scoring **265.514136** at that seed. Its three-seed mean was **108.499136**, rank **10**.

h0's weighted score was **304.660064**, versus the leader's **310.681567**. Mean consistency was **0.403248 versus 0.506885**, and fidelity was **0.937826 versus 0.972800**. The largest factor gap is consistency.

Seed 208 was in h0's candidate window. Seed 481 was in h9's candidate window and the conjunctions' cut search. Seed 809 was outside both groups' candidate windows and the conjunction cut search; the replay above determines actual coverage rather than inferring hits from window membership.

Leader hotkey: `5D4dexxheJrrFaH52KjUCwJ3cr3G4i3af4pKpasZP3HP7oeL`. The leader's uploaded rows, clean sets, HDR bands and per-seed outcomes were unavailable in the inspected data; its aggregate score does not identify those sets.

## Exact uploaded clean sets

Each list covers all seeds 100–999 and includes its uploaded HDR band.

### h0: 81 clean seeds

200, 202, 203, 204, 205, 207, 208, 209, 210, 211, 212, 213, 214, 215, 219, 221, 226, 227, 228, 230, 232, 237, 239, 242, 243, 248, 249, 250, 253, 254, 259, 269, 270, 275, 277, 278, 279, 281, 284, 285, 286, 288, 289, 290, 293, 400, 402, 403, 409, 413, 417, 419, 429, 430, 431, 432, 433, 436, 437, 452, 453, 455, 456, 459, 463, 464, 466, 469, 470, 472, 474, 477, 479, 482, 484, 487, 489, 494, 495, 499, 515.

Cas12a-only clean seeds excluded from the full-row set: 124, 133, 198, 369, 521, 552, 567, 607, 646, 661, 740, 781, 959.

### h1: 74 clean seeds

202, 204, 205, 208, 209, 211, 213, 215, 219, 220, 221, 222, 223, 224, 225, 227, 228, 229, 238, 243, 244, 245, 246, 250, 255, 257, 258, 260, 263, 264, 266, 267, 272, 273, 276, 279, 282, 286, 290, 296, 297, 400, 403, 406, 409, 416, 417, 422, 424, 429, 434, 438, 442, 443, 444, 445, 448, 452, 453, 458, 461, 462, 463, 465, 468, 470, 472, 477, 483, 484, 485, 494, 496, 498.

Cas12a-only clean seeds excluded from the full-row set: 176, 192, 393, 507, 636, 763, 882, 912.

### h2: 77 clean seeds

201, 205, 206, 208, 209, 211, 217, 220, 221, 224, 226, 227, 230, 231, 233, 234, 236, 237, 238, 240, 241, 242, 243, 244, 245, 248, 253, 258, 262, 271, 273, 277, 279, 281, 286, 293, 294, 295, 296, 297, 400, 403, 405, 406, 411, 416, 420, 422, 431, 433, 434, 435, 436, 437, 438, 439, 446, 448, 450, 451, 460, 461, 464, 465, 466, 468, 471, 473, 475, 479, 482, 494, 496, 497, 560, 796, 817.

Cas12a-only clean seeds excluded from the full-row set: 115, 178, 310, 534, 588, 594, 595, 634, 712, 791, 814, 881, 906, 917, 998.

### h3: 83 clean seeds

200, 201, 202, 204, 205, 207, 209, 213, 217, 220, 221, 224, 225, 227, 229, 231, 233, 237, 241, 246, 248, 250, 251, 252, 253, 255, 256, 257, 258, 259, 261, 264, 266, 273, 280, 282, 283, 294, 298, 400, 404, 405, 407, 412, 413, 416, 418, 420, 421, 422, 424, 427, 432, 434, 438, 442, 443, 444, 445, 446, 449, 450, 451, 455, 456, 462, 465, 468, 470, 471, 472, 473, 474, 481, 482, 483, 488, 492, 498, 499, 511, 630, 792.

Cas12a-only clean seeds excluded from the full-row set: 156, 184, 596, 610, 666, 687, 888, 981, 993.

### h4: 79 clean seeds

200, 201, 202, 204, 205, 206, 207, 212, 213, 216, 217, 221, 222, 223, 228, 229, 232, 233, 234, 238, 243, 247, 249, 250, 252, 253, 255, 260, 261, 262, 263, 264, 266, 267, 269, 270, 271, 272, 274, 277, 282, 286, 287, 288, 290, 293, 298, 402, 403, 410, 413, 414, 415, 417, 422, 427, 428, 430, 435, 436, 442, 447, 448, 457, 458, 460, 466, 468, 471, 478, 479, 483, 486, 487, 488, 495, 498, 730, 811.

Cas12a-only clean seeds excluded from the full-row set: 112, 394, 529, 561, 647, 791.

### h5: 81 clean seeds

199, 202, 203, 207, 210, 216, 218, 219, 223, 224, 231, 234, 236, 237, 239, 244, 245, 251, 253, 255, 258, 259, 266, 268, 270, 273, 274, 275, 276, 277, 278, 279, 281, 282, 283, 284, 285, 287, 288, 289, 290, 292, 295, 298, 299, 406, 408, 410, 412, 415, 428, 435, 436, 439, 440, 442, 445, 449, 450, 453, 455, 456, 457, 459, 462, 467, 473, 475, 477, 478, 480, 483, 484, 486, 490, 491, 497, 498, 806, 832, 898.

Cas12a-only clean seeds excluded from the full-row set: 145, 192, 546, 707, 710, 735, 742, 760, 793, 957.

### h6: 73 clean seeds

157, 200, 201, 202, 203, 204, 205, 208, 210, 211, 213, 217, 218, 219, 221, 223, 226, 234, 235, 236, 241, 242, 243, 245, 248, 261, 265, 272, 273, 275, 279, 281, 291, 292, 293, 294, 295, 296, 297, 298, 299, 403, 406, 413, 414, 415, 420, 421, 424, 431, 433, 440, 443, 445, 448, 454, 456, 458, 459, 461, 467, 472, 473, 474, 481, 484, 487, 492, 496, 498, 499, 535, 708.

Cas12a-only clean seeds excluded from the full-row set: 647, 802.

### h7: 19 clean seeds

131, 400, 404, 405, 407, 408, 412, 413, 415, 416, 417, 420, 423, 425, 428, 517, 652, 687, 704.

Cas12a-only clean seeds excluded from the full-row set: 113, 138, 152, 179, 235, 269, 283, 299, 353, 356, 403, 414, 504, 508, 552, 612, 627, 656, 747, 783, 810, 811, 821, 883, 918, 920, 967, 979, 981.

### h8: 25 clean seeds

191, 232, 435, 440, 441, 442, 443, 444, 445, 446, 447, 448, 454, 456, 461, 462, 508, 708, 744, 764, 851, 884, 934, 951, 982.

Cas12a-only clean seeds excluded from the full-row set: 103, 107, 143, 246, 257, 328, 339, 356, 357, 378, 400, 404, 411, 425, 483, 510, 523, 557, 569, 590, 606, 625, 644, 698, 724, 727, 732, 740, 786, 816, 917.

### h9: 23 clean seeds

353, 359, 368, 469, 471, 474, 475, 477, 478, 481, 484, 485, 487, 488, 489, 490, 493, 511, 595, 776, 791, 940, 993.

Cas12a-only clean seeds excluded from the full-row set: 114, 160, 192, 254, 262, 290, 291, 298, 305, 315, 332, 375, 383, 418, 429, 470, 482, 491, 514, 528, 530, 596, 626, 713, 716, 740, 750, 781, 792, 831, 836, 859, 910, 988.

## Verification and provenance

- All ten successful upload records identify this task. Contracts and references match the published task after substituting its scoring seed stamp.
- All 2,500 rows pass validation at HUDEP-2 accessibility 0.82. Every locally reproduced three-seed final score matches its official score within 1e-10.
- Replayed **900 × 250 × 10 = 2,250,000 row/seed combinations**. At each scoring seed, every row's outcome and indel length matches local stage-3 detail.
- Conjunction cut-clean counts match the build logs. All-HDR HDR-clean counts match the recovered bands. All uploads used completed prepared builds.
- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).
- Reproduce from the saved snapshots with `.venv/bin/python reports/e6caf784/validate.py`, then `.venv/bin/python reports/e6caf784/audit.py`, then `.venv/bin/python reports/e6caf784/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.

| Hotkey | Upload UTC | Prepared build completed UTC | Public hotkey | Submission archive |
|---|---|---|---|---|
| h0 | 2026-09-24 04:35:17,320 | 2026-09-24 03:46:20,372 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | [data/inst/niome_hotkey/result/2026-09-24T04:35:47](../../data/inst/niome_hotkey/result/2026-09-24T04:35:47/submission.json) |
| h1 | 2026-09-24 04:24:18,430 | 2026-09-24 03:46:16,024 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | [data/inst/niome_hotkey1/result/2026-09-24T04:24:48](../../data/inst/niome_hotkey1/result/2026-09-24T04:24:48/submission.json) |
| h2 | 2026-09-24 04:36:22,546 | 2026-09-24 03:46:19,112 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | [data/inst/niome_hotkey2/result/2026-09-24T04:36:52](../../data/inst/niome_hotkey2/result/2026-09-24T04:36:52/submission.json) |
| h3 | 2026-09-24 04:31:26,262 | 2026-09-24 03:46:19,852 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | [data/inst/niome_hotkey3/result/2026-09-24T04:31:55](../../data/inst/niome_hotkey3/result/2026-09-24T04:31:55/submission.json) |
| h4 | 2026-09-24 04:12:19,931 | 2026-09-24 03:46:18,598 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | [data/inst/niome_hotkey4/result/2026-09-24T04:12:49](../../data/inst/niome_hotkey4/result/2026-09-24T04:12:49/submission.json) |
| h5 | 2026-09-24 04:27:14,658 | 2026-09-24 03:46:18,251 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | [data/inst/niome_hotkey5/result/2026-09-24T04:27:44](../../data/inst/niome_hotkey5/result/2026-09-24T04:27:44/submission.json) |
| h6 | 2026-09-24 04:36:48,004 | 2026-09-24 03:46:16,254 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | [data/inst/niome_hotkey6/result/2026-09-24T04:37:17](../../data/inst/niome_hotkey6/result/2026-09-24T04:37:17/submission.json) |
| h7 | 2026-09-24 04:21:23,389 | 2026-09-24 03:44:57,075 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | [data/inst/niome_hotkey7/result/2026-09-24T04:21:53](../../data/inst/niome_hotkey7/result/2026-09-24T04:21:53/submission.json) |
| h8 | 2026-09-24 04:03:02,317 | 2026-09-24 03:44:57,018 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | [data/inst/niome_hotkey8/result/2026-09-24T04:03:32](../../data/inst/niome_hotkey8/result/2026-09-24T04:03:32/submission.json) |
| h9 | 2026-09-24 04:29:03,846 | 2026-09-24 03:44:56,004 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | [data/inst/niome_hotkey9/result/2026-09-24T04:29:32](../../data/inst/niome_hotkey9/result/2026-09-24T04:29:32/submission.json) |
