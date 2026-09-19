# Submission audit: 11276911-854f-4809-8610-e7d5c7bc8772

**CD34+_HSPC**, created **2026-09-18T05:25:14.214966 UTC**. Scoring seeds: **448, 931, 493**.

**Best of h0–h6: h3, 32.050277, rank 172/248. Leader: UID 124, 160.154107.** h3 reached 20.01% of the leader, with a gap of 128.103830. The top-10 cutoff was 80.651690, 48.601413 above h3. All seven uploaded successfully, have 250 valid rows, and received zero reward weight.

**h0 and h6 submitted identical ordinary HDR-construction fallbacks.** Their prepared conjunctions finished after their wait deadlines and uploads. h1–h5 submitted their prepared conjunction rows. The prepared-band records for h0/h6 do not describe their actual uploads.

Scores: [task-specific public API](https://niome-api.genomes.io/api/v3/miners/scores?task_id=11276911-854f-4809-8610-e7d5c7bc8772&limit=40000), fetched 2026-09-18T07:57:06.109925+00:00; 248 records, 248 distinct miners, 1 validator. Leader's score timestamp: 2026-09-18T07:49:04.977082 UTC. Audit generated 2026-09-18T08:03:15.510957+00:00.

## Official score comparison

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 124 | 1 | 160.154107 | 235.902612 | 0.694242 | 0.977900 | 0.294000 |
| h0 | 248 | 232 | 23.477629 | 236.500371 | 0.106477 | 0.932327 | 0.000000 |
| h1 | 41 | 184 | 27.840058 | 242.760513 | 0.122424 | 0.936751 | 0.000000 |
| h2 | 163 | 198 | 26.813198 | 241.007092 | 0.118217 | 0.941110 | 0.000000 |
| h3 | 118 | 172 | 32.050277 | 236.876348 | 0.142998 | 0.946195 | 0.000000 |
| h4 | 190 | 213 | 25.758486 | 244.427175 | 0.113305 | 0.930083 | 0.000000 |
| h5 | 234 | 234 | 23.137120 | 242.926591 | 0.101378 | 0.939485 | 0.000000 |
| h6 | 175 | 232 | 23.477629 | 236.500371 | 0.106477 | 0.932327 | 0.000000 |

Ranks are calculated from official scores using competition ranking: exact ties share rank. h0 and h6 both rank **232**; they occupy sorted positions 233 and 232 when API order is retained within the tie. The submitted files have identical SHA-256 hashes.

Leader hotkey: `5Cr2i6SoUWpbYm9DLZktweW5Dw2BbqCfD6JmrMQUDQ5nbYMR`.

Final score is the mean of the three per-seed scores. For these seven submissions, weighted score and fidelity are constant across seeds, so final = weighted score × mean consistency × fidelity.

## Joined windows, clean sets and actual bands

The plan generated **2026-09-18 05:17:11 UTC** selected joined band space **300–399 ∪ 600–799**, or 100-wide windows **300/600/700**. Conjunction cut search covered **100–999 (900 seeds)**. The HDR candidate slices were 150 seeds wide, at offsets **0/50/100/150/200/250/25** for h0–h6. Each of the five submitted conjunctions has **100 Cas12a + 150 Cas9 rows** and an **11-seed HDR band**. The two ordinary fallbacks each have **74 Cas12a + 176 Cas9 rows**.

**Clean set:** every uploaded row cuts at that seed. **HDR band:** every uploaded row returns HDR. Both are measured from the exact uploaded rows over all seeds 100–999. For submitted conjunctions, the full-row clean set equals the Cas12a group clean set and its size agrees with the build log. Their recovered HDR bands exactly match the original band records.

| Hotkey | Uploaded construction | Planned candidate slice | Offset | Clean /900 | Clean inside joined /300 | Clean outside joined |
|---|---|---|---:|---:|---:|---:|
| h0 | ordinary fallback | 300–399, 600–649 (unused) | 0 | 4 | 1 | 3 |
| h1 | conjunction | 350–399, 600–699 | 50 | 142 | 52 | 90 |
| h2 | conjunction | 600–749 | 100 | 141 | 54 | 87 |
| h3 | conjunction | 650–799 | 150 | 138 | 48 | 90 |
| h4 | conjunction | 300–349, 700–799 | 200 | 148 | 54 | 94 |
| h5 | conjunction | 300–399, 750–799 | 250 | 148 | 52 | 96 |
| h6 | ordinary fallback | 325–399, 600–674 (unused) | 25 | 4 | 1 | 3 |

| Hotkey | Actual uploaded HDR band seeds | Clean scoring-seed hits | HDR scoring-seed hits |
|---|---|---|---|
| h0 | none | none | none |
| h1 | 362, 368, 388, 607, 624, 625, 633, 636, 641, 648, 666 | 493 | none |
| h2 | 625, 642, 648, 661, 665, 666, 670, 685, 708, 727, 732 | none | none |
| h3 | 690, 691, 697, 721, 730, 739, 746, 764, 784, 787, 795 | 493 | none |
| h4 | 330, 334, 727, 731, 732, 738, 742, 744, 757, 762, 781 | 931 | none |
| h5 | 309, 313, 339, 341, 350, 374, 385, 391, 753, 761, 792 | none | none |
| h6 | none | none | none |

Candidate slices are reconstructed from the logged joined space and configured offsets. All submitted conjunction bands lie inside the corresponding candidate slice.

Fleet union: **523/900 clean seeds**, **50 distinct HDR-band seeds** (55 band memberships before overlap). Clean scoring-seed hits: **493, 931**; HDR-band hits: **none**.

All three scoring seeds **448/931/493** are outside the predicted joined band space. The table above identifies any all-cut hits recovered by the wider cut search.

## Per-seed scores and outcomes

Each cell is **per-seed final score / no-cut rows**. Zero no-cut rows means an all-cut clean hit.

| Hotkey | Seed 448 | Seed 931 | Seed 493 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 24.040181 / 4 | 24.980909 / 2 | 21.411796 / 5 | 23.477629 |
| h1 | 20.465932 / 6 | 26.209721 / 2 | 36.844521 / 0 | 27.840058 |
| h2 | 28.111269 / 4 | 27.783566 / 2 | 24.544758 / 5 | 26.813198 |
| h3 | 23.126314 / 4 | 22.502303 / 7 | 50.522214 / 0 | 32.050277 |
| h4 | 20.756698 / 9 | 32.177518 / 0 | 24.341242 / 3 | 25.758486 |
| h5 | 22.503410 / 4 | 22.231647 / 5 | 24.676302 / 2 | 23.137120 |
| h6 | 24.040181 / 4 | 24.980909 / 2 | 21.411796 / 5 | 23.477629 |

## Comparison with the leader

h3's weighted score is **236.876348**, compared with the leader's **235.902612**. Fidelity is **0.946195 versus 0.977900**. The main difference is mean consistency: **0.142998 versus 0.694242**. The leader's clean set, HDR band and per-seed outcomes cannot be recovered from its public aggregate score; its submitted rows were not available in the inspected data.

## h0 and h6 fallback provenance

| Hotkey | Prepared wait expired (UTC) | Ordinary upload succeeded (UTC) | Prepared conjunction finished (UTC) | Unused prepared clean count | Uploaded clean count |
|---|---|---|---|---:|---:|
| h0 | 2026-09-18 05:30:11.000 | 2026-09-18 05:30:15.270 | 2026-09-18 05:34:29.573 | 142 | 4 |
| h6 | 2026-09-18 05:31:50.000 | 2026-09-18 05:31:53.526 | 2026-09-18 05:34:24.163 | 150 | 4 |

Both handlers waited about 179–180 seconds for their prepared builds. When those waits expired, the GPU was still occupied, so they used ordinary construction. Both uploads succeeded with about 101 seconds of URL TTL remaining. Their prepared conjunctions completed later, so the completed bands were never submitted.

h5 also received a request before its build finished, but its wait lasted only about 36 seconds; it successfully uploaded the prepared conjunction with 247 seconds of TTL remaining.

The final `window_used.json` records for h0/h6 combine `all_hdr_not_attempted` with `all_hdr_built: true` and bands written by the later prepared builds. Uploaded-row replay and exact official score agreement establish which rows were actually submitted.

| Hotkey | Recorded prepared band — not submitted |
|---|---|
| h0 | 307, 317, 324, 327, 330, 336, 397, 614, 617, 625, 648 |
| h6 | 330, 339, 365, 375, 614, 617, 619, 625, 642, 648, 668 |

## Exact uploaded clean sets

These lists include HDR-band seeds and cover all seeds 100–999. Full arrays are also saved in [analysis.json](analysis.json).

### h0: 4 clean seeds

175, 712, 981, 982.

### h1: 142 clean seeds

101, 103, 110, 119, 124, 127, 135, 146, 151, 157, 161, 165, 168, 173, 184, 194, 200, 213, 219, 222, 225, 234, 242, 257, 261, 269, 270, 273, 276, 282, 284, 294, 296, 305, 311, 326, 331, 334, 340, 349, 353, 362, 368, 372, 388, 390, 397, 402, 409, 417, 420, 422, 425, 428, 429, 431, 445, 447, 451, 456, 461, 476, 488, 489, 490, 493, 505, 512, 513, 517, 519, 550, 553, 554, 562, 563, 578, 581, 584, 589, 590, 607, 611, 613, 615, 617, 622, 624, 625, 628, 629, 631, 633, 634, 636, 641, 647, 648, 655, 656, 663, 664, 666, 669, 683, 686, 689, 690, 693, 699, 701, 704, 712, 722, 739, 741, 745, 767, 796, 800, 803, 806, 813, 818, 823, 832, 844, 851, 856, 858, 869, 880, 911, 938, 943, 946, 959, 970, 978, 979, 987, 999.

### h2: 141 clean seeds

103, 105, 107, 111, 114, 116, 117, 126, 140, 150, 151, 161, 162, 179, 183, 187, 195, 202, 208, 211, 212, 213, 229, 250, 261, 264, 273, 275, 285, 289, 295, 305, 309, 327, 348, 352, 354, 364, 376, 377, 380, 383, 394, 408, 409, 415, 418, 422, 431, 449, 469, 473, 480, 488, 489, 492, 496, 499, 508, 513, 517, 521, 547, 555, 557, 560, 571, 575, 582, 596, 597, 602, 606, 613, 623, 625, 629, 634, 639, 642, 643, 644, 648, 652, 657, 661, 665, 666, 670, 672, 678, 682, 685, 693, 697, 699, 705, 708, 712, 714, 718, 721, 727, 728, 732, 738, 741, 755, 769, 780, 786, 788, 798, 806, 815, 822, 829, 840, 841, 874, 875, 881, 882, 889, 891, 902, 903, 912, 916, 922, 933, 934, 940, 952, 954, 956, 969, 972, 976, 987, 989.

### h3: 138 clean seeds

107, 110, 111, 112, 114, 126, 138, 145, 151, 158, 161, 167, 168, 177, 181, 190, 205, 207, 231, 234, 244, 247, 252, 256, 257, 262, 274, 277, 278, 280, 288, 291, 300, 316, 329, 341, 353, 373, 387, 389, 391, 406, 409, 410, 412, 422, 433, 434, 439, 442, 445, 454, 467, 470, 476, 488, 491, 492, 493, 499, 510, 529, 544, 547, 549, 551, 569, 578, 581, 585, 588, 596, 599, 601, 614, 621, 622, 634, 642, 643, 656, 657, 660, 665, 668, 677, 686, 690, 691, 697, 698, 702, 721, 726, 730, 737, 739, 744, 746, 753, 756, 764, 770, 783, 784, 785, 787, 792, 793, 794, 795, 798, 805, 809, 821, 828, 836, 839, 840, 846, 850, 852, 864, 866, 881, 885, 906, 915, 923, 930, 946, 961, 966, 980, 985, 988, 994, 995.

### h4: 148 clean seeds

101, 107, 122, 126, 127, 128, 144, 155, 156, 160, 173, 174, 178, 180, 191, 194, 206, 217, 221, 223, 230, 233, 234, 238, 239, 248, 251, 253, 256, 265, 268, 270, 274, 277, 290, 295, 311, 314, 323, 325, 330, 332, 334, 339, 344, 349, 361, 369, 390, 416, 419, 433, 439, 454, 455, 459, 470, 472, 475, 485, 486, 491, 494, 506, 533, 537, 541, 560, 567, 568, 571, 574, 581, 600, 603, 607, 629, 631, 633, 637, 655, 678, 680, 683, 685, 687, 693, 695, 698, 705, 714, 718, 719, 721, 722, 727, 731, 732, 738, 742, 744, 745, 747, 749, 755, 757, 759, 761, 762, 776, 777, 779, 780, 781, 804, 806, 808, 813, 828, 831, 842, 847, 852, 857, 862, 869, 883, 884, 896, 897, 907, 908, 913, 919, 921, 924, 925, 931, 939, 942, 943, 956, 961, 965, 979, 985, 990, 996.

### h5: 148 clean seeds

113, 118, 122, 124, 126, 129, 132, 137, 152, 157, 162, 175, 178, 185, 186, 199, 200, 214, 221, 232, 260, 272, 273, 279, 283, 286, 287, 293, 309, 313, 318, 324, 331, 332, 334, 339, 341, 348, 350, 353, 362, 364, 373, 374, 375, 385, 387, 390, 391, 396, 397, 400, 415, 429, 431, 433, 434, 438, 445, 447, 470, 484, 486, 492, 496, 497, 504, 505, 512, 517, 524, 525, 529, 535, 556, 570, 572, 576, 590, 593, 595, 608, 617, 625, 627, 641, 645, 655, 676, 679, 681, 685, 688, 695, 698, 700, 701, 707, 714, 721, 722, 726, 732, 736, 753, 754, 761, 765, 789, 792, 800, 801, 803, 805, 806, 807, 808, 810, 811, 819, 834, 844, 850, 852, 856, 874, 877, 881, 885, 888, 889, 899, 912, 914, 920, 937, 954, 957, 960, 964, 968, 974, 977, 980, 985, 988, 993, 994.

### h6: 4 clean seeds

175, 712, 981, 982.

## Verification and artifacts

- All seven upload records confirm this task; contracts and reference files match the published task after the scoring-seed stamp.
- All seven submissions were scored through calc.py into this report directory. Each has 250/250 valid rows, cell accessibility 0.87, and a three-seed score matching the official API exactly.
- Stage-3 outcomes and indel lengths agree with the full local validation for all 250 rows at all three real scoring seeds.
- Recovered actual clean sets and HDR bands with **900 × 250 × 7 = 1,575,000 row simulations**.
- [Official scores](scores.json), [published task](task.json), [sanitized logs](build_evidence.json), [archive manifest](manifest.json), [summary CSV](summary.csv), [full analysis](analysis.json).
- Reproduce from the repository root with `.venv/bin/python reports/11276911/validate.py`, then `.venv/bin/python reports/11276911/audit.py`, then `.venv/bin/python reports/11276911/render_report.py`.

Uploaded archives:

- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — [data/inst/niome_hotkey/result/2026-09-18T05:30:43](../../data/inst/niome_hotkey/result/2026-09-18T05:30:43/submission.json); SHA-256 `15882312e47caeedcea36a927822bf45ded98130e6562904573d106e5d362296`.
- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-18T06:22:33](../../data/inst/niome_hotkey1/result/2026-09-18T06:22:33/submission.json); SHA-256 `a045293d4a1be4cc2411a2b114b9fef6ebad2b0a63ec9d5153dbe49e8f193513`.
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-18T06:03:34](../../data/inst/niome_hotkey2/result/2026-09-18T06:03:34/submission.json); SHA-256 `5819537fb2d63874110abc037dc92ddb8232f3931dd44388a27f62ab2afb0f46`.
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-18T06:07:33](../../data/inst/niome_hotkey3/result/2026-09-18T06:07:33/submission.json); SHA-256 `a4af3207260b3a857d67f3b3c4db7efcdffb82b0c86335dc57ed2bf8351c0be2`.
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — [data/inst/niome_hotkey4/result/2026-09-18T06:21:49](../../data/inst/niome_hotkey4/result/2026-09-18T06:21:49/submission.json); SHA-256 `d569bc325de53d7e6c682a71848daec0bbb550ee85f424f00352a626153d1dbc`.
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — [data/inst/niome_hotkey5/result/2026-09-18T05:34:51](../../data/inst/niome_hotkey5/result/2026-09-18T05:34:51/submission.json); SHA-256 `1f8631e93d2fc8022cc73b196f81e896ad0e1403018e87ff9482dbb883cc34e6`.
- **h6** `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` — [data/inst/niome_hotkey6/result/2026-09-18T05:32:22](../../data/inst/niome_hotkey6/result/2026-09-18T05:32:22/submission.json); SHA-256 `15882312e47caeedcea36a927822bf45ded98130e6562904573d106e5d362296`.
