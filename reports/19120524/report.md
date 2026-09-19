# Submission audit: 19120524-b0f1-4061-b99b-92e4599e98d3

**HUDEP-2**, created **2026-09-18T17:30:15.846726 UTC**. Scoring seeds: **236, 596, 156**.

**Best of h0–h9: h8, 41.002839, rank 150/248. Leader: UID 104, 182.340777.** h8 achieved 22.49% of the leader, a gap of 141.337938. Top-10 cutoff: 126.048544; h8 was 85.045705 below it.

**Six uploaded conjunctions were found: h1, h2, h3, h7, h8 and h9. Each has 250 valid rows. h0, h4, h5 and h6 have no matching upload archive or received-task/upload event; each has official score 0 and 0 valid rows.** All ten hotkeys have zero reward weight. Prepared builds for the four missing uploads are documented separately below.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=19120524-b0f1-4061-b99b-92e4599e98d3&limit=40000), fetched 2026-09-18T20:56:12.229347+00:00. 248 records, 248 miners, 1 validator. Leader score timestamp: 2026-09-18T19:51:37.056316 UTC. Audit generated 2026-09-18T20:59:07.160429+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 104 | 1 | 182.340777 | 325.825107 | 0.589330 | 0.949600 | 250 | 0.294000 |
| h0 | 248 | 237 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0 | 0.000000 |
| h1 | 41 | 229 | 29.062040 | 285.710877 | 0.106300 | 0.956898 | 250 | 0.000000 |
| h2 | 163 | 222 | 30.052147 | 298.447972 | 0.107290 | 0.938532 | 250 | 0.000000 |
| h3 | 118 | 233 | 26.438157 | 285.001752 | 0.096822 | 0.958098 | 250 | 0.000000 |
| h4 | 190 | 237 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0 | 0.000000 |
| h5 | 234 | 237 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0 | 0.000000 |
| h6 | 175 | 237 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0 | 0.000000 |
| h7 | 65 | 231 | 27.877264 | 301.555429 | 0.100081 | 0.923699 | 250 | 0.000000 |
| h8 | 211 | 150 | 41.002839 | 301.325952 | 0.145272 | 0.936686 | 250 | 0.000000 |
| h9 | 245 | 228 | 29.500824 | 302.686967 | 0.104686 | 0.931003 | 250 | 0.000000 |

Rank is 1 plus the number of strictly higher scores. All four zero-score hotkeys share rank **237** with other zero-score miners; that rank does not imply a valid upload.

Leader hotkey: `5Eh3xn7ZpqaHqSHbB2uruFmRB83WZLwwysGADcEcJAvapT3o`.

## Original joined windows and candidate slices

The original **16:40 plan** selected windows **100/200/300**, joined as **100–399**. All ten hotkeys prepared against that band space. Seeds **236 and 156** are inside it; **596** is outside. The submitted builds used the original ten-hotkey **stride-30** layout, with 150-seed circular candidate slices. Every submitted conjunction used a **900-seed cut search (100–999)**, **100 Cas12a + 150 Cas9 rows**, and an **11-seed HDR band**.

| Hotkey | Upload evidence | Band candidate slice | Offset | Actual clean /900 |
|---|---|---|---:|---:|
| h0 | prepared only; no upload found | 100–249 | 0 | N/A |
| h1 | uploaded conjunction | 130–279 | 30 | 141 |
| h2 | uploaded conjunction | 160–309 | 60 | 131 |
| h3 | uploaded conjunction | 190–339 | 90 | 139 |
| h4 | prepared only; no upload found | 220–369 | 120 | N/A |
| h5 | prepared only; no upload found | 250–399 | 150 | N/A |
| h6 | prepared only; no upload found | 100–129, 280–399 | 180 | N/A |
| h7 | uploaded conjunction | 100–159, 310–399 | 210 | 140 |
| h8 | uploaded conjunction | 100–189, 340–399 | 240 | 139 |
| h9 | uploaded conjunction | 100–219, 370–399 | 270 | 133 |

Candidate slices are reconstructed from the historical plan's offsets and the logged 100–399 band space. Every actual uploaded band lies inside its corresponding original slice. The later split between predicted and complementary windows was not the layout used for these uploads.

## Actual uploaded clean sets and HDR bands

**Clean set:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. Both were recovered over seeds 100–999. The full-row clean sets equal the Cas12a group clean sets, and their sizes match the original conjunction build logs.

| Hotkey | Clean /900 | Clean inside joined /300 | Clean outside joined | Actual uploaded HDR band seeds | Clean hits | HDR hits |
|---|---:|---:|---:|---|---|---|
| h0 | N/A | N/A | N/A | No uploaded rows available | N/A | N/A |
| h1 | 141 | 50 | 91 | 136, 141, 144, 145, 222, 251, 252, 253, 256, 261, 274 | 236 | none |
| h2 | 131 | 51 | 80 | 195, 223, 249, 251, 256, 260, 261, 271, 281, 282, 300 | 156 | none |
| h3 | 139 | 52 | 87 | 210, 244, 251, 256, 261, 272, 281, 300, 301, 307, 318 | none | none |
| h4 | N/A | N/A | N/A | No uploaded rows available | N/A | N/A |
| h5 | N/A | N/A | N/A | No uploaded rows available | N/A | N/A |
| h6 | N/A | N/A | N/A | No uploaded rows available | N/A | N/A |
| h7 | 140 | 51 | 89 | 137, 159, 312, 320, 332, 343, 352, 354, 364, 385, 398 | none | none |
| h8 | 139 | 47 | 92 | 101, 107, 117, 139, 140, 159, 164, 348, 355, 357, 385 | 156 | none |
| h9 | 133 | 50 | 83 | 112, 116, 159, 192, 201, 204, 210, 218, 385, 398, 399 | none | none |

**h1 hit clean seed 236; h2 and h8 hit clean seed 156. No other uploaded clean hits or HDR-band hits occurred.** Seed 596 was outside all six uploaded clean sets. The six-upload fleet union contains **565/900 clean seeds** and **52 distinct HDR-band seeds** (66 band memberships before overlap). Preparation-only records are excluded from these totals.

## Per-seed results

Each entry is **per-seed final score / no-cut rows**. Zero no-cut rows identifies a clean hit.

| Hotkey | Seed 236 | Seed 596 | Seed 156 | Mean final score |
|---|---:|---:|---:|---:|
| h1 | 30.487978 / 0 | 31.297975 / 5 | 25.400166 / 12 | 29.062040 |
| h2 | 25.398668 / 7 | 18.291036 / 8 | 46.466735 / 0 | 30.052147 |
| h3 | 26.310879 / 8 | 27.058583 / 6 | 25.945010 / 7 | 26.438157 |
| h7 | 28.474569 / 4 | 26.910807 / 3 | 28.246416 / 4 | 27.877264 |
| h8 | 27.862710 / 6 | 34.670339 / 3 | 60.475469 / 0 | 41.002839 |
| h9 | 33.012972 / 3 | 27.584191 / 6 | 27.905309 / 2 | 29.500824 |

## Comparison with the leader

h8's weighted score is **301.325952**, versus the leader's **325.825107**. Fidelity is **0.936686 versus 0.949600**. The largest gap is mean consistency: **0.145272 versus 0.589330**. h8's clean hit at 156 scored **60.475469**, but its other seed scores of **27.862710 at 236** and **34.670339 at 596** reduced the three-seed mean to **41.002839**. The leader's submitted rows, clean set, bands and individual seed outcomes were not available in the inspected data; aggregate scores do not identify those sets.

## Build and upload provenance

The logs show preparation attempts around **17:30**, **17:34**, and **17:47 UTC**. The completed 17:47 builds were used for all six successful uploads, which occurred between **17:50:23 and 18:23:34 UTC**, with 264–289 seconds of URL TTL remaining. No ordinary or all-HDR fallback was used by these six submissions.

At **18:47–18:49 UTC**, those same six hotkeys rebuilt the task after their uploads. The newer configuration changed h1–h3's offsets and moved h7–h9 into complementary band space. Their latest `window_used.json` bands therefore differ from all six actual uploaded bands. The report uses uploaded-row replay to avoid crediting these later bands.

| Hotkey | Actual uploaded band | Later recorded band — not submitted |
|---|---|---|
| h1 | 136, 141, 144, 145, 222, 251, 252, 253, 256, 261, 274 | 158, 198, 214, 249, 251, 253, 256, 261, 266, 281, 297 |
| h2 | 195, 223, 249, 251, 256, 260, 261, 271, 281, 282, 300 | 210, 251, 256, 261, 272, 281, 283, 300, 307, 318, 347 |
| h3 | 210, 244, 251, 256, 261, 272, 281, 300, 301, 307, 318 | 251, 256, 261, 271, 272, 281, 295, 296, 300, 318, 381 |
| h7 | 137, 159, 312, 320, 332, 343, 352, 354, 364, 385, 398 | 551, 567, 576, 577, 601, 602, 627, 629, 631, 664, 674 |
| h8 | 101, 107, 117, 139, 140, 159, 164, 348, 355, 357, 385 | 737, 738, 750, 758, 760, 787, 803, 819, 829, 848, 849 |
| h9 | 112, 116, 159, 192, 201, 204, 210, 218, 385, 398, 399 | 850, 863, 880, 893, 905, 916, 932, 933, 959, 963, 982 |

## h0/h4/h5/h6: preparation only

These four hotkeys have completed prepared builds in their logs, but no matching upload archives or received-task/successful-upload events were found. The API gives each 0 valid experiments and score 0. The reason the validator obtained no valid rows is not established by the available records.

The following are **logged prepared-build values**, not recovered submitted sets:

| Hotkey | Logged prepared clean count /900 | Recorded prepared band |
|---|---:|---|
| h0 | 137 | 110, 119, 128, 145, 159, 161, 176, 182, 235, 245, 249 |
| h4 | 141 | 236, 242, 251, 256, 261, 272, 281, 300, 307, 318, 359 |
| h5 | 136 | 251, 256, 261, 271, 272, 281, 295, 296, 300, 318, 381 |
| h6 | 135 | 110, 113, 124, 125, 282, 314, 320, 335, 339, 373, 380 |

h4's prepared-band record includes scoring seed **236**, but there is no uploaded submission available to validate or credit that potential band hit. Its official score remains zero.

## Exact uploaded clean sets

These lists include each uploaded HDR band and cover the full 100–999 seed space. No submitted clean set can be recovered for h0/h4/h5/h6.

### h1: 141 clean seeds

101, 108, 126, 127, 131, 136, 141, 144, 145, 147, 153, 176, 194, 197, 201, 204, 222, 225, 236, 249, 251, 252, 253, 256, 259, 261, 271, 274, 275, 277, 290, 292, 294, 296, 300, 316, 336, 339, 341, 355, 361, 366, 369, 370, 372, 376, 387, 388, 395, 397, 425, 426, 434, 446, 452, 455, 467, 477, 478, 481, 482, 484, 487, 488, 490, 500, 502, 509, 516, 521, 524, 526, 533, 537, 542, 544, 549, 553, 572, 575, 603, 616, 623, 641, 644, 651, 654, 663, 664, 669, 673, 674, 685, 688, 704, 706, 712, 725, 730, 731, 736, 743, 760, 770, 779, 784, 793, 802, 803, 805, 807, 810, 814, 831, 832, 836, 843, 858, 859, 867, 876, 877, 887, 891, 898, 908, 921, 928, 937, 944, 945, 947, 955, 956, 958, 967, 976, 978, 985, 989, 995.

### h2: 131 clean seeds

105, 106, 111, 112, 116, 117, 127, 137, 138, 142, 154, 156, 166, 188, 189, 195, 201, 206, 214, 223, 242, 249, 251, 256, 257, 260, 261, 268, 270, 271, 277, 281, 282, 290, 297, 300, 308, 318, 319, 324, 336, 346, 358, 364, 368, 379, 381, 382, 389, 395, 398, 400, 405, 412, 419, 426, 431, 448, 454, 464, 475, 476, 486, 487, 491, 504, 508, 514, 517, 519, 520, 531, 542, 549, 553, 561, 568, 577, 589, 590, 606, 618, 628, 629, 641, 645, 646, 649, 654, 665, 669, 677, 685, 692, 697, 707, 723, 732, 745, 751, 757, 758, 783, 786, 790, 791, 797, 826, 832, 834, 839, 862, 878, 879, 888, 892, 893, 910, 912, 917, 919, 922, 924, 932, 933, 944, 955, 958, 962, 969, 991.

### h3: 139 clean seeds

103, 105, 117, 125, 138, 149, 155, 158, 163, 165, 172, 173, 178, 179, 180, 186, 195, 197, 202, 205, 210, 215, 219, 221, 224, 229, 239, 244, 249, 251, 256, 261, 262, 272, 281, 300, 301, 307, 316, 318, 324, 335, 339, 340, 364, 367, 373, 376, 377, 381, 384, 385, 402, 404, 406, 408, 419, 422, 423, 426, 430, 431, 435, 446, 449, 452, 474, 475, 481, 485, 497, 501, 503, 506, 518, 531, 537, 545, 546, 549, 566, 567, 571, 576, 579, 584, 589, 597, 612, 623, 625, 636, 653, 657, 658, 664, 666, 667, 681, 684, 702, 703, 705, 709, 714, 715, 747, 750, 754, 766, 777, 796, 800, 802, 807, 809, 822, 824, 834, 842, 850, 851, 856, 863, 875, 884, 894, 910, 920, 927, 934, 935, 942, 944, 948, 960, 970, 989, 994.

### h7: 140 clean seeds

137, 139, 144, 146, 151, 152, 154, 159, 166, 169, 175, 176, 196, 199, 213, 218, 227, 231, 237, 238, 246, 250, 253, 258, 281, 287, 305, 311, 312, 313, 320, 328, 329, 331, 332, 333, 336, 343, 345, 349, 352, 354, 364, 365, 369, 377, 380, 384, 385, 386, 398, 403, 406, 407, 420, 432, 451, 454, 455, 460, 466, 476, 480, 485, 490, 497, 504, 512, 524, 525, 530, 549, 566, 568, 572, 575, 580, 586, 591, 593, 611, 622, 625, 627, 634, 655, 656, 657, 665, 667, 671, 696, 698, 705, 710, 725, 733, 736, 740, 743, 747, 749, 751, 754, 755, 758, 769, 775, 797, 804, 808, 811, 812, 816, 817, 819, 829, 833, 834, 836, 839, 854, 877, 882, 885, 888, 894, 895, 897, 898, 899, 917, 929, 949, 954, 955, 972, 977, 982, 994.

### h8: 139 clean seeds

101, 103, 105, 107, 117, 124, 127, 139, 140, 142, 152, 156, 159, 164, 172, 183, 187, 192, 196, 203, 211, 215, 218, 222, 224, 238, 255, 291, 298, 301, 303, 315, 321, 325, 326, 332, 334, 340, 345, 348, 355, 356, 357, 365, 377, 385, 395, 404, 412, 413, 418, 424, 436, 438, 441, 449, 452, 455, 462, 468, 479, 491, 499, 501, 538, 541, 545, 560, 568, 574, 590, 599, 600, 602, 620, 621, 624, 630, 631, 632, 633, 635, 640, 648, 649, 651, 652, 667, 686, 695, 699, 714, 725, 746, 763, 771, 775, 778, 783, 785, 795, 802, 815, 818, 820, 822, 823, 827, 832, 848, 849, 860, 862, 873, 874, 887, 889, 892, 895, 899, 904, 907, 915, 926, 927, 929, 938, 942, 943, 944, 945, 949, 955, 965, 978, 991, 993, 994, 997.

### h9: 133 clean seeds

104, 105, 107, 110, 112, 114, 116, 117, 129, 131, 133, 139, 146, 147, 159, 184, 192, 199, 201, 204, 210, 214, 218, 235, 255, 256, 263, 267, 270, 278, 289, 307, 317, 322, 324, 326, 328, 333, 334, 342, 346, 349, 368, 370, 382, 385, 389, 394, 398, 399, 412, 423, 432, 443, 447, 449, 451, 457, 471, 476, 481, 486, 506, 511, 517, 519, 522, 524, 534, 543, 566, 568, 575, 578, 590, 594, 595, 624, 626, 634, 648, 649, 656, 669, 672, 678, 685, 686, 694, 695, 696, 701, 705, 711, 726, 752, 764, 765, 783, 788, 791, 795, 797, 810, 814, 816, 826, 831, 839, 849, 860, 863, 864, 870, 891, 898, 903, 908, 909, 911, 912, 918, 921, 922, 929, 938, 941, 946, 948, 959, 980, 982, 984.

## Verification and artifacts

- All six upload archives identify this task, with contracts and references matching the published task after the seed stamp.
- All six submissions have 250 stage-1/2 valid rows and accessibility 0.82. Archived local three-seed scores match the API exactly.
- Replayed outcomes and indel lengths match archived stage-3 detail for every uploaded row at all three scoring seeds.
- Recovered clean sets and bands over **900 × 250 × 6 = 1,350,000 row/seed combinations**.
- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized build/upload logs](build_evidence.json), [historical layout](plan_history.json), [archive manifest](manifest.json), [summary CSV](summary.csv).
- Reproduce with `.venv/bin/python reports/19120524/audit.py` and `.venv/bin/python reports/19120524/render_report.py` from the repository root, using the archived validations.

Uploaded archives:

- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-18T18:05:29](../../data/inst/niome_hotkey1/result/2026-09-18T18:05:29/submission.json); SHA-256 `88d4a8740bf23bb58ff50b4d95efb6e450d30a1c4f33c1dfaee113901f2e6cb3`.
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-18T17:59:01](../../data/inst/niome_hotkey2/result/2026-09-18T17:59:01/submission.json); SHA-256 `82cefac71ee5f53fa27f52315f97f98cdb4372ce15cb366f14b793455edbe41b`.
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-18T17:50:52](../../data/inst/niome_hotkey3/result/2026-09-18T17:50:52/submission.json); SHA-256 `0b06855f005ef0fee1e2007393a7fb303ef99a148349abb5f6c1ba25689b0b6a`.
- **h7** `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` — [data/inst/niome_hotkey7/result/2026-09-18T18:24:08](../../data/inst/niome_hotkey7/result/2026-09-18T18:24:08/submission.json); SHA-256 `48cbcda8dd27fd1ed3f4e28cbd302333530b425a858f6b89b6277f9e02fb4a59`.
- **h8** `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` — [data/inst/niome_hotkey8/result/2026-09-18T18:20:57](../../data/inst/niome_hotkey8/result/2026-09-18T18:20:57/submission.json); SHA-256 `c0dab369bc34150e0f1963b28a395c8748da685a6c661dfc12b51de2c6b0b75d`.
- **h9** `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` — [data/inst/niome_hotkey9/result/2026-09-18T18:02:35](../../data/inst/niome_hotkey9/result/2026-09-18T18:02:35/submission.json); SHA-256 `5b57f466b586b87cab1ea220a79bae7b655c9749d631dc9d6bfe6269b0f5538c`.
- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — no upload archive for this task.
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — no upload archive for this task.
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — no upload archive for this task.
- **h6** `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` — no upload archive for this task.
