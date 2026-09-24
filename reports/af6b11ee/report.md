# Submission audit: af6b11ee-dbc2-45b5-86cb-9ce888dba0c3

**HUDEP-2**, created **2026-09-20T03:28:47.355981 UTC**. Scoring seeds: **353, 218, 143**.

**Best of h0–h9: h8, 93.099061, rank 14/248. Leader: UID 36, 139.646088.** h8 achieved 66.67% of the leader, a gap of 46.547027. Top-10 cutoff: 124.583520; h8 was 31.484459 below it.

**Seven uploads were found: h0 and h4–h9.** Each has 250 valid rows, **100 Cas12a + 150 Cas9**, and an **11-seed HDR band**. All seven have zero reward weight in the official snapshot. **h1–h3 have no task upload archive, received-task/upload event or official score entry.** Their results are unavailable, not zero. Prepared builds for these hotkeys are documented separately.

**h8 hit HDR seed 353: all 250 uploaded rows returned HDR**, scoring **235.918078** at that seed. Its scores at 218 and 143 were **21.320120** and **22.058984**, giving the task mean **93.099061**. h7 hit clean seed **218**, scoring **33.948359** there. No other uploaded clean or HDR hits occurred; 143 was outside every uploaded clean set.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=af6b11ee-dbc2-45b5-86cb-9ce888dba0c3&limit=40000), fetched 2026-09-20T06:02:24.276363+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-20T05:51:10.803709 UTC. Audit generated 2026-09-20T06:04:54.488871+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 36 | 1 | 139.646088 | 248.407129 | 0.605652 | 0.928200 | 250 | 0.294000 |
| h0 | 248 | 235 | 22.965340 | 259.207721 | 0.097626 | 0.907530 | 250 | 0.000000 |
| h1 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h2 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h3 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h4 | 190 | 234 | 23.038672 | 255.380222 | 0.098653 | 0.914445 | 250 | 0.000000 |
| h5 | 234 | 232 | 24.572203 | 257.711746 | 0.104533 | 0.912128 | 250 | 0.000000 |
| h6 | 175 | 233 | 23.548017 | 257.415399 | 0.100775 | 0.907753 | 250 | 0.000000 |
| h7 | 65 | 190 | 27.542749 | 258.873208 | 0.117653 | 0.904309 | 250 | 0.000000 |
| h8 | 211 | 14 | 93.099061 | 259.870737 | 0.394625 | 0.907829 | 250 | 0.000000 |
| h9 | 245 | 223 | 26.062625 | 256.611111 | 0.111379 | 0.911881 | 250 | 0.000000 |

Rank is 1 plus the number of strictly higher scores; ties share rank. Final score and breakdown factors are means over three seeds. Multiplying mean factors is not generally equivalent to averaging per-seed products. Missing hotkeys have no rank or score; they are not assigned a zero-score rank.

Leader hotkey: `5EF3SzQyf2YbWGbhsqddQ1un7uRtzmXJu3fBR16RhsfohX4b`.

## Historical joined windows and candidate slices

The **03:17 UTC plan** selected width-100 windows **100 / 200 / 300**, joined as **100–399**. Actual windows were **300 / 200 / 100** for seeds **353 / 218 / 143**: **all three matched**.

The task's original build logs show **ten disjoint 30-seed candidate slices**, at offsets 0, 30, …, 270 within the predicted 300 seeds. Every prepared conjunction used **k=11, group=100** and a full **100–999 cut search**. All seven uploaded bands lie inside their corresponding slices. This report uses the task's logged layout; the repository's later seven-hotkey width-300 layout was not used for these submissions.

| Hotkey | Upload status | Assigned 30-seed candidate slice | Offset | Actual clean /900 |
|---|---|---|---:|---:|
| h0 | uploaded conjunction | 100–129 | 0 | 124 |
| h1 | prepared only; no upload found | 130–159 | 30 | N/A |
| h2 | prepared only; no upload found | 160–189 | 60 | N/A |
| h3 | prepared only; no upload found | 190–219 | 90 | N/A |
| h4 | uploaded conjunction | 220–249 | 120 | 124 |
| h5 | uploaded conjunction | 250–279 | 150 | 131 |
| h6 | uploaded conjunction | 280–309 | 180 | 122 |
| h7 | uploaded conjunction | 310–339 | 210 | 133 |
| h8 | uploaded conjunction | 340–369 | 240 | 123 |
| h9 | uploaded conjunction | 370–399 | 270 | 136 |

The seven uploaded candidate slices cover **100–129 ∪ 220–399** (210 seeds). The assigned slices **130–219** belong to h1–h3, for which no upload exists. Scoring seeds **143 and 218** fell into those unsubmitted slices; **353** fell into h8's uploaded slice. Clean sets can extend beyond a candidate slice because the cut search spans all 900 seeds.

## Actual uploaded clean sets and HDR bands

**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. Both sets were recovered over 100–999. For all seven uploads, the full-row clean set equals the Cas12a group clean set, its count matches the build log, and the HDR band matches the window record.

| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |
|---|---:|---|---|---|
| h0 | 124 | 101, 104, 106, 109, 113, 114, 115, 119, 124, 125, 126 | none | none |
| h1 | N/A | No uploaded rows available | N/A | N/A |
| h2 | N/A | No uploaded rows available | N/A | N/A |
| h3 | N/A | No uploaded rows available | N/A | N/A |
| h4 | 124 | 221, 225, 227, 229, 231, 239, 240, 242, 244, 246, 247 | none | none |
| h5 | 131 | 255, 257, 259, 263, 264, 268, 270, 273, 276, 277, 278 | none | none |
| h6 | 122 | 280, 285, 287, 288, 290, 298, 304, 305, 306, 308, 309 | none | none |
| h7 | 133 | 310, 311, 312, 318, 319, 321, 322, 328, 331, 332, 333 | 218 | none |
| h8 | 123 | 340, 341, 344, 349, 352, 353, 354, 362, 364, 367, 368 | 353 | 353 |
| h9 | 136 | 376, 377, 380, 384, 386, 387, 388, 389, 391, 392, 398 | none | none |

The uploaded fleet union contains **579/900 clean seeds** and **77 distinct HDR-band seeds**. Its seven bands are disjoint (7 × 11 = 77 seeds). All totals exclude h1–h3's prepared builds.

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.

| Hotkey | Seed 353 | Seed 218 | Seed 143 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 18.035317 / 126 / 10 | 27.876161 / 135 / 3 | 22.984542 / 127 / 5 | 22.965340 |
| h4 | 24.878660 / 123 / 4 | 23.828457 / 135 / 3 | 20.408899 / 131 / 10 | 23.038672 |
| h5 | 25.953162 / 139 / 7 | 24.093441 / 132 / 2 | 23.670006 / 118 / 2 | 24.572203 |
| h6 | 26.772907 / 120 / 5 | 20.810337 / 121 / 8 | 23.060806 / 122 / 3 | 23.548017 |
| h7 | 28.174534 / 108 / 3 | 33.948359 / 120 / 0 | 20.505355 / 129 / 3 | 27.542749 |
| h8 | 235.918078 / 250 / 0 | 21.320120 / 129 / 6 | 22.058984 / 138 / 6 | 93.099061 |
| h9 | 23.972188 / 127 / 5 | 25.090784 / 129 / 2 | 29.124904 / 130 / 5 | 26.062625 |

## Comparison with the leader

h8's weighted score was **259.870737**, versus the leader's **248.407129**. Mean consistency was **0.394625 versus 0.605652**, and fidelity was **0.907829 versus 0.928200**. The largest factor gap is consistency. h8's strong score at its HDR hit was averaged with two low seed scores. The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the inspected data; its aggregate score does not identify those sets.

## h1–h3: prepared builds without uploads

Runtime logs reported these three hotkeys unregistered on netuid 55 on **2026-09-19**, before this task. Their processes still prepared datasets for this task, but no received-task/upload event or matching upload archive was found, and none of their public hotkeys appears in this task's complete 248-record score feed.

The following counts and bands are **prepared-build records, not reconstructed submitted sets**:

| Hotkey | Unregistered log time UTC | Logged prepared clean /900 | Recorded prepared band | Scoring seeds in prepared band |
|---|---|---:|---|---|
| h1 | 2026-09-19 21:11:45,103 | 126 | 130, 133, 135, 139, 142, 149, 153, 154, 155, 156, 159 | none |
| h2 | 2026-09-19 21:39:22,213 | 132 | 160, 161, 167, 168, 171, 174, 176, 182, 187, 188, 189 | none |
| h3 | 2026-09-19 22:02:35,440 | 138 | 193, 195, 196, 197, 202, 209, 210, 212, 213, 218, 219 | 218 |

**h3's prepared band includes 218**, but it was not uploaded. It cannot be counted as an actual HDR hit or assigned a hypothetical official score. Exact submitted clean sets cannot be recovered for h1–h3 because there are no submitted rows to replay.

## Exact uploaded clean sets

Each list includes its HDR band and covers the entire 100–999 seed space.

### h0: 124 clean seeds

101, 104, 106, 109, 113, 114, 115, 119, 122, 124, 125, 126, 129, 132, 139, 148, 154, 156, 161, 166, 185, 192, 206, 208, 223, 228, 235, 237, 274, 276, 277, 285, 288, 296, 298, 303, 305, 306, 311, 313, 317, 319, 322, 323, 326, 328, 333, 343, 374, 381, 382, 430, 452, 453, 458, 474, 482, 483, 484, 485, 492, 507, 510, 516, 518, 538, 540, 545, 548, 562, 565, 567, 586, 596, 599, 614, 616, 634, 652, 671, 676, 682, 692, 693, 701, 708, 715, 718, 732, 738, 743, 746, 747, 749, 752, 753, 755, 764, 774, 789, 798, 799, 809, 825, 830, 835, 843, 850, 860, 861, 864, 866, 868, 877, 878, 889, 901, 902, 917, 939, 941, 973, 990, 997.

### h4: 124 clean seeds

100, 102, 111, 118, 122, 124, 144, 159, 168, 169, 179, 183, 184, 189, 207, 211, 212, 215, 216, 221, 225, 227, 229, 231, 233, 239, 240, 242, 244, 246, 247, 248, 252, 263, 273, 279, 286, 306, 309, 315, 319, 321, 325, 337, 340, 346, 347, 354, 369, 374, 390, 418, 437, 440, 442, 452, 474, 490, 491, 494, 516, 517, 519, 531, 533, 536, 540, 541, 553, 555, 557, 568, 574, 599, 608, 611, 612, 620, 633, 634, 635, 641, 644, 653, 669, 670, 679, 691, 717, 727, 728, 732, 746, 757, 760, 771, 772, 774, 788, 796, 799, 804, 810, 812, 818, 829, 833, 839, 844, 881, 886, 895, 913, 917, 923, 926, 930, 934, 942, 950, 963, 986, 990, 993.

### h5: 131 clean seeds

106, 128, 133, 150, 157, 160, 164, 173, 193, 196, 201, 203, 206, 208, 220, 227, 232, 238, 247, 252, 255, 257, 259, 263, 264, 268, 269, 270, 271, 273, 276, 277, 278, 284, 288, 293, 301, 312, 328, 330, 336, 341, 343, 351, 352, 355, 359, 362, 363, 365, 383, 392, 441, 464, 480, 488, 494, 500, 503, 507, 509, 513, 516, 521, 534, 537, 545, 553, 558, 562, 566, 569, 571, 573, 585, 586, 587, 591, 592, 600, 609, 610, 612, 616, 619, 624, 628, 631, 646, 653, 665, 669, 671, 677, 686, 707, 721, 725, 740, 747, 785, 786, 789, 791, 793, 795, 796, 818, 820, 821, 824, 856, 867, 869, 877, 889, 895, 899, 906, 913, 921, 937, 941, 944, 949, 962, 979, 984, 989, 994, 998.

### h6: 122 clean seeds

105, 106, 109, 116, 135, 139, 147, 170, 185, 189, 204, 207, 209, 212, 215, 226, 231, 237, 246, 266, 280, 285, 287, 288, 290, 292, 298, 304, 305, 306, 308, 309, 313, 315, 331, 337, 339, 348, 351, 362, 368, 369, 370, 391, 392, 404, 409, 410, 414, 417, 451, 455, 456, 458, 463, 471, 473, 482, 486, 487, 511, 517, 524, 527, 529, 533, 534, 535, 543, 566, 589, 604, 607, 616, 625, 628, 630, 632, 635, 636, 657, 671, 673, 683, 687, 690, 695, 704, 706, 715, 719, 722, 725, 728, 737, 793, 808, 820, 821, 822, 823, 831, 835, 841, 842, 844, 850, 854, 871, 881, 883, 884, 889, 892, 894, 930, 932, 935, 937, 949, 962, 980.

### h7: 133 clean seeds

109, 111, 119, 122, 128, 129, 130, 138, 142, 144, 148, 154, 157, 159, 162, 166, 172, 206, 208, 211, 212, 218, 220, 233, 239, 240, 241, 242, 244, 283, 288, 292, 295, 301, 308, 310, 311, 312, 314, 318, 319, 321, 322, 328, 331, 332, 333, 334, 339, 351, 354, 372, 382, 384, 386, 387, 407, 408, 410, 412, 414, 447, 449, 453, 476, 486, 491, 495, 505, 513, 522, 540, 541, 550, 552, 561, 566, 571, 582, 588, 596, 616, 628, 630, 644, 648, 653, 658, 662, 663, 664, 682, 685, 688, 696, 713, 743, 749, 752, 754, 760, 763, 776, 783, 801, 802, 811, 813, 829, 833, 840, 842, 854, 855, 868, 869, 872, 875, 884, 895, 898, 910, 926, 933, 937, 948, 960, 965, 966, 969, 978, 979, 998.

### h8: 123 clean seeds

103, 109, 110, 113, 119, 123, 134, 171, 181, 188, 189, 197, 220, 223, 232, 237, 253, 261, 273, 274, 275, 280, 285, 286, 290, 291, 295, 298, 302, 312, 333, 335, 340, 341, 344, 349, 352, 353, 354, 356, 362, 364, 367, 368, 369, 377, 392, 399, 406, 440, 451, 455, 464, 466, 471, 479, 483, 490, 491, 493, 500, 518, 540, 547, 558, 576, 588, 589, 590, 591, 592, 595, 601, 602, 611, 623, 636, 640, 645, 648, 652, 670, 680, 687, 689, 692, 693, 705, 707, 732, 743, 752, 769, 770, 771, 772, 780, 787, 789, 792, 794, 800, 802, 807, 818, 820, 836, 839, 843, 900, 903, 907, 909, 918, 925, 928, 929, 939, 940, 957, 989, 994, 998.

### h9: 136 clean seeds

100, 104, 110, 114, 130, 131, 136, 139, 140, 169, 170, 171, 195, 197, 199, 200, 207, 215, 221, 229, 233, 240, 250, 259, 261, 265, 270, 277, 285, 287, 294, 295, 303, 312, 319, 324, 327, 332, 337, 339, 355, 359, 376, 377, 380, 381, 384, 386, 387, 388, 389, 391, 392, 394, 396, 397, 398, 401, 410, 420, 425, 426, 434, 441, 442, 455, 457, 477, 479, 483, 487, 489, 492, 507, 522, 531, 532, 533, 539, 542, 570, 575, 581, 591, 592, 593, 598, 601, 618, 620, 634, 635, 644, 660, 662, 664, 668, 677, 711, 719, 726, 727, 738, 747, 753, 754, 795, 805, 819, 820, 838, 842, 843, 845, 848, 851, 857, 859, 862, 875, 881, 884, 904, 905, 913, 915, 937, 941, 946, 956, 962, 967, 981, 982, 985, 988.

## Build provenance and verification

- All seven successful uploads used prepared conjunctions. No ordinary fallback appears in their task logs. h0 and h5 waited for their ongoing builds, then uploaded with 199 and 117 seconds of URL TTL remaining.
- All seven archives identify this task; contracts and references match the published task after the seed stamp.
- All 1,750 submitted rows passed stage 1/2, with accessibility 0.82. Local three-seed scores exactly match all seven official scores.
- Reused h4's archived validation and generated the other six in report-local directories.
- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.
- Recovered sets over **900 × 250 × 7 = 1,575,000 row/seed combinations**.
- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized build/upload and registration evidence](build_evidence.json), [manifest](manifest.json), [CSV](summary.csv).
- Reproduce with `.venv/bin/python reports/af6b11ee/validate.py`, then `.venv/bin/python reports/af6b11ee/audit.py`, then `.venv/bin/python reports/af6b11ee/render_report.py`. The analysis includes validator-code and submission SHA-256 hashes.

| Hotkey | Public hotkey | Upload time UTC | Submission archive |
|---|---|---|---|
| h0 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | 2026-09-20 03:36:28,101 | [data/inst/niome_hotkey/result/2026-09-20T03:37:00](../../data/inst/niome_hotkey/result/2026-09-20T03:37:00/submission.json) |
| h1 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | N/A | No upload archive |
| h2 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | N/A | No upload archive |
| h3 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | N/A | No upload archive |
| h4 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | 2026-09-20 04:12:39,154 | [data/inst/niome_hotkey4/result/2026-09-20T04:13:07](../../data/inst/niome_hotkey4/result/2026-09-20T04:13:07/submission.json) |
| h5 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | 2026-09-20 03:36:32,757 | [data/inst/niome_hotkey5/result/2026-09-20T03:37:01](../../data/inst/niome_hotkey5/result/2026-09-20T03:37:01/submission.json) |
| h6 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | 2026-09-20 03:47:16,843 | [data/inst/niome_hotkey6/result/2026-09-20T03:47:45](../../data/inst/niome_hotkey6/result/2026-09-20T03:47:45/submission.json) |
| h7 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | 2026-09-20 04:14:27,859 | [data/inst/niome_hotkey7/result/2026-09-20T04:14:56](../../data/inst/niome_hotkey7/result/2026-09-20T04:14:56/submission.json) |
| h8 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | 2026-09-20 04:08:30,907 | [data/inst/niome_hotkey8/result/2026-09-20T04:08:59](../../data/inst/niome_hotkey8/result/2026-09-20T04:08:59/submission.json) |
| h9 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | 2026-09-20 03:48:13,428 | [data/inst/niome_hotkey9/result/2026-09-20T03:48:41](../../data/inst/niome_hotkey9/result/2026-09-20T03:48:41/submission.json) |
