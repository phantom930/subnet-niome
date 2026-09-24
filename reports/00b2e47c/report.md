# Submission audit: 00b2e47c-f674-4f13-939d-cb4ad624eacd

**CD34+_HSPC**, created **2026-09-20T20:29:10.506241 UTC**. Scoring seeds: **143, 806, 625**.

**Best of h0–h9: h9, 78.159178, rank 13/248. Leader: UID 41, 151.602520.** h9 achieved 51.56% of the leader, a gap of 73.443342. Top-10 cutoff: 117.295730; h9 was 39.136553 below it.

**Seven uploads were found: h0 and h4–h9.** Each has **250 valid rows, 100 Cas12a + 150 Cas9**, and an **11-seed HDR band**. All seven used prepared conjunctions; no fallback was used. All seven have zero reward weight in this snapshot. h1–h3 have no task build/window/upload records or official score entries; their results are unavailable, not zero.

**h9 hit HDR seed 806: all 250 rows returned HDR**, scoring **197.660227** at that seed. Its other seed scores were **17.062024** at 143 and **19.755283** at 625, giving the mean **78.159178**. Additional clean hits occurred at **143 on h0/h4** and **806 on h7**. Seed **625** was outside every uploaded clean set. No other clean or HDR hits occurred.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=00b2e47c-f674-4f13-939d-cb4ad624eacd&limit=40000), fetched 2026-09-20T23:24:52.078407+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-20T22:51:43.334761 UTC. Audit generated 2026-09-20T23:26:36.611071+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 41 | 1 | 151.602520 | 211.410399 | 0.742417 | 0.965900 | 250 | 0.294000 |
| h0 | 248 | 177 | 24.516270 | 212.060902 | 0.122512 | 0.943660 | 250 | 0.000000 |
| h1 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h2 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h3 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h4 | 190 | 171 | 26.326320 | 207.657460 | 0.133292 | 0.951125 | 250 | 0.000000 |
| h5 | 234 | 230 | 19.948968 | 208.922720 | 0.101818 | 0.937802 | 250 | 0.000000 |
| h6 | 175 | 229 | 20.706760 | 210.906708 | 0.104357 | 0.940803 | 250 | 0.000000 |
| h7 | 65 | 199 | 22.928872 | 209.216241 | 0.116060 | 0.944286 | 250 | 0.000000 |
| h8 | 211 | 231 | 19.559221 | 210.857601 | 0.099216 | 0.934931 | 250 | 0.000000 |
| h9 | 245 | 13 | 78.159178 | 205.630610 | 0.395422 | 0.961239 | 250 | 0.000000 |

Rank is 1 plus the number of strictly higher scores; ties share rank. Final score and breakdown factors are means over three seeds. Multiplying mean factors is not generally equivalent to averaging per-seed products. Missing hotkeys have no score or rank.

**UID 41 is the leader's current UID, not our h1's identity.** Leader hotkey: `5DCxsxCbSALgrYDfRXyJ1vCFZbCBYNGsUvyvxLWskFJpNed2`. Our h1 hotkey: `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb`. All identity matching in this audit uses full public hotkeys.

## Joined windows and actual band placement

The **20:17 UTC plan** recorded width-100 windows **100 / 200 / 300**, joined as **100–399**. Actual windows were **100 / 800 / 600** for scoring seeds **143 / 806 / 625**: one matched that prediction.

For these conjunctions, the active band space was the full **100–999**. Seven **300-seed contiguous candidate windows at stride 100** were assigned to h0, h4, h5, h6, h7, h8 and h9. The plan's predicted classes did not steer this placement. The build log label `predicted space of 900` denotes the full space, not the plan's 300 predicted seeds.

All seven **CD34+_HSPC** builds used **cut 900 seeds (wide), k=11, group=100**, confirmed by the task logs. Band candidates span 300 seeds per hotkey; cut selection spans 100–999. Clean sets therefore extend beyond the candidate windows.

| Hotkey | Assigned 300-seed band candidates | Offset in 100–999 | Cut search | Clean /900 | Clean inside candidate window | Clean outside candidate window |
|---|---|---:|---|---:|---:|---:|
| h0 | 100–399 | 0 | 100–999 | 140 | 58 | 82 |
| h1 | No task build or upload | N/A | N/A | N/A | N/A | N/A |
| h2 | No task build or upload | N/A | N/A | N/A | N/A | N/A |
| h3 | No task build or upload | N/A | N/A | N/A | N/A | N/A |
| h4 | 200–499 | 100 | 100–999 | 149 | 72 | 77 |
| h5 | 300–599 | 200 | 100–999 | 149 | 54 | 95 |
| h6 | 400–699 | 300 | 100–999 | 141 | 53 | 88 |
| h7 | 500–799 | 400 | 100–999 | 150 | 57 | 93 |
| h8 | 600–899 | 500 | 100–999 | 154 | 56 | 98 |
| h9 | 700–999 | 600 | 100–999 | 159 | 71 | 88 |

The seven overlapping candidate windows together cover all 100–999. Seed 143 was in h0's candidate window; 806 was in h8/h9's; 625 was in h6/h7/h8's. These are candidate memberships, not guarantees of clean or HDR hits.

## Actual uploaded clean sets and HDR bands

**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. Both were recovered by replay over all 900 seeds (100–999). Each full-row clean set equals the Cas12a group clean set and matches the build log count. Every recovered band matches its task window record and lies within its candidate window.

| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |
|---|---:|---|---|---|
| h0 | 140 | 101, 124, 126, 151, 207, 238, 241, 279, 350, 376, 381 | 143 | none |
| h1 | N/A | No uploaded rows available | N/A | N/A |
| h2 | N/A | No uploaded rows available | N/A | N/A |
| h3 | N/A | No uploaded rows available | N/A | N/A |
| h4 | 149 | 238, 284, 299, 333, 366, 381, 419, 437, 438, 444, 461 | 143 | none |
| h5 | 149 | 371, 381, 387, 399, 419, 447, 537, 543, 551, 564, 598 | none | none |
| h6 | 141 | 400, 421, 454, 461, 524, 553, 600, 626, 629, 631, 664 | none | none |
| h7 | 150 | 642, 669, 691, 697, 716, 718, 739, 764, 775, 780, 785 | 806 | none |
| h8 | 154 | 642, 645, 651, 668, 691, 709, 718, 780, 785, 837, 899 | none | none |
| h9 | 159 | 704, 706, 716, 771, 788, 806, 823, 838, 862, 931, 946 | 806 | 806 |

The uploaded fleet union contains **622/900 clean seeds** and **66 distinct HDR-band seeds** (77 band memberships before overlaps). h9 at 806 was the only uploaded all-HDR hit.

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.

| Hotkey | Seed 143 | Seed 806 | Seed 625 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 34.143661 / 139 / 0 | 20.118494 / 128 / 5 | 19.286655 / 127 / 5 | 24.516270 |
| h4 | 40.887142 / 131 / 0 | 20.676674 / 127 / 5 | 17.415144 / 133 / 4 | 26.326320 |
| h5 | 19.762264 / 134 / 5 | 21.191422 / 127 / 2 | 18.893219 / 133 / 5 | 19.948968 |
| h6 | 20.944284 / 120 / 2 | 19.804510 / 130 / 5 | 21.371485 / 135 / 2 | 20.706760 |
| h7 | 18.908465 / 136 / 6 | 30.173642 / 131 / 0 | 19.704509 / 134 / 3 | 22.928872 |
| h8 | 14.805593 / 140 / 9 | 21.520469 / 124 / 4 | 22.351601 / 128 / 3 | 19.559221 |
| h9 | 17.062024 / 125 / 9 | 197.660227 / 250 / 0 | 19.755283 / 146 / 5 | 78.159178 |

## Comparison with the leader

h9's weighted score was **205.630610**, versus the leader's **211.410399**. Mean consistency was **0.395422 versus 0.742417**, and fidelity was **0.961239 versus 0.965900**. The largest factor gap is consistency. h9's high score at its HDR hit was averaged with two low seed scores. The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the inspected data; its aggregate score does not identify those sets.

## h1–h3: no task activity

No task preparation, received-task/upload event, window record or upload archive was found for h1–h3. None appears in the complete official score feed. Their earlier runtime logs reported them unregistered on netuid 55 on 2026-09-19. No submitted clean set, band, score or rank is available for these hotkeys on this task.

## Exact uploaded clean sets

Each list covers all 100–999 seeds and includes its HDR band.

### h0: 140 clean seeds

101, 104, 107, 112, 119, 121, 124, 126, 127, 134, 143, 144, 151, 153, 155, 156, 158, 159, 169, 170, 180, 190, 191, 193, 198, 207, 213, 217, 218, 224, 227, 233, 234, 238, 241, 242, 248, 268, 269, 276, 279, 293, 304, 305, 309, 317, 320, 324, 332, 341, 350, 355, 376, 380, 381, 383, 387, 398, 401, 403, 419, 429, 431, 434, 444, 449, 457, 460, 474, 476, 479, 489, 491, 494, 496, 499, 501, 520, 535, 556, 557, 561, 562, 574, 583, 593, 611, 622, 624, 632, 649, 652, 669, 674, 680, 687, 697, 707, 724, 736, 738, 743, 752, 754, 755, 769, 772, 774, 792, 799, 811, 812, 817, 821, 824, 827, 829, 830, 832, 838, 845, 848, 857, 872, 883, 887, 897, 904, 911, 915, 920, 925, 943, 946, 948, 964, 968, 988, 992, 994.

### h4: 149 clean seeds

100, 116, 118, 125, 131, 132, 142, 143, 145, 160, 164, 172, 189, 194, 204, 216, 221, 222, 224, 228, 237, 238, 239, 243, 249, 250, 252, 255, 259, 260, 262, 264, 268, 277, 284, 286, 293, 296, 299, 304, 305, 312, 316, 322, 324, 326, 333, 338, 339, 340, 341, 350, 351, 354, 360, 364, 366, 371, 375, 379, 381, 387, 397, 398, 405, 413, 419, 424, 425, 426, 437, 438, 440, 444, 446, 448, 449, 460, 461, 463, 470, 474, 481, 484, 488, 492, 501, 504, 512, 522, 523, 540, 551, 556, 561, 576, 580, 584, 592, 608, 619, 637, 644, 645, 664, 667, 689, 705, 714, 718, 741, 742, 745, 746, 749, 755, 766, 778, 780, 792, 799, 801, 802, 804, 808, 819, 820, 833, 839, 848, 849, 877, 879, 893, 897, 902, 910, 912, 920, 939, 964, 972, 974, 980, 987, 988, 993, 996, 997.

### h5: 149 clean seeds

113, 125, 148, 150, 151, 169, 171, 176, 183, 188, 197, 203, 213, 217, 221, 226, 229, 238, 245, 249, 252, 267, 268, 277, 278, 280, 284, 289, 299, 308, 343, 345, 347, 349, 351, 352, 355, 356, 363, 370, 371, 379, 381, 387, 391, 399, 412, 413, 418, 419, 423, 428, 434, 440, 447, 454, 455, 466, 467, 489, 493, 494, 498, 513, 532, 536, 537, 538, 540, 543, 546, 551, 552, 554, 555, 559, 560, 564, 565, 567, 578, 591, 598, 601, 602, 614, 615, 617, 630, 631, 635, 638, 659, 672, 684, 688, 690, 693, 713, 717, 718, 720, 722, 723, 731, 732, 743, 745, 746, 758, 769, 774, 775, 781, 782, 792, 794, 795, 800, 804, 805, 810, 812, 813, 833, 834, 838, 842, 848, 857, 859, 862, 866, 871, 872, 876, 881, 888, 899, 910, 911, 922, 931, 953, 958, 962, 970, 976, 981.

### h6: 141 clean seeds

114, 119, 121, 122, 126, 139, 145, 152, 159, 170, 174, 182, 189, 192, 217, 219, 227, 232, 236, 247, 254, 255, 267, 277, 279, 280, 288, 293, 295, 301, 308, 315, 323, 335, 338, 342, 347, 364, 365, 393, 395, 400, 410, 417, 421, 445, 449, 454, 455, 457, 461, 468, 469, 471, 474, 476, 478, 481, 508, 510, 521, 523, 524, 525, 528, 529, 537, 538, 546, 551, 553, 557, 566, 567, 577, 582, 600, 607, 610, 611, 626, 629, 630, 631, 637, 651, 664, 666, 667, 669, 675, 679, 693, 698, 701, 705, 712, 716, 718, 728, 730, 731, 738, 739, 743, 756, 761, 762, 769, 793, 805, 808, 811, 830, 831, 840, 847, 849, 850, 862, 878, 879, 886, 907, 908, 909, 922, 929, 932, 944, 946, 950, 952, 954, 965, 969, 987, 988, 989, 991, 992.

### h7: 150 clean seeds

106, 113, 114, 121, 122, 133, 137, 148, 165, 171, 175, 179, 184, 193, 206, 214, 215, 219, 223, 237, 248, 252, 254, 262, 264, 271, 284, 288, 289, 298, 305, 306, 308, 329, 332, 338, 344, 348, 350, 352, 355, 358, 379, 381, 385, 387, 388, 392, 393, 397, 400, 401, 409, 412, 414, 426, 435, 437, 445, 448, 452, 471, 479, 483, 491, 501, 504, 505, 517, 520, 522, 526, 546, 547, 549, 555, 559, 570, 573, 576, 577, 581, 583, 610, 612, 614, 629, 630, 642, 652, 664, 669, 678, 682, 684, 687, 691, 697, 706, 709, 715, 716, 717, 718, 722, 724, 731, 734, 737, 739, 741, 748, 764, 767, 771, 774, 775, 778, 779, 780, 785, 795, 806, 814, 816, 837, 855, 857, 859, 863, 877, 879, 888, 891, 902, 918, 933, 934, 937, 940, 945, 954, 965, 971, 979, 980, 983, 987, 993, 999.

### h8: 154 clean seeds

110, 114, 116, 118, 120, 122, 124, 125, 133, 145, 158, 167, 187, 188, 193, 198, 200, 228, 232, 237, 245, 248, 250, 254, 257, 260, 263, 268, 274, 276, 281, 292, 293, 294, 295, 296, 300, 308, 322, 333, 334, 337, 339, 362, 363, 364, 371, 378, 379, 380, 381, 387, 393, 396, 399, 401, 413, 427, 428, 433, 434, 435, 436, 443, 461, 466, 467, 468, 472, 473, 484, 485, 488, 489, 498, 499, 506, 511, 533, 535, 539, 543, 548, 564, 567, 571, 576, 577, 580, 593, 597, 605, 607, 621, 631, 638, 641, 642, 645, 649, 651, 652, 655, 659, 663, 666, 668, 670, 674, 691, 693, 703, 704, 705, 706, 709, 714, 718, 728, 729, 732, 736, 737, 738, 739, 740, 744, 746, 763, 771, 780, 785, 809, 821, 835, 837, 843, 846, 847, 859, 864, 870, 877, 885, 893, 897, 899, 927, 948, 972, 980, 986, 989, 996.

### h9: 159 clean seeds

103, 108, 117, 120, 127, 135, 138, 150, 171, 180, 181, 190, 201, 204, 211, 222, 226, 228, 229, 236, 248, 252, 263, 265, 266, 289, 299, 300, 306, 308, 315, 317, 322, 327, 332, 335, 348, 358, 360, 370, 371, 382, 385, 386, 412, 417, 426, 428, 431, 443, 447, 450, 456, 460, 475, 480, 481, 482, 495, 498, 504, 516, 520, 522, 525, 526, 537, 538, 557, 579, 581, 601, 611, 612, 620, 622, 626, 637, 648, 649, 654, 660, 661, 667, 670, 675, 678, 691, 703, 704, 706, 707, 709, 710, 712, 716, 717, 718, 721, 726, 727, 729, 735, 741, 750, 752, 753, 754, 758, 766, 767, 771, 786, 788, 790, 799, 803, 805, 806, 807, 811, 816, 823, 825, 827, 832, 838, 842, 847, 857, 858, 862, 881, 882, 888, 896, 904, 909, 912, 914, 919, 923, 931, 942, 946, 951, 959, 961, 963, 967, 973, 979, 981, 985, 986, 987, 989, 992, 998.

## Verification and artifacts

- All seven uploads used their prepared conjunctions; no ordinary fallback was used. All recorded bands match actual uploads.
- All seven archives identify this task; contracts and references match the published task after the seed stamp.
- All 1,750 submitted rows passed stage 1/2, with accessibility 0.87. Local three-seed scores exactly match all seven official scores.
- Reused h7/h8/h9 archived validations and generated h0/h4/h5/h6 validations in report-local directories.
- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.
- Recovered sets over **900 × 250 × 7 = 1,575,000 row/seed combinations**.
- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized build/upload and registration evidence](build_evidence.json), [layout source evidence](layout_evidence.json), [manifest](manifest.json), [CSV](summary.csv).
- Reproduce with `.venv/bin/python reports/00b2e47c/validate.py`, then `.venv/bin/python reports/00b2e47c/audit.py`, then `.venv/bin/python reports/00b2e47c/render_report.py`. The analysis includes validator-code and submission SHA-256 hashes.

| Hotkey | Public hotkey | Upload time UTC | Submission archive |
|---|---|---|---|
| h0 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | 2026-09-20 21:28:21,380 | [data/inst/niome_hotkey/result/2026-09-20T21:28:51](../../data/inst/niome_hotkey/result/2026-09-20T21:28:51/submission.json) |
| h1 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | N/A | No upload archive |
| h2 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | N/A | No upload archive |
| h3 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | N/A | No upload archive |
| h4 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | 2026-09-20 21:02:07,956 | [data/inst/niome_hotkey4/result/2026-09-20T21:02:37](../../data/inst/niome_hotkey4/result/2026-09-20T21:02:37/submission.json) |
| h5 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | 2026-09-20 21:20:56,049 | [data/inst/niome_hotkey5/result/2026-09-20T21:21:25](../../data/inst/niome_hotkey5/result/2026-09-20T21:21:25/submission.json) |
| h6 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | 2026-09-20 21:22:42,604 | [data/inst/niome_hotkey6/result/2026-09-20T21:23:12](../../data/inst/niome_hotkey6/result/2026-09-20T21:23:12/submission.json) |
| h7 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | 2026-09-20 21:22:08,876 | [data/inst/niome_hotkey7/result/2026-09-20T21:22:39](../../data/inst/niome_hotkey7/result/2026-09-20T21:22:39/submission.json) |
| h8 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | 2026-09-20 20:49:24,604 | [data/inst/niome_hotkey8/result/2026-09-20T20:49:54](../../data/inst/niome_hotkey8/result/2026-09-20T20:49:54/submission.json) |
| h9 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | 2026-09-20 20:37:40,054 | [data/inst/niome_hotkey9/result/2026-09-20T20:38:09](../../data/inst/niome_hotkey9/result/2026-09-20T20:38:09/submission.json) |
