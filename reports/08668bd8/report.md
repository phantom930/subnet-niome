# Submission audit: 08668bd8-9461-41b1-b1e4-083dec3c9f18

**CD34+_HSPC**, created **2026-09-19T00:44:02.247879 UTC**. Scoring seeds: **709, 369, 545**.

**Best of h0–h9: h6, 24.747346, rank 209/248. Leader: UID 247, 151.038996.** h6 achieved 16.38% of the leader, a gap of 126.291649. Top-10 cutoff: 111.604780; h6 was 86.857434 below it.

**All ten uploaded 250 valid rows. Nine uploaded prepared conjunctions with 100 Cas12a + 150 Cas9 rows and 12-seed HDR bands. h4 uploaded an ordinary fallback with 75 Cas12a + 175 Cas9 rows, eight clean seeds and no HDR band.** All ten have zero reward weight in this score snapshot.

**No uploaded clean set or HDR band hit any of the three scoring seeds.** The nine uploaded conjunction bands match their window records. h4's record contains a prepared 12-seed band that completed after its fallback upload and was not submitted.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=08668bd8-9461-41b1-b1e4-083dec3c9f18&limit=40000), fetched 2026-09-19T03:34:55.678958+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-19T03:07:44.213690 UTC. Audit generated 2026-09-19T03:36:39.520694+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 247 | 1 | 151.038996 | 240.878051 | 0.666704 | 0.940500 | 250 | 0.294000 |
| h0 | 248 | 233 | 20.976572 | 231.304787 | 0.097562 | 0.929543 | 250 | 0.000000 |
| h1 | 41 | 229 | 21.178910 | 233.513260 | 0.099404 | 0.912405 | 250 | 0.000000 |
| h2 | 163 | 228 | 21.218087 | 236.092321 | 0.097756 | 0.919347 | 250 | 0.000000 |
| h3 | 118 | 234 | 20.890433 | 237.233853 | 0.097456 | 0.903574 | 250 | 0.000000 |
| h4 | 190 | 231 | 21.158120 | 223.548712 | 0.101594 | 0.931616 | 250 | 0.000000 |
| h5 | 234 | 226 | 22.567550 | 236.431543 | 0.105615 | 0.903758 | 250 | 0.000000 |
| h6 | 175 | 209 | 24.747346 | 235.905393 | 0.115092 | 0.911474 | 250 | 0.000000 |
| h7 | 65 | 227 | 22.378630 | 234.942519 | 0.103306 | 0.922031 | 250 | 0.000000 |
| h8 | 211 | 225 | 23.279031 | 233.721638 | 0.107400 | 0.927385 | 250 | 0.000000 |
| h9 | 245 | 235 | 20.667429 | 234.109033 | 0.096666 | 0.913257 | 250 | 0.000000 |

Rank is 1 plus the number of strictly higher scores; ties share rank. Final score is the mean of the three per-seed scores. Breakdown factors are also means; multiplying the mean factors is not generally equivalent to averaging the per-seed products.

Leader hotkey: `5DAeKBa9KqXpQf4L61zywz5Yzwjz5dz8KZbbZ96rDdRgTU5X`.

## Joined windows and per-hotkey candidates

The **00:17 UTC plan** selected width-100 windows **100 / 200 / 300**, joined as **100–399** (300 seeds). Actual windows were **700 / 300 / 500** for seeds **709 / 369 / 545**: **1 of 3** matched the prediction, through seed **369**.

h0–h5 were assigned 150-seed circular slices of the predicted space at stride 50. h6–h9 were assigned 150-seed slices of its complement, **400–999**, at stride 150. Their four slices partition the 600-seed complement. The nine submitted conjunctions all used a full **100–999 cut search**. h4's assigned candidate slice was used by its prepared build, but its actual upload used the ordinary fallback.

| Hotkey | Actual construction | Band group | Assigned 150-seed candidate slice | Offset | Actual clean /900 |
|---|---|---|---|---:|---:|
| h0 | conjunction | predicted | 100–249 | 0 | 98 |
| h1 | conjunction | predicted | 150–299 | 50 | 109 |
| h2 | conjunction | predicted | 200–349 | 100 | 85 |
| h3 | conjunction | predicted | 250–399 | 150 | 85 |
| h4 | ordinary fallback | predicted | 100–149 ∪ 300–399 (not used by upload) | 200 | 8 |
| h5 | conjunction | predicted | 100–199 ∪ 350–399 | 250 | 100 |
| h6 | conjunction | complement | 400–549 | 0 | 89 |
| h7 | conjunction | complement | 550–699 | 150 | 66 |
| h8 | conjunction | complement | 700–849 | 300 | 92 |
| h9 | conjunction | complement | 850–999 | 450 | 91 |

The exact slices preserve gaps in the joined space. The build log's min/max display can be misleading: h5's `100-399 (150 seeds)` means **100–199 ∪ 350–399**, because its slice wraps.

## Actual uploaded clean sets and HDR bands

**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. Both sets were recovered over 100–999. For each submitted conjunction, the full-row clean set equals the Cas12a group clean set, and its count matches the original 12-band build log. The h4 fallback was replayed separately from its unsubmitted prepared builds.

| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |
|---|---:|---|---|---|
| h0 | 98 | 120, 162, 171, 173, 180, 182, 183, 190, 195, 200, 210, 223 | none | none |
| h1 | 109 | 152, 169, 183, 198, 206, 221, 223, 226, 230, 246, 254, 288 | none | none |
| h2 | 85 | 200, 220, 223, 225, 242, 263, 288, 301, 308, 316, 326, 343 | none | none |
| h3 | 85 | 251, 260, 265, 293, 300, 308, 319, 353, 354, 357, 363, 372 | none | none |
| h4 | 8 | none | none | none |
| h5 | 100 | 108, 125, 132, 162, 172, 173, 174, 188, 189, 191, 199, 384 | none | none |
| h6 | 89 | 401, 414, 416, 429, 434, 435, 481, 485, 507, 524, 538, 549 | none | none |
| h7 | 66 | 557, 564, 579, 581, 603, 604, 621, 625, 634, 648, 651, 696 | none | none |
| h8 | 92 | 701, 738, 749, 765, 767, 778, 792, 795, 804, 815, 824, 848 | none | none |
| h9 | 91 | 865, 886, 889, 894, 914, 945, 953, 954, 964, 972, 986, 996 | none | none |

The uploaded fleet union contains **566/900 clean seeds** and **100 distinct HDR-band seeds** (108 band memberships before overlaps). These totals exclude all unsubmitted prepared/rebuilt rows. A candidate window does not guarantee a clean or HDR hit: 369 lay in the submitted candidate slices of h3 and h5; 545 lay in h6's; 709 lay in h8's. All missed their respective clean sets and HDR bands. All other uploaded clean sets also missed every scoring seed.

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.

| Hotkey | Seed 709 | Seed 369 | Seed 545 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 18.394721 / 123 / 5 | 22.178434 / 119 / 2 | 22.356562 / 134 / 4 | 20.976572 |
| h1 | 17.573897 / 135 / 5 | 24.446710 / 141 / 2 | 21.516124 / 131 / 6 | 21.178910 |
| h2 | 21.423859 / 138 / 4 | 18.644267 / 133 / 9 | 23.586135 / 122 / 3 | 21.218087 |
| h3 | 18.669300 / 120 / 7 | 21.709808 / 141 / 4 | 22.292193 / 137 / 3 | 20.890433 |
| h4 | 26.041321 / 124 / 2 | 19.291172 / 128 / 5 | 18.141868 / 137 / 9 | 21.158120 |
| h5 | 23.739102 / 138 / 3 | 22.520001 / 124 / 5 | 21.443546 / 133 / 3 | 22.567550 |
| h6 | 27.481407 / 146 / 1 | 26.165091 / 129 / 3 | 20.595541 / 141 / 5 | 24.747346 |
| h7 | 20.204662 / 131 / 3 | 20.810703 / 138 / 6 | 26.120525 / 134 / 3 | 22.378630 |
| h8 | 23.334458 / 130 / 4 | 21.882876 / 135 / 4 | 24.619759 / 135 / 5 | 23.279031 |
| h9 | 22.968276 / 122 / 1 | 18.004161 / 132 / 5 | 21.029850 / 120 / 7 | 20.667429 |

## Comparison with the leader

h6's weighted score was **235.905393**, versus the leader's **240.878051**. Mean consistency was **0.115092 versus 0.666704**, and fidelity was **0.911474 versus 0.940500**. The largest factor gap is consistency. h6 scored **27.481407** at 709, **26.165091** at 369, and **20.595541** at 545. It had no clean or all-HDR seed to lift the three-seed mean. The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the inspected data; they cannot be recovered from its aggregate score.

## h4 fallback and prepared-band provenance

h4 began preparation at **00:44:11 UTC**. The validator requested its submission at **00:46:27**. At **00:46:28** it waited up to **179 seconds** for the prepared build. At **00:49:27**, that build was still running (315 seconds elapsed), so the miner generated the ordinary fallback while the prepared build still held the GPU. It uploaded at **00:49:32.265**, with **100 seconds** of URL TTL left. The prepared 12-band conjunction finished at **00:51:56.733**, about **144 seconds after the upload**. Its logged **94 clean seeds** do not describe the fallback's **eight clean seeds**.

h4's recorded prepared band, **not submitted**, is: **121, 124, 130, 300, 309, 324, 335, 357, 362, 378, 380, 396**.

h4's latest window record mixes `source=all_hdr_not_attempted` and `window=null` from its fallback with `all_hdr_built=true` and a band from the completed background preparation. Uploaded-row replay confirms that this recorded band is absent from the submitted fallback. The other nine hotkeys used their original prepared 12-band submissions, and every one of their recorded bands matches the corresponding uploaded rows.


## Exact uploaded clean sets

Each list includes its HDR band and covers the entire 100–999 seed space.

### h0: 98 clean seeds

102, 106, 120, 145, 162, 168, 171, 173, 180, 182, 183, 190, 191, 195, 196, 199, 200, 210, 212, 218, 221, 223, 229, 252, 255, 265, 272, 286, 292, 293, 298, 302, 321, 323, 326, 329, 330, 360, 361, 367, 368, 370, 381, 392, 394, 398, 431, 437, 440, 441, 448, 453, 455, 466, 499, 515, 524, 527, 535, 546, 548, 549, 570, 592, 606, 613, 653, 670, 671, 687, 696, 707, 720, 736, 739, 755, 770, 772, 795, 824, 832, 846, 855, 858, 865, 886, 899, 904, 911, 925, 929, 939, 955, 957, 959, 962, 985, 987.

### h1: 109 clean seeds

110, 116, 143, 149, 152, 162, 164, 165, 169, 170, 176, 183, 189, 198, 203, 204, 206, 211, 221, 222, 223, 226, 228, 230, 231, 246, 248, 254, 270, 288, 299, 305, 320, 322, 324, 329, 331, 334, 345, 363, 382, 388, 390, 407, 427, 433, 435, 464, 468, 485, 494, 501, 510, 512, 514, 519, 526, 530, 538, 544, 550, 570, 574, 575, 586, 608, 613, 617, 624, 625, 631, 640, 641, 650, 657, 670, 699, 701, 723, 739, 743, 746, 752, 764, 771, 783, 786, 797, 804, 820, 828, 841, 845, 858, 886, 908, 912, 924, 928, 929, 931, 948, 951, 960, 963, 976, 978, 986, 998.

### h2: 85 clean seeds

104, 129, 135, 147, 166, 169, 178, 180, 188, 193, 200, 220, 222, 223, 225, 233, 242, 245, 263, 268, 276, 288, 291, 298, 301, 308, 309, 310, 316, 326, 327, 332, 333, 339, 343, 355, 356, 364, 398, 402, 412, 465, 466, 475, 496, 505, 508, 533, 540, 551, 555, 576, 577, 583, 585, 603, 622, 629, 632, 634, 640, 649, 671, 683, 686, 704, 727, 734, 735, 736, 746, 751, 758, 820, 837, 840, 845, 846, 849, 855, 909, 916, 942, 988, 998.

### h3: 85 clean seeds

113, 114, 123, 149, 157, 166, 171, 185, 190, 206, 210, 214, 219, 227, 232, 251, 253, 260, 265, 273, 278, 293, 300, 308, 319, 334, 339, 350, 353, 354, 357, 363, 368, 372, 376, 382, 404, 407, 424, 444, 489, 500, 508, 524, 541, 573, 591, 597, 602, 636, 666, 671, 708, 719, 739, 743, 748, 760, 766, 774, 779, 782, 792, 794, 803, 806, 827, 831, 836, 840, 845, 859, 862, 878, 901, 903, 907, 919, 938, 951, 958, 975, 988, 992, 999.

### h4: 8 clean seeds

250, 333, 373, 575, 780, 865, 904, 962.

### h5: 100 clean seeds

106, 108, 119, 125, 132, 142, 159, 162, 172, 173, 174, 183, 188, 189, 191, 194, 197, 199, 221, 233, 236, 239, 246, 273, 282, 285, 305, 315, 319, 326, 344, 349, 371, 384, 403, 404, 409, 415, 427, 428, 439, 443, 459, 476, 477, 496, 498, 505, 508, 512, 547, 554, 559, 571, 589, 596, 600, 606, 616, 619, 642, 649, 660, 662, 673, 694, 714, 722, 739, 740, 744, 747, 748, 760, 769, 775, 790, 798, 803, 807, 819, 826, 828, 831, 842, 853, 859, 860, 862, 871, 879, 895, 898, 905, 922, 926, 942, 951, 955, 967.

### h6: 89 clean seeds

111, 124, 125, 126, 139, 142, 185, 192, 193, 201, 214, 216, 239, 246, 251, 256, 288, 308, 324, 344, 353, 378, 394, 401, 408, 414, 416, 419, 429, 434, 435, 448, 451, 455, 456, 466, 467, 471, 481, 485, 507, 509, 515, 524, 528, 534, 537, 538, 549, 550, 562, 575, 579, 600, 626, 629, 633, 638, 677, 685, 689, 693, 698, 711, 729, 746, 751, 754, 756, 789, 795, 796, 813, 823, 826, 839, 846, 864, 865, 875, 878, 899, 917, 930, 935, 954, 955, 991, 998.

### h7: 66 clean seeds

125, 134, 136, 140, 141, 165, 181, 186, 200, 201, 292, 295, 305, 310, 330, 333, 363, 365, 390, 397, 398, 421, 476, 519, 544, 557, 564, 579, 581, 582, 587, 590, 595, 600, 603, 604, 614, 621, 625, 634, 638, 647, 648, 651, 669, 670, 685, 687, 696, 723, 775, 784, 796, 817, 846, 879, 880, 883, 885, 896, 909, 957, 977, 987, 991, 997.

### h8: 92 clean seeds

104, 114, 122, 130, 136, 142, 155, 190, 193, 198, 205, 238, 248, 253, 270, 276, 280, 286, 293, 304, 312, 313, 333, 350, 382, 407, 440, 443, 447, 449, 470, 482, 492, 513, 522, 527, 536, 565, 567, 579, 586, 598, 619, 622, 626, 628, 643, 687, 700, 701, 707, 716, 717, 718, 732, 738, 749, 754, 763, 764, 765, 767, 769, 771, 778, 792, 795, 802, 804, 805, 813, 815, 824, 829, 848, 856, 858, 860, 864, 865, 871, 882, 886, 902, 912, 932, 942, 951, 962, 964, 995, 997.

### h9: 91 clean seeds

129, 132, 135, 147, 163, 184, 200, 203, 207, 210, 224, 241, 271, 275, 276, 283, 288, 308, 392, 410, 412, 422, 442, 444, 451, 454, 463, 465, 484, 491, 507, 508, 511, 516, 529, 542, 566, 598, 644, 655, 656, 682, 683, 684, 690, 692, 703, 716, 744, 745, 751, 752, 758, 760, 764, 776, 789, 811, 837, 844, 851, 865, 874, 878, 886, 889, 892, 894, 902, 914, 918, 919, 929, 932, 934, 935, 945, 950, 953, 954, 955, 959, 964, 972, 982, 986, 988, 989, 995, 996, 998.

## Verification and artifacts

- All ten archives confirm successful uploads of this task; contracts and references match the published task after the seed stamp.
- All 2,500 submitted rows passed stage 1/2. Accessibility is 0.87. Local three-seed final scores exactly match all ten official API scores.
- Reused nine archived validations and generated h8's missing validation in this report's h8 directory, without modifying its upload archive.
- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.
- Recovered clean sets and bands over **900 × 250 × 10 = 2,250,000 row/seed combinations**.
- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized build/upload evidence](build_evidence.json), [archive manifest](manifest.json), [summary CSV](summary.csv).
- Reproduce using `.venv/bin/python reports/08668bd8/validate.py`, then `.venv/bin/python reports/08668bd8/audit.py`, then `.venv/bin/python reports/08668bd8/render_report.py`. The analysis includes current validator-code and submission SHA-256 hashes.

| Hotkey | Public hotkey | Upload time UTC | Submission archive |
|---|---|---|---|
| h0 / UID 248 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | 2026-09-19 01:00:30,786 | [data/inst/niome_hotkey/result/2026-09-19T01:01:00](../../data/inst/niome_hotkey/result/2026-09-19T01:01:00/submission.json) |
| h1 / UID 41 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | 2026-09-19 01:01:06,855 | [data/inst/niome_hotkey1/result/2026-09-19T01:01:35](../../data/inst/niome_hotkey1/result/2026-09-19T01:01:35/submission.json) |
| h2 / UID 163 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | 2026-09-19 01:13:50,812 | [data/inst/niome_hotkey2/result/2026-09-19T01:14:19](../../data/inst/niome_hotkey2/result/2026-09-19T01:14:19/submission.json) |
| h3 / UID 118 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | 2026-09-19 01:24:17,033 | [data/inst/niome_hotkey3/result/2026-09-19T01:24:45](../../data/inst/niome_hotkey3/result/2026-09-19T01:24:45/submission.json) |
| h4 / UID 190 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | 2026-09-19 00:49:32,265 | [data/inst/niome_hotkey4/result/2026-09-19T00:49:59](../../data/inst/niome_hotkey4/result/2026-09-19T00:49:59/submission.json) |
| h5 / UID 234 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | 2026-09-19 01:17:34,369 | [data/inst/niome_hotkey5/result/2026-09-19T01:18:03](../../data/inst/niome_hotkey5/result/2026-09-19T01:18:03/submission.json) |
| h6 / UID 175 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | 2026-09-19 01:31:02,550 | [data/inst/niome_hotkey6/result/2026-09-19T01:31:31](../../data/inst/niome_hotkey6/result/2026-09-19T01:31:31/submission.json) |
| h7 / UID 65 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | 2026-09-19 01:25:38,329 | [data/inst/niome_hotkey7/result/2026-09-19T01:26:07](../../data/inst/niome_hotkey7/result/2026-09-19T01:26:07/submission.json) |
| h8 / UID 211 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | 2026-09-19 01:07:21,565 | [data/inst/niome_hotkey8/result/2026-09-19T01:07:50](../../data/inst/niome_hotkey8/result/2026-09-19T01:07:50/submission.json) |
| h9 / UID 245 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | 2026-09-19 01:13:07,807 | [data/inst/niome_hotkey9/result/2026-09-19T01:13:36](../../data/inst/niome_hotkey9/result/2026-09-19T01:13:36/submission.json) |
