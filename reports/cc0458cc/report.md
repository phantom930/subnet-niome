# Submission audit: cc0458cc-f62c-43db-946e-a7d835e5bf69

**CD34+_HSPC**, created **2026-09-18T22:19:29.341104 UTC**. Scoring seeds: **161, 326, 557**.

**Best of h0–h9: h8, 27.797746, rank 89/248. Leader: UID 24, 143.501616.** h8 achieved 19.37% of the leader, a gap of 115.703871. Top-10 cutoff: 89.429841; h8 was 61.632095 below it.

**All ten uploaded 250 valid rows. Nine uploaded prepared conjunctions with 100 Cas12a + 150 Cas9 rows and 11-seed HDR bands. h4 uploaded an ordinary fallback with 76 Cas12a + 174 Cas9 rows, nine clean seeds and no HDR band.** All ten have zero reward weight in this score snapshot.

**No uploaded HDR band hit a scoring seed.** h1 and h6 were clean at **326**; h8 was clean at **161**. Seed **557** was outside every uploaded clean set. All ten later rebuilt this task with 12-seed bands, overwriting their latest window records after uploading; those later bands were not submitted for this task.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=cc0458cc-f62c-43db-946e-a7d835e5bf69&limit=40000), fetched 2026-09-19T02:01:27.204026+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-19T00:42:55.253346 UTC. Audit generated 2026-09-19T02:03:11.738347+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 24 | 1 | 143.501616 | 226.609895 | 0.664694 | 0.952700 | 250 | 0.294000 |
| h0 | 248 | 232 | 18.570534 | 203.755706 | 0.097663 | 0.933225 | 250 | 0.000000 |
| h1 | 41 | 140 | 23.354009 | 202.780643 | 0.122349 | 0.941317 | 250 | 0.000000 |
| h2 | 163 | 227 | 19.419739 | 207.794882 | 0.102142 | 0.914961 | 250 | 0.000000 |
| h3 | 118 | 230 | 19.114511 | 202.054837 | 0.099096 | 0.954632 | 250 | 0.000000 |
| h4 | 190 | 206 | 21.038027 | 197.146791 | 0.110409 | 0.966519 | 250 | 0.000000 |
| h5 | 234 | 204 | 21.183703 | 202.036617 | 0.111498 | 0.940386 | 250 | 0.000000 |
| h6 | 175 | 105 | 24.744912 | 203.042922 | 0.131139 | 0.929319 | 250 | 0.000000 |
| h7 | 65 | 240 | 17.223376 | 203.310113 | 0.090873 | 0.932228 | 250 | 0.000000 |
| h8 | 211 | 89 | 27.797746 | 202.866658 | 0.147170 | 0.931064 | 250 | 0.000000 |
| h9 | 245 | 225 | 19.798691 | 201.534692 | 0.105833 | 0.928249 | 250 | 0.000000 |

Rank is 1 plus the number of strictly higher scores; ties share rank. Final score is the mean of the three per-seed scores. Breakdown factors are also means; multiplying the mean factors is not generally equivalent to averaging the per-seed products.

Leader hotkey: `5EAS7Sn71N787SstNTbPREbfzXCB7Syz8iBWrAGJPhhJcADe`.

## Joined windows and per-hotkey candidates

The **22:17 UTC plan** selected width-100 windows **200 / 400 / 900**, joined as **200–299 ∪ 400–499 ∪ 900–999** (300 seeds). Actual windows were **100 / 300 / 500** for seeds **161 / 326 / 557**: **0 of 3** matched the prediction.

h0–h5 were assigned 150-seed circular slices of the predicted space at stride 50. h6–h9 were assigned 150-seed slices of its complement, **100–199 ∪ 300–399 ∪ 500–899**, at stride 150. Their four slices partition the 600-seed complement. The nine submitted conjunctions all used a full **100–999 cut search**. h4's assigned candidate slice was used by its prepared build, but its actual upload used the ordinary fallback.

| Hotkey | Actual construction | Band group | Assigned 150-seed candidate slice | Offset | Actual clean /900 |
|---|---|---|---|---:|---:|
| h0 | conjunction | predicted | 200–299 ∪ 400–449 | 0 | 141 |
| h1 | conjunction | predicted | 250–299 ∪ 400–499 | 50 | 137 |
| h2 | conjunction | predicted | 400–499 ∪ 900–949 | 100 | 144 |
| h3 | conjunction | predicted | 450–499 ∪ 900–999 | 150 | 139 |
| h4 | ordinary fallback | predicted | 200–249 ∪ 900–999 (not used by upload) | 200 | 9 |
| h5 | conjunction | predicted | 200–299 ∪ 950–999 | 250 | 139 |
| h6 | conjunction | complement | 100–199 ∪ 300–349 | 0 | 149 |
| h7 | conjunction | complement | 350–399 ∪ 500–599 | 150 | 137 |
| h8 | conjunction | complement | 600–749 | 300 | 142 |
| h9 | conjunction | complement | 750–899 | 450 | 132 |

The exact slices preserve gaps in the joined space. The build log's min/max display can be misleading: for example, h6's `100-349 (150 seeds)` means **100–199 ∪ 300–349**.

## Actual uploaded clean sets and HDR bands

**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. Both sets were recovered over 100–999. For each submitted conjunction, the full-row clean set equals the Cas12a group clean set, and its count matches the original 11-band build log. The h4 fallback was replayed separately from its unsubmitted prepared builds.

| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |
|---|---:|---|---|---|
| h0 | 141 | 204, 213, 235, 247, 256, 260, 268, 284, 409, 419, 442 | none | none |
| h1 | 137 | 263, 264, 266, 287, 400, 410, 411, 432, 453, 481, 489 | 326 | none |
| h2 | 144 | 400, 406, 411, 414, 432, 453, 468, 470, 489, 926, 929 | none | none |
| h3 | 139 | 453, 459, 461, 464, 466, 486, 908, 923, 925, 975, 994 | none | none |
| h4 | 9 | none | none | none |
| h5 | 139 | 204, 205, 231, 256, 268, 285, 294, 972, 973, 988, 993 | none | none |
| h6 | 149 | 114, 122, 141, 143, 149, 152, 168, 190, 302, 321, 342 | 326 | none |
| h7 | 137 | 363, 365, 382, 387, 390, 516, 540, 544, 576, 587, 595 | none | none |
| h8 | 142 | 608, 620, 629, 637, 638, 649, 696, 709, 720, 731, 743 | 161 | none |
| h9 | 132 | 760, 769, 795, 812, 856, 857, 862, 865, 881, 887, 890 | none | none |

The uploaded fleet union contains **706/900 clean seeds** and **90 distinct HDR-band seeds** (99 band memberships before overlaps). These totals exclude all unsubmitted prepared/rebuilt rows. A candidate window does not guarantee a clean or HDR hit: h6's candidate slice contained 161 and 326, but only 326 was clean and neither was in its HDR band; h7's candidate slice contained 557, which was not clean.

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.

| Hotkey | Seed 161 | Seed 326 | Seed 557 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 17.875330 / 134 / 4 | 18.468977 / 135 / 6 | 19.367294 / 136 / 6 | 18.570534 |
| h1 | 17.648707 / 116 / 5 | 36.034192 / 141 / 0 | 16.379127 / 142 / 6 | 23.354009 |
| h2 | 19.523752 / 140 / 7 | 19.130469 / 132 / 5 | 19.604996 / 146 / 4 | 19.419739 |
| h3 | 18.809747 / 137 / 3 | 17.800899 / 141 / 6 | 20.732887 / 124 / 4 | 19.114511 |
| h4 | 19.980974 / 136 / 3 | 20.194172 / 122 / 3 | 22.938934 / 140 / 1 | 21.038027 |
| h5 | 23.642354 / 128 / 6 | 20.947099 / 139 / 4 | 18.961657 / 128 / 6 | 21.183703 |
| h6 | 17.092636 / 138 / 6 | 39.133034 / 132 / 0 | 18.009065 / 128 / 3 | 24.744912 |
| h7 | 15.686851 / 133 / 3 | 18.296780 / 116 / 6 | 17.686496 / 126 / 4 | 17.223376 |
| h8 | 39.986504 / 144 / 0 | 21.218207 / 136 / 5 | 22.188526 / 137 / 2 | 27.797746 |
| h9 | 18.011162 / 137 / 6 | 21.585007 / 142 / 4 | 19.799902 / 128 / 4 | 19.798691 |

## Comparison with the leader

h8's weighted score was **202.866658**, versus the leader's **226.609895**. Mean consistency was **0.147170 versus 0.664694**, and fidelity was **0.931064 versus 0.952700**. The largest factor gap is consistency. h8's clean hit at 161 scored **39.986504**, with **21.218207** at 326 and **22.188526** at 557. It had no all-HDR seed to lift the three-seed mean. The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the inspected data; they cannot be recovered from its aggregate score.

## h4 fallback and later rebuilds

h4 began preparation at **22:19:43 UTC**. The validator requested its submission at **22:21:03**. At **22:21:08** it waited up to **174 seconds** for the prepared build. At **22:24:02**, that build was still running (259 seconds elapsed), so the miner generated the ordinary fallback while the prepared build still held the GPU. It uploaded at **22:24:05.683**, with **101 seconds** of URL TTL left. The prepared 11-band conjunction finished at **22:27:53.941**, about **228 seconds after the upload**. Its logged 137 clean seeds do not describe the fallback that was submitted.

The nine other uploads used their original prepared 11-band builds. After each hotkey had uploaded, a later rebuild began at **23:07–23:16 UTC** using **k=12**. All latest window records therefore contain 12-seed bands. The later logs also report different clean counts; these are prepared-build observations, not validated submitted sets. For the nine original conjunctions, the later band record retains the 11 uploaded band seeds and adds one seed. h4's entire later band is unrelated to any uploaded HDR band.

| Hotkey | Uploaded clean /900 | Later logged clean /900 (unsubmitted) | Later recorded band seeds absent from uploaded band |
|---|---:|---:|---|
| h0 | 141 | 88 | 433 |
| h1 | 137 | 86 | 449 |
| h2 | 144 | 92 | 407 |
| h3 | 139 | 101 | 456 |
| h4 | 9 | 79 | 204, 247, 905, 924, 936, 957, 960, 972, 988, 993, 994, 999 |
| h5 | 139 | 93 | 242 |
| h6 | 149 | 77 | 320 |
| h7 | 137 | 104 | 530 |
| h8 | 142 | 92 | 623 |
| h9 | 132 | 86 | 888 |

## Exact uploaded clean sets

Each list includes its HDR band and covers the entire 100–999 seed space.

### h0: 141 clean seeds

110, 128, 130, 147, 155, 159, 160, 162, 167, 169, 172, 182, 197, 204, 213, 217, 218, 219, 221, 222, 225, 233, 235, 240, 241, 247, 256, 260, 261, 264, 267, 268, 284, 311, 320, 323, 328, 330, 332, 337, 342, 346, 351, 355, 356, 359, 364, 374, 376, 380, 390, 393, 395, 409, 412, 419, 423, 435, 436, 442, 443, 444, 448, 449, 453, 461, 478, 483, 505, 511, 519, 528, 536, 542, 543, 550, 571, 581, 583, 589, 597, 599, 602, 619, 622, 623, 629, 636, 642, 643, 657, 665, 671, 678, 679, 692, 703, 721, 732, 738, 742, 749, 758, 761, 783, 795, 801, 802, 818, 824, 825, 831, 832, 833, 840, 842, 847, 867, 871, 872, 879, 883, 894, 899, 918, 920, 924, 929, 936, 946, 962, 963, 964, 965, 969, 971, 979, 982, 992, 994, 998.

### h1: 137 clean seeds

110, 114, 119, 131, 135, 137, 140, 141, 150, 170, 174, 188, 197, 198, 205, 229, 233, 237, 246, 249, 263, 264, 266, 274, 279, 284, 287, 291, 295, 306, 318, 324, 326, 329, 333, 335, 336, 352, 375, 378, 380, 388, 397, 400, 409, 410, 411, 414, 432, 433, 438, 441, 450, 453, 467, 478, 481, 483, 489, 518, 521, 537, 546, 552, 554, 555, 559, 587, 588, 598, 600, 602, 614, 618, 635, 637, 641, 647, 649, 651, 653, 654, 655, 657, 658, 659, 663, 667, 668, 671, 685, 695, 699, 707, 741, 742, 743, 749, 752, 767, 773, 777, 779, 784, 789, 805, 810, 812, 816, 837, 840, 851, 865, 866, 873, 898, 899, 901, 904, 905, 906, 907, 908, 910, 918, 921, 923, 931, 934, 936, 942, 948, 956, 959, 966, 968, 969.

### h2: 144 clean seeds

104, 105, 109, 114, 117, 120, 122, 128, 130, 156, 173, 176, 182, 184, 187, 189, 190, 202, 206, 208, 218, 219, 222, 235, 237, 238, 249, 266, 280, 281, 283, 286, 290, 294, 297, 303, 306, 307, 308, 314, 319, 321, 334, 335, 352, 353, 368, 371, 386, 400, 403, 406, 411, 414, 418, 421, 422, 432, 433, 447, 453, 461, 462, 467, 468, 470, 473, 479, 480, 489, 490, 506, 513, 518, 524, 539, 550, 555, 556, 563, 584, 587, 590, 607, 620, 636, 638, 642, 651, 659, 664, 670, 679, 681, 682, 683, 691, 696, 707, 725, 738, 752, 762, 781, 782, 783, 790, 793, 794, 798, 799, 810, 817, 820, 821, 822, 823, 826, 828, 842, 848, 849, 851, 854, 858, 861, 875, 888, 893, 904, 906, 908, 910, 920, 923, 926, 929, 934, 936, 944, 948, 957, 976, 985.

### h3: 139 clean seeds

107, 113, 128, 139, 141, 145, 149, 156, 165, 173, 178, 185, 190, 201, 217, 224, 228, 238, 246, 252, 266, 272, 274, 275, 287, 290, 313, 317, 328, 330, 338, 351, 359, 371, 375, 378, 380, 382, 384, 385, 387, 399, 406, 407, 416, 417, 426, 440, 442, 447, 449, 452, 453, 459, 461, 463, 464, 466, 482, 483, 486, 490, 495, 510, 513, 530, 534, 548, 550, 553, 559, 573, 588, 603, 608, 611, 626, 637, 638, 644, 658, 659, 664, 667, 680, 685, 694, 697, 708, 714, 719, 726, 727, 731, 736, 742, 746, 747, 758, 761, 779, 782, 783, 793, 816, 819, 821, 833, 842, 845, 849, 856, 861, 864, 868, 876, 886, 895, 908, 915, 918, 919, 923, 925, 934, 947, 949, 953, 955, 975, 979, 980, 982, 990, 993, 994, 996, 997, 999.

### h4: 9 clean seeds

226, 386, 446, 581, 817, 858, 869, 947, 956.

### h5: 139 clean seeds

115, 118, 119, 123, 124, 125, 126, 142, 149, 151, 171, 173, 178, 179, 184, 191, 204, 205, 210, 228, 229, 231, 235, 236, 241, 247, 256, 257, 266, 268, 273, 277, 285, 291, 294, 312, 320, 329, 334, 338, 353, 355, 366, 372, 376, 377, 380, 385, 396, 401, 407, 419, 427, 432, 442, 471, 492, 493, 497, 499, 501, 505, 525, 529, 537, 542, 546, 563, 572, 577, 601, 602, 607, 620, 621, 632, 641, 647, 651, 659, 660, 670, 680, 681, 684, 686, 696, 698, 711, 713, 721, 723, 726, 730, 735, 739, 741, 752, 755, 769, 770, 774, 778, 788, 794, 807, 822, 823, 828, 837, 853, 861, 862, 876, 877, 878, 889, 891, 893, 897, 904, 905, 909, 916, 934, 938, 940, 950, 955, 960, 961, 964, 972, 973, 975, 978, 987, 988, 993.

### h6: 149 clean seeds

110, 114, 122, 124, 134, 141, 143, 149, 150, 152, 155, 165, 166, 168, 176, 188, 190, 205, 215, 218, 224, 225, 232, 233, 246, 251, 252, 273, 284, 300, 302, 304, 316, 317, 321, 324, 326, 333, 342, 350, 354, 358, 363, 367, 370, 373, 383, 387, 390, 404, 407, 411, 413, 436, 440, 446, 449, 452, 456, 459, 460, 464, 467, 470, 474, 482, 484, 486, 495, 501, 514, 515, 517, 532, 533, 534, 549, 553, 554, 556, 559, 562, 564, 569, 575, 577, 579, 585, 586, 590, 601, 602, 603, 606, 607, 621, 623, 630, 635, 637, 654, 667, 670, 678, 679, 680, 682, 687, 700, 704, 710, 714, 719, 739, 741, 756, 763, 774, 796, 801, 812, 814, 819, 825, 830, 831, 833, 835, 837, 838, 850, 854, 880, 883, 886, 887, 891, 912, 923, 931, 943, 946, 953, 965, 966, 975, 993, 995, 999.

### h7: 137 clean seeds

109, 117, 138, 149, 150, 156, 171, 172, 174, 195, 197, 198, 209, 214, 217, 223, 224, 227, 234, 240, 249, 258, 261, 262, 264, 268, 271, 284, 296, 298, 315, 319, 324, 331, 332, 334, 341, 344, 346, 351, 362, 363, 365, 366, 382, 384, 387, 390, 404, 407, 413, 414, 415, 429, 435, 440, 444, 445, 447, 463, 488, 499, 509, 511, 516, 520, 531, 535, 540, 541, 544, 545, 556, 567, 569, 573, 576, 581, 586, 587, 589, 595, 598, 601, 602, 604, 613, 622, 626, 646, 650, 653, 661, 682, 686, 692, 697, 705, 711, 719, 723, 741, 756, 760, 762, 767, 768, 779, 780, 781, 784, 785, 791, 794, 796, 799, 801, 808, 823, 831, 847, 853, 855, 889, 902, 911, 921, 925, 929, 964, 969, 973, 979, 980, 987, 988, 996.

### h8: 142 clean seeds

105, 111, 135, 138, 142, 161, 166, 170, 172, 177, 179, 186, 191, 204, 205, 206, 214, 221, 227, 231, 232, 253, 265, 268, 272, 276, 295, 317, 330, 333, 343, 361, 368, 386, 406, 431, 438, 447, 449, 468, 474, 476, 477, 495, 496, 497, 499, 500, 501, 502, 504, 511, 512, 519, 521, 523, 532, 534, 537, 538, 547, 549, 552, 553, 560, 561, 564, 568, 578, 596, 602, 605, 608, 612, 620, 622, 624, 629, 637, 638, 641, 645, 646, 649, 651, 657, 678, 683, 686, 696, 697, 704, 709, 720, 722, 731, 735, 737, 738, 743, 755, 758, 764, 765, 769, 771, 774, 781, 789, 799, 800, 801, 802, 806, 807, 812, 827, 839, 841, 855, 864, 877, 880, 891, 895, 902, 905, 907, 911, 925, 935, 939, 946, 956, 969, 970, 971, 972, 977, 979, 984, 989.

### h9: 132 clean seeds

101, 105, 120, 126, 129, 132, 158, 183, 184, 187, 189, 194, 195, 216, 223, 228, 236, 252, 269, 278, 279, 289, 292, 308, 318, 327, 329, 331, 334, 337, 341, 348, 350, 351, 352, 353, 362, 364, 367, 368, 391, 406, 417, 423, 428, 429, 430, 432, 436, 443, 447, 453, 455, 456, 458, 470, 474, 476, 480, 490, 496, 497, 501, 554, 563, 589, 591, 594, 599, 606, 608, 627, 631, 634, 638, 648, 674, 677, 681, 683, 685, 699, 700, 702, 703, 712, 717, 728, 747, 757, 760, 766, 767, 768, 769, 790, 795, 797, 801, 807, 812, 826, 829, 833, 835, 848, 856, 857, 862, 865, 868, 869, 874, 877, 880, 881, 885, 887, 890, 895, 901, 910, 913, 929, 932, 935, 936, 942, 950, 961, 965, 967.

## Verification and artifacts

- All ten archives confirm successful uploads of this task; contracts and references match the published task after the seed stamp.
- All 2,500 submitted rows passed stage 1/2. Accessibility is 0.87. Archived three-seed final scores exactly match all ten official API scores.
- Replayed outcomes and indel lengths match archived stage-3 detail for every uploaded row at all three scoring seeds.
- Recovered clean sets and bands over **900 × 250 × 10 = 2,250,000 row/seed combinations**.
- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized build/upload evidence](build_evidence.json), [archive manifest](manifest.json), [summary CSV](summary.csv).
- Reproduce using archived validations: `.venv/bin/python reports/cc0458cc/audit.py`, then `.venv/bin/python reports/cc0458cc/render_report.py`. The analysis includes current validator-code and submission SHA-256 hashes.

| Hotkey | Public hotkey | Upload time UTC | Submission archive |
|---|---|---|---|
| h0 / UID 248 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | 2026-09-18 22:51:20,468 | [data/inst/niome_hotkey/result/2026-09-18T22:51:49](../../data/inst/niome_hotkey/result/2026-09-18T22:51:49/submission.json) |
| h1 / UID 41 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | 2026-09-18 22:35:10,961 | [data/inst/niome_hotkey1/result/2026-09-18T22:35:40](../../data/inst/niome_hotkey1/result/2026-09-18T22:35:40/submission.json) |
| h2 / UID 163 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | 2026-09-18 23:13:55,993 | [data/inst/niome_hotkey2/result/2026-09-18T23:14:25](../../data/inst/niome_hotkey2/result/2026-09-18T23:14:25/submission.json) |
| h3 / UID 118 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | 2026-09-18 23:05:59,645 | [data/inst/niome_hotkey3/result/2026-09-18T23:06:28](../../data/inst/niome_hotkey3/result/2026-09-18T23:06:28/submission.json) |
| h4 / UID 190 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | 2026-09-18 22:24:05,683 | [data/inst/niome_hotkey4/result/2026-09-18T22:24:38](../../data/inst/niome_hotkey4/result/2026-09-18T22:24:38/submission.json) |
| h5 / UID 234 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | 2026-09-18 22:48:03,289 | [data/inst/niome_hotkey5/result/2026-09-18T22:48:34](../../data/inst/niome_hotkey5/result/2026-09-18T22:48:34/submission.json) |
| h6 / UID 175 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | 2026-09-18 22:59:20,783 | [data/inst/niome_hotkey6/result/2026-09-18T22:59:48](../../data/inst/niome_hotkey6/result/2026-09-18T22:59:48/submission.json) |
| h7 / UID 65 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | 2026-09-18 22:31:34,644 | [data/inst/niome_hotkey7/result/2026-09-18T22:32:04](../../data/inst/niome_hotkey7/result/2026-09-18T22:32:04/submission.json) |
| h8 / UID 211 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | 2026-09-18 23:15:57,468 | [data/inst/niome_hotkey8/result/2026-09-18T23:16:00](../../data/inst/niome_hotkey8/result/2026-09-18T23:16:00/submission.json) |
| h9 / UID 245 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | 2026-09-18 22:35:35,401 | [data/inst/niome_hotkey9/result/2026-09-18T22:36:06](../../data/inst/niome_hotkey9/result/2026-09-18T22:36:06/submission.json) |
