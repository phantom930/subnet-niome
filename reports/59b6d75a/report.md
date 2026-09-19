# Submission audit: 59b6d75a-824b-4021-bd7a-3360b43ef0e4

**CD34+_HSPC**, created **2026-09-18T03:03:14.557809 UTC**. Scoring seeds: **733, 312, 818**.

**Best of h0–h6: h4, score 36.160546, rank 198/248. Top miner: UID 132, score 202.741290.** h4 achieved 17.84% of the leader's score, a gap of 166.580744. The top-10 cutoff was 140.548283; h4 was 104.387737 below it. All seven uploads succeeded, contain 250 valid rows, and have zero reward weight.

**h1 submitted an ordinary HDR-construction fallback.** Its prepared conjunction build finished after the wait expired. The other six submitted their prepared conjunction rows. h1's recorded prepared band does not describe its uploaded rows.

Official scores: [task-specific API](https://niome-api.genomes.io/api/v3/miners/scores?task_id=59b6d75a-824b-4021-bd7a-3360b43ef0e4&limit=40000), fetched 2026-09-18T05:44:43.145985+00:00. There are 248 score records, 248 distinct miners, and 1 validator. Leader's score timestamp: 2026-09-18T05:24:15.822336 UTC. Audit generated 2026-09-18T05:46:13.100591+00:00.

## Official score comparison

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 132 | 1 | 202.741290 | 315.569149 | 0.662879 | 0.969200 | 0.294000 |
| h0 | 248 | 241 | 29.228199 | 338.622487 | 0.094453 | 0.913841 | 0.000000 |
| h1 | 41 | 224 | 32.077144 | 323.841011 | 0.108441 | 0.913417 | 0.000000 |
| h2 | 163 | 226 | 31.162057 | 315.772927 | 0.103603 | 0.952532 | 0.000000 |
| h3 | 118 | 240 | 29.839030 | 337.471574 | 0.096093 | 0.920146 | 0.000000 |
| h4 | 190 | 198 | 36.160546 | 332.428168 | 0.117579 | 0.925143 | 0.000000 |
| h5 | 234 | 218 | 33.980081 | 343.165571 | 0.108332 | 0.914039 | 0.000000 |
| h6 | 175 | 243 | 25.606939 | 332.901467 | 0.083069 | 0.925985 | 0.000000 |

Top miner hotkey: `5FEV7ePiuiNVmqP96Fn1A2695EoPFkfbTnopvsWZQp1vjuBJ`.

Final score is the mean over three seeds. For our seven submissions, weighted score and fidelity are constant across seeds, so final = weighted score × mean consistency × fidelity.

## Joined windows, clean sets and bands

The original plan, generated **2026-09-18 02:17:09 UTC**, selected joined band space **100–299 ∪ 900–999** (predicted 100-wide windows **100/200/900**). For the conjunction builds, **cut search covered the full 100–999 space (900 seeds)**. A 150-seed candidate slice within the joined space supplied each HDR band. Offsets were **0/50/100/150/200/250/25** for h0–h6, with wraparound. Each submitted conjunction has **100 Cas12a + 150 Cas9 rows**, and an **11-seed HDR band**. h1's ordinary fallback has **74 Cas12a + 176 Cas9 rows**.

**Clean set:** all 250 uploaded rows cut. **HDR band:** all 250 uploaded rows return HDR. Both are measured by replay over every seed 100–999. The HDR band is a subset of the clean set. For all six submitted conjunction builds, the full-row clean set equals the Cas12a group clean set and matches the build's logged count. h1's uploaded sets are measured independently of its unused prepared build.

| Hotkey | Uploaded construction | Planned candidate slice | Offset | Clean / 900 | Clean inside joined / 300 | Clean outside joined |
|---|---|---|---:|---:|---:|---:|
| h0 | conjunction | 100–249 | 0 | 145 | 62 | 83 |
| h1 | ordinary fallback | 150–299 (unused) | 50 | 5 | 1 | 4 |
| h2 | conjunction | 200–299, 900–949 | 100 | 144 | 56 | 88 |
| h3 | conjunction | 250–299, 900–999 | 150 | 139 | 47 | 92 |
| h4 | conjunction | 100–149, 900–999 | 200 | 142 | 56 | 86 |
| h5 | conjunction | 100–199, 950–999 | 250 | 139 | 46 | 93 |
| h6 | conjunction | 125–274 | 25 | 141 | 61 | 80 |

| Hotkey | Actual uploaded HDR band seeds | Clean scoring-seed hits | HDR scoring-seed hits |
|---|---|---|---|
| h0 | 104, 109, 116, 118, 137, 147, 165, 188, 204, 230, 235 | none | none |
| h1 | none | none | none |
| h2 | 204, 215, 221, 222, 230, 285, 289, 900, 908, 940, 943 | 733 | none |
| h3 | 261, 263, 279, 296, 908, 924, 934, 952, 957, 982, 988 | none | none |
| h4 | 118, 120, 126, 132, 144, 908, 918, 942, 953, 958, 982 | 818 | none |
| h5 | 103, 104, 118, 125, 157, 183, 184, 187, 188, 982, 998 | none | none |
| h6 | 125, 126, 137, 158, 161, 188, 205, 214, 236, 267, 269 | none | none |

Candidate slices are reconstructed from the logged joined space and configured offsets; the actual bands for all six submitted conjunctions agree exactly with both their original band records and the uploaded-row replay.

Fleet union: **578/900 clean seeds** and **52 distinct HDR band seeds** (66 memberships before overlap). Scoring-seed hits in the clean union: **733, 818**. HDR-band hits: **none**.

All three scoring seeds **733/312/818** lie outside the predicted joined band space. The wider cut search still gave **h2 an all-cut hit at 733** and **h4 an all-cut hit at 818**. Seed 312 is outside every submitted clean set. No hotkey hit its HDR band.

## Per-seed scores and outcomes

Each entry is **per-seed final score / no-cut rows**. A clean hit has zero no-cut rows.

| Hotkey | Seed 733 | Seed 312 | Seed 818 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 30.099806 / 4 | 29.985132 / 2 | 27.599660 / 7 | 29.228199 |
| h1 | 28.978379 / 6 | 36.903877 / 3 | 30.349177 / 7 | 32.077144 |
| h2 | 37.362618 / 0 | 33.476818 / 5 | 22.646734 / 5 | 31.162057 |
| h3 | 30.488735 / 2 | 29.307441 / 4 | 29.720914 / 3 | 29.839030 |
| h4 | 31.400497 / 4 | 38.949663 / 6 | 38.131479 / 0 | 36.160546 |
| h5 | 33.126410 / 4 | 37.639079 / 2 | 31.174755 / 4 | 33.980081 |
| h6 | 25.984259 / 3 | 24.722797 / 6 | 26.113760 / 10 | 25.606939 |

## Why the leader scores higher

h4's weighted score is **332.428168**, above the leader's **315.569149**. Its fidelity is somewhat lower (**0.925143 versus 0.969200**), but the largest gap is mean consistency: **0.117579 versus 0.662879**. h4's one clean hit did not produce a competitive three-seed mean. Across h0–h6, mean consistency ranges from **0.083069 to 0.117579**. The public leader record contains aggregate scores; its submission, clean set, HDR band and individual seed outcomes were not available in the inspected data and cannot be inferred from those aggregates.

## h1: prepared build versus actual upload

- **03:07:33.404 UTC:** validator handler finds the prepared build still running and waits up to 179 seconds.
- **03:10:32.981:** wait ends; handler reports the prepared round unusable and starts ordinary construction. The prepared build still holds the GPU, so conjunction/hedge paths are skipped.
- **03:10:34.138:** ordinary construction has generated 250 valid rows; local scoring begins.
- **03:10:38.892:** the prepared conjunction finishes, about six seconds after the fallback started. It reports 140 clean seeds and an 11-seed band, but those are not the rows submitted.
- **03:10:40.447:** ordinary fallback upload succeeds, with 97 seconds of URL TTL remaining.

The unused prepared band is **155, 160, 167, 175, 183, 192, 220, 260, 271, 285, 289**. The actual uploaded fallback has only **5 clean seeds: 189, 517, 742, 765, 810**, and **no all-HDR band**. Its local score reproduces the official **32.077144**. `window_used.json` combines source `all_hdr_not_attempted` with `all_hdr_built: true` and the prepared band; the concurrent prepared build completed after the fallback started. Treating that record's band as the uploaded band would be incorrect.

## Exact uploaded clean sets

These lists cover every seed 100–999 and include the HDR band seeds. Full arrays and per-seed outcomes are also in [analysis.json](analysis.json).

### h0: 145 clean seeds

104, 109, 110, 112, 115, 116, 118, 123, 124, 128, 137, 141, 143, 147, 149, 151, 155, 161, 165, 166, 186, 188, 190, 191, 192, 194, 195, 204, 208, 224, 225, 227, 230, 232, 235, 243, 244, 252, 254, 258, 261, 266, 276, 277, 278, 282, 283, 284, 304, 315, 316, 325, 335, 338, 350, 363, 365, 369, 377, 384, 388, 390, 392, 395, 398, 425, 433, 445, 446, 462, 470, 479, 484, 489, 491, 498, 526, 528, 531, 532, 540, 547, 565, 569, 570, 572, 578, 583, 592, 603, 613, 617, 628, 637, 663, 666, 687, 711, 737, 742, 746, 750, 756, 760, 764, 766, 777, 791, 793, 799, 801, 802, 810, 812, 816, 833, 839, 841, 847, 848, 849, 863, 865, 868, 873, 875, 876, 880, 889, 892, 897, 915, 921, 924, 925, 934, 942, 949, 956, 967, 980, 981, 990, 996, 997.

### h1: 5 clean seeds

189, 517, 742, 765, 810.

### h2: 144 clean seeds

101, 105, 114, 115, 124, 134, 155, 158, 163, 165, 168, 169, 173, 177, 182, 194, 196, 204, 206, 209, 215, 221, 222, 226, 230, 246, 248, 251, 255, 261, 264, 266, 279, 280, 285, 289, 291, 298, 299, 303, 309, 311, 313, 318, 322, 332, 337, 342, 347, 353, 357, 386, 395, 396, 409, 419, 423, 426, 430, 433, 436, 441, 448, 450, 462, 464, 470, 472, 488, 489, 493, 494, 506, 530, 539, 541, 544, 550, 552, 567, 577, 592, 611, 633, 641, 644, 653, 656, 657, 659, 662, 664, 672, 675, 678, 682, 683, 690, 700, 708, 709, 719, 724, 733, 747, 767, 771, 777, 781, 792, 804, 805, 811, 819, 840, 847, 848, 850, 867, 873, 877, 881, 882, 884, 887, 889, 892, 900, 908, 925, 929, 930, 931, 940, 941, 943, 946, 955, 958, 964, 975, 990, 996, 998.

### h3: 139 clean seeds

100, 102, 107, 116, 126, 129, 130, 133, 137, 140, 141, 144, 145, 153, 156, 167, 170, 173, 176, 185, 223, 224, 230, 244, 255, 261, 263, 265, 273, 279, 293, 296, 308, 311, 318, 319, 327, 333, 337, 342, 343, 346, 353, 355, 358, 366, 373, 382, 383, 387, 408, 410, 411, 412, 414, 428, 436, 443, 450, 467, 470, 471, 478, 489, 496, 531, 534, 535, 536, 545, 554, 555, 558, 563, 579, 593, 598, 618, 620, 625, 637, 644, 650, 652, 658, 667, 679, 680, 686, 687, 689, 697, 700, 714, 722, 726, 728, 737, 745, 758, 768, 769, 771, 775, 777, 778, 781, 789, 794, 815, 821, 825, 831, 837, 844, 847, 855, 858, 875, 879, 884, 886, 888, 891, 908, 914, 918, 919, 924, 933, 934, 950, 952, 957, 971, 972, 982, 983, 988.

### h4: 142 clean seeds

103, 108, 118, 120, 122, 126, 132, 134, 136, 137, 138, 144, 155, 158, 167, 169, 176, 203, 205, 231, 241, 242, 252, 261, 262, 263, 266, 275, 285, 294, 297, 298, 299, 302, 310, 321, 326, 328, 333, 343, 348, 350, 353, 382, 389, 405, 416, 424, 426, 430, 431, 439, 440, 447, 452, 453, 465, 466, 468, 470, 473, 495, 497, 499, 501, 506, 514, 515, 523, 526, 530, 534, 535, 539, 541, 544, 548, 550, 577, 584, 586, 601, 608, 611, 615, 620, 633, 639, 652, 653, 656, 668, 670, 675, 677, 699, 712, 716, 731, 741, 746, 780, 782, 787, 794, 798, 812, 815, 818, 823, 826, 835, 837, 844, 847, 885, 892, 894, 896, 900, 908, 913, 914, 916, 918, 934, 939, 941, 942, 943, 945, 950, 953, 954, 955, 958, 974, 982, 990, 991, 992, 997.

### h5: 139 clean seeds

103, 104, 106, 118, 125, 131, 157, 159, 160, 171, 172, 173, 176, 183, 184, 187, 188, 200, 202, 207, 221, 225, 226, 235, 236, 237, 241, 272, 273, 279, 280, 291, 293, 295, 306, 311, 327, 332, 334, 336, 345, 355, 359, 365, 367, 386, 401, 407, 408, 409, 419, 421, 422, 427, 436, 444, 447, 450, 456, 459, 479, 484, 485, 486, 489, 504, 508, 512, 532, 543, 544, 550, 574, 581, 587, 598, 609, 613, 616, 619, 620, 626, 629, 630, 633, 635, 638, 639, 647, 671, 677, 684, 688, 690, 692, 693, 695, 696, 708, 717, 721, 731, 732, 736, 740, 759, 764, 770, 775, 783, 784, 795, 796, 799, 805, 808, 814, 824, 832, 844, 850, 873, 875, 877, 878, 879, 887, 924, 928, 937, 938, 946, 977, 982, 989, 990, 995, 998, 999.

### h6: 141 clean seeds

102, 103, 107, 108, 109, 110, 112, 121, 124, 125, 126, 137, 155, 158, 161, 183, 184, 185, 188, 205, 208, 214, 217, 220, 221, 225, 232, 236, 242, 247, 256, 258, 259, 262, 266, 267, 269, 270, 271, 274, 284, 287, 288, 293, 295, 296, 297, 299, 309, 319, 325, 328, 339, 340, 347, 382, 394, 405, 411, 424, 429, 438, 441, 446, 451, 454, 462, 464, 480, 494, 495, 498, 503, 510, 515, 532, 536, 549, 551, 552, 558, 567, 569, 577, 584, 598, 599, 600, 627, 635, 638, 647, 668, 673, 681, 695, 703, 713, 715, 718, 724, 738, 760, 764, 767, 769, 777, 784, 792, 795, 799, 811, 812, 815, 839, 844, 845, 851, 855, 857, 862, 863, 871, 875, 879, 887, 892, 897, 900, 903, 905, 907, 908, 917, 920, 927, 928, 955, 971, 972, 992.

## Verification and artifacts

- All seven archived upload records confirm submission of this task.
- Archived contracts and references match the published task, allowing only the subsequent scoring-seed stamp.
- All 250 uploaded rows per hotkey pass stages 1/2 and match the valid experiment entries; cell accessibility is 0.87.
- All seven locally calculated three-seed scores match the official API exactly. Existing validation artifacts were reused for h0/h5/h6; the other four were scored with calc.py into this report directory.
- Stage-3 replay matches saved or regenerated outcomes and indel lengths for every row at all three scoring seeds.
- Scanned **900 × 250 × 7 = 1,575,000 row simulations** to recover actual clean sets and HDR bands.
- [Official scores](scores.json), [published task](task.json), [sanitized build/upload logs](build_evidence.json), [archive and validation locations](manifest.json), [summary CSV](summary.csv), [full analysis](analysis.json).
- Reproduce with `.venv/bin/python reports/59b6d75a/audit.py` then `.venv/bin/python reports/59b6d75a/render_report.py` from the repository root. To regenerate a missing validation, run `calc.py --folder <archive> --out-dir reports/59b6d75a/<hotkey> --task reports/59b6d75a/task.json --cell-types reports/59b6d75a/cell_types.json --no-compare --quiet`.

Uploaded archives:

- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — [data/inst/niome_hotkey/result/2026-09-18T03:17:33](../../data/inst/niome_hotkey/result/2026-09-18T03:17:33/submission.json), SHA-256 `b439bb91a5348d6298de517acfdb22ac329cf30184b58375e8576cfaa9e8a46e`.
- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-18T03:11:11](../../data/inst/niome_hotkey1/result/2026-09-18T03:11:11/submission.json), SHA-256 `c07c68b800c14e33f1b9c32c0462a7e547a1757fc7d25db73db16b0fe6244cb9`.
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-18T03:52:42](../../data/inst/niome_hotkey2/result/2026-09-18T03:52:42/submission.json), SHA-256 `acde6e2e2cedb66a2bd336dac242dc590531c603780b81014860e8326acd72f0`.
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-18T03:57:20](../../data/inst/niome_hotkey3/result/2026-09-18T03:57:20/submission.json), SHA-256 `dbfe57f7ec4674975b8ae028c9bc74389cc683ec5eca8571dbd25130a25e2b27`.
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — [data/inst/niome_hotkey4/result/2026-09-18T03:11:06](../../data/inst/niome_hotkey4/result/2026-09-18T03:11:06/submission.json), SHA-256 `f76197e8c843a5d0492a4ba169b1217ddd4036667f0d4160a5387be0131a25b7`.
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — [data/inst/niome_hotkey5/result/2026-09-18T03:40:53](../../data/inst/niome_hotkey5/result/2026-09-18T03:40:53/submission.json), SHA-256 `e0b427a82c2483f2cfc2f2751071a889a7bacbe4fee0f030166a34af25f30d97`.
- **h6** `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` — [data/inst/niome_hotkey6/result/2026-09-18T03:54:07](../../data/inst/niome_hotkey6/result/2026-09-18T03:54:07/submission.json), SHA-256 `62e20a3455521efd77f4aaa5997254b6074f35259ac1ad3a40e9aca29497ed3a`.
