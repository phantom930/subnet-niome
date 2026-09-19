# Submission audit: a31f1224-969b-4195-bb49-8220aa8800e5

**CD34+_HSPC**, created **2026-09-18T15:04:31.395506 UTC**. Scoring seeds: **287, 486, 998**.

**Best of h0–h9: h8, 24.754256, rank 197/248. Leader: UID 250, 131.657702.** h8 achieved 18.80% of the top score, a gap of 106.903447. The top-10 cutoff was 98.820194, 74.065938 above h8.

**Nine successful uploads were found, h0–h8, each with 250 valid rows. h9 has no task upload archive or received-task/upload event in its log; its official score is 0 with 0 valid rows.** All ten hotkeys have zero reward weight. A prepared h9 build exists in the log, but its clean set and band must not be represented as submitted.

h5 uploaded an **ordinary fallback**; h7 uploaded an **on-demand all-HDR** build; h0–h4, h6 and h8 uploaded **conjunctions**. Later rebuilds overwrote several per-task band records. All actual clean sets and bands below come from the uploaded archives.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=a31f1224-969b-4195-bb49-8220aa8800e5&limit=40000), fetched 2026-09-18T19:58:10.909111+00:00. 248 score records, 248 miners, 1 validator. Leader score timestamp: 2026-09-18T17:29:11.146437 UTC. Audit generated 2026-09-18T20:04:02.447552+00:00.

## Official score comparison

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 250 | 1 | 131.657702 | 221.488857 | 0.638682 | 0.930700 | 250 | 0.294000 |
| h0 | 248 | 244 | 20.312937 | 231.753979 | 0.093844 | 0.933987 | 250 | 0.000000 |
| h1 | 41 | 243 | 21.141536 | 232.575571 | 0.097482 | 0.932502 | 250 | 0.000000 |
| h2 | 163 | 234 | 21.860324 | 233.504940 | 0.100325 | 0.933149 | 250 | 0.000000 |
| h3 | 118 | 236 | 21.713498 | 233.861299 | 0.099657 | 0.931672 | 250 | 0.000000 |
| h4 | 190 | 208 | 24.438038 | 234.454268 | 0.112135 | 0.929542 | 250 | 0.000000 |
| h5 | 234 | 231 | 22.427310 | 225.286416 | 0.107107 | 0.929446 | 250 | 0.000000 |
| h6 | 175 | 222 | 23.053654 | 237.351916 | 0.104820 | 0.926622 | 250 | 0.000000 |
| h7 | 65 | 233 | 22.000477 | 230.692318 | 0.100828 | 0.945839 | 250 | 0.000000 |
| h8 | 211 | 197 | 24.754256 | 232.331193 | 0.114101 | 0.933798 | 250 | 0.000000 |
| h9 | 245 | 246 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0 | 0.000000 |

Ranks use competition ranking: 1 plus the count of strictly higher scores. h9 shares the zero-score rank with other zero-score miners; its displayed rank does not imply a valid submission.

Leader hotkey: `5CmEZemqhpTqzNLQbbSFAMwBU1f7DTAap9FdYPSgBffNmE1X`.

## Original joined windows and candidate slices

**h0–h7:** the original 14:17 plan supplied **200–299 ∪ 400–499 ∪ 800–899** (windows **200/400/800**). Seeds **287 and 486** are inside this joined space; 998 is outside. h5 attempted this plan but ultimately uploaded ordinary fallback rows.

**h8/h9:** their logs show the built-in fallback joined layout **100–199 ∪ 300–399 ∪ 800–899**, source `env_pin_no_entry_joined`. All three scoring seeds are outside this space. h8 uploaded its conjunction; h9 only has preparation evidence.

Submitted conjunctions searched for clean seeds over the entire **100–999** range, while their HDR bands came from 150-seed slices of the joined band space. Their rows are **100 Cas12a + 150 Cas9**. h7's all-HDR build used the full joined 300-seed candidate space and the same Cas mix. h5's ordinary fallback is **76 Cas12a + 174 Cas9**.

| Hotkey | Uploaded construction | Joined band space used/planned | Band candidate slice | Offset | Clean /900 |
|---|---|---|---|---:|---:|
| h0 | conjunction | 200–299, 400–499, 800–899 | 200–299, 400–449 | 0 | 137 |
| h1 | conjunction | 200–299, 400–499, 800–899 | 250–299, 400–499 | 50 | 134 |
| h2 | conjunction | 200–299, 400–499, 800–899 | 400–499, 800–849 | 100 | 150 |
| h3 | conjunction | 200–299, 400–499, 800–899 | 450–499, 800–899 | 150 | 139 |
| h4 | conjunction | 200–299, 400–499, 800–899 | 200–249, 800–899 | 200 | 141 |
| h5 | ordinary fallback | 200–299, 400–499, 800–899 | 200–299, 850–899 (unused) | 250 | 9 |
| h6 | conjunction | 200–299, 400–499, 800–899 | 225–299, 400–474 | 25 | 128 |
| h7 | all-HDR | 200–299, 400–499, 800–899 | 200–299, 400–499, 800–899 | full 300 | 15 |
| h8 | conjunction | 100–199, 300–399, 800–899 | 100–189, 840–899 (reconstructed) | 240 | 131 |
| h9 | no upload archive | 100–199, 300–399, 800–899 | 100–199, 300–319, 870–899 (prepared only; reconstructed) | 270 | N/A |

h0–h5's original offsets were **0/50/100/150/200/250**, with **25 for h6**, as supported by the historical stride-50 plan logs. h7's uploaded all-HDR build used all 300 candidates. h8/h9's 240/270 offsets are reconstructed under the new ten-hotkey stride-30 layout; their exact candidate slices were not directly logged. h8's recovered band is consistent with that slice. The exact uploaded bands below are independently established by replay, regardless of that reconstruction.

## Actual clean sets and HDR bands

**Clean set:** all 250 uploaded rows cut. **HDR band:** all 250 uploaded rows return HDR. Both sets were scanned over seeds 100–999. An HDR-band seed is also clean. All seven uploaded conjunctions have exactly 11 HDR seeds; their full-row clean sets equal their Cas12a clean sets and match the original build-log counts.

| Hotkey | Clean /900 | Clean inside joined /300 | Clean outside joined | Actual uploaded HDR band seeds | Clean hits | HDR hits |
|---|---:|---:|---:|---|---|---|
| h0 | 137 | 52 | 85 | 213, 221, 235, 247, 269, 273, 298, 411, 431, 437, 449 | none | none |
| h1 | 134 | 54 | 80 | 253, 295, 401, 405, 411, 431, 437, 442, 461, 465, 470 | none | none |
| h2 | 150 | 61 | 89 | 403, 431, 435, 437, 445, 497, 499, 800, 802, 808, 830 | none | none |
| h3 | 139 | 59 | 80 | 450, 487, 800, 806, 809, 822, 824, 856, 859, 860, 874 | none | none |
| h4 | 141 | 52 | 89 | 209, 218, 231, 248, 249, 800, 817, 859, 860, 874, 876 | 998 | none |
| h5 | 9 | 3 | 6 | none | none | none |
| h6 | 128 | 44 | 84 | 227, 249, 270, 276, 285, 297, 411, 431, 437, 465, 470 | none | none |
| h7 | 15 | 13 | 2 | 212, 262, 282, 292, 299, 405, 412, 434, 457, 807, 816 | none | none |
| h8 | 131 | 56 | 75 | 116, 120, 143, 157, 168, 184, 864, 868, 885, 892, 895 | 486 | none |
| h9 | N/A | N/A | N/A | No uploaded rows available | N/A | N/A |

**h4 hits clean seed 998; h8 hits clean seed 486. No other clean hits and no HDR-band hits occur.** Both clean hits are outside those hotkeys' joined band spaces, demonstrating that the wider cut search provided coverage outside the band prediction. Seed 287 is absent from every uploaded clean set.

Across the **nine actual uploads**, the clean union contains **610/900** seeds and the HDR-band union contains **71 distinct seeds** (88 band memberships before overlap). h9's prepared build is excluded from these unions.

## Per-seed results

Each cell is **per-seed final score / no-cut rows**. Zero no-cut rows denotes a clean hit.

| Hotkey | Seed 287 | Seed 486 | Seed 998 | Mean score |
|---|---:|---:|---:|---:|
| h0 | 19.210719 / 4 | 19.309969 / 9 | 22.418124 / 6 | 20.312937 |
| h1 | 20.527325 / 4 | 22.452762 / 6 | 20.444521 / 5 | 21.141536 |
| h2 | 25.608935 / 11 | 18.080324 / 11 | 21.891713 / 7 | 21.860324 |
| h3 | 21.245614 / 8 | 20.374624 / 7 | 23.520255 / 7 | 21.713498 |
| h4 | 22.477221 / 4 | 20.038201 / 5 | 30.798690 / 0 | 24.438038 |
| h5 | 20.889495 / 3 | 22.249356 / 5 | 24.143077 / 3 | 22.427310 |
| h6 | 21.613544 / 7 | 20.200416 / 2 | 27.347003 / 2 | 23.053654 |
| h7 | 18.967313 / 8 | 25.308619 / 6 | 21.725497 / 5 | 22.000477 |
| h8 | 21.305381 / 4 | 25.672447 / 0 | 27.284938 / 5 | 24.754256 |

## Comparison with the leader

h8's weighted score (**232.331193**) and fidelity (**0.933798**) are above the leader's **221.488857** and **0.930700**. The decisive gap is mean consistency: **0.114101 versus 0.638682**. The leader's submitted rows, clean set, bands and individual seed outcomes were not available in the inspected data and cannot be reconstructed from its aggregate score.

## Build and upload provenance

- **h0–h4 and h6:** prepared conjunctions started around 15:09 and finished around 15:21 UTC, taking 711–719 seconds. Later validator requests used those prepared rows.
- **h5:** the validator handler waited from 15:16:29 to 15:19:28, then started ordinary construction because the prepared build still held the GPU. It uploaded the fallback at **15:19:32**, with 101 seconds of TTL remaining. Its prepared conjunction finished at 15:21:05 and was not submitted.
- **h7:** at 15:06:45 its prepared cache was for another task, so it built this task on demand. It completed an **all-HDR** build and uploaded at **15:09:45**, with 103 seconds of TTL remaining. Its later prepared conjunction was not the upload.
- **h8:** it prepared against the fallback joined layout at 15:11:30, finished at 15:21:17, and uploaded that prepared conjunction at 15:44:40.
- **h9:** its log records a prepared conjunction from 15:11:30 to 15:21:14, but no received-task or successful upload event was found for this task. There is no matching upload archive. The API records **0 valid rows and score 0**. The reason the validator received no valid rows is not established by these records.

### Later rebuilds and overwritten band records

At approximately **15:50–15:56 UTC**, h0–h7 rebuilt the same task under the newer stride-30 layout. This happened after every archived upload. The current `window_used.json` entries for these hotkeys therefore refer to later builds; h0 retained the same band because its offset remained zero. h1–h7's recorded bands differ from their actual uploaded bands. h8's record still matches its upload.

| Hotkey | Latest task-record time (UTC) | Latest recorded band matches actual upload? |
|---|---|---|
| h0 | 2026-09-18T15:50:23.518010+00:00 | yes |
| h1 | 2026-09-18T15:50:25.513669+00:00 | NO |
| h2 | 2026-09-18T15:50:26.942282+00:00 | NO |
| h3 | 2026-09-18T15:50:29.066958+00:00 | NO |
| h4 | 2026-09-18T15:50:30.933801+00:00 | NO |
| h5 | 2026-09-18T15:50:32.015575+00:00 | NO |
| h6 | 2026-09-18T15:50:33.288500+00:00 | NO |
| h7 | 2026-09-18T15:50:35.243741+00:00 | NO |
| h8 | 2026-09-18T15:11:30.310311+00:00 | yes |

### h9 preparation only — not submitted evidence

h9's prepared build logged **134/900 clean seeds** and the recorded band **116, 123, 143, 157, 158, 163, 168, 182, 885, 892, 895**. These describe preparation only. No uploaded h9 rows are available to recover its actual submitted clean set or band, so the main tables use N/A.

## Exact uploaded clean sets

The lists include HDR-band seeds and cover the full 100–999 seed space. There is no reconstructable submitted clean set for h9.

### h0: 137 clean seeds

101, 109, 111, 115, 126, 130, 138, 143, 146, 149, 151, 152, 159, 160, 192, 195, 196, 200, 201, 204, 213, 221, 222, 223, 228, 235, 247, 258, 259, 269, 273, 274, 275, 281, 292, 298, 312, 315, 325, 328, 329, 332, 336, 353, 359, 362, 370, 375, 376, 377, 385, 389, 394, 395, 400, 404, 410, 411, 424, 427, 431, 433, 437, 441, 446, 449, 456, 457, 469, 488, 489, 491, 492, 505, 519, 521, 522, 541, 544, 545, 558, 560, 572, 574, 579, 625, 656, 664, 672, 673, 693, 696, 714, 725, 726, 731, 735, 751, 755, 758, 759, 760, 767, 771, 772, 790, 797, 803, 808, 829, 830, 856, 859, 860, 861, 862, 864, 871, 875, 880, 899, 903, 906, 911, 918, 923, 924, 941, 957, 967, 969, 972, 974, 979, 982, 989, 996.

### h1: 134 clean seeds

126, 127, 135, 140, 153, 164, 165, 167, 193, 198, 201, 202, 206, 209, 220, 232, 233, 237, 242, 253, 257, 295, 303, 310, 320, 324, 327, 331, 332, 342, 353, 356, 373, 377, 401, 402, 404, 405, 407, 411, 417, 427, 428, 431, 437, 439, 441, 442, 448, 461, 464, 465, 467, 470, 475, 478, 483, 491, 493, 498, 501, 506, 508, 510, 519, 527, 530, 532, 542, 550, 553, 556, 557, 558, 566, 571, 579, 582, 595, 598, 599, 600, 609, 611, 632, 639, 652, 658, 670, 674, 679, 682, 683, 693, 701, 709, 724, 727, 735, 739, 756, 762, 779, 784, 792, 807, 812, 821, 826, 834, 847, 855, 863, 865, 866, 881, 885, 886, 889, 893, 896, 910, 917, 926, 931, 933, 943, 952, 963, 965, 967, 968, 974, 979.

### h2: 150 clean seeds

103, 107, 116, 133, 136, 150, 156, 173, 175, 187, 200, 209, 214, 221, 226, 227, 239, 242, 252, 254, 255, 265, 266, 274, 277, 278, 279, 291, 301, 303, 308, 314, 317, 328, 332, 335, 352, 356, 358, 361, 369, 374, 376, 403, 405, 411, 414, 427, 431, 434, 435, 437, 442, 445, 451, 453, 460, 463, 465, 470, 485, 497, 499, 512, 513, 523, 530, 534, 543, 548, 556, 560, 571, 572, 574, 584, 585, 590, 600, 601, 609, 619, 620, 634, 636, 642, 643, 654, 656, 659, 670, 679, 685, 688, 696, 703, 710, 711, 714, 725, 727, 731, 735, 740, 743, 749, 753, 754, 755, 766, 770, 771, 780, 782, 800, 802, 808, 811, 815, 817, 821, 822, 824, 830, 835, 841, 852, 856, 859, 864, 866, 872, 875, 882, 884, 893, 896, 903, 904, 910, 919, 939, 940, 959, 961, 968, 979, 982, 984, 993.

### h3: 139 clean seeds

102, 104, 106, 109, 111, 123, 125, 128, 131, 134, 137, 139, 142, 144, 157, 160, 163, 166, 172, 176, 194, 198, 199, 203, 204, 209, 212, 213, 214, 218, 222, 231, 232, 253, 255, 258, 259, 260, 265, 278, 283, 291, 299, 309, 318, 320, 325, 328, 338, 340, 367, 368, 381, 396, 399, 405, 424, 431, 434, 441, 443, 445, 450, 455, 456, 468, 472, 480, 487, 488, 489, 494, 530, 533, 538, 545, 548, 550, 578, 600, 604, 606, 609, 613, 617, 625, 652, 658, 666, 668, 670, 681, 686, 696, 700, 701, 702, 705, 708, 727, 742, 749, 763, 766, 778, 781, 792, 800, 806, 809, 810, 815, 822, 824, 825, 848, 849, 853, 855, 856, 859, 860, 868, 870, 874, 891, 893, 894, 895, 909, 911, 934, 935, 941, 946, 961, 971, 979, 987.

### h4: 141 clean seeds

100, 107, 111, 115, 118, 132, 135, 141, 142, 151, 157, 158, 170, 171, 180, 183, 185, 193, 201, 209, 215, 218, 225, 230, 231, 238, 248, 249, 251, 254, 263, 270, 273, 278, 279, 285, 286, 293, 303, 304, 308, 309, 310, 317, 320, 330, 334, 348, 354, 361, 369, 371, 372, 384, 385, 399, 400, 401, 405, 407, 417, 422, 427, 431, 434, 451, 466, 478, 492, 493, 523, 540, 557, 558, 569, 571, 577, 579, 581, 591, 595, 597, 599, 604, 607, 620, 626, 630, 635, 662, 672, 675, 679, 689, 696, 697, 704, 708, 714, 717, 718, 725, 731, 739, 752, 755, 760, 766, 772, 777, 781, 785, 788, 789, 790, 795, 800, 803, 810, 817, 825, 828, 830, 835, 836, 847, 852, 853, 859, 860, 864, 874, 876, 887, 915, 919, 930, 934, 939, 943, 998.

### h5: 9 clean seeds

320, 417, 447, 499, 503, 646, 657, 768, 788.

### h6: 128 clean seeds

101, 102, 108, 117, 119, 135, 173, 178, 182, 185, 190, 191, 194, 196, 200, 222, 225, 227, 230, 232, 239, 247, 249, 251, 252, 270, 274, 276, 280, 282, 284, 285, 293, 297, 298, 305, 308, 309, 313, 329, 343, 345, 350, 353, 359, 360, 369, 374, 377, 381, 393, 405, 409, 411, 419, 421, 422, 431, 435, 437, 465, 470, 483, 485, 487, 496, 509, 513, 520, 530, 537, 541, 543, 547, 556, 558, 559, 561, 563, 566, 568, 576, 596, 602, 606, 613, 621, 627, 630, 649, 675, 682, 683, 713, 714, 720, 731, 739, 745, 751, 752, 753, 762, 763, 781, 782, 789, 792, 794, 796, 803, 811, 813, 831, 864, 878, 879, 884, 951, 954, 956, 959, 960, 962, 973, 990, 992, 994.

### h7: 15 clean seeds

212, 262, 282, 292, 299, 405, 412, 425, 434, 457, 577, 807, 816, 835, 911.

### h8: 131 clean seeds

106, 116, 120, 121, 127, 134, 138, 139, 143, 145, 156, 157, 159, 163, 168, 173, 175, 177, 179, 180, 181, 183, 184, 192, 207, 213, 231, 234, 245, 246, 256, 259, 265, 284, 286, 297, 306, 316, 318, 319, 337, 341, 348, 353, 366, 368, 372, 373, 393, 397, 398, 411, 412, 426, 439, 478, 479, 484, 486, 507, 513, 524, 530, 540, 569, 570, 581, 593, 598, 599, 601, 604, 616, 623, 626, 628, 632, 639, 643, 653, 663, 668, 673, 679, 681, 682, 683, 686, 688, 694, 696, 701, 716, 733, 738, 748, 761, 766, 775, 786, 793, 796, 799, 806, 808, 810, 816, 821, 823, 827, 864, 865, 868, 878, 885, 886, 888, 889, 892, 895, 906, 925, 927, 929, 939, 949, 951, 977, 979, 992, 995.

## Verification and artifacts

- All nine archived uploads identify this task. Their contracts and references match the public task after accounting for its later seed stamp.
- All nine submissions contain 250 stage-1/2 valid rows, with cell accessibility 0.87, and their local three-seed scores match the API exactly.
- Saved validations were reused for h0–h3, h7 and h8; h4–h6 were validated with calc.py into this report directory.
- Stage-3 outcomes and indel lengths match saved/regenerated detail at all three scoring seeds for every uploaded row.
- Recovered sets by evaluating **900 × 250 × 9 = 2,025,000 row/seed combinations**; h9 was not simulated without uploaded rows.
- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized miner logs](build_evidence.json), [historical layout logs](plan_history.json), [archive manifest](manifest.json), [summary CSV](summary.csv).
- Reproduce from the repository root with `.venv/bin/python reports/a31f1224/validate.py`, `.venv/bin/python reports/a31f1224/audit.py`, and `.venv/bin/python reports/a31f1224/render_report.py`.

Uploaded archives:

- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — [data/inst/niome_hotkey/result/2026-09-18T15:44:53](../../data/inst/niome_hotkey/result/2026-09-18T15:44:53/submission.json); SHA-256 `b168e2c4bc83efda522384ccafacdc6e8902588529dc4ba599a4b1db729bcac7`.
- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-18T15:35:02](../../data/inst/niome_hotkey1/result/2026-09-18T15:35:02/submission.json); SHA-256 `8390ea4fbac721929b3ae00d8114df3d88a59f304af94eb38da797c8e47f28cd`.
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-18T15:49:31](../../data/inst/niome_hotkey2/result/2026-09-18T15:49:31/submission.json); SHA-256 `e255fb6d178377e0d71919bb321907b78f42bb3c6c4242c34b16d7220f271bf4`.
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-18T15:29:45](../../data/inst/niome_hotkey3/result/2026-09-18T15:29:45/submission.json); SHA-256 `ce18eb3fda285bc89f5db4870579b5143e26b7cf52950fba941ab407ba84ccfb`.
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — [data/inst/niome_hotkey4/result/2026-09-18T15:22:48](../../data/inst/niome_hotkey4/result/2026-09-18T15:22:48/submission.json); SHA-256 `4445f3ef3da200999ff91ad2fbca8235ae77c75e62a7aa30c233d77c81978f5e`.
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — [data/inst/niome_hotkey5/result/2026-09-18T15:20:01](../../data/inst/niome_hotkey5/result/2026-09-18T15:20:01/submission.json); SHA-256 `dd5e321af51557f04055fc67562b9b1fb449bb8d12ce6418c9cc1a47f2af7329`.
- **h6** `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` — [data/inst/niome_hotkey6/result/2026-09-18T15:43:56](../../data/inst/niome_hotkey6/result/2026-09-18T15:43:56/submission.json); SHA-256 `d06da2f02e801d3eac9274a2b9dfb2ad6c3a2c4a1bfee475b8df73623731087b`.
- **h7** `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` — [data/inst/niome_hotkey7/result/2026-09-18T15:10:14](../../data/inst/niome_hotkey7/result/2026-09-18T15:10:14/submission.json); SHA-256 `45282bf5e919b4bc8a9eece4a2d82dcccf57e7856252e11075d22aa59eafd01d`.
- **h8** `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` — [data/inst/niome_hotkey8/result/2026-09-18T15:45:09](../../data/inst/niome_hotkey8/result/2026-09-18T15:45:09/submission.json); SHA-256 `9acfb836fc4f6b6bbc8fc50e92a5e146854a39fdb9cc9c953c265f59433aa20b`.
- **h9** `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` — no upload archive for this task.
