# Submission audit: cef21622-627f-4171-bb83-53bfc44b4f91

CD34+_HSPC · task opened 2026-09-17T22:08:27.490979 UTC · scoring seeds **115, 977, 940**.

**Best of h0–h5: h0, 29.413795, rank 156/248. Top miner: UID 199, 151.991292.** h0 reached 19.35% of the top score, a gap of 122.577497 points. It missed the top-10 cutoff (102.205275) by **72.791480** points. All six submissions were uploaded successfully, contained 250 valid rows, and received zero reward weight.

Scores are from the [task-specific public API](https://niome-api.genomes.io/api/v3/miners/scores?task_id=cef21622-627f-4171-bb83-53bfc44b4f91&limit=40000); 248 records, 248 distinct miners, 1 validator. Score timestamp: 2026-09-18T00:34:29.201881 UTC. Audit generated: 2026-09-18T04:13:37.700347+00:00.

## Official score comparison

| Miner | UID this round | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 199 | 1 | 151.991292 | 256.752820 | 0.631912 | 0.936800 | 0.294000 |
| h0 | 248 | 156 | 29.413795 | 249.058314 | 0.126023 | 0.937133 | 0.000000 |
| h1 | 41 | 214 | 24.413403 | 250.176829 | 0.103971 | 0.938576 | 0.000000 |
| h2 | 163 | 216 | 24.195833 | 242.272932 | 0.105066 | 0.950549 | 0.000000 |
| h3 | 118 | 241 | 19.795969 | 253.146155 | 0.083607 | 0.935324 | 0.000000 |
| h4 | 190 | 232 | 23.024480 | 247.052635 | 0.098987 | 0.941503 | 0.000000 |
| h5 | 234 | 211 | 24.843243 | 248.976213 | 0.106005 | 0.941291 | 0.000000 |

Top miner hotkey: `5DFXpWVqdiz1WSYiAHgryVPLacmJCozh6RvyPpfKPz5zBVM2`.

Final score is the mean of the per-seed scores. For these submissions weighted score and fidelity are constant across seeds, so final = weighted score × mean consistency × fidelity.

## Joined seed windows and recovered bands

All six original builds used joined band space **300–499 ∪ 800–899** (300 seeds), from the plan generated at 2026-09-17 21:20 UTC. The separate cut search covered **100–999** (900 seeds). Band candidate slices were 150 seeds wide, at offsets 0/50/100/150/200/250 for h0–h5. Every submission has 100 Cas12a and 150 Cas9 rows and an 11-seed HDR band.

**Clean set:** every submitted row cuts at that seed. **HDR band:** every submitted row returns HDR. The band is a subset of the clean set. Both sets were recovered by simulating the exact uploaded rows at every seed 100–999; the recovered bands match the original recorded bands exactly. The full-row clean sets equal the Cas12a group clean sets, and their sizes match the build logs.

| Hotkey | Candidate slice within joined space | Clean / 900 | Actual HDR band seeds |
|---|---|---:|---|
| h0 | 300–449 | 143 | 344, 348, 350, 364, 372, 397, 400, 403, 412, 427, 428 |
| h1 | 350–499 | 137 | 352, 403, 412, 427, 447, 449, 450, 453, 471, 486, 499 |
| h2 | 400–499, 800–849 | 140 | 403, 406, 412, 417, 427, 450, 453, 456, 470, 471, 804 |
| h3 | 450–499, 800–899 | 142 | 450, 461, 465, 806, 816, 817, 818, 832, 852, 880, 894 |
| h4 | 300–349, 800–899 | 134 | 300, 310, 316, 333, 824, 829, 841, 852, 885, 886, 890 |
| h5 | 300–399, 850–899 | 133 | 314, 377, 385, 392, 858, 860, 867, 873, 874, 876, 898 |

Candidate slices are reconstructed from the original logged joined space and the six-hotkey offset layout; every recovered band lies inside its corresponding slice.

Fleet unions: **564 distinct clean seeds** and **55 distinct HDR band seeds** (66 band memberships before overlap). The scoring seeds covered by the fleet's clean union are 115; the band union hits none.

## Hits and per-seed scores

| Hotkey | Clean seed hits (includes HDR hits) | HDR band hits |
|---|---|---|
| h0 | 115 | none |
| h1 | none | none |
| h2 | none | none |
| h3 | none | none |
| h4 | none | none |
| h5 | none | none |

All three scoring seeds, **115, 977, and 940**, are outside the joined band space. None hits any hotkey's HDR band. The wider 900-seed cut search supplies exactly one clean hit among these six hotkeys: **h0 at seed 115**. Neither 977 nor 940 is clean for any hotkey.

Each cell below is **per-seed final score / no-cut rows**.

| Hotkey | Seed 115 | Seed 977 | Seed 940 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 42.100869 / 0 | 22.283709 / 10 | 23.856806 / 4 | 29.413795 |
| h1 | 26.148313 / 5 | 20.853366 / 6 | 26.238531 / 5 | 24.413403 |
| h2 | 21.102905 / 9 | 24.798287 / 3 | 26.686307 / 3 | 24.195833 |
| h3 | 20.990043 / 4 | 18.875818 / 8 | 19.522045 / 7 | 19.795969 |
| h4 | 22.381399 / 7 | 25.966480 / 2 | 20.725563 / 3 | 23.024480 |
| h5 | 23.929520 / 6 | 24.275411 / 7 | 26.324799 / 6 | 24.843243 |

h0's clean hit at seed 115 scores **42.100869**, but it is not an HDR band hit. Its scores at 977 and 940 are **22.283709** and **23.856806**, respectively, giving a mean of **29.413795**. All remaining hotkey/seed combinations contain no-cut rows.

## Comparison with the leader

h0's fidelity is essentially equal to the leader's (0.937133 versus 0.936800), and its weighted score is only modestly lower (**249.058314 versus 256.752820**). The dominant gap is mean consistency: **0.126023 versus 0.631912**. Across all six hotkeys, consistency ranges from 0.083607 to 0.126023. The observed clean or band hits explain the fleet's per-seed outcomes on this task; one task does not establish long-run strategy performance. The leader's clean sets, bands, and per-seed outcomes are not available from its aggregate score and were not reconstructed.

## Exact clean sets

These lists include the HDR band seeds and cover the full 100–999 seed space. Machine-readable lists are also in [analysis.json](analysis.json).

### h0: 143 clean seeds

100, 115, 127, 128, 133, 140, 143, 144, 150, 170, 173, 178, 179, 180, 192, 194, 199, 202, 205, 213, 215, 219, 221, 223, 224, 226, 231, 239, 242, 245, 251, 252, 254, 256, 258, 264, 278, 281, 292, 293, 295, 298, 302, 336, 339, 344, 347, 348, 350, 358, 364, 370, 372, 377, 379, 385, 397, 399, 400, 403, 412, 417, 427, 428, 433, 434, 439, 444, 449, 456, 464, 473, 478, 483, 492, 495, 521, 522, 529, 535, 556, 557, 561, 569, 576, 577, 588, 592, 599, 617, 620, 625, 642, 654, 655, 657, 658, 660, 669, 676, 683, 686, 691, 693, 695, 698, 710, 714, 724, 730, 733, 739, 744, 745, 748, 749, 755, 757, 760, 777, 783, 786, 787, 806, 808, 830, 841, 842, 859, 868, 875, 876, 893, 910, 913, 924, 925, 936, 941, 947, 950, 965, 976.

### h1: 137 clean seeds

104, 116, 118, 119, 140, 146, 147, 149, 154, 160, 166, 169, 171, 172, 184, 195, 196, 202, 203, 246, 252, 262, 264, 265, 266, 277, 280, 281, 282, 287, 291, 295, 307, 308, 314, 315, 317, 318, 322, 323, 324, 326, 332, 335, 338, 341, 343, 352, 354, 358, 374, 377, 403, 408, 412, 424, 427, 447, 449, 450, 453, 456, 458, 464, 466, 471, 476, 478, 482, 484, 486, 499, 500, 502, 507, 510, 513, 515, 539, 560, 561, 563, 578, 586, 595, 610, 620, 629, 660, 662, 665, 667, 670, 695, 703, 711, 721, 729, 735, 737, 739, 747, 763, 765, 767, 768, 777, 786, 787, 797, 806, 809, 813, 828, 840, 842, 843, 844, 865, 879, 880, 891, 894, 896, 897, 899, 906, 908, 913, 935, 946, 947, 948, 957, 988, 996, 997.

### h2: 140 clean seeds

105, 112, 121, 129, 138, 141, 152, 171, 180, 184, 191, 193, 200, 207, 209, 210, 211, 220, 225, 237, 262, 264, 276, 285, 301, 306, 307, 322, 326, 339, 342, 346, 371, 395, 401, 403, 406, 407, 410, 412, 417, 427, 448, 450, 452, 453, 456, 459, 467, 469, 470, 471, 474, 475, 482, 494, 501, 511, 520, 522, 524, 536, 541, 543, 553, 555, 559, 560, 562, 564, 577, 579, 583, 584, 590, 591, 597, 598, 600, 601, 611, 627, 628, 632, 637, 644, 646, 649, 670, 679, 683, 689, 697, 726, 743, 745, 754, 764, 773, 776, 781, 792, 794, 796, 804, 809, 814, 819, 826, 831, 841, 845, 850, 855, 856, 862, 871, 873, 874, 875, 883, 885, 916, 931, 937, 943, 944, 952, 955, 961, 964, 971, 973, 979, 982, 983, 985, 989, 990, 994.

### h3: 142 clean seeds

127, 128, 131, 132, 140, 142, 144, 168, 171, 174, 183, 194, 201, 219, 220, 227, 228, 229, 234, 239, 246, 254, 256, 264, 270, 275, 281, 285, 288, 298, 302, 308, 315, 326, 330, 337, 340, 343, 360, 370, 375, 376, 381, 401, 407, 412, 414, 415, 421, 422, 435, 440, 442, 443, 448, 450, 461, 465, 478, 490, 498, 502, 503, 505, 507, 509, 519, 541, 542, 544, 545, 554, 560, 564, 566, 568, 573, 579, 581, 584, 595, 614, 617, 637, 656, 664, 667, 669, 670, 673, 675, 680, 681, 682, 703, 705, 713, 714, 731, 740, 766, 773, 775, 783, 787, 801, 806, 816, 817, 818, 819, 823, 829, 832, 834, 847, 852, 858, 860, 861, 880, 881, 882, 894, 905, 906, 911, 912, 915, 917, 923, 928, 934, 949, 953, 966, 969, 972, 980, 988, 990, 991.

### h4: 134 clean seeds

109, 112, 126, 129, 130, 131, 132, 147, 149, 155, 157, 160, 166, 175, 183, 193, 204, 205, 207, 208, 218, 228, 233, 234, 241, 250, 264, 267, 269, 274, 300, 310, 312, 316, 317, 318, 329, 333, 339, 341, 342, 363, 366, 372, 378, 380, 381, 395, 397, 398, 413, 426, 427, 438, 444, 449, 463, 470, 478, 488, 502, 503, 504, 510, 522, 523, 539, 544, 547, 552, 555, 557, 566, 569, 578, 580, 585, 590, 595, 625, 626, 632, 646, 653, 658, 660, 682, 683, 700, 729, 739, 747, 748, 750, 759, 768, 783, 816, 819, 824, 829, 830, 839, 841, 852, 858, 859, 865, 867, 872, 873, 885, 886, 887, 890, 894, 901, 908, 909, 910, 913, 919, 923, 931, 932, 936, 946, 951, 952, 966, 981, 984, 990, 992.

### h5: 133 clean seeds

102, 109, 116, 120, 131, 134, 138, 145, 153, 157, 160, 164, 168, 179, 198, 201, 206, 210, 211, 214, 215, 229, 234, 272, 275, 289, 292, 314, 316, 324, 326, 336, 348, 358, 362, 377, 379, 380, 385, 391, 392, 403, 404, 407, 413, 417, 418, 419, 421, 422, 431, 437, 439, 452, 456, 457, 463, 476, 486, 494, 496, 502, 504, 521, 530, 531, 532, 536, 540, 552, 554, 562, 564, 568, 580, 583, 586, 590, 595, 597, 600, 619, 620, 622, 629, 631, 635, 641, 652, 657, 659, 668, 674, 699, 707, 714, 720, 721, 744, 772, 783, 789, 791, 795, 800, 812, 815, 817, 821, 824, 846, 858, 860, 861, 866, 867, 873, 874, 876, 885, 898, 906, 923, 928, 939, 943, 951, 970, 972, 986, 987, 989, 996.

## Submission provenance and verification

The six builds started at 22:08:37 UTC and finished between 22:10:56 and 22:11:06 UTC (139–149 seconds). All six hotkeys served the validator's requests from prepared submissions. Uploads succeeded between 22:13:18 and 22:57:58 UTC, with 266–288 seconds of upload TTL remaining. Each hotkey has one recorded build and one recorded successful upload for this task. The recorded bands agree with the archived rows' full-seed replay.

- All archived contracts and reference files match the published task after accounting for the seed changing from 0 at submission time to the three stamped scoring seeds.
- All 250 archived stage-1/2 valid-row entries match each submission exactly; accessibility is 0.87.
- All six archived three-seed final scores match the official API exactly.
- Replayed stage-3 outcomes and indel lengths match the archived detail for every row at all three real seeds.
- Scanned 900 seeds × 250 rows × 6 hotkeys = 1,350,000 row simulations. No miner configuration, running process, or uploaded submission was changed.
- [Full analysis](analysis.json), [official score snapshot](scores.json), [published task](task.json), [sanitized build/upload evidence](build_evidence.json).
- Reproduce with `.venv/bin/python reports/cef21622/audit.py` and `.venv/bin/python reports/cef21622/render_report.py` from the repository root.

Uploaded archives:

- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — [data/inst/niome_hotkey/result/2026-09-17T22:22:52](../../data/inst/niome_hotkey/result/2026-09-17T22:22:52/submission.json)
- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-17T22:13:48](../../data/inst/niome_hotkey1/result/2026-09-17T22:13:48/submission.json)
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-17T22:51:31](../../data/inst/niome_hotkey2/result/2026-09-17T22:51:31/submission.json)
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-17T22:34:45](../../data/inst/niome_hotkey3/result/2026-09-17T22:34:45/submission.json)
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — [data/inst/niome_hotkey4/result/2026-09-17T22:55:42](../../data/inst/niome_hotkey4/result/2026-09-17T22:55:42/submission.json)
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — [data/inst/niome_hotkey5/result/2026-09-17T22:58:26](../../data/inst/niome_hotkey5/result/2026-09-17T22:58:26/submission.json)
