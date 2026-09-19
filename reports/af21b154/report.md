# Submission audit: af21b154-5ee4-4c34-98cb-ebd58f04cdc2

K562 · task opened 2026-09-17T19:44:49.986698 UTC · scoring seeds **391, 292, 785**.

**Best of h0–h5: h5, 123.602543, rank 11/248. Top miner: UID 157, 191.740106.** h5 reached 64.46% of the top score, a gap of 68.137563 points. It missed the top-10 cutoff (131.974875) by **8.372332** points. All six submissions were uploaded successfully, contained 250 valid rows, and received zero reward weight.

Scores are from the [task-specific public API](https://niome-api.genomes.io/api/v3/miners/scores?task_id=af21b154-5ee4-4c34-98cb-ebd58f04cdc2&limit=40000); 248 records, 248 distinct miners, 1 validator. Score timestamp: 2026-09-17T22:07:19.678314 UTC. Audit generated: 2026-09-18T04:08:12.493498+00:00.

## Official score comparison

| Miner | UID this round | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 157 | 1 | 191.740106 | 355.007471 | 0.578825 | 0.933100 | 0.294000 |
| h0 | 248 | 160 | 40.336034 | 324.205304 | 0.133145 | 0.934433 | 0.000000 |
| h1 | 41 | 159 | 48.000694 | 322.194357 | 0.159074 | 0.936547 | 0.000000 |
| h2 | 163 | 215 | 32.049747 | 322.411575 | 0.105892 | 0.938748 | 0.000000 |
| h3 | 118 | 207 | 33.031917 | 326.391945 | 0.108434 | 0.933317 | 0.000000 |
| h4 | 190 | 87 | 62.939083 | 327.924882 | 0.206180 | 0.930894 | 0.000000 |
| h5 | 234 | 11 | 123.602543 | 322.083998 | 0.410924 | 0.933892 | 0.000000 |

Top miner hotkey: `5EvPaqeUKZXX3NC36o6AC9tJPn9Q25nipmxEzyCck2D45s77`.

Final score is the mean of the per-seed scores. For these submissions weighted score and fidelity are constant across seeds, so final = weighted score × mean consistency × fidelity.

## Joined seed windows and recovered bands

All six original builds used joined band space **200–299 ∪ 800–999** (300 seeds), from the plan generated at 2026-09-17 19:17 UTC. The separate cut search covered **100–999** (900 seeds). Band candidate slices were 150 seeds wide, at offsets 0/50/100/150/200/250 for h0–h5. Every submission has 100 Cas12a and 150 Cas9 rows and an 11-seed HDR band.

**Clean set:** every submitted row cuts at that seed. **HDR band:** every submitted row returns HDR. The band is a subset of the clean set. Both sets were recovered by simulating the exact uploaded rows at every seed 100–999; the recovered bands match the original recorded bands exactly. The full-row clean sets equal the Cas12a group clean sets, and their sizes match the build logs.

| Hotkey | Candidate slice within joined space | Clean / 900 | Actual HDR band seeds |
|---|---|---:|---|
| h0 | 200–299, 800–849 | 138 | 218, 233, 245, 274, 288, 298, 806, 809, 815, 839, 845 |
| h1 | 250–299, 800–899 | 135 | 280, 285, 288, 295, 800, 826, 832, 857, 877, 879, 880 |
| h2 | 800–949 | 135 | 808, 827, 844, 846, 848, 859, 870, 873, 874, 919, 924 |
| h3 | 850–999 | 135 | 875, 882, 893, 916, 924, 925, 927, 930, 931, 956, 974 |
| h4 | 200–249, 900–999 | 137 | 205, 215, 226, 230, 233, 902, 911, 924, 936, 938, 996 |
| h5 | 200–299, 950–999 | 139 | 214, 222, 233, 264, 292, 953, 960, 978, 988, 991, 995 |

Candidate slices are reconstructed from the original logged joined space and the six-hotkey offset layout; every recovered band lies inside its corresponding slice.

Fleet unions: **579 distinct clean seeds** and **61 distinct HDR band seeds** (66 band memberships before overlap). The scoring seeds covered by the fleet's clean union are 292, 391, 785; the band union hits 292.

## Hits and per-seed scores

| Hotkey | Clean seed hits (includes HDR hits) | HDR band hits |
|---|---|---|
| h0 | 292 | none |
| h1 | 292, 785 | none |
| h2 | none | none |
| h3 | none | none |
| h4 | 292, 391 | none |
| h5 | 292 | 292 |

Seed **292** falls in the candidate windows of h0/h1/h5; only **h5** pins it as an HDR band seed. Seeds **391 and 785** are outside the joined band space. The wider cut search nevertheless supplies a clean hit at 391 for h4 and at 785 for h1.

Each cell below is **per-seed final score / no-cut rows**.

| Hotkey | Seed 391 | Seed 292 | Seed 785 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 35.766292 / 3 | 59.907309 / 0 | 25.334502 / 7 | 40.336034 |
| h1 | 31.949246 / 3 | 53.554270 / 0 | 58.498567 / 0 | 48.000694 |
| h2 | 33.046770 / 5 | 29.974345 / 5 | 33.128126 / 6 | 32.049747 |
| h3 | 26.600275 / 5 | 37.440237 / 2 | 35.055240 / 3 | 33.031917 |
| h4 | 86.514238 / 0 | 73.372271 / 0 | 28.930739 / 6 | 62.939083 |
| h5 | 34.119542 / 2 | 300.791678 / 0 | 35.896408 / 6 | 123.602543 |

h5's seed-292 outcome is **250/250 HDR rows**, giving consistency 1.0 and score **300.791678** on that seed. Its scores on 391 and 785 are **34.119542** and **35.896408**, so the three-seed mean is **123.602543**. The band hit worked; the other two seeds kept the round score below the top-10 cutoff.

## Comparison with the leader

h5's fidelity is essentially equal to the leader's (0.933892 versus 0.933100). The gap comes from lower mean consistency (**0.410924 versus 0.578825**) and lower weighted score (**322.083998 versus 355.007471**). The observed clean or band hits explain the fleet's per-seed outcomes on this task; one task does not establish long-run strategy performance. The leader's clean sets, bands, and per-seed outcomes are not available from its aggregate score and were not reconstructed.

## Exact clean sets

These lists include the HDR band seeds and cover the full 100–999 seed space. Machine-readable lists are also in [analysis.json](analysis.json).

### h0: 138 clean seeds

105, 107, 109, 110, 112, 118, 119, 127, 145, 153, 172, 191, 196, 197, 213, 214, 218, 230, 233, 234, 245, 250, 265, 274, 278, 286, 287, 288, 292, 297, 298, 310, 314, 316, 317, 323, 331, 335, 347, 350, 365, 367, 382, 384, 386, 388, 389, 396, 404, 408, 409, 451, 456, 462, 482, 484, 485, 486, 487, 497, 524, 527, 539, 542, 547, 548, 555, 561, 565, 573, 586, 591, 601, 622, 624, 625, 631, 633, 649, 650, 662, 668, 671, 679, 686, 687, 697, 715, 723, 739, 744, 761, 766, 775, 777, 789, 790, 796, 806, 809, 812, 813, 814, 815, 821, 839, 840, 845, 848, 852, 853, 858, 861, 864, 868, 874, 878, 879, 895, 897, 909, 910, 918, 921, 930, 933, 936, 947, 948, 949, 958, 959, 960, 964, 987, 993, 994, 999.

### h1: 135 clean seeds

100, 109, 117, 118, 123, 124, 146, 147, 150, 169, 179, 196, 201, 211, 223, 225, 234, 235, 246, 248, 256, 260, 266, 280, 282, 285, 288, 290, 292, 293, 294, 295, 296, 301, 313, 318, 319, 320, 333, 353, 358, 359, 385, 390, 400, 413, 421, 423, 431, 436, 442, 443, 463, 470, 497, 505, 516, 519, 521, 524, 534, 549, 553, 568, 577, 585, 589, 596, 604, 618, 627, 634, 635, 644, 645, 646, 649, 657, 678, 694, 695, 702, 704, 706, 729, 733, 734, 736, 738, 744, 745, 755, 757, 761, 764, 765, 769, 770, 772, 785, 795, 799, 800, 804, 814, 815, 817, 826, 832, 839, 843, 846, 856, 857, 862, 869, 872, 874, 877, 879, 880, 895, 901, 911, 917, 919, 920, 923, 925, 952, 956, 969, 971, 990, 999.

### h2: 135 clean seeds

104, 115, 118, 122, 148, 149, 153, 156, 169, 199, 231, 235, 236, 238, 247, 252, 257, 269, 275, 289, 290, 291, 298, 300, 304, 305, 308, 309, 311, 315, 336, 339, 343, 344, 392, 400, 401, 407, 412, 427, 428, 455, 460, 462, 463, 467, 472, 481, 482, 486, 490, 495, 498, 504, 506, 509, 513, 514, 516, 523, 531, 538, 545, 546, 547, 557, 560, 564, 568, 572, 574, 581, 590, 596, 604, 605, 618, 623, 637, 641, 661, 663, 664, 668, 669, 678, 684, 686, 688, 690, 709, 723, 731, 754, 765, 768, 770, 782, 784, 794, 795, 797, 799, 808, 810, 816, 825, 827, 833, 844, 846, 848, 859, 861, 863, 865, 870, 873, 874, 875, 888, 889, 891, 911, 919, 920, 924, 929, 935, 946, 948, 981, 992, 995, 998.

### h3: 135 clean seeds

101, 106, 111, 118, 121, 130, 134, 154, 156, 159, 162, 175, 186, 206, 209, 237, 240, 245, 249, 265, 267, 270, 275, 276, 291, 293, 297, 300, 302, 318, 319, 322, 324, 333, 345, 348, 349, 352, 373, 378, 380, 410, 417, 427, 433, 439, 440, 446, 449, 450, 452, 456, 463, 464, 477, 479, 490, 491, 522, 528, 536, 538, 541, 544, 545, 547, 554, 565, 570, 589, 591, 599, 602, 603, 604, 608, 620, 640, 642, 650, 662, 668, 672, 682, 683, 684, 685, 689, 702, 703, 716, 720, 736, 766, 769, 773, 780, 782, 783, 787, 796, 797, 801, 807, 816, 826, 834, 839, 846, 863, 867, 872, 875, 882, 893, 899, 903, 916, 924, 925, 927, 930, 931, 940, 947, 951, 956, 960, 962, 966, 974, 982, 984, 991, 998.

### h4: 137 clean seeds

102, 105, 108, 119, 121, 122, 133, 135, 136, 146, 160, 176, 194, 201, 205, 207, 211, 212, 215, 219, 220, 226, 230, 231, 233, 249, 250, 256, 283, 292, 306, 309, 316, 328, 331, 338, 346, 347, 355, 357, 373, 391, 396, 398, 400, 404, 408, 423, 438, 440, 441, 449, 451, 452, 467, 469, 471, 475, 482, 497, 501, 505, 520, 522, 523, 544, 552, 583, 588, 598, 599, 604, 607, 611, 638, 643, 649, 652, 653, 658, 669, 670, 672, 702, 709, 710, 713, 718, 719, 721, 736, 747, 774, 775, 777, 786, 794, 802, 804, 811, 823, 824, 833, 838, 845, 847, 851, 865, 868, 876, 881, 886, 887, 892, 893, 902, 911, 921, 924, 931, 936, 938, 939, 941, 942, 947, 952, 953, 960, 967, 973, 977, 978, 985, 992, 994, 996.

### h5: 139 clean seeds

102, 106, 115, 127, 136, 137, 143, 150, 154, 167, 181, 186, 214, 215, 216, 222, 226, 227, 233, 234, 245, 255, 261, 263, 264, 265, 267, 278, 281, 283, 286, 288, 292, 293, 295, 298, 318, 324, 326, 329, 337, 354, 361, 362, 387, 388, 396, 402, 403, 413, 415, 424, 442, 444, 449, 455, 477, 478, 493, 497, 503, 513, 535, 544, 551, 553, 560, 561, 567, 583, 592, 594, 595, 598, 607, 612, 617, 637, 638, 648, 664, 671, 676, 677, 690, 714, 717, 723, 730, 733, 745, 747, 759, 766, 776, 780, 783, 790, 792, 804, 808, 811, 839, 842, 860, 861, 864, 865, 866, 876, 881, 890, 894, 897, 900, 901, 912, 915, 921, 925, 926, 927, 931, 933, 934, 941, 944, 950, 953, 955, 958, 960, 968, 976, 977, 978, 988, 991, 995.

## Submission provenance and verification

The six builds started at 19:44:57 UTC and finished between 19:53:20 and 19:53:25 UTC (503–508 seconds). All six hotkeys served the validator's requests from prepared submissions. Uploads succeeded between 19:54:53 and 20:44:47 UTC, with 282–288 seconds of upload TTL remaining. Each hotkey has one recorded build and one recorded successful upload for this task. The recorded bands agree with the archived rows' full-seed replay.

- All archived contracts and reference files match the published task after accounting for the seed changing from 0 at submission time to the three stamped scoring seeds.
- All 250 archived stage-1/2 valid-row entries match each submission exactly; accessibility is 0.77.
- All six archived three-seed final scores match the official API exactly.
- Replayed stage-3 outcomes and indel lengths match the archived detail for every row at all three real seeds.
- Scanned 900 seeds × 250 rows × 6 hotkeys = 1,350,000 row simulations. No miner configuration, running process, or uploaded submission was changed.
- [Full analysis](analysis.json), [official score snapshot](scores.json), [published task](task.json), [sanitized build/upload evidence](build_evidence.json).
- Reproduce with `.venv/bin/python reports/af21b154/audit.py` and `.venv/bin/python reports/af21b154/render_report.py` from the repository root.

Uploaded archives:

- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — [data/inst/niome_hotkey/result/2026-09-17T19:55:23](../../data/inst/niome_hotkey/result/2026-09-17T19:55:23/submission.json)
- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-17T20:04:48](../../data/inst/niome_hotkey1/result/2026-09-17T20:04:48/submission.json)
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-17T20:21:49](../../data/inst/niome_hotkey2/result/2026-09-17T20:21:49/submission.json)
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-17T20:19:42](../../data/inst/niome_hotkey3/result/2026-09-17T20:19:42/submission.json)
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — [data/inst/niome_hotkey4/result/2026-09-17T20:45:20](../../data/inst/niome_hotkey4/result/2026-09-17T20:45:20/submission.json)
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — [data/inst/niome_hotkey5/result/2026-09-17T20:36:51](../../data/inst/niome_hotkey5/result/2026-09-17T20:36:51/submission.json)
