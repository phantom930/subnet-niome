"""Render the saved task audit into a Markdown report."""
import json
import re
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
evidence = json.loads((OUT / "build_evidence.json").read_text())
top = d["top_miner"]
ours = d["hotkeys"]
best_h = max(ours, key=lambda h: ours[h]["score"]["final_score"])
best = ours[best_h]


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {score['weight']:.6f} |")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"HEK293 · task opened {d['task_created_at']} UTC · scoring seeds **117, 842, 989**.", "",
    f"**Best of h0–h5: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Top miner: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} reached {best['percent_of_top']:.2f}% of the top score, a gap of {best['gap_to_top']:.6f} points. "
    "All six submissions were uploaded successfully, contained 250 valid rows, and received zero reward weight.", "",
    f"Scores are from the [task-specific public API]({d['source_api']}); "
    f"{d['score_records']} records, {d['unique_miners']} distinct miners, one validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated: {d['generated_at']}.", "",
    "## Official score comparison", "",
    "| Miner | UID this round | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", f"Top miner hotkey: `{top['miner_hotkey']}`. "
          f"The rank-10 cutoff was **{d['top10_cutoff']:.6f}**.", "",
          "Final score is the mean of the per-seed scores. For these submissions the weighted score "
          "and fidelity do not change across seeds, so final = weighted score × mean consistency × fidelity.", "",
          "## Joined window, clean sets, and bands", "",
          "The original 17:19 builds used joined band space **200–399 ∪ 900–999** (300 seeds), "
          "with a separate **100–999 cut search** (900 seeds). "
          "Band candidate slices were 150 seeds wide, at offsets 0/50/100/150/200/250 for h0–h5. "
          "All six shipped the conjunction build with 80 Cas12a and 170 Cas9 rows, pinning eight HDR seeds.", "",
          "**Clean** means all 250 submitted rows cut at that seed. "
          "**Band** means all 250 return HDR; the HDR band is a subset of the clean set. "
          "Both sets below were recovered by replaying the exact uploaded rows on every seed 100–999. "
          "The full-row clean sets equal the Cas12a group clean sets and their sizes match the original build logs.", "",
          "| Hotkey | Candidate slice within joined space | Clean / 900 | HDR band seeds | Hits among 117, 842, 989 |",
          "|---|---|---:|---|---|"]
for h, r in ours.items():
    first_build = next(s for s in evidence[h]["events"] if "Build: conjunction " in s)
    match = re.search(r"clean (\d+)/(\d+)", first_build)
    assert (int(match[1]), int(match[2])) == (len(r["clean_seeds"]), 900)
    assert r["clean_seeds"] == r["cas12a_clean_seeds"]
    assert r["cas_mix"] == {"Cas12a": 80, "Cas9": 170}
    lines.append(f"| {h} | {r['candidate_window_ranges']} | {len(r['clean_seeds'])} | "
                 f"{', '.join(map(str, r['band_seeds']))} | 0 clean; 0 band |")
lines += ["", "Candidate slices are reconstructed from the original logged joined space and the six-hotkey "
          "offset layout; all recovered bands lie inside their corresponding slices.", "",
          "The fleet covers **99 distinct clean seeds** and **41 distinct HDR band seeds** "
          "(48 band memberships before overlap). **None of 117, 842, or 989 lies in either union.** "
          "Seeds 117 and 842 are outside the joined band space. Seed 989 is inside the candidate slices "
          "of h3/h4/h5 but was not selected into any band and is not clean for any hotkey.", "",
          "Exact clean sets (including band seeds):", ""]
for h, r in ours.items():
    lines.append(f"- **{h} ({len(r['clean_seeds'])}):** " + ", ".join(map(str, r["clean_seeds"])) + ".")
lines += ["", "## What happened at the three real seeds", "",
          "Each cell below shows **per-seed final score / no-cut rows**. "
          "Every hotkey had at least one no-cut row on every drawn seed, losing the all-cut consistency benefit.", "",
          "| Hotkey | Seed 117 | Seed 842 | Seed 989 | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    values = [f"{r['per_seed'][str(s)]['final_score']:.4f} / {r['per_seed'][str(s)]['outcomes'].get('no_cut',0)}"
              for s in d["seeds"]]
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
lines += ["", "## Why the scores trail the leader", "",
          "All six weighted scores (306.42–311.65) exceed the top miner's 292.74. "
          "Fidelity is somewhat lower (0.9146–0.9219 versus 0.9599), but the dominant difference is "
          "consistency: **0.0702–0.0862 versus 0.561675**. "
          "Missing both the clean sets and HDR bands leaves all six on the low-consistency outcomes for this round. "
          "This establishes what happened on this task; it does not estimate long-run strategy performance. "
          "The top miner's bands and clean sets cannot be inferred from its aggregate score, and were not reconstructed.", "",
          "## Submission history caveat", "",
          "The uploaded rows came from the **17:19–17:23 UTC builds** and were submitted between "
          "17:23:59 and 17:55:56 UTC. All six used prepared submissions with 283–286 seconds of upload TTL remaining. "
          "Later, at **18:54 UTC**, each hotkey rebuilt this same task after HEK293's cut search was narrowed "
          "to 300 seeds. Those later builds overwrite the task's `window_used.json` band entry and have different "
          "bands and clean counts. They were **not the uploaded submissions**: each hotkey has one recorded upload "
          "for this task, before the rebuild. This report therefore uses the uploaded archives and original logs, "
          "not the later overwritten band entries. h0's UID for this round was **248**, as shown by its upload "
          "target and the official score record.", "",
          "## Verification and artifacts", "",
          "- All six archived contracts match the published task, allowing for the subsequently stamped seed; "
          "their reference files also match after that same seed update.",
          "- Each archived stage-1/2 valid-row list matches the submitted designs, with accessibility 0.35.",
          "- All six archived three-seed final scores match the official API exactly.",
          "- Stage-3 replay matches every row's archived outcome and indel length at all three real seeds.",
          "- Scanned all 900 seeds for each hotkey: 1,350,000 row simulations; no new build or upload was performed.",
          "- [Machine-readable analysis and full seed sets](analysis.json), [official score snapshot](scores.json), "
          "[published task snapshot](task.json), and [sanitized build/upload evidence](build_evidence.json).",
          "- Reproduce with `.venv/bin/python reports/7bca65ce/audit.py` then "
          "`.venv/bin/python reports/7bca65ce/render_report.py` from the repository root.", "",
          "Uploaded archives:", ""]
for h, r in ours.items():
    lines.append(f"- **{h}** `{r['hotkey']}` — [{r['archive']}](../../{r['archive']}/submission.json)")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
print(OUT / "report.md")
