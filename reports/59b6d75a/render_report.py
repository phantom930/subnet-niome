"""Render the task audit from official scores, uploaded rows and full-seed replay."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours = d["hotkeys"]
top = d["top_miner"]
seeds = d["seeds"]
best_h = max(ours, key=lambda h: ours[h]["score"]["final_score"])
best = ours[best_h]


def nums(values):
    return ", ".join(map(str, values)) or "none"


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {score['weight']:.6f} |")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"**{d['cell_type']}**, created **{d['task_created_at']} UTC**. Scoring seeds: **{nums(seeds)}**.", "",
    f"**Best of h0–h6: {best_h}, score {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Top miner: UID {top['miner_uid']}, score {top['final_score']:.6f}.** "
    f"{best_h} achieved {best['percent_of_top']:.2f}% of the leader's score, a gap of {best['gap_to_top']:.6f}. "
    f"The top-10 cutoff was {d['top10_cutoff']:.6f}; {best_h} was {d['top10_cutoff'] - best['score']['final_score']:.6f} below it. "
    "All seven uploads succeeded, contain 250 valid rows, and have zero reward weight.", "",
    "**h1 submitted an ordinary HDR-construction fallback.** Its prepared conjunction build finished after the wait expired. "
    "The other six submitted their prepared conjunction rows. h1's recorded prepared band does not describe its uploaded rows.", "",
    f"Official scores: [task-specific API]({d['source_api']}), fetched {d['fetched_at']}. "
    f"There are {d['score_records']} score records, {d['unique_miners']} distinct miners, "
    f"and {len(d['validators'])} validator. Leader's score timestamp: {top['created_at']} UTC. "
    f"Audit generated {d['generated_at']}.", "",
    "## Official score comparison", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", f"Top miner hotkey: `{top['miner_hotkey']}`.", "",
          "Final score is the mean over three seeds. For our seven submissions, weighted score and fidelity "
          "are constant across seeds, so final = weighted score × mean consistency × fidelity.", "",
          "## Joined windows, clean sets and bands", "",
          "The original plan, generated **2026-09-18 02:17:09 UTC**, selected joined band space "
          "**100–299 ∪ 900–999** (predicted 100-wide windows **100/200/900**). "
          "For the conjunction builds, **cut search covered the full 100–999 space (900 seeds)**. "
          "A 150-seed candidate slice within the joined space supplied each HDR band. "
          "Offsets were **0/50/100/150/200/250/25** for h0–h6, with wraparound. "
          "Each submitted conjunction has **100 Cas12a + 150 Cas9 rows**, and an **11-seed HDR band**. "
          "h1's ordinary fallback has **74 Cas12a + 176 Cas9 rows**.", "",
          "**Clean set:** all 250 uploaded rows cut. **HDR band:** all 250 uploaded rows return HDR. "
          "Both are measured by replay over every seed 100–999. The HDR band is a subset of the clean set. "
          "For all six submitted conjunction builds, the full-row clean set equals the Cas12a group clean set "
          "and matches the build's logged count. h1's uploaded sets are measured independently of its unused prepared build.", "",
          "| Hotkey | Uploaded construction | Planned candidate slice | Offset | Clean / 900 | Clean inside joined / 300 | Clean outside joined |",
          "|---|---|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    construction = "ordinary fallback" if h == "h1" else "conjunction"
    slice_text = r["candidate_window_ranges"] + (" (unused)" if h == "h1" else "")
    lines.append(f"| {h} | {construction} | {slice_text} | {r['candidate_offset']} | {len(r['clean_seeds'])} | "
                 f"{len(r['clean_seeds_in_joined_space'])} | {len(r['clean_seeds_outside_joined_space'])} |")
lines += ["", "| Hotkey | Actual uploaded HDR band seeds | Clean scoring-seed hits | HDR scoring-seed hits |",
          "|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", "Candidate slices are reconstructed from the logged joined space and configured offsets; "
          "the actual bands for all six submitted conjunctions agree exactly with both their original band records and the uploaded-row replay.", "",
          f"Fleet union: **{len(d['fleet_clean_union'])}/900 clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR band seeds** "
          f"({sum(len(r['band_seeds']) for r in ours.values())} memberships before overlap). "
          f"Scoring-seed hits in the clean union: **{nums(d['fleet_clean_hits'])}**. HDR-band hits: **{nums(d['fleet_band_hits'])}**.", "",
          "All three scoring seeds **733/312/818** lie outside the predicted joined band space. "
          "The wider cut search still gave **h2 an all-cut hit at 733** and **h4 an all-cut hit at 818**. "
          "Seed 312 is outside every submitted clean set. No hotkey hit its HDR band.", "",
          "## Per-seed scores and outcomes", "",
          "Each entry is **per-seed final score / no-cut rows**. A clean hit has zero no-cut rows.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    values = [f"{r['per_seed'][str(s)]['final_score']:.6f} / {r['per_seed'][str(s)]['outcomes'].get('no_cut', 0)}"
              for s in seeds]
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
lines += ["", "## Why the leader scores higher", "",
          "h4's weighted score is **332.428168**, above the leader's **315.569149**. "
          "Its fidelity is somewhat lower (**0.925143 versus 0.969200**), but the largest gap is mean consistency: "
          "**0.117579 versus 0.662879**. h4's one clean hit did not produce a competitive three-seed mean. "
          "Across h0–h6, mean consistency ranges from **0.083069 to 0.117579**. "
          "The public leader record contains aggregate scores; its submission, clean set, HDR band and individual seed outcomes "
          "were not available in the inspected data and cannot be inferred from those aggregates.", "",
          "## h1: prepared build versus actual upload", "",
          "- **03:07:33.404 UTC:** validator handler finds the prepared build still running and waits up to 179 seconds.",
          "- **03:10:32.981:** wait ends; handler reports the prepared round unusable and starts ordinary construction. "
          "The prepared build still holds the GPU, so conjunction/hedge paths are skipped.",
          "- **03:10:34.138:** ordinary construction has generated 250 valid rows; local scoring begins.",
          "- **03:10:38.892:** the prepared conjunction finishes, about six seconds after the fallback started. "
          "It reports 140 clean seeds and an 11-seed band, but those are not the rows submitted.",
          "- **03:10:40.447:** ordinary fallback upload succeeds, with 97 seconds of URL TTL remaining.", "",
          f"The unused prepared band is **{nums(ours['h1']['recorded_prepared_band'])}**. "
          f"The actual uploaded fallback has only **{len(ours['h1']['clean_seeds'])} clean seeds: {nums(ours['h1']['clean_seeds'])}**, "
          "and **no all-HDR band**. Its local score reproduces the official **32.077144**. "
          "`window_used.json` combines source `all_hdr_not_attempted` with `all_hdr_built: true` and the prepared band; "
          "the concurrent prepared build completed after the fallback started. Treating that record's band as the uploaded band would be incorrect.", "",
          "## Exact uploaded clean sets", "",
          "These lists cover every seed 100–999 and include the HDR band seeds. "
          "Full arrays and per-seed outcomes are also in [analysis.json](analysis.json).", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Verification and artifacts", "",
          "- All seven archived upload records confirm submission of this task.",
          "- Archived contracts and references match the published task, allowing only the subsequent scoring-seed stamp.",
          "- All 250 uploaded rows per hotkey pass stages 1/2 and match the valid experiment entries; cell accessibility is 0.87.",
          "- All seven locally calculated three-seed scores match the official API exactly. Existing validation artifacts were reused "
          "for h0/h5/h6; the other four were scored with calc.py into this report directory.",
          "- Stage-3 replay matches saved or regenerated outcomes and indel lengths for every row at all three scoring seeds.",
          "- Scanned **900 × 250 × 7 = 1,575,000 row simulations** to recover actual clean sets and HDR bands.",
          "- [Official scores](scores.json), [published task](task.json), [sanitized build/upload logs](build_evidence.json), "
          "[archive and validation locations](manifest.json), [summary CSV](summary.csv), [full analysis](analysis.json).",
          "- Reproduce with `.venv/bin/python reports/59b6d75a/audit.py` then "
          "`.venv/bin/python reports/59b6d75a/render_report.py` from the repository root. "
          "To regenerate a missing validation, run `calc.py --folder <archive> --out-dir reports/59b6d75a/<hotkey> "
          "--task reports/59b6d75a/task.json --cell-types reports/59b6d75a/cell_types.json --no-compare --quiet`.", "",
          "Uploaded archives:", ""]
for h, r in ours.items():
    lines.append(f"- **{h}** `{r['hotkey']}` — [{r['archive']}](../../{r['archive']}/submission.json), "
                 f"SHA-256 `{r['submission_sha256']}`.")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Score", "Construction", "Planned candidate slice", "Clean count /900",
                     "Clean hits", "HDR band", "HDR hits", "Percent of top", "Gap to top"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["construction"],
                         r["candidate_window_ranges"], len(r["clean_seeds"]), nums(r["clean_hits"]), nums(r["band_seeds"]),
                         nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"]])
print(OUT / "report.md")
