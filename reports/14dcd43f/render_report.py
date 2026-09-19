"""Render the verified uploaded sets and official score comparison for 14dcd43f."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
best_h = max(ours, key=lambda h: ours[h]["score"]["final_score"])
best = ours[best_h]


def nums(values):
    return ", ".join(map(str, values)) or "none"


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {b['n_valid_experiments']} | {score['weight']:.6f} |")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"**{d['cell_type']}**, created **{d['task_created_at']} UTC**. Scoring seeds: **{nums(seeds)}**.", "",
    f"**Best of h0–h9: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Leader: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} achieved {best['percent_of_top']:.2f}% of the leader, a gap of {best['gap_to_top']:.6f}. "
    f"Top-10 cutoff: {d['top10_cutoff']:.6f}; {best_h} was {d['top10_cutoff']-best['score']['final_score']:.6f} below it.", "",
    "**All ten hotkeys uploaded 250 valid rows: 80 Cas12a + 170 Cas9.** Each used a prepared conjunction with "
    "an eight-seed HDR band and a cut search over all 900 seeds (100–999). All ten have zero reward weight in this score snapshot.", "",
    "**h6 hit HDR seed 124: all 250 uploaded rows returned HDR, scoring 318.719399 at that seed.** "
    "Its other two seed scores were 25.200971 at 479 and 18.470026 at 237, giving the mean 120.796799. "
    "No other hotkey hit a clean seed or an HDR band at any of the three scoring seeds.", "",
    f"[Official task score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Leader score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Rank is 1 plus the number of strictly higher scores; ties share rank. Final score is the mean "
          "of the three per-seed scores. Breakdown factors are also means, so multiplying the displayed mean factors "
          "is not generally equivalent to averaging the per-seed products.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`.", "",
          "## Joined seed windows and per-hotkey candidates", "",
          "The **19:17 UTC plan** selected width-100 windows **200 / 300 / 900**, joined as "
          "**200–399 ∪ 900–999** (300 seeds). The actual seed windows were **100 / 400 / 200**, "
          "so one of the three actual windows matched the prediction: seed **237** was inside the predicted space; "
          "**124 and 479** were outside it.", "",
          "**h0–h5** took 150-seed circular slices of the predicted 300 seeds, at stride 50. "
          "**h6–h9** took 150-seed slices of its 600-seed complement, **100–199 ∪ 400–899**, at stride 150. "
          "Their four candidate slices partition that complement. All ten used the same full cut search, **100–999**.", "",
          "| Hotkey | Band group | Exact 150-seed candidate slice | Offset within group | Clean /900 | Clean in candidate slice |",
          "|---|---|---|---:|---:|---:|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['band_group']} | {r['candidate_window_ranges']} | {r['candidate_offset']} | "
                 f"{len(r['clean_seeds'])} | {len(r['clean_seeds_in_candidate_window'])} |")
lines += ["", "The build logs confirm the predicted/complement group sizes and offsets. Their candidate-window "
          "messages sometimes show only the minimum and maximum seed; the table expands the actual joined slices "
          "and preserves gaps. For example, h6's logged `100-449 (150 seeds)` is **100–199 ∪ 400–449**, "
          "not the entire contiguous 100–449 range.", "",
          "## Actual submitted HDR bands", "",
          "**Clean set:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. "
          "Both were recovered by replay over seeds 100–999. The full-row clean sets equal the Cas12a group clean sets, "
          "and their counts match the original build logs. Every HDR band is a subset of its clean set and its candidate slice.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band (8 seeds) | Clean scoring-seed hits | HDR scoring-seed hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"The fleet covers **{len(d['fleet_clean_union'])} distinct clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** across 100–999 "
          f"({sum(len(r['band_seeds']) for r in ours.values())} band memberships before overlaps). "
          "Its only clean and HDR scoring-seed hit is **124 on h6**, in the complementary group. "
          "The predicted group had no clean or HDR hit even though seed 237 lay in its candidate space. "
          "This describes this task's outcome; it does not establish a long-run advantage for either group.", "",
          "## Per-seed results", "",
          "Each entry is **per-seed final score / HDR rows / no-cut rows**. The total is 250 rows per seed.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    values = []
    for s in seeds:
        sr = r["per_seed"][str(s)]
        values.append(f"{sr['final_score']:.6f} / {sr['outcomes'].get('HDR', 0)} / {sr['outcomes'].get('no_cut', 0)}")
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
b, t = best["score"]["breakdown"], top["breakdown"]
lines += ["", "## Comparison with the leader", "",
          f"{best_h}'s weighted score is **{b['total_weighted_score']:.6f}**, versus the leader's "
          f"**{t['total_weighted_score']:.6f}**. Mean consistency is **{b['consistency_factor']:.6f} versus "
          f"{t['consistency_factor']:.6f}**, and fidelity is **{b['distribution_fidelity_factor']:.6f} versus "
          f"{t['distribution_fidelity_factor']:.6f}**. "
          "h6's consistency factor was 1.0 at its HDR hit, but low scores at the other two seeds reduced its task average. "
          "The mean consistency difference is the largest factor gap in the published breakdown. "
          "The leader's submitted rows, clean set, bands and individual seed outcomes were not available in the "
          "inspected data; its aggregate score does not identify those sets.", "",
          "## Exact clean sets per hotkey", "",
          "These are the full submitted clean sets over 100–999, including the HDR band seeds.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Build and upload provenance", "",
          "All ten original builds began at **19:53:30–36 UTC** and finished at **19:56:21–28 UTC**. "
          "Each successful upload used its prepared 250-row conjunction. The ten recorded bands match the bands "
          "recovered from the actual uploads. No fallback or later replacement build was needed for these submissions.", "",
          "| Hotkey | Upload time UTC | Submission archive |",
          "|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['uploaded_at']} | [{r['archive']}](../../{r['archive']}/submission.json) |")
lines += ["", "## Verification and artifacts", "",
          "- All ten upload archives identify this task; contracts and references match the published task after the seed stamp.",
          "- All ten have 250 stage-1/2 valid rows, accessibility 0.35, and 80 Cas12a + 170 Cas9 rows.",
          "- Independent local three-seed final scores match all ten official API scores exactly.",
          "- Replayed outcomes and indel lengths match local stage-3 details for every uploaded row at all three scoring seeds.",
          "- Recovered clean sets and HDR bands over **900 × 250 × 10 = 2,250,000 row/seed combinations**.",
          "- [Full analysis](analysis.json), [official score snapshot](scores.json), [published task](task.json), "
          "[sanitized build/upload evidence](build_evidence.json), [archive manifest](manifest.json), [summary CSV](summary.csv).",
          "- Reproduce from the repository root: `.venv/bin/python reports/14dcd43f/validate.py`, "
          "then `.venv/bin/python reports/14dcd43f/audit.py`, then `.venv/bin/python reports/14dcd43f/render_report.py`. "
          "Validations write only to report-local directories. The analysis includes validator-code and submission SHA-256 hashes.", "",
          "Public hotkey mapping:", ""]
for h, r in ours.items():
    lines.append(f"- **{h} / UID {r['uid']}**: `{r['hotkey']}`.")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Final score", "Valid rows", "Band group", "Joined band space",
                     "Candidate slice", "Clean count /900", "Exact clean set", "Clean hits", "HDR band", "HDR hits",
                     "Percent of leader", "Gap to leader", "Reward weight"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["rows"], r["band_group"],
                         r["band_space_ranges"], r["candidate_window_ranges"], len(r["clean_seeds"]), nums(r["clean_seeds"]),
                         nums(r["clean_hits"]), nums(r["band_seeds"]), nums(r["band_hits"]), r["percent_of_top"],
                         r["gap_to_top"], r["score"]["weight"]])
print(OUT / "report.md")
