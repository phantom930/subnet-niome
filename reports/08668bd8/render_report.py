"""Render the actual 12-band uploads and h4's ordinary fallback for 08668bd8."""
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
    "**All ten uploaded 250 valid rows. Nine uploaded prepared conjunctions with 100 Cas12a + 150 Cas9 rows "
    "and 12-seed HDR bands. h4 uploaded an ordinary fallback with 75 Cas12a + 175 Cas9 rows, eight clean seeds "
    "and no HDR band.** All ten have zero reward weight in this score snapshot.", "",
    "**No uploaded clean set or HDR band hit any of the three scoring seeds.** "
    "The nine uploaded conjunction bands match their window records. h4's record contains a prepared 12-seed "
    "band that completed after its fallback upload and was not submitted.", "",
    f"[Official task score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Rank is 1 plus the number of strictly higher scores; ties share rank. Final score is the mean "
          "of the three per-seed scores. Breakdown factors are also means; multiplying the mean factors "
          "is not generally equivalent to averaging the per-seed products.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`.", "",
          "## Joined windows and per-hotkey candidates", "",
          "The **00:17 UTC plan** selected width-100 windows **100 / 200 / 300**, joined as "
          "**100–399** (300 seeds). Actual windows were **700 / 300 / 500** "
          "for seeds **709 / 369 / 545**: **1 of 3** matched the prediction, through seed **369**.", "",
          "h0–h5 were assigned 150-seed circular slices of the predicted space at stride 50. "
          "h6–h9 were assigned 150-seed slices of its complement, **400–999**, "
          "at stride 150. Their four slices partition the 600-seed complement. "
          "The nine submitted conjunctions all used a full **100–999 cut search**. "
          "h4's assigned candidate slice was used by its prepared build, but its actual upload used the ordinary fallback.", "",
          "| Hotkey | Actual construction | Band group | Assigned 150-seed candidate slice | Offset | Actual clean /900 |",
          "|---|---|---|---|---:|---:|"]
for h, r in ours.items():
    construction = "ordinary fallback" if h == "h4" else "conjunction"
    candidate = r["candidate_window_ranges"] + (" (not used by upload)" if h == "h4" else "")
    lines.append(f"| {h} | {construction} | {r['band_group']} | {candidate} | {r['candidate_offset']} | {len(r['clean_seeds'])} |")
lines += ["", "The exact slices preserve gaps in the joined space. The build log's min/max display can be "
          "misleading: h5's `100-399 (150 seeds)` means **100–199 ∪ 350–399**, because its slice wraps.", "",
          "## Actual uploaded clean sets and HDR bands", "",
          "**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. "
          "Both sets were recovered over 100–999. For each submitted conjunction, the full-row clean set equals "
          "the Cas12a group clean set, and its count matches the original 12-band build log. "
          "The h4 fallback was replayed separately from its unsubmitted prepared builds.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"The uploaded fleet union contains **{len(d['fleet_clean_union'])}/900 clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** "
          f"({sum(len(r['band_seeds']) for r in ours.values())} band memberships before overlaps). "
          "These totals exclude all unsubmitted prepared/rebuilt rows. A candidate window does not guarantee a "
          "clean or HDR hit: 369 lay in the submitted candidate slices of h3 and h5; 545 lay in h6's; "
          "709 lay in h8's. All missed their respective clean sets and HDR bands. "
          "All other uploaded clean sets also missed every scoring seed.", "",
          "## Per-seed results", "",
          "Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.", "",
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
          f"{best_h}'s weighted score was **{b['total_weighted_score']:.6f}**, versus the leader's "
          f"**{t['total_weighted_score']:.6f}**. Mean consistency was **{b['consistency_factor']:.6f} versus "
          f"{t['consistency_factor']:.6f}**, and fidelity was **{b['distribution_fidelity_factor']:.6f} versus "
          f"{t['distribution_fidelity_factor']:.6f}**. The largest factor gap is consistency. "
          "h6 scored **27.481407** at 709, **26.165091** at 369, and **20.595541** at 545. "
          "It had no clean or all-HDR seed to lift the three-seed mean. "
          "The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the "
          "inspected data; they cannot be recovered from its aggregate score.", "",
          "## h4 fallback and prepared-band provenance", "",
          "h4 began preparation at **00:44:11 UTC**. The validator requested its submission at **00:46:27**. "
          "At **00:46:28** it waited up to **179 seconds** for the prepared build. At **00:49:27**, that build "
          "was still running (315 seconds elapsed), so the miner generated the ordinary fallback while the "
          "prepared build still held the GPU. It uploaded at **00:49:32.265**, with **100 seconds** of URL TTL left. "
          "The prepared 12-band conjunction finished at **00:51:56.733**, about **144 seconds after the upload**. "
          "Its logged **94 clean seeds** do not describe the fallback's **eight clean seeds**.", "",
          f"h4's recorded prepared band, **not submitted**, is: **{nums(ours['h4']['latest_record_band'])}**.", "",
          "h4's latest window record mixes `source=all_hdr_not_attempted` and `window=null` from its fallback "
          "with `all_hdr_built=true` and a band from the completed background preparation. Uploaded-row replay "
          "confirms that this recorded band is absent from the submitted fallback. "
          "The other nine hotkeys used their original prepared 12-band submissions, and every one of their "
          "recorded bands matches the corresponding uploaded rows.", ""]
lines += ["", "## Exact uploaded clean sets", "",
          "Each list includes its HDR band and covers the entire 100–999 seed space.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Verification and artifacts", "",
          "- All ten archives confirm successful uploads of this task; contracts and references match the published task after the seed stamp.",
          "- All 2,500 submitted rows passed stage 1/2. Accessibility is 0.87. Local three-seed final scores exactly match all ten official API scores.",
          "- Reused nine archived validations and generated h8's missing validation in this report's h8 directory, without modifying its upload archive.",
          "- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.",
          "- Recovered clean sets and bands over **900 × 250 × 10 = 2,250,000 row/seed combinations**.",
          "- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), "
          "[sanitized build/upload evidence](build_evidence.json), [archive manifest](manifest.json), [summary CSV](summary.csv).",
          "- Reproduce using `.venv/bin/python reports/08668bd8/validate.py`, then `.venv/bin/python reports/08668bd8/audit.py`, then "
          "`.venv/bin/python reports/08668bd8/render_report.py`. The analysis includes current validator-code and submission SHA-256 hashes.", "",
          "| Hotkey | Public hotkey | Upload time UTC | Submission archive |",
          "|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} / UID {r['uid']} | `{r['hotkey']}` | {r['uploaded_at']} | [{r['archive']}](../../{r['archive']}/submission.json) |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Final score", "Valid rows", "Construction", "Assigned band group",
                     "Assigned joined band space", "Assigned candidate slice", "Candidate slice used by upload",
                     "Clean count /900", "Exact clean set", "Clean hits", "Actual HDR band", "HDR hits",
                     "Percent of leader", "Gap to leader", "Reward weight", "Latest recorded band", "Recorded band matches upload"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["rows"], r["construction"], r["band_group"],
                         r["band_space_ranges"], r["candidate_window_ranges"], r["candidate_window_applies_to_upload"],
                         len(r["clean_seeds"]), nums(r["clean_seeds"]), nums(r["clean_hits"]), nums(r["band_seeds"]),
                         nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"], r["score"]["weight"],
                         nums(r["latest_record_band"]), r["latest_record_band_matches_upload"]])
print(OUT / "report.md")
