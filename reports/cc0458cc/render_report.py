"""Render the actual uploads, fallback and later unsubmitted rebuilds for cc0458cc."""
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
    "and 11-seed HDR bands. h4 uploaded an ordinary fallback with 76 Cas12a + 174 Cas9 rows, nine clean seeds "
    "and no HDR band.** All ten have zero reward weight in this score snapshot.", "",
    "**No uploaded HDR band hit a scoring seed.** h1 and h6 were clean at **326**; h8 was clean at **161**. "
    "Seed **557** was outside every uploaded clean set. All ten later rebuilt this task with 12-seed bands, "
    "overwriting their latest window records after uploading; those later bands were not submitted for this task.", "",
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
          "The **22:17 UTC plan** selected width-100 windows **200 / 400 / 900**, joined as "
          "**200–299 ∪ 400–499 ∪ 900–999** (300 seeds). Actual windows were **100 / 300 / 500** "
          "for seeds **161 / 326 / 557**: **0 of 3** matched the prediction.", "",
          "h0–h5 were assigned 150-seed circular slices of the predicted space at stride 50. "
          "h6–h9 were assigned 150-seed slices of its complement, **100–199 ∪ 300–399 ∪ 500–899**, "
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
          "misleading: for example, h6's `100-349 (150 seeds)` means **100–199 ∪ 300–349**.", "",
          "## Actual uploaded clean sets and HDR bands", "",
          "**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. "
          "Both sets were recovered over 100–999. For each submitted conjunction, the full-row clean set equals "
          "the Cas12a group clean set, and its count matches the original 11-band build log. "
          "The h4 fallback was replayed separately from its unsubmitted prepared builds.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"The uploaded fleet union contains **{len(d['fleet_clean_union'])}/900 clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** "
          f"({sum(len(r['band_seeds']) for r in ours.values())} band memberships before overlaps). "
          "These totals exclude all unsubmitted prepared/rebuilt rows. A candidate window does not guarantee a "
          "clean or HDR hit: h6's candidate slice contained 161 and 326, but only 326 was clean and neither was "
          "in its HDR band; h7's candidate slice contained 557, which was not clean.", "",
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
          "h8's clean hit at 161 scored **39.986504**, with **21.218207** at 326 and **22.188526** at 557. "
          "It had no all-HDR seed to lift the three-seed mean. "
          "The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the "
          "inspected data; they cannot be recovered from its aggregate score.", "",
          "## h4 fallback and later rebuilds", "",
          "h4 began preparation at **22:19:43 UTC**. The validator requested its submission at **22:21:03**. "
          "At **22:21:08** it waited up to **174 seconds** for the prepared build. At **22:24:02**, that build "
          "was still running (259 seconds elapsed), so the miner generated the ordinary fallback while the "
          "prepared build still held the GPU. It uploaded at **22:24:05.683**, with **101 seconds** of URL TTL left. "
          "The prepared 11-band conjunction finished at **22:27:53.941**, about **228 seconds after the upload**. "
          "Its logged 137 clean seeds do not describe the fallback that was submitted.", "",
          "The nine other uploads used their original prepared 11-band builds. After each hotkey had uploaded, "
          "a later rebuild began at **23:07–23:16 UTC** using **k=12**. All latest window records therefore contain "
          "12-seed bands. The later logs also report different clean counts; these are prepared-build observations, "
          "not validated submitted sets. For the nine original conjunctions, the later band record retains the "
          "11 uploaded band seeds and adds one seed. h4's entire later band is unrelated to any uploaded HDR band.", "",
          "| Hotkey | Uploaded clean /900 | Later logged clean /900 (unsubmitted) | Later recorded band seeds absent from uploaded band |",
          "|---|---:|---:|---|"]
for h, r in ours.items():
    extra = sorted(set(r["latest_record_band"]) - set(r["band_seeds"]))
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {r['later_logged_clean_count']} | {nums(extra)} |")
lines += ["", "## Exact uploaded clean sets", "",
          "Each list includes its HDR band and covers the entire 100–999 seed space.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Verification and artifacts", "",
          "- All ten archives confirm successful uploads of this task; contracts and references match the published task after the seed stamp.",
          "- All 2,500 submitted rows passed stage 1/2. Accessibility is 0.87. Archived three-seed final scores exactly match all ten official API scores.",
          "- Replayed outcomes and indel lengths match archived stage-3 detail for every uploaded row at all three scoring seeds.",
          "- Recovered clean sets and bands over **900 × 250 × 10 = 2,250,000 row/seed combinations**.",
          "- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), "
          "[sanitized build/upload evidence](build_evidence.json), [archive manifest](manifest.json), [summary CSV](summary.csv).",
          "- Reproduce using archived validations: `.venv/bin/python reports/cc0458cc/audit.py`, then "
          "`.venv/bin/python reports/cc0458cc/render_report.py`. The analysis includes current validator-code and submission SHA-256 hashes.", "",
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
                     "Percent of leader", "Gap to leader", "Reward weight", "Later unsubmitted band"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["rows"], r["construction"], r["band_group"],
                         r["band_space_ranges"], r["candidate_window_ranges"], r["candidate_window_applies_to_upload"],
                         len(r["clean_seeds"]), nums(r["clean_seeds"]), nums(r["clean_hits"]), nums(r["band_seeds"]),
                         nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"], r["score"]["weight"], nums(r["latest_record_band"])])
print(OUT / "report.md")
