"""Render the mixed-generation submission audit without reading live miner state."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
assert len(ours) == 10 and "fleet_clean_union" in d
best_h = d["best_our_hotkey"]
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
    f"{best_h} reached {best['percent_of_top']:.2f}% of the leader, a gap of {best['gap_to_top']:.6f}. "
    f"Top-10 cutoff: {d['top10_cutoff']:.6f}.", "",
    "All ten hotkeys uploaded prepared conjunctions with **250 valid rows: 80 Cas12a + 170 Cas9**. "
    "Each uploaded an **8-seed HDR band**. No fallback was used. "
    "**None of the three scoring seeds was clean or all-HDR on any hotkey.**", "",
    "**This task spans a configuration change. h0–h3 and h9 uploaded the original narrow-cut builds; "
    "h4–h8 uploaded the later shared-window, wide-cut builds.** The latest `window_used.json` entries "
    "describe the rebuild for all ten and therefore contain the wrong uploaded bands for h0–h3 and h9. "
    "This audit derives bands and clean sets from the actual archived submission rows.", "",
    f"[Official score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Ranks count strictly higher scores, so ties share rank. All identities were matched by full "
          "public hotkeys. h1–h3 appear at UIDs **189, 45 and 136** for this task; their older UIDs are not used. "
          "Final scores and breakdown factors are means across three seeds. The product of mean factors need "
          "not equal the mean final score.", "",
          "## Joined seed windows and cut search", "",
          "The original round plan (logged **06:04 UTC**) predicted width-100 windows **100 / 200 / 300**, "
          "joined as **100–399**. The rebuild's plan (generated **06:17:51 UTC**) records the same prediction. "
          "Actual scoring windows were **400 / 900 / 100** for seeds **423 / 999 / 189**: **one matched window**.", "",
          "The later shared layout instead fixes its band candidates to **100–299 ∪ 500–599** "
          "(classes **100 / 200 / 500**, 300 seeds), independently of the plan. It likewise covers one "
          "of the three actual windows. The latest record's `window: [100,399]` describes the plan, "
          "not the rebuild's actual shared band candidates. h4–h8 submitted shared loops **5–9 of 10**.", "",
          "Original h1–h3 each use a 100-seed band slice, while their cut search spans the full predicted "
          "300 seeds. Original h0/h9 each use 300 band candidates. Rebuilt h4–h8 use all 900 seeds "
          "for cut selection. Band candidates and cut-search seeds are separate quantities.", "",
          "| Hotkey | Uploaded build | Joined band candidates | Candidate count | Shared loop | Cut search | Clean in cut search | Clean over 100–999 |",
          "|---|---|---|---:|---:|---|---:|---:|"]
for h, r in ours.items():
    generation = "original / narrow" if r["uploaded_generation"] == "original" else "rebuilt / wide"
    lines.append(f"| {h} | {generation} | {r['candidate_window_ranges']} | {len(r['candidate_window'])} | "
                 f"{r['band_loop'] or '—'} | {r['cut_search_ranges']} | "
                 f"{len(r['clean_seeds_in_cut_space'])}/{len(r['cut_search_space'])} | {len(r['clean_seeds'])} |")
lines += ["", "## Uploaded HDR bands", "",
          "**Clean** means all 250 rows cut. **HDR band** means all 250 rows return HDR. "
          "Sets were recovered by replay over every seed from 100 to 999. "
          "Every band is contained in its uploaded candidate window and clean set. "
          "The build's Cas12a clean count is checked within its cut-search space, separately from "
          "the full-row clean definition.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits | Latest band record matches upload |",
          "|---|---:|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | "
                 f"{nums(r['band_hits'])} | {'yes' if r['latest_record_band_matches_upload'] else 'NO — later rebuild'} |")
lines += ["", f"Fleet union: **{len(d['fleet_clean_union'])} distinct clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds**, from 80 band memberships. "
          "The five uploaded shared-loop bands are mutually disjoint, but bands across the mixed fleet "
          "overlap. Overlapping seeds: " + "; ".join(f"{s} ({'/'.join(hs)})" for s, hs in sorted(d["duplicate_band_seeds"].items(), key=lambda x: int(x[0]))) + ".", "",
          "## Per-seed results", "",
          "Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.", "",
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
          "Our best weighted score exceeds the leader's, but its much lower consistency accompanies a "
          "much lower final score. Every uploaded hotkey missed all three seeds with its full clean set "
          "and HDR band. This single mixed-configuration task does not establish which layout performs "
          "better across tasks.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`. Its submitted rows, clean sets, HDR bands and per-seed "
          "outcomes were unavailable in the inspected data. Its aggregate score does not identify those sets.", "",
          "## Exact uploaded clean sets", "",
          "These lists cover the full 100–999 range and include the corresponding HDR bands.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
    extra = sorted(set(r["cas12a_clean_seeds"]) - set(r["clean_seeds"]))
    if extra:
        lines += [f"Cas12a-only clean seeds excluded from the full-row set: {nums(extra)}.", ""]
lines += ["## Build and upload provenance", "",
          "Original builds completed around 06:12–06:13 UTC. Rebuilds began at 06:39 and completed "
          "around 06:42. Each log explicitly reports using the prepared 250-row submission. "
          "The timestamps below select the build actually available at upload time, and replay confirms "
          "the uploaded band and clean count.", "",
          "| Hotkey | Uploaded build completed UTC | Uploaded UTC | Latest window record UTC |",
          "|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['build_completed_at']} | {r['uploaded_at']} | {r['latest_record_at']} |")
lines += ["", "## Verification and artifacts", "",
          "- All 10 archive upload records identify this task and a successful upload. Contracts and references match the published task after substituting its seed stamp.",
          "- All 2,500 rows passed validation with HEK293 accessibility 0.35; all 10 local three-seed final scores match official scores within 1e-10.",
          "- Reused h7's archived validation; generated the other nine validations under this report directory.",
          "- Replayed every row at all 900 seeds: **2,250,000 row/seed combinations**. At the three scoring seeds, every outcome and indel length matches local stage-3 detail.",
          "- Clean counts within cut-search windows match all selected build logs. Recovered bands match the latest records for h4–h8 and differ for h0–h3/h9, as their upload timing predicts.",
          "- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).",
          "- Reproduce with `.venv/bin/python reports/09b7b842/validate.py`, then `.venv/bin/python reports/09b7b842/audit.py`, then `.venv/bin/python reports/09b7b842/render_report.py`. Validator and submission SHA-256 hashes are included in the analysis.", "",
          "| Hotkey | Public hotkey | Submission archive |", "|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | `{r['hotkey']}` | [{r['archive']}](../../{r['archive']}/submission.json) |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Uploaded generation", "Valid rows", "Band candidate window", "Shared loop",
                     "Cut window", "Clean count /900", "Clean count in cut space", "Exact clean set", "Clean hits", "HDR band", "HDR hits",
                     "Latest record matches uploaded band", "Percent of leader", "Gap to leader", "Reward weight"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["uploaded_generation"], r["rows"],
                         r["candidate_window_ranges"], r["band_loop"], r["cut_search_ranges"], len(r["clean_seeds"]),
                         len(r["clean_seeds_in_cut_space"]), nums(r["clean_seeds"]), nums(r["clean_hits"]), nums(r["band_seeds"]),
                         nums(r["band_hits"]), r["latest_record_band_matches_upload"], r["percent_of_top"], r["gap_to_top"], r["score"]["weight"]])
print(OUT / "report.md")
