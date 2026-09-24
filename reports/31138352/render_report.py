"""Render six narrow-cut conjunctions, h6's fallback, and three absent hotkeys."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
uploaded = {h: r for h, r in ours.items() if r["submitted"]}
missing = {h: r for h, r in ours.items() if not r["submitted"]}
best_h = max(uploaded, key=lambda h: uploaded[h]["score"]["final_score"])
best = uploaded[best_h]


def nums(values):
    return "N/A" if values is None else (", ".join(map(str, values)) or "none")


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
    "**Seven uploads were found: h0 and h4–h9**, each with **250 valid rows**. Six uploaded prepared "
    "conjunctions with **80 Cas12a + 170 Cas9** and **eight-seed HDR bands**. **h6 uploaded an ordinary "
    "fallback with 75 Cas12a + 175 Cas9, zero clean seeds and no HDR band across 100–999.** "
    "All seven have zero reward weight in this snapshot. h1–h3 have no task build/window/upload records "
    "or official score entries; their results are unavailable, not zero.", "",
    "**h4 hit HDR seed 431: all 250 rows returned HDR**, scoring **307.934233** at that seed. "
    "Its other seed scores were **24.383393** at 189 and **24.192852** at 566, giving the mean **118.836826**. "
    "This was the only uploaded clean or HDR hit. Seeds 189 and 566 were outside every uploaded clean set.", "",
    f"[Official task score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
for h, r in ours.items():
    lines.append(score_line(h, r["rank"], r["score"]) if r["submitted"] else f"| {h} | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |")
lines += ["", "Rank is 1 plus the number of strictly higher scores; ties share rank. Final score and breakdown "
          "factors are means over three seeds. Multiplying mean factors is not generally equivalent to averaging "
          "per-seed products. Missing hotkeys have no score or rank.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`.", "",
          "## Joined windows, actual band space and cut windows", "",
          "The **12:24 UTC plan** recorded width-100 windows **100 / 200 / 300**, joined as **100–399**. "
          "Actual windows were **400 / 100 / 500** for scoring seeds **431 / 189 / 566**: one matched that prediction.", "",
          "For these conjunctions, the plan's predicted space did **not** determine band placement. "
          "The active layout used the full **100–999** space, with **300-seed contiguous windows at stride 100** "
          "for h0, h4, h5, h6, h7, h8 and h9. The source's `band_space()` ignores the predicted argument for "
          "these hotkeys. The logs retain the label `predicted space of 900`; that denotes this full space, "
          "not the plan's 300 predicted seeds.", "",
          "For **HEK293**, `cut mode=union` resolves to each hotkey's **own 300-seed band candidate window**. "
          "The original build logs confirm `cut 300 seeds (union) | k=8`; this task did not use a 900-seed "
          "cut search. h6 prepared that configuration but uploaded an ordinary fallback instead.", "",
          "| Hotkey | Actual construction | Assigned band candidates / conjunction cut window | Offset in 100–999 | Actual clean /900 | Clean inside assigned window |",
          "|---|---|---|---:|---:|---:|"]
for h, r in ours.items():
    if not r["submitted"]:
        lines.append(f"| {h} | No task build or upload | N/A | N/A | N/A | N/A |")
        continue
    fallback = r["construction"] == "ordinary_fallback"
    label = "ordinary fallback" if fallback else "conjunction"
    window = r["candidate_window_ranges"] + (" (prepared only)" if fallback else "")
    lines.append(f"| {h} | {label} | {window} | {r['candidate_offset']} | {len(r['clean_seeds'])} | {len(r['clean_seeds_in_candidate_window'])} |")
lines += ["", "The six uploaded conjunction candidate windows together still cover all **100–999**, despite "
          "h6's fallback. They overlap. All recovered full-row clean seeds lie within their corresponding "
          "conjunction cut windows; replay found no additional full-row clean seeds outside them.", "",
          "## Actual uploaded HDR bands", "",
          "**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. "
          "Both were recovered by scanning all 900 seeds (100–999), independently of the 300-seed build search.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    if not r["submitted"]:
        lines.append(f"| {h} | N/A | No uploaded rows available | N/A | N/A |")
    else:
        lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"The uploaded fleet union contains **{len(d['fleet_clean_union'])}/900 clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** "
          f"({sum(len(r['band_seeds']) for r in uploaded.values())} band memberships before overlaps). "
          "h6's prepared band is excluded from these totals.", "",
          "## Per-seed results", "",
          "Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in uploaded.items():
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
          "h4's strong score at its HDR hit was averaged with two low seed scores. "
          "The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the "
          "inspected data; its aggregate score does not identify those sets.", "",
          "## h6: uploaded fallback versus prepared band", "",
          "h6 began preparation at **13:13:36 UTC**. The validator requested its submission at **13:15:06**. "
          "At **13:15:09**, it waited up to **177 seconds** for the prepared build. At **13:18:06**, that build "
          "was still running (270 seconds elapsed), so the miner used ordinary construction while preparation "
          "still held the GPU. It uploaded at **13:18:22.734**, with **88 seconds** of URL TTL left. "
          "The prepared conjunction finished at **13:19:43.059**, about **80 seconds after the upload**.", "",
          f"That unsubmitted preparation logged **26 clean seeds within 400–699** and recorded band "
          f"**{nums(ours['h6']['latest_record_band'])}**. The actual fallback has **zero clean seeds and no HDR band "
          "across 100–999**. Its window record combines `source=all_hdr_not_attempted` / `window=null` with "
          "a band from the completed background build; the band does not describe uploaded rows. "
          "All six other uploaded bands match their task window records.", "",
          "## h1–h3: no task activity", "",
          "No task preparation, received-task/upload event, window record or upload archive was found for "
          "h1–h3. None appears in the complete official score feed. Their earlier runtime logs reported them "
          "unregistered on netuid 55 on 2026-09-19. No submitted clean set, band, score or rank is available "
          "for these hotkeys on this task.", "",
          "## Exact uploaded clean sets", "",
          "Each list covers all 100–999 seeds; band seeds are included. h6's set is empty.", ""]
for h, r in uploaded.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
ca_extra = sorted(set(ours["h7"]["cas12a_clean_seeds"]) - set(ours["h7"]["clean_seeds"]))
lines += ["## Verification and artifacts", "",
          "- All seven archives identify this task; contracts and references match the published task after the seed stamp.",
          "- All 1,750 submitted rows passed stage 1/2, with accessibility 0.35. Local three-seed scores exactly match all seven official scores.",
          "- Reused h6/h8/h9 archived validations and generated h0/h4/h5/h7 validations in report-local directories.",
          "- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.",
          "- Full-row clean counts inside each conjunction's cut window match the original logged Cas12a clean counts.",
          f"- h7 has a Cas12a-only clean seed outside its cut window, **{nums(ca_extra)}**, where the full 250-row submission is not clean. It is excluded from the reported full-row clean set.",
          "- Recovered sets over **900 × 250 × 7 = 1,575,000 row/seed combinations**.",
          "- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), "
          "[sanitized build/upload and registration evidence](build_evidence.json), [layout source evidence](layout_evidence.json), "
          "[manifest](manifest.json), [CSV](summary.csv).",
          "- Reproduce with `.venv/bin/python reports/31138352/validate.py`, then "
          "`.venv/bin/python reports/31138352/audit.py`, then `.venv/bin/python reports/31138352/render_report.py`. "
          "The analysis includes validator-code and submission SHA-256 hashes.", "",
          "| Hotkey | Public hotkey | Upload time UTC | Submission archive |",
          "|---|---|---|---|"]
for h, r in ours.items():
    archive = f"[{r['archive']}](../../{r['archive']}/submission.json)" if r["submitted"] else "No upload archive"
    lines.append(f"| {h} | `{r['hotkey']}` | {r.get('uploaded_at', 'N/A')} | {archive} |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Uploaded", "Valid rows", "Construction", "Assigned band/cut window",
                     "Window applied to upload", "Clean count /900", "Exact submitted clean set", "Clean hits", "Uploaded HDR band",
                     "HDR hits", "Percent of leader", "Gap to leader", "Unsubmitted prepared band"])
    for h, r in ours.items():
        yes = r["submitted"]
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"] if yes else "N/A", yes, r["rows"],
                         r["construction"], r["candidate_window_ranges"], r.get("candidate_window_applies_to_upload"),
                         len(r["clean_seeds"]) if yes else "N/A", nums(r["clean_seeds"]), nums(r["clean_hits"]),
                         nums(r["band_seeds"]), nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"],
                         nums(r["latest_record_band"]) if h == "h6" else ""])
print(OUT / "report.md")
