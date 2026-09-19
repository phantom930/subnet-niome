"""Render the cef21622 audit from its saved evidence and replay results."""
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
evidence = json.loads((OUT / "build_evidence.json").read_text())
top = d["top_miner"]
ours = d["hotkeys"]
seeds = d["seeds"]
best_h = max(ours, key=lambda h: ours[h]["score"]["final_score"])
best = ours[best_h]
gap10 = d["top10_cutoff"] - best["score"]["final_score"]


def nums(values):
    return ", ".join(map(str, values)) or "none"


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {score['weight']:.6f} |")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"{d['cell_type']} · task opened {d['task_created_at']} UTC · scoring seeds **{nums(seeds)}**.", "",
    f"**Best of h0–h5: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Top miner: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} reached {best['percent_of_top']:.2f}% of the top score, a gap of {best['gap_to_top']:.6f} points. "
    f"It missed the top-10 cutoff ({d['top10_cutoff']:.6f}) by **{gap10:.6f}** points. "
    "All six submissions were uploaded successfully, contained 250 valid rows, and received zero reward weight.", "",
    f"Scores are from the [task-specific public API]({d['source_api']}); "
    f"{d['score_records']} records, {d['unique_miners']} distinct miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated: {d['generated_at']}.", "",
    "## Official score comparison", "",
    f"| Miner | UID this round | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", f"Top miner hotkey: `{top['miner_hotkey']}`.", "",
          "Final score is the mean of the per-seed scores. For these submissions weighted score and fidelity "
          "are constant across seeds, so final = weighted score × mean consistency × fidelity.", "",
          "## Joined seed windows and recovered bands", "",
          "All six original builds used joined band space **300–499 ∪ 800–899** (300 seeds), "
          "from the plan generated at 2026-09-17 21:20 UTC. "
          "The separate cut search covered **100–999** (900 seeds). "
          "Band candidate slices were 150 seeds wide, at offsets 0/50/100/150/200/250 for h0–h5. "
          "Every submission has 100 Cas12a and 150 Cas9 rows and an 11-seed HDR band.", "",
          "**Clean set:** every submitted row cuts at that seed. **HDR band:** every submitted row returns HDR. "
          "The band is a subset of the clean set. Both sets were recovered by simulating the exact uploaded rows "
          "at every seed 100–999; the recovered bands match the original recorded bands exactly. "
          "The full-row clean sets equal the Cas12a group clean sets, and their sizes match the build logs.", "",
          "| Hotkey | Candidate slice within joined space | Clean / 900 | Actual HDR band seeds |",
          "|---|---|---:|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['candidate_window_ranges']} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} |")
lines += ["", "Candidate slices are reconstructed from the original logged joined space and the six-hotkey "
          "offset layout; every recovered band lies inside its corresponding slice.", "",
          f"Fleet unions: **{len(d['fleet_clean_union'])} distinct clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR band seeds** "
          f"({sum(len(r['band_seeds']) for r in ours.values())} band memberships before overlap). "
          f"The scoring seeds covered by the fleet's clean union are {nums(d['fleet_clean_hits'])}; "
          f"the band union hits {nums(d['fleet_band_hits'])}.", "",
          "## Hits and per-seed scores", "",
          "| Hotkey | Clean seed hits (includes HDR hits) | HDR band hits |",
          "|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", "All three scoring seeds, **115, 977, and 940**, are outside the joined band space. "
          "None hits any hotkey's HDR band. The wider 900-seed cut search supplies exactly one clean hit "
          "among these six hotkeys: **h0 at seed 115**. Neither 977 nor 940 is clean for any hotkey.", "",
          "Each cell below is **per-seed final score / no-cut rows**.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    values = [f"{r['per_seed'][str(s)]['final_score']:.6f} / {r['per_seed'][str(s)]['outcomes'].get('no_cut', 0)}"
              for s in seeds]
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
lines += ["", "h0's clean hit at seed 115 scores **42.100869**, but it is not an HDR band hit. "
          "Its scores at 977 and 940 are **22.283709** and **23.856806**, respectively, giving a mean of "
          "**29.413795**. All remaining hotkey/seed combinations contain no-cut rows.", "",
          "## Comparison with the leader", "",
          "h0's fidelity is essentially equal to the leader's (0.937133 versus 0.936800), and its weighted "
          "score is only modestly lower (**249.058314 versus 256.752820**). The dominant gap is mean "
          "consistency: **0.126023 versus 0.631912**. Across all six hotkeys, consistency ranges from "
          "0.083607 to 0.126023. The observed clean or band hits explain the fleet's "
          "per-seed outcomes on this task; one task does not establish long-run strategy performance. "
          "The leader's clean sets, bands, and per-seed outcomes are not available from its aggregate score "
          "and were not reconstructed.", "",
          "## Exact clean sets", "",
          "These lists include the HDR band seeds and cover the full 100–999 seed space. "
          "Machine-readable lists are also in [analysis.json](analysis.json).", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Submission provenance and verification", "",
          "The six builds started at 22:08:37 UTC and finished between 22:10:56 and 22:11:06 UTC "
          "(139–149 seconds). All six hotkeys served the validator's requests from prepared submissions. "
          "Uploads succeeded between 22:13:18 and 22:57:58 UTC, with 266–288 seconds of upload TTL remaining. "
          "Each hotkey has one recorded build and one recorded successful upload for this task. "
          "The recorded bands agree with the archived rows' full-seed replay.", "",
          "- All archived contracts and reference files match the published task after accounting for the "
          "seed changing from 0 at submission time to the three stamped scoring seeds.",
          "- All 250 archived stage-1/2 valid-row entries match each submission exactly; accessibility is 0.87.",
          "- All six archived three-seed final scores match the official API exactly.",
          "- Replayed stage-3 outcomes and indel lengths match the archived detail for every row at all three real seeds.",
          "- Scanned 900 seeds × 250 rows × 6 hotkeys = 1,350,000 row simulations. "
          "No miner configuration, running process, or uploaded submission was changed.",
          "- [Full analysis](analysis.json), [official score snapshot](scores.json), "
          "[published task](task.json), [sanitized build/upload evidence](build_evidence.json).",
          "- Reproduce with `.venv/bin/python reports/cef21622/audit.py` and "
          "`.venv/bin/python reports/cef21622/render_report.py` from the repository root.", "",
          "Uploaded archives:", ""]
for h, r in ours.items():
    lines.append(f"- **{h}** `{r['hotkey']}` — [{r['archive']}](../../{r['archive']}/submission.json)")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
print(OUT / "report.md")
