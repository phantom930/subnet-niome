"""Replay seven conjunctions and three all-HDR uploads for task e6caf784."""
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TASK_ID = "e6caf784-f2af-48cc-938a-1377c441a1dd"
sys.path.insert(0, str(ROOT))


def read(path):
    return json.loads(Path(path).read_text())


def ranges(values):
    groups = []
    for value in sorted(set(values)):
        if groups and value == groups[-1][-1] + 1:
            groups[-1].append(value)
        else:
            groups.append([value])
    return " ∪ ".join(str(g[0]) if len(g) == 1 else f"{g[0]}–{g[-1]}" for g in groups)


def main():
    spec = importlib.util.spec_from_file_location("audit_stage3", ROOT / "niome_subnet/genomics/validation/stage3.py")
    stage3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage3)
    task, manifest, evidence = (read(OUT / f"{f}.json") for f in ("task", "manifest", "build_evidence"))
    stamp = task["content"]["contract"]["seed"]
    seeds = [int(x) for x in str(stamp).split(",")]
    assert seeds == [481, 208, 809]
    validation_seeds = seeds
    snapshot = read(OUT / "scores.json")
    scores = [r for r in snapshot["items"] if r["task_id"] == TASK_ID]
    assert len(scores) == snapshot["pagination"]["total"]
    best = {}
    for row in scores:
        hk = row["miner_hotkey"]
        if hk not in best or row["final_score"] > best[hk]["final_score"]:
            best[hk] = row
    ranked = sorted(best.values(), key=lambda r: -r["final_score"])
    ranks = {hk: 1 + sum(x["final_score"] > r["final_score"] for x in ranked) for hk, r in best.items()}
    conjunction_space = list(range(200, 300))
    hdr_space = list(range(400, 500))
    conjunction_cut = conjunction_space + hdr_space
    result = {
        "task_id": TASK_ID, "task_created_at": task["created_at"],
        "cell_type": task["content"]["contract"]["cell_type"], "seeds": seeds,
        "api_seed_stamp": stamp, "scoring_seeds_available": seeds is not None,
        "generated_at": datetime.now(timezone.utc).isoformat(), "fetched_at": manifest["fetched_at"],
        "source_api": f"https://niome-api.genomes.io/api/v3/miners/scores?task_id={TASK_ID}&limit=40000",
        "score_records": len(scores), "unique_miners": len(ranked),
        "validators": sorted({r["validator_hotkey"] for r in scores}), "top_miner": ranked[0],
        "top10_cutoff": ranked[9]["final_score"], "top10": ranked[:10],
        "rank_definition": "1 plus the number of strictly higher scores; ties share rank.",
        "conjunction_band_joined": conjunction_space, "all_hdr_band_joined": hdr_space,
        "shared_cut_joined": conjunction_cut, "layout_source": "hdr_fixed_layout",
        "plan_generated_at": None,
        "code_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in
                        ("calc.py", "niome_subnet/genomics/validation/stage12.py", "niome_subnet/genomics/validation/stage3.py",
                         "niome_subnet/genomics/validation/stage4.py", "niome_subnet/genomics/validation/stage5.py")},
        "definitions": {
            "clean_seeds": "All 250 uploaded rows cut at this seed, replayed over 100–999.",
            "band_seeds": "All 250 uploaded rows return HDR at this seed, replayed over 100–999.",
            "candidate_window": "h0–h6: width-15 windows at stride 15 over 200–299, loop 1, cut 200–299 ∪ 400–499. h7–h9: all-HDR search over width-34 windows at stride 34 over 400–499. h6 and h9 wrap.",
            "clean_in_cut_space": "Fully clean seeds restricted to the construction's cut search space. Build logs report Cas12a cleanliness; both are checked separately.",
            "all_hdr_log_clean": "For h7–h9, log clean count/34 is the Cas12a all-HDR intersection in the 34-seed candidate window, not the full-row no-cut-free set. Cas9 rows are then selected to return HDR on those HDR seeds.",
        },
        "hotkeys": {},
    }
    for i in range(10):
        h = f"h{i}"
        entry, ev = manifest["hotkeys"][h], evidence[h]
        score = best[entry["hotkey"]]
        is_conjunction = i < 7
        assert entry["construction"] == ("conjunction" if is_conjunction else "all_hdr")
        joined = conjunction_space if is_conjunction else hdr_space
        offset, width = (15 * i, 15) if is_conjunction else (34 * (i - 7), 34)
        candidates = sorted({joined[(offset + j) % len(joined)] for j in range(width)})
        cut = conjunction_cut if is_conjunction else []
        events = ev["events"]
        window_lines = [e for e in events if "Build: window " in e]
        assert len(window_lines) == 1
        assert f"Build: window {ranges(candidates).replace('–', '-').replace(' ∪ ', ',')} from the fixed all-HDR layout" in window_lines[0]
        candidate_lines = [e for e in events if "Build: band candidates " in e]
        if is_conjunction:
            assert len(candidate_lines) == 1
            assert "15 seeds (shared window, loop 1 of 7) | cut 200 seeds (narrow) | k=11" in candidate_lines[0]
        else:
            assert not candidate_lines
        assert any("using the prepared submission (250 rows" in e for e in events)
        uploaded_events = [e for e in events if f"Submitted 250 rows for task {TASK_ID}" in e]
        assert len(uploaded_events) == 1
        folder, validation_dir = ROOT / entry["archive"], ROOT / entry["validation_dir"]
        rows = read(folder / "submission.json")
        valid = read(validation_dir / "validation/valid_experiments.json")
        validation = read(validation_dir / "validation.json")
        upload = read(folder / "last_upload.json")
        assert upload["submitted"] and upload["task_id"] == TASK_ID and len(rows) == 250
        assert rows == [e["experiment"] for e in valid]
        assert {e["features"]["cell_type_accessibility"] for e in valid} == {0.82}
        assert validation["scored_seeds"] == validation_seeds and validation["contract_pairing"]["verified"]
        contract = read(folder / "contract.json")
        assert {k: v for k, v in contract.items() if k != "seed"} == {k: v for k, v in task["content"]["contract"].items() if k != "seed"}
        reference = read(folder / "hbb_reference.json")
        reference["challenge"]["seed"] = task["content"]["hbb_reference"]["challenge"]["seed"]
        assert reference == task["content"]["hbb_reference"]
        for key in ("total_weighted_score", "distribution_fidelity_factor"):
            assert abs(validation["score"][key] - score["breakdown"][key]) < 1e-9, (h, key)
        if seeds is not None:
            assert abs(validation["score"]["final_score"] - score["final_score"]) < 1e-10
            assert abs(sum(validation["per_seed_final_score"].values()) / len(seeds) - score["final_score"]) < 1e-10
        clean, band, ca_clean, ca_hdr, seed_results = [], [], [], [], {}
        for seed in range(100, 1000):
            simulated = [stage3.simulate(e, seed) for e in valid]
            counts = Counter(r["outcome"] for r in simulated)
            if counts["no_cut"] == 0:
                clean.append(seed)
            if counts["HDR"] == len(rows):
                band.append(seed)
            if all(r["outcome"] != "no_cut" for r in simulated if r["cas"] == "Cas12a"):
                ca_clean.append(seed)
            if all(r["outcome"] == "HDR" for r in simulated if r["cas"] == "Cas12a"):
                ca_hdr.append(seed)
            if seed in validation_seeds:
                detail_path = validation_dir / (f"stage3_detail.seed{seed}.json" if len(validation_seeds) > 1 else "stage3_detail.json")
                detail = read(detail_path)
                old = {r["experiment_id"]: r for r in detail["experiments"]}
                assert all(r["outcome"] == old[r["experiment_id"]]["outcome"] and r["indel_length"] == old[r["experiment_id"]]["indel_length"] for r in simulated)
                seed_results[str(seed)] = {
                    "outcomes": dict(counts),
                    "outcomes_by_cas": {cas: dict(Counter(r["outcome"] for r in simulated if r["cas"] == cas)) for cas in ("Cas9", "Cas12a")},
                    "in_band_space": seed in joined,
                    "in_candidate_window": seed in candidates, "in_cut_space": seed in cut,
                    "clean": seed in clean, "band": seed in band,
                    "final_score": validation["per_seed_final_score"][str(seed)] if len(validation_seeds) > 1 else validation["score"]["final_score"],
                }
        assert set(band) <= set(clean) <= set(ca_clean)
        cas_mix = dict(Counter(r["cas_system"] for r in rows))
        builds = [e for e in events if ("Build: conjunction " if is_conjunction else "Build: all-HDR ") in e]
        assert len(builds) == 1
        clean_in_cut = sorted(set(clean) & set(cut))
        ca_clean_in_cut = sorted(set(ca_clean) & set(cut))
        if is_conjunction:
            logged = tuple(map(int, re.search(r"band (\d+) k=(\d+) group (\d+) clean (\d+)/(\d+)", builds[0]).groups()))
            assert logged == (11, 11, 100, len(ca_clean_in_cut), 200)
            assert clean_in_cut == ca_clean_in_cut
            logged_count, logged_span = logged[3:5]
            record_band = sorted(json.loads(ev["current_window_record"]["band"]))
            assert record_band == band
        else:
            logged = tuple(map(int, re.search(r"group (\d+) clean (\d+)/(\d+)", builds[0]).groups()))
            assert logged == (cas_mix['Cas12a'], len(set(ca_hdr) & set(candidates)), 34)
            assert sorted(set(ca_hdr) & set(candidates)) == band
            logged_count, logged_span = logged[1:3]
            label = ranges(candidates).replace('–', '-').replace(' ∪ ', ',')
            expected_record = label if i < 9 else "joined 34 seeds " + label
            assert ev["current_window_record"]["band"] == expected_record
            record_band = None
        assert uploaded_events[0][:23] > builds[0][:23]
        assert len(band) == (11 if is_conjunction else logged_count) and set(band) <= set(candidates)
        assert sum(cas_mix.values()) == 250 and cas_mix['Cas12a'] == (100 if is_conjunction else logged[0])
        rec = {
            "hotkey": entry["hotkey"], "uid": score["miner_uid"], "construction": entry["construction"],
            "archive": entry["archive"], "validation_dir": entry["validation_dir"],
            "submitted": True, "uploaded_at": uploaded_events[0][:23],
            "build_completed_at": builds[0][:23],
            "prepared_build_completed_at": builds[0][:23], "prepared_build_uploaded": True,
            "prepared_build_evidence": {"window": window_lines[0], "candidates": candidate_lines[0] if candidate_lines else None, "completed": builds[0]},
            "rows": len(rows), "cas_mix": cas_mix, "score": score, "rank": ranks[entry["hotkey"]],
            "gap_to_top": ranked[0]["final_score"] - score["final_score"],
            "percent_of_top": 100 * score["final_score"] / ranked[0]["final_score"],
            "band_group": "conjunction_15_seed_window" if is_conjunction else "all_hdr_34_seed_window", "band_space": joined, "band_space_ranges": ranges(joined),
            "candidate_offset": offset, "band_loop": 1 if is_conjunction else None,
            "candidate_window": candidates, "candidate_window_ranges": ranges(candidates),
            "cut_search_space": cut, "cut_search_ranges": ranges(cut),
            "prepared_logged_clean_count": logged_count, "logged_clean_denominator": logged_span,
            "logged_clean_definition": "Cas12a no-cut-free inside cut window" if is_conjunction else "Cas12a all-HDR inside candidate window",
            "clean_seeds_in_cut_space": clean_in_cut,
            "clean_seeds": clean, "cas12a_clean_seeds": ca_clean, "band_seeds": band,
            "cas12a_hdr_seeds": ca_hdr,
            "clean_seeds_in_band_space": sorted(set(clean) & set(joined)),
            "clean_seeds_in_candidate_window": sorted(set(clean) & set(candidates)),
            "clean_hits": sorted(set(clean) & set(seeds)) if seeds is not None else None,
            "band_hits": sorted(set(band) & set(seeds)) if seeds is not None else None,
            "latest_record_at": ev["current_window_record"]["at"], "latest_record_band_matches_upload": record_band == band if is_conjunction else None,
            "latest_record_band": record_band,
            "latest_record_band_field": ev["current_window_record"]["band"],
            "latest_record_band_semantics": "actual HDR seeds" if is_conjunction else "candidate window; not the actual HDR seeds",
            "submission_sha256": hashlib.sha256((folder / "submission.json").read_bytes()).hexdigest(),
            "per_seed": seed_results if seeds is not None else None,
            "diagnostic_seed_results": seed_results if seeds is None else None,
            "local_score_matches_api": True if seeds is not None else None,
            "seed_independent_factors_match_api": True,
            "stage3_replay_matches_local_detail": True,
        }
        result["hotkeys"][h] = rec
        (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f"{h}: rank={rec['rank']} final={score['final_score']:.8f} clean={len(clean)} band={band} clean_hits={rec['clean_hits']} HDR_hits={rec['band_hits']}", flush=True)
    records = result["hotkeys"].values()
    result["fleet_clean_union"] = sorted(set().union(*(set(r["clean_seeds"]) for r in records)))
    result["fleet_band_union"] = sorted(set().union(*(set(r["band_seeds"]) for r in records)))
    result["fleet_clean_hits"] = sorted(set(seeds) & set(result["fleet_clean_union"])) if seeds is not None else None
    result["fleet_band_hits"] = sorted(set(seeds) & set(result["fleet_band_union"])) if seeds is not None else None
    result["best_our_hotkey"] = max(result["hotkeys"], key=lambda h: result["hotkeys"][h]["score"]["final_score"])
    result["duplicate_band_seeds"] = {str(s): [h for h, r in result["hotkeys"].items() if s in r["band_seeds"]]
                                      for s, n in Counter(s for r in records for s in r["band_seeds"]).items() if n > 1}
    assert set(result["fleet_band_union"]) <= set(conjunction_cut)
    (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Fleet unions:", len(result["fleet_clean_union"]), "clean;", len(result["fleet_band_union"]), "HDR band", flush=True)


if __name__ == "__main__":
    main()
