"""Replay the ten actual uploads, including both generations used for task 09b7b842."""
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
TASK_ID = "09b7b842-ccff-4ea5-ab63-7486437e1b59"
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
    seeds = [int(x) for x in task["content"]["contract"]["seed"].split(",")]
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
    predicted = list(range(100, 400))
    shared = list(range(100, 300)) + list(range(500, 600))
    result = {
        "task_id": TASK_ID, "task_created_at": task["created_at"],
        "cell_type": task["content"]["contract"]["cell_type"], "seeds": seeds,
        "generated_at": datetime.now(timezone.utc).isoformat(), "fetched_at": manifest["fetched_at"],
        "source_api": f"https://niome-api.genomes.io/api/v3/miners/scores?task_id={TASK_ID}&limit=40000",
        "score_records": len(scores), "unique_miners": len(ranked),
        "validators": sorted({r["validator_hotkey"] for r in scores}), "top_miner": ranked[0],
        "top10_cutoff": ranked[9]["final_score"], "top10": ranked[:10],
        "rank_definition": "1 plus the number of strictly higher scores; ties share rank.",
        "predicted_windows": [100, 200, 300], "predicted_joined": predicted,
        "shared_joined": shared, "shared_windows": [100, 200, 500],
        "original_plan_generated_at_log_precision": "2026-09-21T06:04",
        "rebuild_plan_generated_at": "2026-09-21T06:17:51.904738+00:00",
        "code_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in
                        ("calc.py", "niome_subnet/genomics/validation/stage12.py", "niome_subnet/genomics/validation/stage3.py",
                         "niome_subnet/genomics/validation/stage4.py", "niome_subnet/genomics/validation/stage5.py")},
        "definitions": {
            "clean_seeds": "All 250 uploaded rows cut at this seed, replayed over 100–999.",
            "band_seeds": "All 250 uploaded rows return HDR at this seed, replayed over 100–999.",
            "candidate_window": "Actual uploaded build's band candidates, reconstructed from build logs and checked against replayed bands. Original: h0 100–399; h1 100–199; h2 200–299; h3 300–399; h9 700–999. Rebuilt h4–h8: fixed shared 100–299 ∪ 500–599, loops 5–9.",
            "clean_in_cut_space": "Fully clean seeds restricted to the construction's cut search space. Build logs report Cas12a cleanliness; both are checked separately.",
        },
        "hotkeys": {},
    }
    for i in range(10):
        h = f"h{i}"
        entry, ev = manifest["hotkeys"][h], evidence[h]
        score = best[entry["hotkey"]]
        original = entry["uploaded_generation"] == "original"
        assert original == (i in (0, 1, 2, 3, 9))
        joined = (predicted if i in (1, 2, 3) else list(range(100, 1000))) if original else shared
        group = ("predicted" if i in (1, 2, 3) else "full900") if original else "shared"
        offset = {0: 0, 1: 0, 2: 100, 3: 200, 9: 600}.get(i) if original else None
        candidates = (list(range(i * 100, (i + 1) * 100)) if i in (1, 2, 3)
                      else list(range(700, 1000)) if i == 9 else predicted) if original else shared
        cut = (list(range(700, 1000)) if i == 9 else predicted) if original else list(range(100, 1000))
        events = ev["events"]
        selected_events = [e for e in events if (e[:19] < "2026-09-21 06:30:00") == original]
        window_lines = [e for e in selected_events if "Build: window " in e]
        assert len(window_lines) == 1 and "window 100-399" in window_lines[0]
        candidate_lines = [e for e in selected_events if "Build: band candidates " in e]
        assert len(candidate_lines) == 1
        if original:
            assert f"Build: band candidates {min(candidates)}-{max(candidates)} ({len(candidates)} seeds," in candidate_lines[0]
            assert f"space of {len(joined)}, offset {offset}) | cut 300 seeds (narrow) | k=8" in candidate_lines[0]
        else:
            assert f"300 seeds (shared window, loop {i + 1} of 10) | cut 900 seeds (wide) | k=8" in candidate_lines[0]
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
        assert {e["features"]["cell_type_accessibility"] for e in valid} == {0.35}
        assert validation["scored_seeds"] == seeds and validation["contract_pairing"]["verified"]
        contract = read(folder / "contract.json")
        assert {k: v for k, v in contract.items() if k != "seed"} == {k: v for k, v in task["content"]["contract"].items() if k != "seed"}
        reference = read(folder / "hbb_reference.json")
        reference["challenge"]["seed"] = task["content"]["hbb_reference"]["challenge"]["seed"]
        assert reference == task["content"]["hbb_reference"]
        assert abs(validation["score"]["final_score"] - score["final_score"]) < 1e-10
        assert abs(sum(validation["per_seed_final_score"].values()) / len(seeds) - score["final_score"]) < 1e-10
        clean, band, ca_clean, seed_results = [], [], [], {}
        for seed in range(100, 1000):
            simulated = [stage3.simulate(e, seed) for e in valid]
            counts = Counter(r["outcome"] for r in simulated)
            if counts["no_cut"] == 0:
                clean.append(seed)
            if counts["HDR"] == len(rows):
                band.append(seed)
            if all(r["outcome"] != "no_cut" for r in simulated if r["cas"] == "Cas12a"):
                ca_clean.append(seed)
            if seed in seeds:
                detail = read(validation_dir / f"stage3_detail.seed{seed}.json")
                old = {r["experiment_id"]: r for r in detail["experiments"]}
                assert all(r["outcome"] == old[r["experiment_id"]]["outcome"] and r["indel_length"] == old[r["experiment_id"]]["indel_length"] for r in simulated)
                seed_results[str(seed)] = {
                    "outcomes": dict(counts),
                    "outcomes_by_cas": {cas: dict(Counter(r["outcome"] for r in simulated if r["cas"] == cas)) for cas in ("Cas9", "Cas12a")},
                    "in_predicted_space": seed in predicted, "in_band_space": seed in joined,
                    "in_candidate_window": seed in candidates, "in_cut_space": seed in cut,
                    "clean": seed in clean, "band": seed in band,
                    "final_score": validation["per_seed_final_score"][str(seed)],
                }
        assert set(band) <= set(clean) <= set(ca_clean)
        cas_mix = dict(Counter(r["cas_system"] for r in rows))
        builds = [e for e in selected_events if "Build: conjunction " in e]
        assert len(builds) == 1
        logged = tuple(map(int, re.search(r"band (\d+) k=(\d+) group (\d+) clean (\d+)/(\d+)", builds[0]).groups()))
        clean_in_cut = sorted(set(clean) & set(cut))
        ca_clean_in_cut = sorted(set(ca_clean) & set(cut))
        assert logged == (8, 8, 80, len(ca_clean_in_cut), len(cut))
        assert clean_in_cut == ca_clean_in_cut
        assert uploaded_events[0][:23] > builds[0][:23]
        if original:
            assert uploaded_events[0][:23] < ev["current_window_record"]["at"].replace("T", " ").replace(".", ",")[:23]
        assert len(band) == 8 and set(band) <= set(candidates)
        assert cas_mix == {"Cas12a": 80, "Cas9": 170}
        record_band = sorted(json.loads(ev["current_window_record"]["band"]))
        assert (record_band == band) == (not original)
        rec = {
            "hotkey": entry["hotkey"], "uid": score["miner_uid"], "construction": "conjunction",
            "archive": entry["archive"], "validation_dir": entry["validation_dir"],
            "submitted": True, "uploaded_at": uploaded_events[0][:23],
            "uploaded_generation": entry["uploaded_generation"], "build_completed_at": builds[0][:23],
            "build_evidence": {"window": window_lines[0], "candidates": candidate_lines[0], "completed": builds[0]},
            "rows": len(rows), "cas_mix": cas_mix, "score": score, "rank": ranks[entry["hotkey"]],
            "gap_to_top": ranked[0]["final_score"] - score["final_score"],
            "percent_of_top": 100 * score["final_score"] / ranked[0]["final_score"],
            "band_group": group, "band_space": joined, "band_space_ranges": ranges(joined),
            "candidate_offset": offset, "band_loop": None if original else i + 1,
            "candidate_window": candidates, "candidate_window_ranges": ranges(candidates),
            "cut_search_space": cut, "cut_search_ranges": ranges(cut),
            "logged_clean_count": logged[3], "clean_seeds_in_cut_space": clean_in_cut,
            "clean_seeds": clean, "cas12a_clean_seeds": ca_clean, "band_seeds": band,
            "clean_seeds_in_predicted_space": sorted(set(clean) & set(predicted)),
            "clean_seeds_in_band_space": sorted(set(clean) & set(joined)),
            "clean_seeds_in_candidate_window": sorted(set(clean) & set(candidates)),
            "clean_hits": sorted(set(clean) & set(seeds)), "band_hits": sorted(set(band) & set(seeds)),
            "latest_record_at": ev["current_window_record"]["at"], "latest_record_band_matches_upload": record_band == band,
            "latest_record_band": record_band,
            "submission_sha256": hashlib.sha256((folder / "submission.json").read_bytes()).hexdigest(),
            "per_seed": seed_results, "local_score_matches_api": True, "stage3_replay_matches_local_detail": True,
        }
        result["hotkeys"][h] = rec
        (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f"{h}: rank={rec['rank']} final={score['final_score']:.8f} clean={len(clean)} band={band} clean_hits={rec['clean_hits']} HDR_hits={rec['band_hits']}", flush=True)
    records = result["hotkeys"].values()
    result["fleet_clean_union"] = sorted(set().union(*(set(r["clean_seeds"]) for r in records)))
    result["fleet_band_union"] = sorted(set().union(*(set(r["band_seeds"]) for r in records)))
    result["fleet_clean_hits"] = sorted(set(seeds) & set(result["fleet_clean_union"]))
    result["fleet_band_hits"] = sorted(set(seeds) & set(result["fleet_band_union"]))
    result["best_our_hotkey"] = max(result["hotkeys"], key=lambda h: result["hotkeys"][h]["score"]["final_score"])
    result["duplicate_band_seeds"] = {str(s): [h for h, r in result["hotkeys"].items() if s in r["band_seeds"]]
                                      for s, n in Counter(s for r in records for s in r["band_seeds"]).items() if n > 1}
    (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Fleet unions:", len(result["fleet_clean_union"]), "clean;", len(result["fleet_band_union"]), "HDR band", flush=True)


if __name__ == "__main__":
    main()
