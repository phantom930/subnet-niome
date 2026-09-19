"""Audit the uploaded h0-h6 archives for task 11276911, without changing miner data.

Run from the repository root with .venv/bin/python reports/11276911/audit.py.
Uses saved API responses and the validator's stage-3 simulator. No builds or uploads.
"""
import hashlib
import importlib.util
import json
import re
from pathlib import Path
import sys
import time
from collections import Counter
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
OUT = Path(__file__).resolve().parent
TASK_ID = "11276911-854f-4809-8610-e7d5c7bc8772"
HOTKEYS = [
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW",
    "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU",
    "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2",
    "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",
    "5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3",
]


def read(path):
    return json.loads(Path(path).read_text())


def ranges(seeds):
    groups = []
    for seed in sorted(set(seeds)):
        if groups and seed == groups[-1][-1] + 1:
            groups[-1].append(seed)
        else:
            groups.append([seed])
    return ", ".join(str(g[0]) if len(g) == 1 else f"{g[0]}–{g[-1]}" for g in groups)


def main():
    spec = importlib.util.spec_from_file_location(
        "audit_stage3", ROOT / "niome_subnet/genomics/validation/stage3.py")
    stage3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage3)
    task = read(OUT / "task.json")
    evidence = read(OUT / "build_evidence.json")
    manifest = read(OUT / "manifest.json")
    seeds = [int(s) for s in task["content"]["contract"]["seed"].split(",")]
    raw_scores = read(OUT / "scores.json")
    scores = [r for r in raw_scores["items"] if r["task_id"] == TASK_ID]
    best = {}
    for row in scores:
        hk = row["miner_hotkey"]
        if hk not in best or row["final_score"] > best[hk]["final_score"]:
            best[hk] = row
    ranked = sorted(best.values(), key=lambda r: -r["final_score"])
    positions = {r["miner_hotkey"]: i for i, r in enumerate(ranked, 1)}
    ranks = {r["miner_hotkey"]: 1 + sum(x["final_score"] > r["final_score"] for x in ranked)
             for r in ranked}
    joined = list(range(300, 400)) + list(range(600, 800))
    cut_space = list(range(100, 1000))
    result = {
        "task_id": TASK_ID, "task_created_at": task["created_at"],
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cell_type": task["content"]["contract"]["cell_type"], "seeds": seeds,
        "source_api": f"https://niome-api.genomes.io/api/v3/miners/scores?task_id={TASK_ID}&limit=40000",
        "score_records": len(scores), "unique_miners": len(ranked),
        "validators": sorted({r["validator_hotkey"] for r in scores}),
        "rank_definition": "Competition rank: 1 plus number of miners with a strictly higher score. Exact ties share rank.",
        "code_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in
                        ("calc.py", "niome_subnet/genomics/validation/stage12.py",
                         "niome_subnet/genomics/validation/stage3.py", "niome_subnet/genomics/validation/stage4.py",
                         "niome_subnet/genomics/validation/stage5.py", "joined_window.py")},
        "top_miner": ranked[0], "top10_cutoff": ranked[9]["final_score"],
        "fetched_at": manifest["fetched_at"],
        "joined_band_space": joined, "cut_space": cut_space, "scan_space": [100, 999],
        "definitions": {
            "clean_seeds": "All 250 submitted rows cut; HDR band is a subset.",
            "band_seeds": "All 250 submitted rows return HDR, scanned over every seed 100-999.",
            "cas12a_clean_seeds": "All submitted Cas12a rows cut over 100-999 (100 rows for conjunction; 74 for h0/h6 fallbacks).",
            "clean_seeds_in_joined_space": "Full-row clean seeds inside the 300-seed joined band space.",
            "clean_seeds_outside_joined_space": "Full-row clean seeds outside the band space, within the 900-seed cut search space.",
            "candidate_window": "Planned 150-seed slice reconstructed from logged joined space and offsets 0/50/100/150/200/250/25. h0/h6 submitted ordinary fallbacks instead.",
        },
        "hotkeys": {},
    }
    for i, hk in enumerate(HOTKEYS):
        name = "niome_hotkey" + (str(i) if i else "")
        archives = [p.parent for p in (ROOT / f"data/inst/{name}/result").glob("2026-09-18T*/last_upload.json")
                    if read(p).get("task_id") == TASK_ID]
        assert len(archives) == 1, (name, archives)
        folder = archives[0]
        validation_dir = ROOT / manifest["hotkeys"][f"h{i}"]["validation_dir"]
        assert folder == ROOT / manifest["hotkeys"][f"h{i}"]["archive"]
        rows = read(folder / "submission.json")
        valid = read(validation_dir / "validation/valid_experiments.json")
        validation = read(validation_dir / "validation.json")
        upload = read(folder / "last_upload.json")
        assert upload["submitted"] and len(rows) == 250
        assert rows == [e["experiment"] for e in valid]
        assert {e["features"]["cell_type_accessibility"] for e in valid} == {0.87}
        assert validation["scored_seeds"] == seeds
        assert validation["contract_pairing"]["verified"]
        contract = read(folder / "contract.json")
        live_contract = task["content"]["contract"]
        assert {k:v for k,v in contract.items() if k != "seed"} == {k:v for k,v in live_contract.items() if k != "seed"}
        reference = read(folder / "hbb_reference.json")
        reference["challenge"]["seed"] = task["content"]["hbb_reference"]["challenge"]["seed"]
        assert reference == task["content"]["hbb_reference"]
        score = best[hk]
        assert abs(validation["score"]["final_score"] - score["final_score"]) < 1e-10
        offset = 25 if i == 6 else 50*i
        candidate_window = [joined[(offset+j) % 300] for j in range(150)]
        clean, band, ca_clean = [], [], []
        seed_results = {}
        started = time.monotonic()
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
                assert all(r["outcome"] == old[r["experiment_id"]]["outcome"] and
                           r["indel_length"] == old[r["experiment_id"]]["indel_length"]
                           for r in simulated)
                seed_results[str(seed)] = {
                    "outcomes": dict(counts),
                    "outcomes_by_cas": {cas: dict(Counter(r["outcome"] for r in simulated if r["cas"] == cas))
                                        for cas in ("Cas9", "Cas12a")},
                    "in_joined_space": seed in joined, "in_candidate_window": seed in candidate_window,
                    "clean": seed in clean, "band": seed in band,
                    "final_score": validation["per_seed_final_score"][str(seed)],
                }
        assert set(band) <= set(clean) <= set(ca_clean)
        ev = evidence[f"h{i}"]
        recorded_band = sorted(json.loads(ev["window_record"]["band"]))
        build_events = [line for line in ev["events"] if "Build: conjunction " in line]
        assert len(build_events) == 1
        match = re.search(r"band (\d+) k=(\d+) group (\d+) clean (\d+)/(\d+)", build_events[0])
        assert match is not None
        logged = tuple(map(int, match.groups()))
        cas_mix = dict(Counter(r["cas_system"] for r in rows))
        if i not in (0, 6):
            assert len(band) == 11 and set(band) <= set(candidate_window)
            assert recorded_band == band
            assert logged == (11, 11, 100, len(ca_clean), 900)
            assert clean == ca_clean
            assert cas_mix == {"Cas12a": 100, "Cas9": 150}
        else:
            assert cas_mix == {"Cas9": 176, "Cas12a": 74}
            assert any("prepared round unusable" in e for e in ev["events"])
            assert any("hedges skipped" in e for e in ev["events"])
            assert recorded_band != band
        rec = {
            "hotkey": hk, "uid": score["miner_uid"], "archive": str(folder.relative_to(ROOT)),
            "submission_sha256": hashlib.sha256((folder / "submission.json").read_bytes()).hexdigest(),
            "rows": len(rows), "cas_mix": cas_mix,
            "construction": "ordinary_hdr_fallback" if i in (0, 6) else "conjunction",
            "validation_dir": str(validation_dir.relative_to(ROOT)),
            "score": score, "rank": ranks[hk], "sorted_api_position": positions[hk],
            "gap_to_top": ranked[0]["final_score"] - score["final_score"],
            "percent_of_top": 100 * score["final_score"] / ranked[0]["final_score"],
            "candidate_offset": offset, "candidate_window": candidate_window,
            "candidate_window_applies_to_submission": i not in (0, 6),
            "candidate_window_ranges": ranges(candidate_window),
            "clean_seeds": clean, "band_seeds": band, "cas12a_clean_seeds": ca_clean,
            "clean_seeds_in_joined_space": sorted(set(clean) & set(joined)),
            "clean_seeds_outside_joined_space": sorted(set(clean) - set(joined)),
            "clean_hits": sorted(set(clean) & set(seeds)),
            "band_hits": sorted(set(band) & set(seeds)),
            "recorded_band_matches_replay": recorded_band == band,
            "recorded_prepared_band": recorded_band,
            "logged_prepared_clean_count": logged[3],
            "per_seed": seed_results, "local_score_matches_api": True,
            "stage3_replay_matches_local_detail": True,
            "scan_seconds": round(time.monotonic()-started, 3),
        }
        result["hotkeys"][f"h{i}"] = rec
        print(f"h{i}: rank={ranks[hk]} final={score['final_score']:.8f} clean={len(clean)} band={band} clean_hits={rec['clean_hits']} band_hits={rec['band_hits']} ({rec['scan_seconds']}s)", flush=True)
        (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    result["fleet_clean_union"] = sorted(set().union(*(set(r["clean_seeds"]) for r in result["hotkeys"].values())))
    result["fleet_band_union"] = sorted(set().union(*(set(r["band_seeds"]) for r in result["hotkeys"].values())))
    result["fleet_clean_hits"] = sorted(set(seeds) & set(result["fleet_clean_union"]))
    result["fleet_band_hits"] = sorted(set(seeds) & set(result["fleet_band_union"]))
    assert result["hotkeys"]["h0"]["submission_sha256"] == result["hotkeys"]["h6"]["submission_sha256"]
    assert result["hotkeys"]["h0"]["clean_seeds"] == result["hotkeys"]["h6"]["clean_seeds"]
    (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Fleet unions: clean", len(result["fleet_clean_union"]), "band", len(result["fleet_band_union"]), flush=True)


if __name__ == "__main__":
    main()
