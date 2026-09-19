"""Audit the uploaded h0-h5 archives for task af21b154, without changing miner data.

Run from the repository root with .venv/bin/python reports/af21b154/audit.py.
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
TASK_ID = "af21b154-5ee4-4c34-98cb-ebd58f04cdc2"
HOTKEYS = [
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW",
    "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU",
    "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2",
    "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",
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
    seeds = [int(s) for s in task["content"]["contract"]["seed"].split(",")]
    raw_scores = read(OUT / "scores.json")
    scores = [r for r in raw_scores["items"] if r["task_id"] == TASK_ID]
    best = {}
    for row in scores:
        hk = row["miner_hotkey"]
        if hk not in best or row["final_score"] > best[hk]["final_score"]:
            best[hk] = row
    ranked = sorted(best.values(), key=lambda r: -r["final_score"])
    ranks = {r["miner_hotkey"]: i for i, r in enumerate(ranked, 1)}
    joined = list(range(200, 300)) + list(range(800, 1000))
    result = {
        "task_id": TASK_ID, "task_created_at": task["created_at"],
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cell_type": task["content"]["contract"]["cell_type"], "seeds": seeds,
        "source_api": f"https://niome-api.genomes.io/api/v3/miners/scores?task_id={TASK_ID}&limit=40000",
        "score_records": len(scores), "unique_miners": len(ranked),
        "validators": sorted({r["validator_hotkey"] for r in scores}),
        "top_miner": ranked[0], "top10_cutoff": ranked[9]["final_score"],
        "joined_band_space": joined, "cut_space": [100, 999],
        "definitions": {
            "clean_seeds": "All 250 submitted rows cut; HDR band is a subset.",
            "band_seeds": "All 250 submitted rows return HDR, scanned over every seed 100-999.",
            "cas12a_clean_seeds": "All 100 submitted Cas12a rows cut; checked against build log clean count.",
            "candidate_window": "150-seed slice reconstructed from original logged joined space and h0-h5 offsets 0/50/100/150/200/250, checked against recovered band.",
        },
        "hotkeys": {},
    }
    for i, hk in enumerate(HOTKEYS):
        name = "niome_hotkey" + (str(i) if i else "")
        archives = [p.parent for p in (ROOT / f"data/inst/{name}/result").glob("2026-09-17T*/last_upload.json")
                    if read(p).get("task_id") == TASK_ID]
        assert len(archives) == 1, (name, archives)
        folder = archives[0]
        rows = read(folder / "submission.json")
        valid = read(folder / "validation/valid_experiments.json")
        validation = read(folder / "validation.json")
        upload = read(folder / "last_upload.json")
        assert upload["submitted"] and len(rows) == 250
        assert rows == [e["experiment"] for e in valid]
        assert {e["features"]["cell_type_accessibility"] for e in valid} == {0.77}
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
        candidate_window = [joined[(50*i+j) % 300] for j in range(150)]
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
                detail = read(folder / f"stage3_detail.seed{seed}.json")
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
        assert len(band) == 11 and set(band) <= set(candidate_window)
        ev = evidence[f"h{i}"]
        assert sorted(json.loads(ev["window_record"]["band"])) == band
        build_events = [line for line in ev["events"] if "Build: conjunction " in line]
        assert len(build_events) == 1
        match = re.search(r"band (\d+) k=(\d+) group (\d+) clean (\d+)/(\d+)", build_events[0])
        assert tuple(map(int, match.groups())) == (11, 11, 100, len(clean), 900)
        assert clean == ca_clean
        assert Counter(r["cas_system"] for r in rows) == {"Cas12a": 100, "Cas9": 150}
        rec = {
            "hotkey": hk, "uid": score["miner_uid"], "archive": str(folder.relative_to(ROOT)),
            "submission_sha256": hashlib.sha256((folder / "submission.json").read_bytes()).hexdigest(),
            "rows": len(rows), "cas_mix": dict(Counter(r["cas_system"] for r in rows)),
            "score": score, "rank": ranks[hk],
            "gap_to_top": ranked[0]["final_score"] - score["final_score"],
            "percent_of_top": 100 * score["final_score"] / ranked[0]["final_score"],
            "candidate_offset": 50*i, "candidate_window": candidate_window,
            "candidate_window_ranges": ranges(candidate_window),
            "clean_seeds": clean, "band_seeds": band, "cas12a_clean_seeds": ca_clean,
            "clean_hits": sorted(set(clean) & set(seeds)),
            "band_hits": sorted(set(band) & set(seeds)),
            "recorded_band_matches_replay": True,
            "per_seed": seed_results, "archived_score_matches_api": True,
            "stage3_replay_matches_archived_detail": True,
            "scan_seconds": round(time.monotonic()-started, 3),
        }
        result["hotkeys"][f"h{i}"] = rec
        print(f"h{i}: rank={ranks[hk]} final={score['final_score']:.8f} clean={len(clean)} band={band} clean_hits={rec['clean_hits']} band_hits={rec['band_hits']} ({rec['scan_seconds']}s)", flush=True)
        (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    result["fleet_clean_union"] = sorted(set().union(*(set(r["clean_seeds"]) for r in result["hotkeys"].values())))
    result["fleet_band_union"] = sorted(set().union(*(set(r["band_seeds"]) for r in result["hotkeys"].values())))
    result["fleet_clean_hits"] = sorted(set(seeds) & set(result["fleet_clean_union"]))
    result["fleet_band_hits"] = sorted(set(seeds) & set(result["fleet_band_union"]))
    (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Fleet unions: clean", len(result["fleet_clean_union"]), "band", len(result["fleet_band_union"]), flush=True)


if __name__ == "__main__":
    main()
