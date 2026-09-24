"""Validate rows at diagnostic seed 100 while the scoring seed stamp is unavailable.
These local final scores are probes, not reproductions of the official round scores."""
import concurrent.futures
import json
import os
from pathlib import Path
import subprocess

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
manifest = json.loads((OUT / "manifest.json").read_text())["hotkeys"]


def validate(item):
    name, entry = item
    if not entry["archive"]:
        print(name, "no uploaded archive", flush=True)
        return 0
    if entry["existing_validation"]:
        print(name, "using archived validation", flush=True)
        return 0
    dest = ROOT / entry["validation_dir"]
    dest.mkdir(exist_ok=True)
    cmd = [str(ROOT / ".venv/bin/python"), str(ROOT / "calc.py"),
           "--folder", str(ROOT / entry["archive"]), "--out-dir", str(dest),
           "--task", str(OUT / "task.json"), "--cell-types", str(OUT / "cell_types.json"),
           "--no-compare", "--quiet", "--seed", "100"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (dest / "calc.log").open("w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=ROOT)
    print(name, "validation exit", result.returncode, flush=True)
    if result.returncode:
        print((dest / "calc.log").read_text()[-3000:], flush=True)
    return result.returncode


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        codes = list(pool.map(validate, manifest.items()))
    assert not any(codes), codes
