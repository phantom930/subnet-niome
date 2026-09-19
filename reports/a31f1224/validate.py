"""Score uploaded archives into report-local directories; leave miner data intact."""
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
    dest = ROOT / entry["validation_dir"]
    dest.mkdir(exist_ok=True)
    cmd = [str(ROOT / ".venv/bin/python"), str(ROOT / "calc.py"),
           "--folder", str(ROOT / entry["archive"]), "--out-dir", str(dest),
           "--task", str(OUT / "task.json"), "--cell-types", str(OUT / "cell_types.json"),
           "--no-compare", "--quiet"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (dest / "calc.log").open("w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=ROOT)
    print(name, "validation exit", result.returncode, flush=True)
    if result.returncode:
        print((dest / "calc.log").read_text()[-3000:], flush=True)
    return result.returncode


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        codes = list(pool.map(validate, [(h, e) for h, e in manifest.items() if e["archive"] and not e["existing_validation"]]))
    assert not any(codes), codes
