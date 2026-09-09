"""What does the next task of each cell type draw?

Loads the per-cell-type checkpoint, feeds it that cell type's most recent
history, and writes the nine-class distribution plus the three classes it
would bet on. The calibration line in the output is not decoration: if the
walk-forward report says the model does not beat chance, these three classes
are three names with a ~1-in-3 hit rate each, and the file says so.
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch

from .data import (CELL_TYPES, N_CLASSES, SEEDS_PER_TASK, load_rounds,
                   predict_input, refresh_seeds)
from .model import SeedFormer

ROOT = Path(__file__).resolve().parents[1]
CKPT_DIR = ROOT / "seed_model" / "checkpoints"


def class_name(i):
    return f"{(i + 1) * 100}-{(i + 1) * 100 + 99}"


def ckpt_path(cell_type):
    return CKPT_DIR / f"{cell_type.replace('+', 'plus').replace('/', '_')}.pt"


@torch.no_grad()
def predict_cell(cell_type, rows, verdict=None):
    path = ckpt_path(cell_type)
    if not path.exists():
        return {"cell_type": cell_type, "error": f"no checkpoint at {path.name}"}
    ck = torch.load(path, weights_only=True)
    model = SeedFormer(**ck["config"])
    model.load_state_dict(ck["state"])
    model.eval()
    x, mask = predict_input(rows, ck["context"])
    logits = model(torch.from_numpy(x), torch.from_numpy(mask),
                   torch.tensor([ck["cell_index"]]))
    p = torch.softmax(logits, dim=-1).numpy()[0]
    order = np.argsort(-p, kind="stable")
    return {
        "cell_type": cell_type,
        "rounds_seen": len(rows),
        "last_round": {"round": rows[-1]["round"],
                       "created_at": rows[-1]["created_at"],
                       "seeds": rows[-1]["seeds"],
                       "classes": rows[-1]["seed_classes"]},
        "cumulative_counts": {class_name(i): rows[-1]["cumulative_counts"][i]
                              for i in range(N_CLASSES)},
        "distribution": {class_name(i): round(float(p[i]), 5)
                         for i in range(N_CLASSES)},
        "predicted_classes": [class_name(i) for i in order[:SEEDS_PER_TASK]],
        "predicted_mass": round(float(p[order[:SEEDS_PER_TASK]].sum()), 5),
        "uniform_mass_for_three": round(SEEDS_PER_TASK / N_CLASSES, 5),
        "max_deviation_from_uniform": round(
            float(np.abs(p - 1.0 / N_CLASSES).max()), 5),
        "trained_at": ck.get("trained_at"),
        "n_rounds_trained": ck.get("n_rounds"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--start", default="2026-08-27T00:00:00")
    ap.add_argument("--end", default="now")
    ap.add_argument("--out", default=str(ROOT / "seed_model" / "next_prediction.json"))
    ap.add_argument("--report", default=str(ROOT / "seed_model" / "report.json"))
    args = ap.parse_args()

    rounds, meta = (refresh_seeds(args.start, args.end) if args.refresh
                    else load_rounds())
    verdict = None
    if Path(args.report).exists():
        verdict = json.loads(Path(args.report).read_text()).get("verdict")

    preds = [predict_cell(ct, rounds[ct], verdict)
             for ct in CELL_TYPES if ct in rounds]
    beats = (verdict or {}).get("vs_best_reference", {})
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "window": meta["window"],
        "calibration": {
            "walk_forward_model_hits": (verdict or {}).get("model_mean_hits"),
            "walk_forward_best_reference": (verdict or {}).get("best_reference"),
            "walk_forward_reference_hits": (verdict or {}).get("reference_mean_hits"),
            "chance_hits": 1.0,
            "p_value_vs_best_reference": beats.get("p_value"),
            "reading": (
                "the model beat the best reference strategy on held-out rounds"
                if beats.get("p_value") is not None and beats["p_value"] < 0.05
                else "no evidence the model beats chance; treat these three "
                     "classes as an arbitrary pick with a ~1-in-3 hit rate each"),
        },
        "predictions": preds,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n")

    print(f"next-round prediction, window ending {meta['window']['end'][:19]}")
    for pr in preds:
        if "error" in pr:
            print(f"  {pr['cell_type']:<11} {pr['error']}")
            continue
        print(f"  {pr['cell_type']:<11} round {pr['rounds_seen'] + 1:>3}  "
              f"-> {', '.join(pr['predicted_classes'])}  "
              f"(mass {pr['predicted_mass']:.3f} vs uniform "
              f"{pr['uniform_mass_for_three']:.3f}, max dev "
              f"{pr['max_deviation_from_uniform']:.4f})")
    print(f"\n  {payload['calibration']['reading']}")
    print(f"  -> {args.out}")


if __name__ == "__main__":
    main()
