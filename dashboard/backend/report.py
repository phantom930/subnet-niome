"""Parsing scripts/bench_task.py's printed report.

The harness prints a human-readable report and has no JSON mode, so the
numbers are read back out of stdout. Best-effort by design: the raw text is
always returned alongside, and a field the parser does not recognise is simply
absent rather than guessed. The popup shows the raw report, so nothing is lost
if the format changes.

Parsing is section-aware because ``total_weighted_score`` is printed twice,
once under miner and once under validator, with different meanings.
"""

from __future__ import annotations

import re
from typing import Any

# "  consistency_factor             0.1068"
FIELD = re.compile(r"^ {2}(?P<key>[a-z_]+[a-z])\s{2,}(?P<value>\S.*?)\s*$")

# "task 12be08f9-...  (2026-09-08T19:06:32.768413)"
TASK_LINE = re.compile(r"^task\s+(?P<id>[0-9a-f-]{8,})\s+\((?P<created_at>[^)]*)\)")

# "validator — seeds [929, 5001]  (the task's own, recorded after the round closed)". The
# parenthetical says where the seeds came from, and is optional so an older report still parses.
SEEDS_LINE = re.compile(
    r"^validator.*?seeds\s*\[(?P<seeds>[^\]]*)\](?:\s*\((?P<source>[^)]*)\))?"
)

# "  cell_type       CD34+_HSPC  (accessibility 0.87)"
CELL_TYPE = re.compile(r"^(?P<cell_type>\S+)\s+\(accessibility\s+(?P<accessibility>[\d.]+)\)")

# "  recorded seed   239,996,208  (withheld from the miner)". Its own pattern
# because the label contains a space, which FIELD's key does not allow.
RECORDED_SEED = re.compile(r"^ {2}recorded seed\s+(?P<seed>\S+)")

# "   929      263.062       0.1068     0.9225    25.9223"
PER_SEED_ROW = re.compile(
    r"^\s{2,}(?P<seed>\d+)\s+(?P<weighted>[\d.]+)\s+(?P<consistency>[\d.]+)"
    r"\s+(?P<fidelity>[\d.]+)\s+(?P<final>[\d.]+)\s*$"
)

MINER_NUMBERS = {"total_weighted_score", "distinct_feature_vectors"}
VALIDATOR_NUMBERS = {
    "n_valid_experiments",
    "total_weighted_score",
    "consistency_score",
    "consistency_factor",
    "distribution_fidelity_score",
    "distribution_fidelity_factor",
    "final_score",
}


def parse(output: str) -> dict[str, Any]:
    """Pull what is recognisable out of a report. Never raises."""
    task: dict[str, Any] = {}
    miner: dict[str, Any] = {}
    validator: dict[str, Any] = {}
    per_seed: list[dict[str, float]] = []
    seeds: list[int] = []
    seed_source: str | None = None
    spread: float | int | str | None = None
    section = "task"

    for line in output.splitlines():
        if not line.strip():
            continue

        task_match = TASK_LINE.match(line)
        if task_match:
            section = "task"
            task["id"] = task_match.group("id")
            task["created_at"] = task_match.group("created_at")
            continue

        if line.startswith("miner"):
            section = "miner"
            continue

        if line.startswith("validator"):
            section = "validator"
            seeds_match = SEEDS_LINE.match(line)
            if seeds_match:
                seeds = [
                    int(part.strip())
                    for part in seeds_match.group("seeds").split(",")
                    if part.strip().isdigit()
                ]
                seed_source = seeds_match.group("source")
            continue

        if section == "validator" and per_seed is not None:
            row = PER_SEED_ROW.match(line)
            # Only after the per-seed table has started, which the header line
            # ('seed  weighted  ...') does not match because it is not numeric.
            if row:
                per_seed.append(
                    {
                        "seed": int(row.group("seed")),
                        "total_weighted_score": float(row.group("weighted")),
                        "consistency_factor": float(row.group("consistency")),
                        "distribution_fidelity_factor": float(row.group("fidelity")),
                        "final_score": float(row.group("final")),
                    }
                )
                continue

        if section == "task":
            recorded = RECORDED_SEED.match(line)
            if recorded:
                task["recorded_seed"] = recorded.group("seed")
                continue

        field = FIELD.match(line)
        if not field:
            continue
        key, raw = field.group("key"), field.group("value")

        # Printed as the last row of the per-seed table, so it describes the
        # spread across seeds rather than being one of the validator's fields.
        if section == "validator" and key == "spread":
            spread = _number(raw)
            continue

        if section == "task":
            _read_task_field(task, key, raw)
        elif section == "miner":
            miner[key] = _number(raw) if key in MINER_NUMBERS else raw
        elif section == "validator":
            validator[key] = _number(raw) if key in VALIDATOR_NUMBERS else raw

    return {
        "task": task,
        "miner": miner,
        "validator": validator,
        "seeds": seeds,
        # Whether the score is the one the round actually paid, or a draw the
        # round never played. Absent from a report the harness printed before
        # the source was named.
        "seed_source": seed_source,
        "per_seed": per_seed,
        "spread": spread,
    }


def _read_task_field(task: dict[str, Any], key: str, raw: str) -> None:
    if key == "cell_type":
        match = CELL_TYPE.match(raw)
        if match:
            task["cell_type"] = match.group("cell_type")
            task["accessibility"] = _number(match.group("accessibility"))
        else:
            task["cell_type"] = raw
    elif key == "mutations":
        task["mutations"] = [part.strip() for part in raw.split(",") if part.strip()]
    else:
        task[key] = raw


def _number(raw: str) -> float | int | str:
    """A float, an int when it is whole, or the text if it is not a number."""
    try:
        value = float(raw)
    except ValueError:
        return raw
    return int(value) if value.is_integer() and "." not in raw else value
