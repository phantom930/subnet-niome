"""Reading the task snapshot that scripts/bench_task.py writes.

Read-only. The harness owns the file: it fetches, merges and writes, and this
module only reads what is there and diffs two reads to describe what a refresh
changed.
"""

from __future__ import annotations

import json
from typing import Any

import config


class SnapshotMissing(RuntimeError):
    """No snapshot yet. A refresh creates one."""


def parse_seeds(seed: Any) -> list[int]:
    """The round's stamped seeds, or [] when it never closed.

    The comma in this field separates seeds; it is not digit grouping. A round
    has carried three seeds since 2026-08-27 and one before that, and the
    backend joins them into a single string ('263,486,269'), so stripping the
    commas and parsing one number turns three seeds into the nine-digit
    263486269. That reading is why a seed of 1000 used to surface as
    '3,283,711,000'.

    Zero is the 'not stamped yet' placeholder rather than a seed of zero, so it
    is dropped — which also makes an all-zero field parse to [], the same
    answer an unstamped round gives.
    """
    if seed is None or seed == "":
        return []
    seeds = []
    for part in str(seed).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            value = int(float(part))
        except (TypeError, ValueError):
            return []
        if value:
            seeds.append(value)
    return seeds


def read() -> dict[str, Any]:
    """The snapshot as the harness wrote it."""
    if not config.SNAPSHOT_PATH.exists():
        raise SnapshotMissing(
            f"no snapshot at {config.SNAPSHOT_PATH}. Refresh to fetch one, or run "
            "scripts/bench_task.py --fetch"
        )
    with open(config.SNAPSHOT_PATH) as snapshot_file:
        return json.load(snapshot_file)


def read_cell_types() -> dict[str, Any]:
    """The accessibility table, written by the same --fetch run."""
    if not config.CELL_TYPES_PATH.exists():
        return {}
    try:
        with open(config.CELL_TYPES_PATH) as cell_types_file:
            return json.load(cell_types_file)
    except (OSError, ValueError):
        return {}


def read_seed_occurrence() -> dict[str, Any]:
    """The drawn-seed ledger: how often each seed has been stamped.

    Served rather than recomputed from the snapshot so the Tasks page colours
    by the same counts ``design.draw_seeds`` bets against. The two are written
    by one --fetch run, so a second count here could only drift from it.

    ``available`` is false when the ledger is missing or unreadable, and the
    page then shows seeds without occurrence colouring rather than colouring
    them all as though they were unique. ``since`` is the cutoff the ledger was
    built with: rounds older than it were never counted, so their seeds have no
    occurrence rather than an occurrence of zero — a distinction the page has
    to make or a whole era of tasks reads as never-stamped.
    """
    if not config.DRAWN_SEEDS_PATH.exists():
        return {"available": False, "counts": {}, "reason": f"no ledger at {config.DRAWN_SEEDS_PATH}"}
    try:
        with open(config.DRAWN_SEEDS_PATH) as ledger_file:
            document = json.load(ledger_file)
    except (OSError, ValueError) as error:
        return {"available": False, "counts": {}, "reason": f"cannot read the ledger: {error}"}

    counts = document.get("counts") or {}
    if not isinstance(counts, dict):
        return {"available": False, "counts": {}, "reason": "the ledger's counts are not an object"}

    return {
        "available": True,
        "counts": {str(seed): int(times) for seed, times in counts.items()},
        "since": document.get("since") or None,
        "tasks_recorded": document.get("tasks_recorded"),
        "latest_task_at": document.get("latest_task_at"),
        "undrawn": document.get("undrawn"),
        "updated_at": document.get("updated_at"),
    }


def seeds_by_id() -> dict[str, tuple[int, ...]]:
    """Each task's stamped seeds, for diffing across a refresh.

    Parsed rather than raw so a round whose three seeds arrive in a different
    string form is not reported as restamped when nothing about it changed.
    """
    try:
        snapshot = read()
    except (SnapshotMissing, ValueError):
        return {}
    return {
        task["id"]: tuple(parse_seeds(task.get("content", {}).get("contract", {}).get("seed")))
        for task in snapshot.get("tasks", [])
        if "id" in task
    }


def summary() -> dict[str, Any]:
    """Counts without shipping the whole task array."""
    # The cell-type table is written by the same --fetch run but is a separate
    # file, so it is reported whether or not the task snapshot is there.
    cell_types = len(read_cell_types())
    try:
        snapshot = read()
    except (SnapshotMissing, ValueError):
        return {"count": 0, "unstamped": 0, "fetched_at": None, "cell_types": cell_types}
    tasks = snapshot.get("tasks", [])
    return {
        "count": snapshot.get("count", len(tasks)),
        "unstamped": snapshot.get(
            "unstamped",
            sum(
                1
                for task in tasks
                if not parse_seeds(task.get("content", {}).get("contract", {}).get("seed"))
            ),
        ),
        "fetched_at": snapshot.get("fetched_at"),
        "cell_types": cell_types,
    }


def describe_change(before: dict[str, tuple[int, ...]]) -> dict[str, Any]:
    """What changed between a snapshot read earlier and the one on disk now.

    Computed from the files rather than parsed out of the harness's stdout, so
    a change to its output format cannot silently break these numbers.
    """
    after = seeds_by_id()
    added = [task_id for task_id in after if task_id not in before]
    restamped = [
        task_id
        for task_id, seed in after.items()
        if task_id in before and before[task_id] != seed
    ]
    removed = [task_id for task_id in before if task_id not in after]

    return {
        "added": len(added),
        "restamped": len(restamped),
        "removed": len(removed),
        **summary(),
    }
