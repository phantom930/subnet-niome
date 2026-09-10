#!/usr/bin/env python3
"""Run the miner against a recorded task and score it with the validator's own stages.

    scripts/bench_task.py                      # newest task, under the seeds it closed under
    scripts/bench_task.py --task 3 --random-seeds --seeds 5
    scripts/bench_task.py --list

The point of the harness is the asymmetry that the live subnet has and a naive test does not:

* **the miner is run blind.** The contract handed to ``design.build`` has its ``seed`` forced to 0,
  which is what the backend actually broadcasts — the round seed is stamped after the task goes out.
  ``testing/task.json`` is a *closed-round* snapshot and its contracts do carry the real seed, so
  reading it straight into the miner would quietly hand the design an oracle it never has.
* **the validator is run with the seeds the round closed under, which the miner never saw.** A
  recorded task normally carries them, so scoring under them is what a validator actually paid for
  those rows. A task whose ``seed`` is still 0 — the backend's placeholder for "not stamped yet" —
  has nothing to score against, so seeds are drawn at random from ``design.SEED_SUPPORT`` instead;
  ``--random-seeds`` forces that draw for a stamped task too. Either way they go into the ``seed``
  field comma-joined, the way the backend joins a multi-seed round, so ``benchmark_submission``
  averages several draws exactly as it does in production.

Scoring goes through ``benchmark_submission`` itself rather than a reimplementation, so the number
printed is the number a validator would compute. That entry point communicates through the fixed
filenames under ``data/``, which a co-located validator is also using — so this script redirects
every stage's path constant into a temporary directory and leaves ``data/`` untouched. Only
``chr11.fa`` is read from where it already lives.

The rows themselves are kept, at ``testing/submission.json``, since the work directory is thrown
away with everything staged in it. That copy is written straight from the build and is deliberately
not the one handed to the stages: ``truncate_submission`` rewrites the scoring copy in place, so
only a file outside the work directory is guaranteed to be the array the miner would upload.
"""

import argparse
import json
import logging
import os
import random
import shutil
import sys
import tempfile
import time

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import niome_subnet.utils.settings as settings  # noqa: E402

from niome_subnet.genomics import design, validation  # noqa: E402
from niome_subnet.genomics.validation import stage12, stage3, stage4, stage5  # noqa: E402

TASK_SNAPSHOT = PROJECT_ROOT / "testing" / "task.json"
CELL_TYPE_SNAPSHOT = PROJECT_ROOT / "testing" / "cell_types.json"
SUBMISSION_OUTPUT = PROJECT_ROOT / "testing" / "submission.json"
DEFAULT_CELL_TYPES = {"HEK293": {"accessibility": 0.35}}
# How many seeds are drawn when there are none to score under. Only reached for an unstamped task,
# or under --random-seeds.
DEFAULT_SEED_COUNT = 3

# The closed-round history, and the only task endpoint that needs no hotkey signature: /current is
# what the validator signs for, and it answers "Missing required headers" to an unsigned GET. Both
# of these are plain public reads — nothing identifying is sent.
TASK_HISTORY_URL = f"{settings.BASE_URL}/api/v3/tasks"

# Every stage constant that names a file under data/, and the basename it keeps inside the work
# directory. The stages do `from ...settings import CONTRACT_PATH`, so the name is bound in each
# module's own globals and patching settings would have no effect — the modules are patched instead.
REDIRECTED_PATHS = {
    "CONTRACT_PATH": "contract.json",
    "HBB_REFERENCE_PATH": "hbb_reference.json",
    "MINER_SUBMISSION_PATH": "submission.json",
    "VALID_EXPERIMENTS_PATH": "valid_experiments.json",
    "INVALID_EXPERIMENTS_PATH": "invalid_experiments.json",
    "STAGE3_DATASET": "stage3_dataset.json",
    "STAGE3_SUMMARY_PATH": "stage3_summary.json",
    "FINAL_REWARD_PATH": "final_reward.json",
    "DISTRIBUTION_FIDELITY_PATH": "distribution_fidelity_summary.json",
    "KMER_CACHE_DIR": "kmer_cache",
}
STAGE_MODULES = (validation, stage12, stage3, stage4, stage5)


# ---------------------------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------------------------

def redirect_stage_paths(work_dir: Path) -> None:
    """Point the validation stages at ``work_dir`` so ``data/`` is neither read nor written."""
    for module in STAGE_MODULES:
        for name, filename in REDIRECTED_PATHS.items():
            if hasattr(module, name):
                setattr(module, name, str(work_dir / filename))
        # The genome is 130 MB, identical for both roles and read-only, so it stays put.
        if hasattr(module, "CHR11_PATH"):
            setattr(module, "CHR11_PATH", settings.MINER_GENOME_PATH)


def memoise_chr11() -> None:
    """``run_stage12`` re-parses the whole FASTA on every call; once per process is enough here."""
    original_load = stage12.load_chr11
    cache: dict[str, str] = {}

    def cached_load_chr11(path):
        if path not in cache:
            cache[path] = original_load(path)
        return cache[path]

    stage12.load_chr11 = cached_load_chr11


def load_cell_types() -> dict:
    """The accessibility table, from whichever local copy exists.

    Accessibility is the largest term in stage 3's energy, so both roles must be handed the same
    table or the comparison is between two different simulations. ``--fetch`` writes the backend's
    own table next to the task snapshot, which is preferred over the miner's single-cell-type copy
    because a fetched task can name any cell type the backend issues.
    """
    for path in (CELL_TYPE_SNAPSHOT, settings.MINER_CELL_TYPES_PATH, "data/cell_types.json"):
        if os.path.exists(path):
            with open(path) as cell_types_file:
                return json.load(cell_types_file)
    print(f"! no cell-types table on disk; assuming {DEFAULT_CELL_TYPES}. "
          f"Run --fetch to download the backend's own table.")
    return DEFAULT_CELL_TYPES


def check_cell_type(contract: dict, cell_types: dict) -> None:
    """Refuse to score a cell type the table does not cover.

    ``stage2`` defaults a missing accessibility to 1.0 rather than raising, and accessibility is the
    largest term in stage 3's energy — so an absent entry does not fail, it silently scores the task
    as though its chromatin were fully open. That is the one wrong answer worth stopping for.
    """
    cell_type = contract.get("cell_type")
    if cell_type in cell_types:
        return
    raise SystemExit(
        f"the cell-types table has no entry for {cell_type!r} (it has: "
        f"{', '.join(sorted(cell_types)) or 'nothing'}).\n"
        "stage 2 would default its accessibility to 1.0 and the score would describe a different "
        "simulation than a validator's. Run --fetch to download the backend's table."
    )


def load_tasks() -> list[dict]:
    if not TASK_SNAPSHOT.exists():
        raise SystemExit(
            f"{TASK_SNAPSHOT} is missing — it is the recorded task history to test against. "
            "Run --fetch to download it."
        )
    with open(TASK_SNAPSHOT) as snapshot_file:
        return json.load(snapshot_file)["tasks"]


# ---------------------------------------------------------------------------------------------
# Fetching the task history from the backend
# ---------------------------------------------------------------------------------------------

def get_json(url: str, params: dict | None = None) -> dict:
    """A retried, unsigned GET. Raises rather than returning a partial answer."""
    import requests

    last_error = None
    for attempt in range(1, settings.MAX_TASK_RETRIES + 1):
        try:
            response = requests.get(url, params=params,
                                    timeout=settings.TASK_REQUEST_TIMEOUT)
            response.raise_for_status()
            return response.json()
        except Exception as error:
            last_error = error
            print(f"  ! {url} attempt {attempt}/{settings.MAX_TASK_RETRIES}: {error}")
            if attempt < settings.MAX_TASK_RETRIES:
                time.sleep(settings.BASE_DELAY_SECONDS * attempt)
    raise SystemExit(f"could not fetch {url}: {last_error}")


def fetch_task_history() -> list[dict]:
    """Every closed round the backend still lists, newest first.

    The endpoint paginates as ``{items, page, pages, per_page, total}``. Today it answers the whole
    history in one page, but the loop follows ``pages`` so a growing history does not silently
    truncate to the first page.
    """
    tasks: list[dict] = []
    page = 1
    while True:
        payload = get_json(TASK_HISTORY_URL, {"page": page})
        items = payload.get("items")
        if items is None:
            raise SystemExit(
                f"{TASK_HISTORY_URL} did not return an 'items' list "
                f"(got keys: {sorted(payload)}) — the endpoint's shape has changed"
            )
        tasks.extend(items)
        total_pages = int(payload.get("pages") or 1)
        print(f"  page {page}/{total_pages}: {len(items)} tasks "
              f"({len(tasks)} of {payload.get('total', '?')})")
        if page >= total_pages:
            return tasks
        page += 1


def refresh_snapshot(replace: bool = False) -> None:
    """Merge the backend's task history and cell-type table into ``testing/``.

    Merging rather than overwriting is deliberate: the backend's history is a window, so a task that
    has aged out of it stays available to test against instead of disappearing from the snapshot.
    ``--replace`` takes the backend's answer as the whole truth.
    """
    print(f"fetching {TASK_HISTORY_URL}")
    fetched = fetch_task_history()

    existing = []
    if TASK_SNAPSHOT.exists() and not replace:
        with open(TASK_SNAPSHOT) as snapshot_file:
            existing = json.load(snapshot_file).get("tasks", [])

    # Fetched entries win on a shared id: a round that was unstamped when it was last recorded
    # carries its real seed now.
    by_id = {task["id"]: task for task in existing}
    added = sum(1 for task in fetched if task["id"] not in by_id)
    restamped = sum(
        1 for task in fetched
        if task["id"] in by_id
        and by_id[task["id"]]["content"]["contract"].get("seed")
        != task["content"]["contract"].get("seed")
    )
    by_id.update({task["id"]: task for task in fetched})

    tasks = sorted(by_id.values(), key=lambda task: str(task.get("created_at", "")), reverse=True)
    unstamped = [task for task in tasks if not task["content"]["contract"].get("seed")]

    print(f"fetching {settings.CELL_TYPES_URL}")
    cell_types = get_json(settings.CELL_TYPES_URL)

    TASK_SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
    write_json(TASK_SNAPSHOT, {
        "source": TASK_HISTORY_URL,
        "leaderboard": "https://niome-leaderboard.genomes.io/tasks",
        "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "order": "newest first",
        "note": (
            "Closed-round snapshot: each task carries the contract a validator broadcast plus the "
            "seed stamped after that round closed, which a live task does not. Merged across "
            "fetches, so tasks that have aged out of the backend's window are kept."
        ),
        "count": len(tasks),
        "unstamped": len(unstamped),
        "tasks": tasks,
    })
    write_json(CELL_TYPE_SNAPSHOT, cell_types)

    print()
    print(f"  {TASK_SNAPSHOT.relative_to(PROJECT_ROOT)}: {len(tasks)} tasks "
          f"(+{added} new, {restamped} newly stamped, {len(unstamped)} still unstamped)")
    table = ", ".join(
        f"{name} {body.get('accessibility')}" for name, body in sorted(cell_types.items())
    )
    print(f"  {CELL_TYPE_SNAPSHOT.relative_to(PROJECT_ROOT)}: {table}")

    used_cell_types = {task["content"]["contract"]["cell_type"] for task in tasks}
    missing = sorted(used_cell_types - set(cell_types))
    if missing:
        print(f"  ! {len(missing)} cell type(s) in the history are absent from the table: {missing}")


def select_task(tasks: list[dict], selector: str) -> dict:
    """A task by list position (0 is newest) or by id."""
    if selector.lstrip("-").isdigit():
        return tasks[int(selector)]
    for task in tasks:
        if task["id"] == selector or task["id"].startswith(selector):
            return task
    raise SystemExit(f"no task matches {selector!r}; run --list to see them")


# ---------------------------------------------------------------------------------------------
# Which seeds the validator half gets
# ---------------------------------------------------------------------------------------------

def recorded_seeds(contract: dict) -> list[int]:
    """The round seeds stamped on a task, in the order the backend recorded them, or ``[]``.

    Parsed exactly as ``Context.seeds()`` parses it, zeros dropped: 0 is the backend's placeholder
    for "not stamped yet", not a round seed. ``"0"`` therefore has to be filtered *before* the
    emptiness test — a list holding only zeros is still a non-empty list, and treating it as seeds
    would score every row under a seed no validator ever held.
    """
    try:
        parsed = [
            int(seed_text) for seed_text in str(contract.get("seed", "")).split(",")
            if seed_text.strip()
        ]
    except (TypeError, ValueError):
        return []
    return [seed for seed in parsed if seed]


def draw_seeds(count: int, rng_seed: int | None) -> list[int]:
    """``count`` seeds the miner never saw, drawn from every seed the backend has been observed to
    stamp. ``SEED_SUPPORT`` is ``range(100, 1000)``."""
    rng = random.Random(rng_seed)
    return sorted(rng.sample(design.SEED_SUPPORT, min(count, len(design.SEED_SUPPORT))))


def choose_seeds(task: dict, args) -> tuple[list[int], str]:
    """The seeds to score under, and one line saying where they came from.

    The task's own seeds are the default because they are the ones the round closed under: the
    resulting score is the one a validator computed, not an estimate of it. Random draws are the
    fallback for a task that has no seeds — and, under ``--random-seeds``, a way to ask how the same
    rows hold up against seeds that were never played.
    """
    contract = task["content"]["contract"]
    recorded = recorded_seeds(contract)
    seed_count = DEFAULT_SEED_COUNT if args.seeds is None else args.seeds

    if args.random_seeds:
        return draw_seeds(seed_count, args.rng), "drawn at random, not the task's own"

    if recorded:
        if args.seeds is not None or args.rng is not None:
            print("! --seeds/--rng describe a random draw, and this task is stamped, so they are "
                  "ignored. Add --random-seeds to draw seeds instead of using the task's own.")
        return recorded, "the task's own, recorded after the round closed"

    if args.task_seed:
        raise SystemExit(
            f"task {task['id']} is unstamped (seed {contract.get('seed')!r}) — it has no recorded "
            "seed to score under. Drop --task-seed to fall back to random seeds."
        )
    return draw_seeds(seed_count, args.rng), "drawn at random — the task is unstamped"


# ---------------------------------------------------------------------------------------------
# The two halves
# ---------------------------------------------------------------------------------------------

def run_miner(contract: dict, reference: dict, cell_types: dict) -> tuple[list[dict], dict]:
    """Build a submission the way ``Miner._build`` does, from a contract with no round seed."""
    blind_contract = dict(contract, seed=0)
    context = design.build_context(blind_contract, reference, cell_types)
    assert not context.seeds(), "the miner must not be able to see a round seed"

    started = time.time()
    rows, _entries, diagnostics = design.build(context)
    diagnostics["build_seconds"] = time.time() - started
    return rows, diagnostics


def score_with_validator(rows: list[dict], contract: dict, reference: dict, cell_types: dict,
                         seeds: list[int], work_dir: Path, uid: int) -> dict:
    """What a validator holding ``seeds`` would pay for ``rows``.

    The seeds go into the contract's ``seed`` field comma-joined, which is the format
    ``benchmark_submission`` parses and the format the backend uses for a multi-seed round.
    """
    stamped_contract = dict(contract, seed=",".join(str(seed) for seed in seeds))
    write_json(work_dir / "contract.json", stamped_contract)
    write_json(work_dir / "hbb_reference.json", reference)
    write_json(work_dir / "submission.json", rows)

    score = validation.benchmark_submission(cell_types, uid)
    return dict(score.breakdown, final_score=score.final_score)


def write_json(path: Path, document) -> None:
    with open(path, "w") as output_file:
        json.dump(document, output_file)


# ---------------------------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------------------------

def print_report(task: dict, rows: list[dict], diagnostics: dict, seeds: list[int],
                 seed_source: str, breakdown: dict, per_seed: list[tuple[int, dict]],
                 cell_types: dict) -> None:
    contract = task["content"]["contract"]
    rules = contract["rules"]
    accessibility = cell_types.get(contract["cell_type"], {}).get("accessibility")

    print()
    print(f"task {task['id']}  ({task.get('created_at', '?')})")
    # Accessibility scales stage 3's energy, so it decides the cut probability and with it how much
    # of consistency_factor is even reachable. It is the first thing to read off a task.
    print(f"  cell_type       {contract['cell_type']}  (accessibility {accessibility})")
    print(f"  mutations       {', '.join(contract['active_mutations'])}")
    print(f"  weights         {contract.get('mutation_weights', {})}")
    print(f"  rules           cas={rules.get('cas_systems')} max_experiments={rules.get('max_experiments')} "
          f"max_mismatches={rules.get('max_mismatches')} base_padding={rules.get('base_padding')} "
          f"proximity_gate={rules.get('proximity_gate')}")
    print(f"  recorded seed   {contract.get('seed')}  (withheld from the miner)")

    print()
    print(f"miner  — blind, seed 0                                      {diagnostics['build_seconds']:.1f}s")
    print(f"  {'rows':30s} {diagnostics['rows']}/{diagnostics['rows_wanted']}")
    # Named to match the validator's field below: the design computes it from stage 2 directly, so
    # the two lines must agree exactly. A gap means rows were dropped between build and scoring.
    print(f"  {'total_weighted_score':30s} {diagnostics['total_weighted_score']:.3f}")
    print(f"  {'offtarget_factors':30s} {diagnostics['offtarget_factors']}")
    print(f"  {'distinct_feature_vectors':30s} {diagnostics['distinct_feature_vectors']}")
    print(f"  {'allocation':30s} {diagnostics['allocation']}")
    if diagnostics["empty_cells"]:
        print(f"  ! EMPTY CELLS   {diagnostics['empty_cells']}  (stage 5 geometric mean → ~1e-9)")

    print()
    # "seeds [...]" is what the dashboard's report parser looks for, so the source goes after the
    # bracket rather than into the label.
    print(f"validator — seeds {seeds}  ({seed_source})")
    # Every field benchmark_submission puts in MinerScore.breakdown, under its own name. The two
    # *_score fields are the raw stage outputs; the *_factor fields are those clamped to [0, 1],
    # and only the factors enter the product.
    for key, fmt in (
        ("n_valid_experiments",          "d"),
        ("total_weighted_score",         ".3f"),
        ("consistency_score",            ".4f"),
        ("consistency_factor",           ".4f"),
        ("distribution_fidelity_score",  ".4f"),
        ("distribution_fidelity_factor", ".4f"),
    ):
        print(f"  {key:30s} {format(breakdown[key], fmt)}")
    print(f"  {'final_score':30s} {breakdown['final_score']:.4f}")
    print(f"  {'':30s} = {breakdown['total_weighted_score']:.3f}"
          f" x {breakdown['consistency_factor']:.4f}"
          f" x {breakdown['distribution_fidelity_factor']:.4f}")

    if per_seed:
        print()
        print(f"  {'seed':>6s}  {'weighted':>10s}  {'consistency':>11s}  {'fidelity':>9s}  {'final':>9s}")
        for seed, one in per_seed:
            print(f"  {seed:6d}  {one['total_weighted_score']:10.3f}  "
                  f"{one['consistency_factor']:11.4f}  {one['distribution_fidelity_factor']:9.4f}  "
                  f"{one['final_score']:9.4f}")
        finals = [one["final_score"] for _seed, one in per_seed]
        spread = max(finals) - min(finals)
        print(f"  {'spread':>6s}  {'':>10s}  {'':>11s}  {'':>9s}  {spread:9.4f}")


# ---------------------------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the miner blind against a recorded task, then score it with the validator's stages.",
    )
    parser.add_argument("--task", default="0",
                        help="task list position (0 = newest) or task id. Default: 0")
    parser.add_argument("--seeds", type=int, default=None,
                        help=f"how many seeds to draw when they are drawn at random — for an "
                             f"unstamped task, or under --random-seeds. "
                             f"Default: {DEFAULT_SEED_COUNT}")
    parser.add_argument("--rng", type=int, default=None,
                        help="seed the RNG that *picks* random round seeds, to repeat a run exactly")
    # Both of these only say what to do about the task's recorded seed, so asking for both at once
    # is a contradiction rather than a precedence question.
    seed_source = parser.add_mutually_exclusive_group()
    seed_source.add_argument("--random-seeds", action="store_true",
                             help="score under random seeds even though the task carries its own, "
                                  "to see how the same rows hold up against seeds never played")
    seed_source.add_argument("--task-seed", action="store_true",
                             help="require the task's recorded seed: fail on an unstamped task "
                                  "instead of falling back to random seeds")
    parser.add_argument("--per-seed", action="store_true",
                        help="also score each seed on its own, to show the spread")
    parser.add_argument("--uid", type=int, default=0, help="uid to report. Default: 0")
    parser.add_argument("--work-dir", default=None,
                        help="where the staged stage files go. Default: a temporary directory")
    parser.add_argument("--keep", action="store_true", help="keep the work directory")
    parser.add_argument("--list", action="store_true", help="list the recorded tasks and exit")
    parser.add_argument("--fetch", action="store_true",
                        help="download the closed-round history and the cell-type table into "
                             "testing/, merging with what is already there, then exit")
    parser.add_argument("--replace", action="store_true",
                        help="with --fetch, discard the local snapshot instead of merging into it")
    parser.add_argument("--verbose", action="store_true", help="show the stages' own logging")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if args.fetch:
        refresh_snapshot(replace=args.replace)
        return 0
    if args.replace:
        raise SystemExit("--replace only means something with --fetch")

    tasks = load_tasks()
    if args.list:
        for position, task in enumerate(tasks):
            contract = task["content"]["contract"]
            print(f"{position:4d}  {task['id']}  {str(task.get('created_at'))[:19]}  "
                  f"seed={str(contract.get('seed')):12s} "
                  f"cell={contract.get('cell_type', '?'):11s} "
                  f"mismatches={contract['rules'].get('max_mismatches')} "
                  f"rows={contract['rules'].get('max_experiments')}")
        return 0

    task = select_task(tasks, args.task)
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell_types = load_cell_types()
    check_cell_type(contract, cell_types)

    seeds, seed_source = choose_seeds(task, args)

    work_dir = Path(args.work_dir) if args.work_dir else Path(tempfile.mkdtemp(prefix="niome-bench-"))
    work_dir.mkdir(parents=True, exist_ok=True)
    redirect_stage_paths(work_dir)
    memoise_chr11()

    try:
        rows, diagnostics = run_miner(contract, reference, cell_types)
        if not rows:
            raise SystemExit(f"the miner built no rows: {diagnostics.get('error', 'unknown reason')}")

        # Kept before scoring and outside the work directory: the stages rewrite their own copy in
        # place, and row order is part of what is being scored, so this is saved exactly as built.
        write_json(SUBMISSION_OUTPUT, rows)

        breakdown = score_with_validator(
            rows, contract, reference, cell_types, seeds, work_dir, args.uid
        )
        per_seed = []
        if args.per_seed and len(seeds) > 1:
            per_seed = [
                (seed, score_with_validator(
                    rows, contract, reference, cell_types, [seed], work_dir, args.uid
                ))
                for seed in seeds
            ]

        print_report(task, rows, diagnostics, seeds, seed_source, breakdown, per_seed, cell_types)
        print(f"\nsubmission: {SUBMISSION_OUTPUT.relative_to(PROJECT_ROOT)}  "
              f"({len(rows)} rows, as built)")
        if args.keep or args.work_dir:
            print(f"stage files: {work_dir}")
        return 0
    finally:
        if not (args.keep or args.work_dir):
            shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
