"""Invoking scripts/bench_task.py from the backend.

A subprocess rather than an import, on purpose. bench_task pulls in
niome_subnet, the validation stages and a 130 MB chromosome, and it mutates
module-level state: ``redirect_stage_paths`` rewrites the stages' path
constants with setattr. Importing it into a long-lived web process would drag
all of that in and make two concurrent operations share that global state. A
subprocess gets its own, and exits when it is done.
"""

from __future__ import annotations

import asyncio

import config


class FetchError(RuntimeError):
    """The harness could not be run, or it failed."""


class FetchTimeout(FetchError):
    """The harness was still running when the timeout expired."""


# Every harness invocation goes through this, because they share files under
# testing/: a --fetch run rewrites task.json while a benchmark is reading it,
# and every benchmark writes the same testing/submission.json, so two at once
# would clobber each other's rows.
harness_lock = asyncio.Lock()


async def fetch_snapshot(replace: bool = False) -> str:
    """Run ``bench_task.py --fetch`` and return what it printed.

    The script fetches the closed-round history and the cell-type table,
    merges them into testing/ and prints a summary. Its exit code decides
    whether this raises.
    """
    command = [str(config.PYTHON), str(config.BENCH_SCRIPT), "--fetch"]
    if replace:
        command.append("--replace")

    return await _run(command, config.FETCH_TIMEOUT_SECONDS)


async def run_benchmark(
    task: str,
    seeds: int = 3,
    rng: int | None = None,
    task_seed: bool = False,
    per_seed: bool = False,
    uid: int = 0,
) -> str:
    """Run the miner blind against one task and score it. Returns the report.

    Slower than a fetch: the miner builds a submission and the validator puts
    it through five stages once per seed.
    """
    command = [str(config.PYTHON), str(config.BENCH_SCRIPT), "--task", str(task)]
    if task_seed:
        command.append("--task-seed")
    else:
        command += ["--seeds", str(seeds)]
        if rng is not None:
            command += ["--rng", str(rng)]
    if per_seed:
        command.append("--per-seed")
    if uid:
        command += ["--uid", str(uid)]

    return await _run(command, config.BENCHMARK_TIMEOUT_SECONDS)


async def _run(command: list[str], timeout: float) -> str:
    """Run the harness and return its stdout, raising on failure."""
    if not config.BENCH_SCRIPT.exists():
        raise FetchError(f"no harness at {config.BENCH_SCRIPT}")
    if not config.PYTHON.exists():
        raise FetchError(
            f"no interpreter at {config.PYTHON}. Set DASHBOARD_PYTHON to the one "
            "that has the project's dependencies."
        )

    try:
        process = await asyncio.create_subprocess_exec(
            *command,
            cwd=str(config.REPO_ROOT),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
    except OSError as error:
        raise FetchError(f"could not start the harness: {error}") from error

    try:
        raw_stdout, raw_stderr = await asyncio.wait_for(
            process.communicate(), timeout=timeout
        )
    except asyncio.TimeoutError as error:
        process.kill()
        await process.wait()
        raise FetchTimeout(f"the harness was still running after {timeout:.0f}s") from error

    stdout = raw_stdout.decode(errors="replace").strip()
    stderr = raw_stderr.decode(errors="replace").strip()

    if process.returncode != 0:
        detail = stderr or stdout or "no output"
        raise FetchError(f"the harness exited {process.returncode}: {_last_lines(detail)}")

    return stdout or stderr


def _last_lines(text: str, count: int = 4) -> str:
    lines = [line for line in text.splitlines() if line.strip()]
    return " / ".join(lines[-count:])
