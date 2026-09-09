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


async def fetch_snapshot(replace: bool = False) -> str:
    """Run ``bench_task.py --fetch`` and return what it printed.

    The script fetches the closed-round history and the cell-type table,
    merges them into testing/ and prints a summary. Its exit code decides
    whether this raises.
    """
    if not config.BENCH_SCRIPT.exists():
        raise FetchError(f"no harness at {config.BENCH_SCRIPT}")
    if not config.PYTHON.exists():
        raise FetchError(
            f"no interpreter at {config.PYTHON}. Set DASHBOARD_PYTHON to the one "
            "that has the project's dependencies."
        )

    command = [str(config.PYTHON), str(config.BENCH_SCRIPT), "--fetch"]
    if replace:
        command.append("--replace")

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
            process.communicate(), timeout=config.FETCH_TIMEOUT_SECONDS
        )
    except asyncio.TimeoutError as error:
        process.kill()
        await process.wait()
        raise FetchTimeout(
            f"the harness was still running after {config.FETCH_TIMEOUT_SECONDS:.0f}s"
        ) from error

    stdout = raw_stdout.decode(errors="replace").strip()
    stderr = raw_stderr.decode(errors="replace").strip()

    if process.returncode != 0:
        # bench_task raises SystemExit with a message for the failures worth
        # reading, so stderr is the useful half.
        detail = stderr or stdout or "no output"
        raise FetchError(f"the harness exited {process.returncode}: {_last_lines(detail)}")

    return stdout or stderr


def _last_lines(text: str, count: int = 4) -> str:
    lines = [line for line in text.splitlines() if line.strip()]
    return " / ".join(lines[-count:])
