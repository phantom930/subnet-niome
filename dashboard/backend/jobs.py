"""In-process job store for benchmark runs.

A benchmark takes roughly 7 seconds for one seed and longer for more, so the
request starts a job and the client polls for the result rather than holding a
connection open.

State lives in memory, so a restart forgets in-flight and finished jobs. That
is a deliberate trade for not introducing a database: the authoritative output
of a run is the report text, which the client already has once it polls.
"""

from __future__ import annotations

import asyncio
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Literal

import report
import runner

Status = Literal["queued", "running", "done", "failed"]

_MAX_JOBS = 50
_jobs: dict[str, "Job"] = {}
_order: list[str] = []

# Strong references to in-flight tasks, so they are not garbage collected
# mid-run. asyncio only holds a weak reference to a running task.
_running: set["asyncio.Task[None]"] = set()


@dataclass
class Request:
    """What the caller asked for. Mirrors bench_task.py's own flags."""

    task: str
    seeds: int = 3
    rng: int | None = None
    task_seed: bool = False
    per_seed: bool = False
    uid: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "task": self.task,
            "seeds": self.seeds,
            "rng": self.rng,
            "task_seed": self.task_seed,
            "per_seed": self.per_seed,
            "uid": self.uid,
        }


@dataclass
class Job:
    id: str
    request: Request
    status: Status = "queued"
    created_at: float = field(default_factory=time.time)
    started_at: float | None = None
    finished_at: float | None = None
    output: str | None = None
    parsed: dict[str, Any] | None = None
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "status": self.status,
            "request": self.request.to_dict(),
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "duration_seconds": (
                round(self.finished_at - self.started_at, 2)
                if self.started_at and self.finished_at
                else None
            ),
            # The report text is the source of truth; `parsed` is best-effort.
            "output": self.output,
            "result": self.parsed,
            "error": self.error,
        }


def submit(request: Request) -> Job:
    """Queue a benchmark and return it straight away.

    Must be called from the event loop, so the endpoint that calls it has to be
    an async def. A sync FastAPI endpoint runs in a threadpool where there is
    no running loop and create_task fails.
    """
    job = Job(id=uuid.uuid4().hex[:12], request=request)
    _jobs[job.id] = job
    _order.append(job.id)
    while len(_order) > _MAX_JOBS:
        _jobs.pop(_order.pop(0), None)

    job_task = asyncio.create_task(_run(job))
    _running.add(job_task)
    job_task.add_done_callback(_running.discard)
    return job


def get(job_id: str) -> Job | None:
    return _jobs.get(job_id)


def recent(limit: int = 20) -> list[Job]:
    return [_jobs[i] for i in reversed(_order[-limit:]) if i in _jobs]


def active_count() -> int:
    return sum(1 for job in _jobs.values() if job.status in ("queued", "running"))


async def _run(job: Job) -> None:
    # The lock is shared with refresh: harness runs write the same files under
    # testing/, so they queue rather than overlap. A job can therefore sit in
    # 'queued' for a while, which is why the status is polled.
    async with runner.harness_lock:
        job.status = "running"
        job.started_at = time.time()
        try:
            output = await runner.run_benchmark(
                task=job.request.task,
                seeds=job.request.seeds,
                rng=job.request.rng,
                task_seed=job.request.task_seed,
                per_seed=job.request.per_seed,
                uid=job.request.uid,
            )
            job.output = output
            job.parsed = report.parse(output)
            job.status = "done"
        except runner.FetchError as error:
            job.status = "failed"
            job.error = str(error)
        except Exception as error:  # noqa: BLE001 - surfaced to the caller as text
            job.status = "failed"
            job.error = f"{type(error).__name__}: {error}"
        finally:
            job.finished_at = time.time()
