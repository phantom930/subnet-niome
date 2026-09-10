"""HTTP API for the dashboard frontend.

Serves the closed-round task snapshot and refreshes it by running
``scripts/bench_task.py --fetch`` as a subprocess. The harness stays the one
implementation of fetching and merging, and the snapshot in testing/ stays a
single source of truth shared with the CLI.

Run it from this directory:

    ../../.venv/bin/uvicorn main:app --reload --port 8000
"""

from __future__ import annotations

import asyncio
from typing import Any

from fastapi import Body, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

import config
import jobs
import runner
import store

app = FastAPI(
    title="NIOME dashboard API",
    version="0.2.0",
    description="Serves the task snapshot and refreshes it via scripts/bench_task.py --fetch.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=config.CORS_ORIGINS,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

class RefreshBody(BaseModel):
    replace: bool = Field(
        False,
        description="Pass --replace, discarding the local snapshot instead of merging into it",
    )


class BenchmarkBody(BaseModel):
    """Mirrors bench_task.py's flags."""

    task: str = Field(..., description="Task id, or list position with 0 the newest")
    seeds: int = Field(
        runner.DEFAULT_SEEDS, ge=1, le=20,
        description="How many seeds to draw when they are drawn at random: for an unstamped "
                    "task, or with random_seeds",
    )
    rng: int | None = Field(None, description="Seed the RNG that picks the seeds, to repeat a run")
    random_seeds: bool = Field(
        False,
        description="Score under random seeds even though the task carries its own. The default "
                    "is the seeds the round closed under",
    )
    per_seed: bool = Field(False, description="Also score each seed alone, to show the spread")
    uid: int = Field(0, ge=0, description="uid to report")


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {
        "status": "ok",
        "snapshot": store.summary(),
        "snapshot_path": str(config.SNAPSHOT_PATH),
        "harness": str(config.BENCH_SCRIPT),
        "harness_busy": runner.harness_lock.locked(),
        "active_jobs": jobs.active_count(),
    }


@app.get("/api/tasks")
def get_tasks() -> dict[str, Any]:
    """The closed-round snapshot, newest first."""
    try:
        return store.read()
    except store.SnapshotMissing as error:
        raise HTTPException(status_code=404, detail=str(error))
    except (OSError, ValueError) as error:
        raise HTTPException(status_code=500, detail=f"cannot read the snapshot: {error}")


@app.get("/api/cell-types")
def get_cell_types() -> dict[str, Any]:
    """The accessibility table, written by the same --fetch run."""
    return store.read_cell_types()


@app.post("/api/tasks/refresh")
async def refresh_tasks(body: RefreshBody = Body(default_factory=RefreshBody)) -> dict[str, Any]:
    """Run the harness's --fetch, then report what it changed."""
    if runner.harness_lock.locked():
        raise HTTPException(
            status_code=409,
            detail="the harness is busy with another run. Try again in a moment.",
        )

    async with runner.harness_lock:
        before = store.seeds_by_id()

        try:
            output = await runner.fetch_snapshot(replace=body.replace)
        except runner.FetchTimeout as error:
            raise HTTPException(status_code=504, detail=str(error))
        except runner.FetchError as error:
            raise HTTPException(status_code=502, detail=str(error))

        try:
            change = store.describe_change(before)
        except (OSError, ValueError) as error:
            raise HTTPException(
                status_code=500,
                detail=f"the harness ran but its snapshot cannot be read: {error}",
            )

        return {**change, "replaced": body.replace, "output": output}


@app.post("/api/benchmarks", status_code=202)
async def start_benchmark(body: BenchmarkBody) -> dict[str, Any]:
    """Queue a benchmark for one task and return its job immediately.

    A run takes seconds, and it queues behind any other harness run, so the
    outcome is polled from /api/benchmarks/{id} rather than awaited here.
    """
    job = jobs.submit(
        jobs.Request(
            task=body.task,
            seeds=body.seeds,
            rng=body.rng,
            random_seeds=body.random_seeds,
            per_seed=body.per_seed,
            uid=body.uid,
        )
    )
    return job.to_dict()


@app.get("/api/benchmarks")
def list_benchmarks(limit: int = 20) -> dict[str, Any]:
    return {"jobs": [job.to_dict() for job in jobs.recent(limit)]}


@app.get("/api/benchmarks/{job_id}")
def get_benchmark(job_id: str) -> dict[str, Any]:
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(
            status_code=404,
            detail=f"no job {job_id}. Jobs are held in memory, so a restart forgets them.",
        )
    return job.to_dict()
