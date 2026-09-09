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

# One refresh at a time. Two harness runs would both write the snapshot, and
# the second merge would be built from a file the first was mid-way through
# replacing.
_refresh_lock = asyncio.Lock()


class RefreshBody(BaseModel):
    replace: bool = Field(
        False,
        description="Pass --replace, discarding the local snapshot instead of merging into it",
    )


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {
        "status": "ok",
        "snapshot": store.summary(),
        "snapshot_path": str(config.SNAPSHOT_PATH),
        "harness": str(config.BENCH_SCRIPT),
        "refreshing": _refresh_lock.locked(),
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
    if _refresh_lock.locked():
        raise HTTPException(status_code=409, detail="a refresh is already running")

    async with _refresh_lock:
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
