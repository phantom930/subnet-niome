"""HTTP API for the dashboard frontend.

Standalone: it serves the task snapshot and refreshes it from the public
NIOME task history. It does not import niome_subnet or the scripts/ harness,
and it never writes to the subnet's data/ or testing/ directories.

Run it from this directory:

    uvicorn main:app --reload --port 8000
"""

from __future__ import annotations

import asyncio
from typing import Any

from fastapi import Body, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

import config
import store
import upstream

app = FastAPI(
    title="NIOME dashboard API",
    version="0.1.0",
    description="Serves the task snapshot and refreshes it from the public task history.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=config.CORS_ORIGINS,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

# One refresh at a time. Concurrent refreshes would each read the stored
# snapshot, merge into their own copy, and the last write would drop the
# other's additions.
_refresh_lock = asyncio.Lock()


class RefreshBody(BaseModel):
    replace: bool = Field(
        False,
        description="Discard the stored snapshot instead of merging the fetch into it",
    )


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {
        "status": "ok",
        "snapshot": store.summary(),
        "upstream": config.TASK_HISTORY_URL,
        "refreshing": _refresh_lock.locked(),
    }


@app.get("/api/tasks")
def get_tasks() -> dict[str, Any]:
    """The closed-round snapshot, newest first."""
    try:
        return store.read()
    except store.SnapshotMissing as error:
        raise HTTPException(status_code=404, detail=str(error))
    except ValueError as error:
        raise HTTPException(status_code=500, detail=f"the snapshot is not valid JSON: {error}")


@app.post("/api/tasks/refresh")
async def refresh_tasks(body: RefreshBody = Body(default_factory=RefreshBody)) -> dict[str, Any]:
    """Fetch the upstream history and merge it into the stored snapshot."""
    if _refresh_lock.locked():
        raise HTTPException(status_code=409, detail="a refresh is already running")

    async with _refresh_lock:
        try:
            fetched = await upstream.fetch_task_history()
        except upstream.UpstreamError as error:
            raise HTTPException(status_code=502, detail=str(error))

        try:
            # Merging touches the filesystem, so it goes to a worker thread
            # rather than blocking the event loop.
            return await asyncio.to_thread(store.merge, fetched, body.replace)
        except OSError as error:
            raise HTTPException(status_code=500, detail=f"could not write the snapshot: {error}")
