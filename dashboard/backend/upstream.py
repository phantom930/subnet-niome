"""Fetching the closed-round task history from the NIOME API.

Its own httpx call rather than a call into the subnet harness, so the
dashboard has no import-time dependency on niome_subnet.
"""

from __future__ import annotations

import httpx

import config


class UpstreamError(RuntimeError):
    """The upstream was unreachable, slow, or answered in an unexpected shape."""


async def fetch_task_history() -> list[dict]:
    """Every closed round the upstream still lists.

    The endpoint paginates as {items, page, pages, per_page, total}. It
    currently answers the whole history in one page, but the loop follows
    `pages` so a growing history is not silently truncated to the first page.
    """
    tasks: list[dict] = []
    page = 1

    async with httpx.AsyncClient(timeout=config.UPSTREAM_TIMEOUT_SECONDS) as client:
        while page <= config.UPSTREAM_MAX_PAGES:
            payload = await _get_page(client, page)
            items = payload.get("items")
            if items is None:
                raise UpstreamError(
                    f"{config.TASK_HISTORY_URL} returned no 'items' list "
                    f"(keys: {sorted(payload)}). The endpoint's shape has changed."
                )
            tasks.extend(items)

            total_pages = _as_int(payload.get("pages"), default=1)
            if page >= total_pages:
                return tasks
            page += 1

    raise UpstreamError(
        f"stopped after {config.UPSTREAM_MAX_PAGES} pages; raise "
        "DASHBOARD_UPSTREAM_MAX_PAGES if the history is really this long"
    )


async def _get_page(client: httpx.AsyncClient, page: int) -> dict:
    try:
        response = await client.get(config.TASK_HISTORY_URL, params={"page": page})
        response.raise_for_status()
    except httpx.TimeoutException as error:
        raise UpstreamError(
            f"{config.TASK_HISTORY_URL} timed out after "
            f"{config.UPSTREAM_TIMEOUT_SECONDS:.0f}s"
        ) from error
    except httpx.HTTPStatusError as error:
        raise UpstreamError(
            f"{config.TASK_HISTORY_URL} answered {error.response.status_code}"
        ) from error
    except httpx.HTTPError as error:
        raise UpstreamError(f"could not reach {config.TASK_HISTORY_URL}: {error}") from error

    try:
        payload = response.json()
    except ValueError as error:
        raise UpstreamError(f"{config.TASK_HISTORY_URL} did not return JSON") from error

    if not isinstance(payload, dict):
        raise UpstreamError(
            f"{config.TASK_HISTORY_URL} returned {type(payload).__name__}, expected an object"
        )
    return payload


def _as_int(value: object, default: int) -> int:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default
