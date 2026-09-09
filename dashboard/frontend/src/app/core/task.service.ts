import { HttpClient, HttpErrorResponse } from '@angular/common/http';
import { Injectable, inject } from '@angular/core';
import { Observable, catchError, map, throwError } from 'rxjs';

import {
  MutationEntry,
  RawTask,
  RefreshResult,
  SnapshotMeta,
  TaskRow,
  TaskRules,
  TaskSnapshot,
} from './task.models';

export interface TaskData {
  meta: SnapshotMeta;
  rows: TaskRow[];
}

/**
 * Turns the backend's mixed-type seed into a number, or null when the round
 * never closed. Seeds arrive both as raw numbers and as comma-grouped strings,
 * so the commas come out before parsing. Zero means unstamped.
 */
function normalizeSeed(seed: string | number | null | undefined): number | null {
  if (seed === null || seed === undefined || seed === '') return null;
  const value = typeof seed === 'number' ? seed : Number(String(seed).replace(/,/g, ''));
  if (!Number.isFinite(value) || value === 0) return null;
  return value;
}

/**
 * Talks to the dashboard backend in dashboard/backend.
 *
 * The dev server proxies /api to it, see proxy.conf.json. Start the backend
 * first or these calls fail with a connection error, which the Tasks page
 * shows rather than swallowing.
 */
@Injectable({ providedIn: 'root' })
export class TaskService {
  private readonly http = inject(HttpClient);

  /** Not cached: the refresh button needs the next load to see new data. */
  load(): Observable<TaskData> {
    return this.http.get<TaskSnapshot>('/api/tasks').pipe(
      map((snapshot) => this.toTaskData(snapshot)),
      catchError((error) => throwError(() => new Error(describeError(error)))),
    );
  }

  /** Pull the upstream task history and merge it into the stored snapshot. */
  refresh(replace = false): Observable<RefreshResult> {
    return this.http
      .post<RefreshResult>('/api/tasks/refresh', { replace })
      .pipe(catchError((error) => throwError(() => new Error(describeError(error)))));
  }

  private toTaskData(snapshot: TaskSnapshot): TaskData {
    const rows = snapshot.tasks.map((task) => this.toRow(task));
    return {
      meta: {
        source: snapshot.source,
        fetchedAt: new Date(snapshot.fetched_at),
        count: snapshot.count ?? rows.length,
        unstamped: snapshot.unstamped ?? rows.filter((r) => r.seed === null).length,
        sharedRules: this.findSharedRules(snapshot.tasks),
      },
      rows,
    };
  }

  private toRow(task: RawTask): TaskRow {
    const contract = task.content.contract;
    const mutations: MutationEntry[] = contract.active_mutations.map((name) => ({
      name,
      region: contract.mutation_regions?.[name] ?? null,
      weight: contract.mutation_weights?.[name] ?? null,
    }));

    const seed = normalizeSeed(contract.seed);

    return {
      id: task.id,
      shortId: task.id.split('-')[0],
      createdAt: new Date(task.created_at),
      cellType: contract.cell_type,
      mutations,
      seed,
      seedRaw: seed === null ? null : String(contract.seed),
      version: contract.version,
    };
  }

  /**
   * Every task in the snapshot carries the same rules, so they are worth
   * showing once above the table rather than as five identical columns.
   * Returns null if that ever stops being true, and the UI hides the panel.
   */
  private findSharedRules(tasks: RawTask[]): TaskRules | null {
    if (!tasks.length) return null;
    const first = tasks[0].content.contract.rules;
    const fingerprint = JSON.stringify(first);
    const uniform = tasks.every((t) => JSON.stringify(t.content.contract.rules) === fingerprint);
    return uniform ? first : null;
  }
}

/** A message worth showing the user, rather than 'Http failure response'. */
function describeError(error: unknown): string {
  if (!(error instanceof HttpErrorResponse)) {
    return error instanceof Error ? error.message : 'Unknown error';
  }
  if (error.status === 0) {
    return 'Cannot reach the dashboard backend. Start it with "uvicorn main:app --port 8000" in dashboard/backend.';
  }
  const detail = (error.error as { detail?: string } | null)?.detail;
  if (detail) return detail;
  return `Backend answered ${error.status} ${error.statusText}`;
}
