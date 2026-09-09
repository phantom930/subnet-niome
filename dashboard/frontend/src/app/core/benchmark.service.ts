import { HttpClient, HttpErrorResponse } from '@angular/common/http';
import { Injectable, inject } from '@angular/core';
import { EMPTY, Observable, catchError, expand, switchMap, throwError, timer } from 'rxjs';

import { BenchmarkJob, BenchmarkRequest } from './task.models';

/** How often to ask the backend whether a run has finished. */
const POLL_INTERVAL_MS = 1500;

/**
 * Runs benchmarks through the dashboard backend, which shells out to
 * scripts/bench_task.py.
 *
 * A run takes seconds and queues behind any other harness run, so this starts
 * a job and then polls it. `run` emits every state the job passes through,
 * queued and running included, so the caller can show progress, and completes
 * once it settles.
 */
@Injectable({ providedIn: 'root' })
export class BenchmarkService {
  private readonly http = inject(HttpClient);

  run(request: BenchmarkRequest): Observable<BenchmarkJob> {
    return this.start(request).pipe(
      // Emits the first job, then each poll result, stopping on the terminal
      // state. EMPTY ends the recursion.
      expand((job) =>
        isSettled(job) ? EMPTY : timer(POLL_INTERVAL_MS).pipe(switchMap(() => this.get(job.id))),
      ),
      catchError((error) => throwError(() => new Error(describeError(error)))),
    );
  }

  private start(request: BenchmarkRequest): Observable<BenchmarkJob> {
    return this.http.post<BenchmarkJob>('/api/benchmarks', {
      seeds: 3,
      per_seed: true,
      ...request,
    });
  }

  get(jobId: string): Observable<BenchmarkJob> {
    return this.http.get<BenchmarkJob>(`/api/benchmarks/${jobId}`);
  }
}

function isSettled(job: BenchmarkJob): boolean {
  return job.status === 'done' || job.status === 'failed';
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
