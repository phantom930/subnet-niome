import { Component, computed, effect, inject, input, output, signal } from '@angular/core';
import { DecimalPipe } from '@angular/common';
import { ButtonModule } from 'primeng/button';
import { DialogModule } from 'primeng/dialog';
import { MessageModule } from 'primeng/message';
import { ProgressBarModule } from 'primeng/progressbar';
import { TableModule } from 'primeng/table';
import { TagModule } from 'primeng/tag';

import { BenchmarkService } from '../../core/benchmark.service';
import { BenchmarkJob, TaskRow } from '../../core/task.models';

@Component({
  selector: 'app-benchmark-dialog',
  imports: [
    DecimalPipe,
    ButtonModule,
    DialogModule,
    MessageModule,
    ProgressBarModule,
    TableModule,
    TagModule,
  ],
  templateUrl: './benchmark-dialog.html',
  styleUrl: './benchmark-dialog.scss',
})
export class BenchmarkDialog {
  private readonly benchmarks = inject(BenchmarkService);

  /** The row being benchmarked. Null closes the dialog. */
  readonly row = input<TaskRow | null>(null);
  readonly closed = output<void>();

  protected readonly job = signal<BenchmarkJob | null>(null);
  protected readonly error = signal<string | null>(null);
  protected readonly showRaw = signal(false);

  protected readonly running = computed(() => {
    const status = this.job()?.status;
    return status === 'queued' || status === 'running';
  });

  protected readonly result = computed(() => this.job()?.result ?? null);

  protected readonly validatorRows = computed(() => {
    const validator = this.result()?.validator;
    if (!validator) return [];
    // Ordered as the harness prints them: raw stage outputs, then the factors
    // clamped to [0, 1], and only the factors enter the product.
    const order: Array<[keyof typeof validator, string, number]> = [
      ['n_valid_experiments', 'Valid experiments', 0],
      ['total_weighted_score', 'Total weighted score', 3],
      ['consistency_score', 'Consistency score', 4],
      ['consistency_factor', 'Consistency factor', 4],
      ['distribution_fidelity_score', 'Fidelity score', 4],
      ['distribution_fidelity_factor', 'Fidelity factor', 4],
    ];
    return order
      .filter(([key]) => validator[key] !== undefined)
      .map(([key, label, digits]) => ({ label, value: validator[key] as number, digits }));
  });

  /** final_score = weighted x consistency_factor x fidelity_factor. */
  protected readonly formula = computed(() => {
    const v = this.result()?.validator;
    if (
      !v ||
      v.total_weighted_score === undefined ||
      v.consistency_factor === undefined ||
      v.distribution_fidelity_factor === undefined
    ) {
      return null;
    }
    return {
      weighted: v.total_weighted_score,
      consistency: v.consistency_factor,
      fidelity: v.distribution_fidelity_factor,
    };
  });

  /** The row a run has already been started for, so opening starts once. */
  private startedFor: string | null = null;

  constructor() {
    // Starting from the input rather than a method the parent calls means
    // there is no ordering to get right between binding the row and running.
    effect(() => {
      const row = this.row();
      if (!row) {
        this.startedFor = null;
        return;
      }
      if (row.id !== this.startedFor) {
        this.startedFor = row.id;
        this.start();
      }
    });
  }

  /** Runs the benchmark. Also the Run again button. */
  start(): void {
    const row = this.row();
    if (!row) return;

    this.job.set(null);
    this.error.set(null);
    this.showRaw.set(false);

    this.benchmarks.run({ task: row.id }).subscribe({
      next: (job) => this.job.set(job),
      error: (error: Error) => this.error.set(error.message),
    });
  }

  protected close(): void {
    this.closed.emit();
  }

  protected toggleRaw(): void {
    this.showRaw.update((shown) => !shown);
  }

  protected statusLabel(): string {
    const job = this.job();
    if (!job) return 'Starting';
    if (job.status === 'queued') return 'Queued behind another run';
    if (job.status === 'running') return 'Running the miner and five validation stages';
    return job.status;
  }
}
