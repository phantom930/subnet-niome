import { Component, computed, effect, inject, input, output, signal } from '@angular/core';
import { DecimalPipe } from '@angular/common';
import { ButtonModule } from 'primeng/button';
import { DialogModule } from 'primeng/dialog';
import { MessageModule } from 'primeng/message';
import { ProgressBarModule } from 'primeng/progressbar';
import { TableModule } from 'primeng/table';
import { TagModule } from 'primeng/tag';

import { BenchmarkService } from '../../core/benchmark.service';
import { BenchmarkJob, BenchmarkResult, TaskRow } from '../../core/task.models';

/** The hue a number is drawn in, so a value and its label read as one thing. */
type Tone = 'count' | 'weighted' | 'consistency' | 'fidelity';

/**
 * The validator's fields, ordered as the harness prints them: each raw stage
 * output, then the factor clamped to [0, 1] beside it.
 *
 * `tone` groups a raw output with its factor and with the matching term in the
 * formula, so the three colours are enough to read the product off the table.
 * `inProduct` marks the three fields that are actually multiplied — the two
 * raw scores are not.
 */
const VALIDATOR_FIELDS: ReadonlyArray<{
  key: keyof BenchmarkResult['validator'];
  label: string;
  digits: number;
  tone: Tone;
  inProduct: boolean;
}> = [
  {
    key: 'n_valid_experiments',
    label: 'Valid experiments',
    digits: 0,
    tone: 'count',
    inProduct: false,
  },
  {
    key: 'total_weighted_score',
    label: 'Total weighted score',
    digits: 3,
    tone: 'weighted',
    inProduct: true,
  },
  {
    key: 'consistency_score',
    label: 'Consistency score',
    digits: 4,
    tone: 'consistency',
    inProduct: false,
  },
  {
    key: 'consistency_factor',
    label: 'Consistency factor',
    digits: 4,
    tone: 'consistency',
    inProduct: true,
  },
  {
    key: 'distribution_fidelity_score',
    label: 'Fidelity score',
    digits: 4,
    tone: 'fidelity',
    inProduct: false,
  },
  {
    key: 'distribution_fidelity_factor',
    label: 'Fidelity factor',
    digits: 4,
    tone: 'fidelity',
    inProduct: true,
  },
];

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
    return VALIDATOR_FIELDS.filter((field) => validator[field.key] !== undefined).map((field) => ({
      ...field,
      value: validator[field.key] as number,
    }));
  });

  /**
   * The per-seed rows, with the ends of the spread marked.
   *
   * Which seed a round draws is the one thing a miner cannot design for, so the
   * best and worst draws are the numbers worth finding first. A single seed, or
   * seeds that all scored the same, has no spread to point at — marking a row
   * there would be colour that means nothing.
   */
  protected readonly perSeedRows = computed(() => {
    const rows = this.result()?.per_seed ?? [];
    const finals = rows.map((row) => row.final_score);
    const best = Math.max(...finals);
    const worst = Math.min(...finals);
    const hasSpread = rows.length > 1 && best !== worst;
    return rows.map((row) => ({
      ...row,
      best: hasSpread && row.final_score === best,
      worst: hasSpread && row.final_score === worst,
    }));
  });

  /**
   * Where the seeds came from, as a tag.
   *
   * A score under the seeds the round closed under is the one that round paid;
   * a score under a random draw is a hypothetical. That difference changes what
   * the number means, so it reads as a caveat rather than as a detail — amber
   * for a draw nobody played, the score's own emerald for the real thing.
   *
   * The harness says it in a full phrase, which is too long to sit beside the
   * score. Only the one distinction matters here, so the phrase becomes the
   * tooltip and the tag carries a couple of words.
   */
  protected readonly seedTag = computed(() => {
    const source = this.result()?.seed_source;
    if (!source) return null;
    const random = source.includes('random');
    return {
      label: random ? 'random draw' : "the round's own",
      severity: random ? ('warn' as const) : ('success' as const),
      detail: source,
    };
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
