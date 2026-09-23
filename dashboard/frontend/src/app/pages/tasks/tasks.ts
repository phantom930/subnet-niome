import { Component, computed, inject, signal } from '@angular/core';
import { DatePipe, DecimalPipe, SlicePipe } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ButtonModule } from 'primeng/button';
import { CardModule } from 'primeng/card';
import { IconFieldModule } from 'primeng/iconfield';
import { InputIconModule } from 'primeng/inputicon';
import { InputTextModule } from 'primeng/inputtext';
import { MessageModule } from 'primeng/message';
import { SelectModule } from 'primeng/select';
import { SkeletonModule } from 'primeng/skeleton';
import { TableModule } from 'primeng/table';
import { TagModule } from 'primeng/tag';
import { ToastModule } from 'primeng/toast';
import { TooltipModule } from 'primeng/tooltip';
import { MessageService } from 'primeng/api';

import { BenchmarkDialog } from './benchmark-dialog';
import { TaskJsonDialog } from './task-json-dialog';

import { TaskService } from '../../core/task.service';
import { RefreshResult, SeedOccurrence, SnapshotMeta, TaskRow } from '../../core/task.models';

type SeedFilter = 'stamped' | 'unstamped' | 'repeated' | null;

/**
 * One seed of a round, with how often the backend has stamped it.
 *
 * `level` drives the colour and is capped at 4: the tail past four is a
 * handful of seeds out of ~490, so giving each count its own step would spend
 * the ramp's resolution where almost nothing lands. The exact count is in the
 * chip's label and tooltip, which is where a 5 or a 6 is read. `level` 0 means
 * "not counted": the round predates the ledger's window, or the ledger is
 * missing. That is deliberately not the same as a count of zero, which no
 * stamped seed can have.
 */
export interface SeedChip {
  seed: number;
  /** Times stamped, or null when this round was never counted. */
  count: number | null;
  level: 0 | 1 | 2 | 3 | 4;
  /** The tooltip: the count in words, or why there is none. */
  title: string;
}

/** A task row with its seeds resolved against the occurrence ledger. */
export interface TaskRowView extends TaskRow {
  chips: SeedChip[];
  /** The round's highest seed occurrence, for the 'repeated' filter. */
  maxOccurrence: number;
}

@Component({
  selector: 'app-tasks',
  providers: [MessageService],
  imports: [
    DatePipe,
    DecimalPipe,
    SlicePipe,
    FormsModule,
    ButtonModule,
    CardModule,
    IconFieldModule,
    InputIconModule,
    InputTextModule,
    MessageModule,
    SelectModule,
    SkeletonModule,
    TableModule,
    TagModule,
    ToastModule,
    TooltipModule,
    BenchmarkDialog,
    TaskJsonDialog,
  ],
  templateUrl: './tasks.html',
  styleUrl: './tasks.scss',
})
export class Tasks {
  private readonly tasks = inject(TaskService);
  private readonly toast = inject(MessageService);

  protected readonly loading = signal(true);
  protected readonly refreshing = signal(false);
  protected readonly error = signal<string | null>(null);
  protected readonly meta = signal<SnapshotMeta | undefined>(undefined);

  private readonly allRows = signal<TaskRow[]>([]);
  protected readonly allRowCount = computed(() => this.allRows().length);

  protected readonly search = signal('');
  protected readonly cellType = signal<string | null>(null);
  protected readonly seedFilter = signal<SeedFilter>(null);

  /** Page offset, reset whenever a filter shrinks the result set. */
  protected readonly first = signal(0);
  protected readonly rowsPerPage = 25;

  constructor() {
    this.reload();
  }

  protected readonly cellTypeOptions = computed(() => {
    const seen = new Set(this.allRows().map((r) => r.cellType));
    return [...seen].sort().map((value) => ({ label: value, value }));
  });

  protected readonly seedOptions: Array<{ label: string; value: SeedFilter }> = [
    { label: 'Stamped', value: 'stamped' },
    { label: 'Unstamped', value: 'unstamped' },
    { label: 'Has a repeat', value: 'repeated' },
  ];

  /** The occurrence ledger. Absent until it loads, and after it fails. */
  protected readonly occurrence = signal<SeedOccurrence | null>(null);

  protected readonly occurrenceNote = computed(() => {
    const ledger = this.occurrence();
    if (!ledger) return null;
    if (!ledger.available) {
      return `Seed occurrence unavailable (${ledger.reason ?? 'no ledger'}), so seeds are shown
        uncoloured. A refresh writes miner_data/drawn_seeds.json.`;
    }
    return null;
  });

  /**
   * Every row with its seeds resolved against the ledger.
   *
   * Built once per snapshot rather than per render: the filters below read
   * `maxOccurrence` off it, and the template reads the chips, so resolving
   * inside the filter's own computed would redo the lookup on every keystroke.
   */
  private readonly resolvedRows = computed<TaskRowView[]>(() => {
    const ledger = this.occurrence();
    const counts = ledger?.available ? ledger.counts : null;
    // The ledger counts only rounds at or after its cutover, so a seed from an
    // older round has no occurrence — not an occurrence of zero. Colouring
    // those would claim the whole single-seed era was never stamped.
    const since = ledger?.since ? new Date(ledger.since) : null;

    return this.allRows().map((row) => {
      const counted = counts !== null && (since === null || row.createdAt >= since);
      const chips: SeedChip[] = row.seeds.map((seed) => {
        if (!counted) {
          return {
            seed,
            count: null,
            level: 0 as const,
            title:
              counts === null
                ? `${seed} — occurrence unavailable`
                : `${seed} — this round predates the occurrence ledger, so it is not counted`,
          };
        }
        const count = counts[String(seed)] ?? 0;
        return { seed, count, level: levelFor(count), title: describeOccurrence(seed, count) };
      });

      return {
        ...row,
        chips,
        maxOccurrence: chips.reduce((worst, chip) => Math.max(worst, chip.count ?? 0), 0),
      };
    });
  });

  protected readonly rows = computed<TaskRowView[]>(() => {
    const term = this.search().trim().toLowerCase();
    const cellType = this.cellType();
    const seed = this.seedFilter();

    return this.resolvedRows().filter((row) => {
      if (cellType && row.cellType !== cellType) return false;
      if (seed === 'stamped' && row.seeds.length === 0) return false;
      if (seed === 'unstamped' && row.seeds.length > 0) return false;
      if (seed === 'repeated' && row.maxOccurrence < 2) return false;
      if (!term) return true;

      // Matched per seed, so '486' finds the round that stamped it rather than
      // only the rounds whose joined seed string happens to contain those
      // digits — '48' used to match '1,486,203' across two seeds.
      return (
        row.id.toLowerCase().includes(term) ||
        row.cellType.toLowerCase().includes(term) ||
        row.seeds.some((value) => String(value).includes(term)) ||
        row.mutations.some(
          (m) =>
            m.name.toLowerCase().includes(term) ||
            (m.region?.toLowerCase().includes(term) ?? false),
        )
      );
    });
  });

  protected readonly isFiltered = computed(
    () => this.search().trim().length > 0 || this.cellType() !== null || this.seedFilter() !== null,
  );

  /** Fetch the snapshot from the backend. */
  protected reload(): void {
    this.loading.set(true);
    this.error.set(null);
    // Independent of the table's own load: the ledger never rejects, and the
    // rows render uncoloured until it arrives.
    this.tasks.loadSeedOccurrence().subscribe((ledger) => this.occurrence.set(ledger));
    this.tasks.load().subscribe({
      next: (data) => {
        this.meta.set(data.meta);
        this.allRows.set(data.rows);
        this.first.set(0);
        this.loading.set(false);
      },
      error: (error: Error) => {
        this.error.set(error.message);
        this.allRows.set([]);
        this.meta.set(undefined);
        this.loading.set(false);
      },
    });
  }

  /** Ask the backend to pull the upstream history, then show the result. */
  protected refresh(): void {
    if (this.refreshing()) return;
    this.refreshing.set(true);

    this.tasks.refresh().subscribe({
      next: (result) => {
        this.refreshing.set(false);
        this.toast.add({
          severity: this.changed(result) ? 'success' : 'info',
          summary: this.refreshSummary(result),
          detail: `Now holding ${result.count} tasks, ${result.unstamped} unstamped`,
          life: 5000,
        });
        this.reload();
      },
      error: (error: Error) => {
        this.refreshing.set(false);
        this.toast.add({
          severity: 'error',
          summary: 'Refresh failed',
          detail: error.message,
          life: 9000,
        });
      },
    });
  }

  /** The row whose benchmark dialog is open, or null when it is closed. */
  protected readonly benchmarkRow = signal<TaskRow | null>(null);

  protected openBenchmark(row: TaskRow): void {
    // Setting the row both opens the dialog and starts the run: the dialog
    // watches this input, so there is no start() call to sequence against it.
    this.benchmarkRow.set(row);
  }

  protected closeBenchmark(): void {
    this.benchmarkRow.set(null);
  }

  /** The row whose JSON dialog is open, or null when it is closed. */
  protected readonly jsonRow = signal<TaskRow | null>(null);

  protected openJson(row: TaskRow): void {
    this.jsonRow.set(row);
  }

  protected closeJson(): void {
    this.jsonRow.set(null);
  }

  private changed(result: RefreshResult): boolean {
    return result.added > 0 || result.restamped > 0 || result.removed > 0;
  }

  /**
   * A refresh that adds nothing can still have done something: a task that was
   * unstamped last time carries its real seed now. Both are worth saying.
   */
  private refreshSummary(result: RefreshResult): string {
    const parts: string[] = [];
    if (result.added > 0) parts.push(`${result.added} new`);
    if (result.restamped > 0) parts.push(`${result.restamped} newly stamped`);
    if (result.removed > 0) parts.push(`${result.removed} dropped`);
    return parts.length ? parts.join(', ') : 'Already up to date';
  }

  protected setSearch(value: string): void {
    this.search.set(value);
    this.first.set(0);
  }

  protected setCellType(value: string | null): void {
    this.cellType.set(value);
    this.first.set(0);
  }

  protected setSeedFilter(value: SeedFilter): void {
    this.seedFilter.set(value);
    this.first.set(0);
  }

  protected clearFilters(): void {
    this.search.set('');
    this.cellType.set(null);
    this.seedFilter.set(null);
    this.first.set(0);
  }

  /**
   * The legend. Always rendered when the ledger is available, because the
   * colour is an ordinal ramp and a ramp without a key is decoration.
   */
  protected readonly legend: Array<{ level: 0 | 1 | 2 | 3 | 4; label: string }> = [
    { level: 1, label: 'once' },
    { level: 2, label: 'twice' },
    { level: 3, label: '3x' },
    { level: 4, label: '4x or more' },
  ];
}

/**
 * Times stamped to a ramp step. Everything from four up shares the top step —
 * the legend says "4x or more" and the chip prints the real count.
 */
function levelFor(count: number): 0 | 1 | 2 | 3 | 4 {
  if (count <= 0) return 0;
  if (count >= 4) return 4;
  return count as 1 | 2 | 3;
}

/** The chip's tooltip, in words rather than a bare number. */
function describeOccurrence(seed: number, count: number): string {
  if (count <= 0) return `${seed} — not stamped on any counted round`;
  if (count === 1) return `${seed} — stamped once, this round only so far`;
  if (count === 2) return `${seed} — stamped twice, counting this round`;
  return `${seed} — stamped ${count} times, counting this round`;
}
