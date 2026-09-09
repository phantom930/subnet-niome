import { Component, computed, inject, signal } from '@angular/core';
import { DatePipe, DecimalPipe } from '@angular/common';
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
import { MessageService } from 'primeng/api';

import { TaskService } from '../../core/task.service';
import { RefreshResult, SnapshotMeta, TaskRow } from '../../core/task.models';

type SeedFilter = 'stamped' | 'unstamped' | null;

@Component({
  selector: 'app-tasks',
  providers: [MessageService],
  imports: [
    DatePipe,
    DecimalPipe,
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
  ];

  protected readonly rows = computed<TaskRow[]>(() => {
    const term = this.search().trim().toLowerCase();
    const cellType = this.cellType();
    const seed = this.seedFilter();

    return this.allRows().filter((row) => {
      if (cellType && row.cellType !== cellType) return false;
      if (seed === 'stamped' && row.seed === null) return false;
      if (seed === 'unstamped' && row.seed !== null) return false;
      if (!term) return true;

      const seedTerm = term.replace(/,/g, '');
      return (
        row.id.toLowerCase().includes(term) ||
        row.cellType.toLowerCase().includes(term) ||
        (seedTerm.length > 0 && String(row.seed ?? '').includes(seedTerm)) ||
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
}
