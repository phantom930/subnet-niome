import { Component, computed, input, output, signal } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { ButtonModule } from 'primeng/button';
import { DialogModule } from 'primeng/dialog';
import { SelectButtonModule } from 'primeng/selectbutton';

import { TaskRow } from '../../core/task.models';

type Scope = 'task' | 'contract';

@Component({
  selector: 'app-task-json-dialog',
  imports: [FormsModule, ButtonModule, DialogModule, SelectButtonModule],
  templateUrl: './task-json-dialog.html',
  styleUrl: './task-json-dialog.scss',
})
export class TaskJsonDialog {
  /** The row whose JSON to show. Null closes the dialog. */
  readonly row = input<TaskRow | null>(null);
  readonly closed = output<void>();

  protected readonly scope = signal<Scope>('task');
  protected readonly copied = signal(false);
  protected readonly copyFailed = signal(false);

  protected readonly scopeOptions: Array<{ label: string; value: Scope }> = [
    { label: 'Whole task', value: 'task' },
    { label: 'Contract only', value: 'contract' },
  ];

  /**
   * Pretty-printed JSON, rendered as text rather than markup. Highlighting it
   * would mean building HTML from backend data, which is not worth the
   * injection risk for a viewer.
   */
  protected readonly json = computed(() => {
    const row = this.row();
    if (!row) return '';
    const subject = this.scope() === 'contract' ? row.raw.content.contract : row.raw;
    return JSON.stringify(subject, null, 2);
  });

  protected readonly lineCount = computed(() => this.json().split('\n').length);
  protected readonly byteCount = computed(() => new Blob([this.json()]).size);

  protected setScope(scope: Scope): void {
    this.scope.set(scope);
    this.copied.set(false);
    this.copyFailed.set(false);
  }

  protected async copy(): Promise<void> {
    this.copyFailed.set(false);
    try {
      await navigator.clipboard.writeText(this.json());
      this.copied.set(true);
      setTimeout(() => this.copied.set(false), 2000);
    } catch {
      // Blocked clipboard, or a context the browser does not consider secure.
      this.copyFailed.set(true);
    }
  }

  protected close(): void {
    this.copied.set(false);
    this.copyFailed.set(false);
    this.closed.emit();
  }
}
