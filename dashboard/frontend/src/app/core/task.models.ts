/** Shapes of testing/task.json, plus the flattened row the table renders. */

/** Rules block of a contract. Identical across every task in the snapshot. */
export interface TaskRules {
  base_padding: number;
  cas_systems: string[];
  max_experiments: number;
  max_mismatches: number;
  proximity_gate: boolean;
}

export interface TaskContract {
  active_mutations: string[];
  cell_type: string;
  /** Mutation to region name, or null when the mutation is not in a named region. */
  mutation_regions: Record<string, string | null>;
  mutation_weights: Record<string, number>;
  rules: TaskRules;
  /**
   * Stamped after the round closes; 0 means never stamped.
   *
   * Mixed type in practice: the backend sends most seeds as raw numbers and
   * some as comma-grouped strings, and the grouping is not always correct
   * (one snapshot value reads '328,371,1000'). TaskService normalizes it.
   */
  seed: string | number;
  version: string;
}

export interface RawTask {
  id: string;
  created_at: string;
  content: {
    contract: TaskContract;
    hbb_reference: unknown;
  };
}

export interface TaskSnapshot {
  source: string;
  leaderboard: string;
  fetched_at: string;
  order: string;
  note: string;
  count: number;
  unstamped: number;
  tasks: RawTask[];
}

/** One mutation with the region and weight it carries in its task. */
export interface MutationEntry {
  name: string;
  region: string | null;
  weight: number | null;
}

/** Flattened task, one per table row. */
export interface TaskRow {
  id: string;
  /** First segment of the id, enough to identify a row at a glance. */
  shortId: string;
  createdAt: Date;
  cellType: string;
  mutations: MutationEntry[];
  /** Normalized numeric seed, or null when the task was never stamped. */
  seed: number | null;
  /** Seed exactly as the backend sent it, kept because its grouping can differ. */
  seedRaw: string | null;
  version: string;
}

/** What POST /api/tasks/refresh reports back. */
export interface RefreshResult {
  added: number;
  restamped: number;
  fetched: number;
  count: number;
  unstamped: number;
  fetched_at: string;
  replaced: boolean;
}

/** Rule values shared by every task, shown once instead of as columns. */
export interface SnapshotMeta {
  source: string;
  fetchedAt: Date;
  count: number;
  unstamped: number;
  /** Rules common to all tasks, or null if they are not in fact uniform. */
  sharedRules: TaskRules | null;
}
