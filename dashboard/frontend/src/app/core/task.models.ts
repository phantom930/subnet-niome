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
   * A round carries *several* seeds and the backend joins them with commas
   * ('263,486,269'), so this is a list, not digit grouping — a round has had
   * three seeds since 2026-08-27 and one before that. A single-seed round
   * arrives as a raw number. TaskService parses both into `TaskRow.seeds`.
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
  /** The round's stamped seeds, in the order recorded. Empty when unstamped. */
  seeds: number[];
  /** Sort key for the seed column: the round's first seed, null when unstamped. */
  seedSort: number | null;
  /** The seed field exactly as the backend sent it, for the 'as sent' tooltip. */
  seedRaw: string | null;
  version: string;
  /**
   * The task exactly as the backend sent it, for the JSON view.
   *
   * Kept rather than refetched because it is already loaded: a task is about
   * 1.2 kB, so all of them together are the half-megabyte the snapshot
   * response already carried.
   */
  raw: RawTask;
}

/**
 * How often the backend has stamped each seed — GET /api/seed-occurrence,
 * served straight from miner_data/drawn_seeds.json.
 *
 * The same ledger design.draw_seeds bets against, so what the Tasks page
 * colours by is what the miner acts on rather than a second count beside it.
 */
export interface SeedOccurrence {
  /** False when the ledger is missing; the page then shows seeds uncoloured. */
  available: boolean;
  /** Seed (as a string key) to the number of rounds that stamped it. */
  counts: Record<string, number>;
  /** Why it is unavailable, when it is. */
  reason?: string;
  /**
   * Rounds before this were never counted, so their seeds have no occurrence
   * rather than an occurrence of zero. The ledger starts at the 3-seed
   * cutover because the generator before it drew from a different range.
   */
  since?: string | null;
  tasks_recorded?: number;
  latest_task_at?: string | null;
  /** Seeds in [100,1000] the backend has never stamped. */
  undrawn?: number;
  updated_at?: string | null;
}

/**
 * One row of the accessibility table the backend serves.
 *
 * Accessibility is the largest term in stage 3's energy, so it sets the cut
 * probability. A cell type missing from this table is scored as though its
 * chromatin were fully open, because the lookup defaults to 1.0.
 */
export interface CellType {
  accessibility: number;
  /** 'sourced' or 'estimated', per the upstream table. */
  basis?: string;
  note?: string;
}

/**
 * What POST /api/tasks/refresh reports back.
 *
 * The backend runs scripts/bench_task.py --fetch and then diffs the snapshot
 * on disk, so these counts come from the files rather than the script's
 * stdout. `output` is what the script printed, for the details.
 */
export interface RefreshResult {
  added: number;
  restamped: number;
  removed: number;
  count: number;
  unstamped: number;
  fetched_at: string | null;
  cell_types: number;
  replaced: boolean;
  output: string;
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

/** ---- Benchmark jobs ---- */

export type JobStatus = 'queued' | 'running' | 'done' | 'failed';

export interface BenchmarkRequest {
  task: string;
  /** How many seeds to draw when they are drawn: an unstamped task, or `random_seeds`. */
  seeds?: number;
  rng?: number | null;
  /**
   * Score under random seeds even though the task carries its own. The
   * default is the seeds the round actually closed under.
   */
  random_seeds?: boolean;
  per_seed?: boolean;
  uid?: number;
}

export interface PerSeedScore {
  seed: number;
  total_weighted_score: number;
  consistency_factor: number;
  distribution_fidelity_factor: number;
  final_score: number;
}

/**
 * Numbers read back out of the harness's printed report.
 *
 * Best-effort: bench_task.py has no JSON mode, so the backend parses its
 * stdout. Any field can be absent if the format changes, which is why the job
 * also carries the raw `output` and the dialog shows it.
 */
export interface BenchmarkResult {
  task: {
    id?: string;
    created_at?: string;
    cell_type?: string;
    accessibility?: number;
    mutations?: string[];
    recorded_seed?: string;
    weights?: string;
    rules?: string;
  };
  miner: Record<string, string | number>;
  validator: {
    n_valid_experiments?: number;
    total_weighted_score?: number;
    consistency_score?: number;
    consistency_factor?: number;
    distribution_fidelity_score?: number;
    distribution_fidelity_factor?: number;
    final_score?: number;
  };
  seeds: number[];
  /**
   * Where those seeds came from: the task's own recorded ones, or a random
   * draw. Absent on an older backend, or if the report's format changes.
   */
  seed_source?: string | null;
  per_seed: PerSeedScore[];
  /** Null for a single-seed run; absent on older backends. */
  spread?: number | null;
}

export interface BenchmarkJob {
  id: string;
  status: JobStatus;
  request: Required<BenchmarkRequest>;
  created_at: number;
  started_at: number | null;
  finished_at: number | null;
  duration_seconds: number | null;
  /** The report as the harness printed it. The source of truth. */
  output: string | null;
  result: BenchmarkResult | null;
  error: string | null;
}
