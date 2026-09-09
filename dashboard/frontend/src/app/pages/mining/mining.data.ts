/**
 * Content for the Mining page.
 *
 * Sourced from the repo's own documentation and code, so it can be checked:
 *   docs/miner_flow.md                            the request path and the build
 *   docs/validation_pipeline.md                   what each factor measures
 *   niome_subnet/genomics/validation/stage3.py    the energy and cut formulas
 *   niome_subnet/genomics/validation/stage5.py    the six ratios and their floor
 *   niome_subnet/genomics/design.py               the build itself
 *   neurons/miner.py                              request handling and upload
 *   niome_subnet/validator/forward.py             the S3 round trip
 *
 * Kept in one file rather than spread through the template so that when the
 * docs change there is a single place to update, and every step names the file
 * it came from.
 */

export interface Step {
  /** Label the miner logs, e.g. '1/5', or a word for the unnumbered stages. */
  tag: string;
  title: string;
  /** Where in the codebase this happens. */
  source: string;
  /** One or two sentences on what it does. */
  summary: string;
  /** The detail worth knowing, each a self-contained point. */
  details: string[];
  icon: string;
}

/** The path a broadcast takes through the miner, in the order it happens. */
export const REQUEST_STEPS: Step[] = [
  {
    tag: 'verify',
    title: 'Signature and blacklist',
    source: 'niome_subnet/base/miner.py, Miner.blacklist',
    summary:
      'The miner is a FastAPI server with one route, POST /forward, not a bt.Axon. The validator signs the request body with its hotkey, and verification runs before any of the miner’s own code.',
    details: [
      'The server listens on --axon.port, default 8091, and publishes its address on chain via ServeAxon. The dataset never travels in the HTTP response.',
      'Two verify() defaults must be overridden or every broadcast is rejected. require_receiver has to be false, because the validator signs without a receiver_ss58 and the default raises a 401 before the signature is even checked.',
      'max_age has to be raised to about 50 seconds. The window is spent on host clock skew, not latency, so a miner whose clock runs ten seconds fast rejects every task it is sent and looks unreachable.',
      '58 seconds is the hard ceiling on max_age, not 60, because verify() refuses max_age plus allowed_skew above the nonce store’s 60 second retention. Exceeding it surfaces as a 500 on valid traffic rather than a 401.',
      'The blacklist rejects hotkeys that are not in the metagraph. Both of its flags default to false, so out of the box any registered hotkey is accepted whether or not it holds a validator permit. Set force_validator_permit to require one.',
    ],
    icon: 'pi pi-shield',
  },
  {
    tag: '1/5',
    title: 'Parse the task and record the target',
    source: 'Miner.forward, _upload_deadline, _record_upload',
    summary:
      'The synapse is parsed, the upload deadline is worked out, and the target is written to disk. Then forward returns an empty acknowledgement so the validator stops waiting.',
    details: [
      'The presigned URL states its own expiry, Expires on SigV2 or X-Amz-Date plus X-Amz-Expires on SigV4. Reading it beats assuming the URL was minted the moment it arrived, and it is clamped to the local 300 second budget so a slow clock cannot invent a deadline that has already passed.',
      'The target is recorded before any work starts, so a round lost later can still be uploaded by hand while the URL lives.',
      'A task with no presigned_url is abandoned and logged. There is nowhere to upload, so the round cannot be scored.',
    ],
    icon: 'pi pi-inbox',
  },
  {
    tag: '2/5',
    title: 'Fetch the contract and the HBB reference',
    source: 'Miner._fetch_artifacts',
    summary:
      'Two plain unsigned GETs, since the presigning is already in the URL. Three retries with linear backoff.',
    details: [
      'Both documents are persisted under miner_data/, never data/, which belongs to the validation pipeline. That is what lets a restart prewarm its caches and an operator re-score offline.',
      'Everything from here runs in worker threads via asyncio.to_thread, so /forward stays answerable during a build.',
    ],
    icon: 'pi pi-download',
  },
  {
    tag: '3/5',
    title: 'Fetch the cell-type table',
    source: 'Miner._fetch_cell_types',
    summary:
      'The accessibility table for the contract’s cell type. Rows stay valid without it; what is at risk is the local prediction and the row allocation.',
    details: [
      'Accessibility is the largest single term in stage 3’s energy, so it sets both the cut probability and the repair mix.',
      'If the fetch fails the cached table is used. With no cache either, accessibility defaults to 1.0, which makes both the prediction and the pricing optimistic.',
    ],
    icon: 'pi pi-table',
  },
  {
    tag: '4/5',
    title: 'Build the submission',
    source: 'Miner._build, design.build_context, design.build',
    summary:
      'The rows are generated, checked against the invariants, and persisted. This is the part worth optimising.',
    details: [
      'Held under a build lock and memoised on sha256 of the task id plus the canonical contract. Every validator broadcasts the same task with its own URL, and the rows are a deterministic function of the contract, so later broadcasts reuse the first build and repeat only the upload.',
      'A contract that changes under one task id rebuilds, rather than re-uploading rows designed against the old rules.',
    ],
    icon: 'pi pi-cog',
  },
  {
    tag: 'optional',
    title: 'Score locally',
    source: 'Miner._score_locally, design.score_rows',
    summary:
      'Validators compute the score independently and never send it back, so this replica is the only signal available before the next task.',
    details: [
      'It calls stage 4’s and stage 5’s own functions, so the number it reports is the number they would report.',
      'The contract a miner receives carries seed 0, so with no real seed the report is a mean over sampled seeds.',
      'Skipped when the remaining time is too short, and a raised exception is logged and swallowed. A failed prediction must never cost the upload.',
    ],
    icon: 'pi pi-chart-line',
  },
  {
    tag: '5/5',
    title: 'Upload the array',
    source: 'Miner._upload, _upload_headers',
    summary:
      'A PUT of the bare JSON array to the presigned URL, which the validator minted for the S3 key niome/{uid}.json. Three retries, each bounded by the remaining time.',
    details: [
      'The body is the array of row objects with no envelope around it.',
      'The URL goes out exactly as received. A re-encoded query string is a SignatureDoesNotMatch.',
      'A SigV2 URL signs Content-Type as the empty string, so sending Content-Type: application/json makes S3 hash a different string and reject the upload. V2 URLs get no headers at all. A SigV4 URL covers only the headers named in X-Amz-SignedHeaders, so there Content-Type is required exactly when it was signed and forbidden otherwise.',
      'If the deadline passes, the round is lost. A manual PUT to the recorded URL is the only recovery.',
    ],
    icon: 'pi pi-upload',
  },
];

/** The design pipeline inside step 4, in order. */
export const BUILD_STEPS: Step[] = [
  {
    tag: '1',
    title: 'build_context',
    source: 'design.build_context',
    summary:
      'Loads chr11 from the process cache and builds a 12-mer index over the gene region plus or minus 50 kb.',
    details: [
      'The same window and the same k the validator’s off-target check hard-codes.',
      'Built in memory rather than through the pipeline’s pickle cache, which would write into data/.',
    ],
    icon: 'pi pi-database',
  },
  {
    tag: '2',
    title: 'enumerate_coordinates',
    source: 'design.enumerate_coordinates',
    summary:
      'Finds every position where a PAM can exist, then confirms each one through the real gate.',
    details: [
      'check_pam reads a fixed motif at a fixed offset, so the candidate positions are exactly the occurrences of a short literal: GG or CC for Cas9’s NGG, TTT or AAA for Cas12a’s TTTV once the minus strand is read off the reverse complement.',
      'Every hit is then confirmed through check_pam itself, so a change to the gate cannot leave a stale coordinate behind.',
    ],
    icon: 'pi pi-search',
  },
  {
    tag: '3',
    title: 'One coordinate per cell',
    source: 'design._cell_coordinates',
    summary:
      'Each of the eight coverage cells gets exactly one coordinate, filled scarcest group first.',
    details: [
      'A coordinate serves exactly one cell because stage 1 dedups on (cas, start, strand, guide) and the mutation is not in that key. Two mutations sharing a coordinate would collide on any guide they both used.',
      'Cas12a is filled first, since its TTTV PAM is several times rarer than Cas9’s NGG and it would otherwise be left with whatever Cas9 declined.',
    ],
    icon: 'pi pi-th-large',
  },
  {
    tag: '4',
    title: 'enumerate_guides and gate',
    source: 'design.enumerate_guides, design.gate_and_score',
    summary:
      'Spells guide variants within the contract’s mismatch budget, then scores them through the real stage 1 and stage 2.',
    details: [
      'The mismatch budget is the design’s only free lever and does three jobs: a class flip moves a base between {G,C} and {A,T} to pull GC to exactly 50 percent where gc_score peaks at 1.0; a within-class swap leaves the count alone and exists only to spell another distinct guide; and any variant whose 12-mer seed still appears in the index is discarded, which is what takes offtarget_factor from 0.7 to 1.0.',
      'Every guide returned sits at the same coordinate at the same GC count, so stages 1, 2 and 5 see one feature vector for the whole set.',
    ],
    icon: 'pi pi-filter',
  },
  {
    tag: '5',
    title: 'Grow until the rows are there',
    source: 'design.build, the growth loop',
    summary:
      'If the pools cannot fill the row count, the number of sites per cell quadruples and the pools are rebuilt.',
    details: [
      'One coordinate per cell only fills the submission if the mismatch budget can spell that many distinct guides on it. Three free substitutions on a 20-mer give over a thousand.',
      'A contract with max_mismatches 0 allows exactly one guide per coordinate, so a cell contributes a single row. Measured on that case, growing took 8 rows to 250 and total_weighted_score from 7.0 to 233.4.',
    ],
    icon: 'pi pi-arrows-alt',
  },
  {
    tag: '6',
    title: 'allocate_rows',
    source: 'design.allocate_rows',
    summary:
      'Distributes the rows across the eight cells by hill climbing on the objective, starting from as even a split as capacity allows.',
    details: [
      'Neither the allocation nor the fidelity factor reads the seed, so the objective is exact and needs no simulation.',
      'It prices each cell’s cut probability, because is_cut’s normalised error falls as the cut rate rises. That leans the optimum toward Cas9 and pays the coverage entropy it costs. Worth 8.8 percent over a coverage-only allocation.',
      'No cell is ever emptied, but the reason is not the 1e-9 floor that design.py’s comments name. Seven occupied cells out of eight still leave the joint coverage ratio at 0.94, so one empty cell costs about 2 percent. The floor is reached only when a whole dimension collapses, and occupying every cell rules that out as a side effect. See the occupancy figure.',
    ],
    icon: 'pi pi-sliders-h',
  },
  {
    tag: '7',
    title: 'select_for_diversity',
    source: 'design.select_for_diversity',
    summary:
      'Picks each guide where its 12-mer windows are least common, greedily against a running census and round-robin across cells.',
    details: [
      'Collapsing a cell onto one coordinate buys the flat feature matrix, but it also means a cell’s guides differ in at most three positions and share most of their windows.',
      'That showed up as kmer_diversity_entropy_ratio falling from 0.97 to 0.84, a 2.4 percent haircut on the whole score through the sixth root. Levelling the multiplicities recovered it to 0.985.',
    ],
    icon: 'pi pi-share-alt',
  },
  {
    tag: '8',
    title: 'Emit',
    source: 'design.build, final ordering',
    summary:
      'Sorts strongest first, dedups, assigns experiment_id, and re-gates every row in its final position.',
    details: [
      'experiment_id is the key stage 4 merges on and the field truncate_submission dedups, so it has to be unique in exactly the array that gets sent. That is why it is assigned last, after the dedup on (cas, start, strand, guide) that stage 1 would otherwise apply silently.',
      'The array is ordered strongest first so that anything a cap cuts is the cheapest row rather than an arbitrary one.',
    ],
    icon: 'pi pi-check-circle',
  },
];

export interface OccupancyScenario {
  /** Panel title in the occupancy figure. Short: it is drawn inside an SVG. */
  label: string;
  /**
   * The eight coverage cells. Index order matches stage 5's joint support:
   * mutation is `i >> 2`, Cas system is `(i >> 1) & 1`, strand is `i & 1`.
   * True where the cell holds rows. The figure draws these as a 2 x 4 grid,
   * one row per mutation, and the fidelity factor is derived from the pattern.
   */
  cells: boolean[];
  /** Which ratio moves, so the derived number can be traced back. */
  hint: string;
}

/**
 * Three occupancy patterns, drawn side by side to separate the cheap loss from
 * the expensive one.
 *
 * Worth stating plainly, because design.py's own docstring gets it wrong and
 * the error has been copied into docs/miner_flow.md: an empty cell is not the
 * 1e-9 cliff. stage5.coverage_entropy_ratio measures entropy over the declared
 * support, so seven of eight cells still return 0.94, and stage5.geometric_mean
 * floors each ratio at 1e-9 before taking a sixth root. A ratio that actually
 * reaches zero therefore caps the factor at 1e-9^(1/6), which is 0.032, and
 * only a collapsed dimension does that: all rows on one Cas system, one strand
 * or one mutation.
 */
export const OCCUPANCY_SCENARIOS: OccupancyScenario[] = [
  {
    label: 'all eight occupied',
    cells: [true, true, true, true, true, true, true, true],
    hint: 'all four ratios 1.0',
  },
  {
    label: 'one cell empty',
    cells: [true, true, true, true, true, true, true, false],
    hint: 'joint ratio 0.94',
  },
  {
    label: 'one Cas system only',
    cells: [true, true, false, false, true, true, false, false],
    hint: 'cas ratio 0 → floored',
  },
];

export interface Factor {
  name: string;
  question: string;
  range: string;
  summary: string;
  points: string[];
}

/** The three multiplicative factors, from docs/validation_pipeline.md. */
export const FACTORS: Factor[] = [
  {
    name: 'total_weighted_score',
    question: 'How many good designs did you submit?',
    range: '0 to max_experiments × max(mutation_weight)',
    summary:
      'Sum over rows of (0.625·gc_score + 0.375·dist_score) · offtarget_factor · mutation_weight. Deterministic, so it can be maximised exactly rather than estimated.',
    points: [
      'GC pinned to 50 percent puts gc_score at 1.0.',
      'Choosing the nearest usable PAM to each mutation puts dist_score near 1.0.',
      'Pushing the 12-mer seed out of the index puts offtarget_factor at 1.0 on every row, a flat 1.43 times that nothing else in the pipeline charges for.',
    ],
  },
  {
    name: 'consistency_factor',
    question: 'Are the simulated outcomes learnable from your features?',
    range: '0 to 1',
    summary:
      '0.7·max(avg_r2, 0) + 0.3·(1 − avg_nmae) over a random forest fitted to is_cut, is_hdr and indel_length.',
    points: [
      'avg_r2 is negative for a seed-blind design, so the 0.7 term normally contributes nothing and the reachable part is the 0.3·(1 − avg_nmae) term.',
      'The 0.7 term turns positive only through a degeneracy: if no row draws no_cut then is_cut is a constant column, r2_score returns 1.0 and normalized_mae short-circuits to 0. Whether that is in reach depends on the cell type, because the cut probability rises with chromatin accessibility. See the accessibility table below.',
      'The flat feature matrix took measured r² on is_hdr and indel_length from −0.25 and −0.31 to −0.12 and −0.11.',
      'This is the only factor that varies with the round seed, so comparing two designs on a single seed is mostly noise.',
    ],
  },
  {
    name: 'distribution_fidelity_factor',
    question: 'Is the design space covered evenly and non-redundantly?',
    range: '0 to 1',
    summary:
      'Geometric mean of six coverage ratios, clipped to [0, 1]. The geometric mean is the point: one collapsed dimension drags the whole factor down in a way an arithmetic mean would hide.',
    points: [
      'The six are mutation, Cas system, strand and joint coverage entropy, plus k-mer diversity entropy and the distinct guide ratio.',
      'Every one is a veto and the sixth root makes small losses cheap but zeros fatal.',
      'There is a real tension with mutation_weights: piling rows onto the high-weight mutation raises total_weighted_score but lowers the mutation coverage ratio. The optimum is not the uniform split.',
      'Rows are scored in upload order, because stage 4 shuffles its cross-validation folds from the round seed and applies that shuffle in file order.',
    ],
  },
];

export interface Pitfall {
  violation: string;
  consequence: string;
}

/** Checked by Miner._check_invariants before every upload. */
export const PITFALLS: Pitfall[] = [
  {
    violation: 'Blank or non-string experiment_id',
    consequence: 'truncate_submission drops the row.',
  },
  {
    violation: 'Duplicate experiment_id',
    consequence:
      'Dropped. Before the dedup existed it fanned a stage 4 merge out instead, which is why very high historical scores are unlikely to be reproducible.',
  },
  {
    violation: 'Duplicate (cas, start, strand, guide)',
    consequence: 'Stage 1 keeps the first and silently discards the rest.',
  },
  {
    violation: 'More rows than max_experiments',
    consequence:
      'Everything past the cap is cut, which is why the array is ordered strongest first.',
  },
];

export interface HarnessCommand {
  command: string;
  purpose: string;
}

export const HARNESS_COMMANDS: HarnessCommand[] = [
  {
    command: 'scripts/bench_task.py --fetch',
    purpose: 'Refresh the task history from the backend.',
  },
  { command: 'scripts/bench_task.py --list', purpose: 'List the recorded task history.' },
  { command: 'scripts/bench_task.py', purpose: 'Newest task, three random seeds.' },
  {
    command: 'scripts/bench_task.py --task de19c2e0 --task-seed',
    purpose: 'Reproduce a closed round exactly, under its own recorded seed.',
  },
  {
    command: 'scripts/bench_task.py --seeds 5 --per-seed',
    purpose: 'Show the spread across seeds.',
  },
];

export interface RowField {
  field: string;
  meaning: string;
}

/**
 * The whole miner-supplied surface, from design.make_row. Its docstring calls
 * these "the entire miner-supplied surface": every biological outcome is
 * rolled by the validator afterwards.
 */
export const ROW_FIELDS: RowField[] = [
  {
    field: 'experiment_id',
    meaning:
      'Formatted exp-00000 and assigned last, in final upload order. Stage 4 merges on it and truncate_submission dedups on it.',
  },
  {
    field: 'guideRNA',
    meaning: 'The guide sequence. Length must be 20 or 23; nothing else passes stage 1.',
  },
  { field: 'target_alignment_start', meaning: 'Coordinate on chr11 where the guide aligns.' },
  {
    field: 'target_alignment_end',
    meaning: 'Must equal start plus the guide length, or stage 1 rejects the row.',
  },
  { field: 'strand', meaning: 'Either + or -.' },
  { field: 'mutation', meaning: 'Must be one of the contract’s active_mutations.' },
  { field: 'cas_system', meaning: 'Must be one of the contract’s cas_systems.' },
  {
    field: 'cell_type',
    meaning:
      'Copied verbatim from the contract. Stage 1 compares it and zeroes the row on a mismatch.',
  },
];

export interface ContractRule {
  field: string;
  gate: string;
  honoured: string;
}

/** What the contract obliges, and where the validator enforces it. */
export const CONTRACT_RULES: ContractRule[] = [
  {
    field: 'active_mutations',
    gate: 'Stage 1 rejects any row naming a mutation outside the list.',
    honoured: 'One axis of every coverage cell. Rows can only carry these names.',
  },
  {
    field: 'max_experiments',
    gate: 'truncate_submission cuts the file at this count before scoring starts.',
    honoured:
      'Sets the row target, caps the allocation, and is re-checked before upload. The array is ordered strongest first so the cut takes the cheapest rows.',
  },
  {
    field: 'max_mismatches',
    gate: 'Stage 1 counts mismatches against the reference target and rejects anything over the budget.',
    honoured:
      'The design’s only free lever. Spent on pulling GC to 50 percent, escaping the off-target index, and spelling distinct guides. A budget of 0 allows one guide per coordinate, which is what triggers the growth loop.',
  },
  {
    field: 'base_padding',
    gate: 'Stage 2 computes dist_score as exp(−distance / base_padding), and drops consistency from 1.0 to 0.3 beyond it.',
    honoured:
      'Used when ranking coordinates, and as the hard distance limit when proximity_gate is on.',
  },
  {
    field: 'cas_systems',
    gate: 'Decides the PAM motif, the off-target seed window, the cut probability base, and one of stage 5’s coverage supports.',
    honoured: 'One axis of every coverage cell, with the PAM motif chosen per system.',
  },
  {
    field: 'proximity_gate',
    gate: 'When true, stage 1 rejects a row whose coordinate is further than base_padding from the mutation.',
    honoured: 'Switches the coordinate search radius from the wide PAM flank down to base_padding.',
  },
  {
    field: 'mutation_weights',
    gate: 'Stage 2 multiplies each row’s structural score by its mutation’s weight. Stage 4 also uses it as a sample weight.',
    honoured: 'Feeds the allocation objective, pulling rows toward the higher-weighted mutation.',
  },
  {
    field: 'cell_type',
    gate: 'Stage 1 compares the row’s cell_type and zeroes it on a mismatch.',
    honoured:
      'Copied onto every row, and used to look up accessibility. See the accessibility section.',
  },
  {
    field: 'seed',
    gate: 'Stamped by the backend after the broadcast and hashed into each experiment’s RNG.',
    honoured:
      'Always 0 to a miner, and zero is dropped as a placeholder rather than used, so the design never sees a round seed.',
  },
];

/** Values the validator hard-codes, which no contract can change. */
export const FIXED_CONSTRAINTS: RowField[] = [
  {
    field: 'Guide length 20 or 23',
    meaning:
      'Nothing else passes stage 1. Length 20 is preferred because an even length can hit exactly 50 percent GC, where gc_score reaches 1.0; 23 tops out at 0.957.',
  },
  {
    field: 'k = 12, off-target flank 50 kb',
    meaning:
      'The off-target index is a 12-mer index over the gene region plus or minus 50 kb. Both are hard-coded at the validator’s call site, not contract fields.',
  },
  {
    field: 'Off-target tiers',
    meaning:
      'offtarget_factor is 1.0 for no seed hits, 0.7 for 5 or fewer, 0.4 for 20 or fewer, and 0.1 beyond that. Escaping the index entirely is a flat 1.43 times over the 0.7 tier.',
  },
  {
    field: 'Row cap and dedup',
    meaning:
      'Rows past max_experiments are cut, and stage 1 dedups on (cas_system, target_alignment_start, strand, guideRNA), keeping the first.',
  },
];
