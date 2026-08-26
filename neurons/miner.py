# The MIT License (MIT)
# Copyright © 2023 Yuma Rao
# Copyright © 2025 Genomes.io
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of
# the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


import asyncio
import hashlib
import json
import logging
import os
import requests
import shutil
import subprocess
import sys
import threading
import time
import traceback

from collections import Counter
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qs, urlparse

# Add project root to Python path. This has to happen before the niome_subnet imports below —
# running this file as a script puts neurons/ on sys.path, not the repo root — and settings has to
# be the first of them, because it fixes BT_NO_PARSE_CLI_ARGS before bittensor is imported.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import niome_subnet.utils.settings as settings  # noqa: E402

# genExp is the common generator and validator-replica layer. Most cell types mirror
# submission.py's flow directly; HEK293 reuses its context, site, validation and scoring primitives
# through a dedicated clustered builder below. Note the module chdir()s to the repo root on import —
# harmless here (every settings.py path is relative to it and neurons must run from there anyway)
# but it is an import-time side effect, so it stays below settings and above nothing that cares
# about cwd.
import genExp as G  # noqa: E402

from niome_subnet.base.miner import BaseMinerNeuron  # noqa: E402
from niome_subnet.genomics.hek293_generation import (  # noqa: E402
    generate_seed_agnostic_clustered,
)
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.model import Task  # noqa: E402
from niome_subnet.protocol import GenomicsTaskSynapse  # noqa: E402

logger = logging.getLogger(__name__)


class Miner(BaseMinerNeuron):
    """
    Miner neuron. Receives genomics tasks from validators via HTTP and processes them.

    The reply to a validator is an empty ack; the dataset travels out of band as a PUT to the
    presigned S3 URL that came with the task, and must land before that URL expires
    (``SUBMISSION_TIMEOUT``, 300 s from the moment the validator minted it). Generation itself is a
    second or so once the caches are warm, so the budget is not the constraint — see :mod:`genExp`
    for what is built and why it scores.
    """

    MAX_RETRIES = 3

    # Generation knobs, the same ones submission.py builds every task's dataset with. They are
    # class attributes rather than CLI flags because they are design decisions measured over the
    # backend's whole task history, not per-run operator settings.
    #
    # "pure" forces every row's outcome to satisfy CONSTRUCTION and drops rows that will not; the
    # "shaped" strategy instead fits outcomes to a distance ramp and takes whatever the draw gives.
    # Pure is what reaches consistency_factor 1.0.
    STRATEGY = "pure"
    # "packed" ranks each cell's sites by the stage-2 score they can actually reach (proximity is
    # pure total_weighted_score under the pure strategy); "stratified" spreads them over distance
    # bands to give stage 4 a feature axis, which pure does not need.
    SELECTION = "nearest"
    # Treat any guide this close to 50% GC as equally good so distance breaks the tie. Read only by
    # the "nearest" selection. Without the tier the ranking chases the last thousandth of gc_score
    # and reaches for sites hundreds of bp further out; gc_score is flat to first order around 50%,
    # dist_score is not.
    GC_TOLERANCE = 0.03
    # Which construction every row's simulated outcome is forced to satisfy. "hdr" makes every row
    # repair by HDR: all three stage-4 targets go constant, r2_score hits its zero/zero branch and
    # normalized_mae short-circuits on std < 1e-9, so consistency_factor reaches 1.0 without the
    # forest learning anything. It is also easier to hit than "mh" (P(HDR|variant) ~ 0.5 against
    # ~0.37 for the compound rule), so fewer near-mutation sites are abandoned to the reserve —
    # measured +2.3 on term 1 across every cas mix. The cost is identifiability: three degenerate
    # target columns and a cas-shift diagnostic of exactly zero. "mh" gives up ~0.8% to avoid that.
    # See G.CONSTRUCTIONS for the alternatives.
    CONSTRUCTION = "hdr"
    # Hold every Cas system to its own cut_p ceiling — 0.99 for Cas9, 0.96 for Cas12a. These differ
    # because stage 3's base rate does (0.86 vs 0.78) and energy is clamped at 1.0, so one shared
    # floor above 0.96 is satisfiable by Cas9 alone: it empties every Cas12a cell, takes
    # cas_system_coverage_entropy_ratio to 0, and the 1e-9 clip turns that into a 0.0316x
    # multiplier. Measured 8.81 against 289.45 for the per-system form.
    CUT_P_CEILING = True
    # Row share per Cas system, in rules.cas_systems order. Under a cut_p floor this acts as a cap
    # on Cas12a rather than a target: only ~168 Cas12a sites can reach 0.96, so whatever the quota
    # cannot fill is topped up by BACKFILL_CAS.
    CAS_MIX = "70/30"
    # Top the submission up to max_experiments from this Cas system, highest cut_p first, when a
    # filter empties a cell or the construction drops a row. Only Cas9 has spare sites at its own
    # ceiling. Leaving the cap short costs whole rows off term 1, which is a sum — measured +29 at
    # a 50/50 mix. None disables it.
    BACKFILL_CAS = "Cas9"
    # Guide variants searched per site. Each is a distinct stage-3 draw at an identical feature
    # vector, so this is the budget for hitting the construction without paying anything in stage 2.
    VARIANTS = 24
    # Site-enumeration window either side of gene_region. Wider only helps if cells run short of
    # coordinates, and costs a longer enumeration.
    FLANK = 3000
    # Guide lengths enumerated. Stage 1 accepts 20 and 23 only.
    LENGTHS = (20, 23)
    # Predict the validator's score after each build (~5 s on top of a ~1 s build). Worth it:
    # validators never report back, so this is the only signal a miner gets before the next task.
    SCORE_LOCALLY = True

    # ---- Post-upload archive and scoring ---------------------------------------------------
    # Where each round's artifacts are kept, one folder per upload. A submission can only be
    # scored *after* its own round: the contract is broadcast with seed 0 and the backend stamps
    # the real seed once the task closes, so stages 3 and 4 have nothing to run against until
    # then. Each folder is therefore scored one round later, when calc.py can fetch that seed
    # from /api/v3/tasks — which is what makes the archive worth keeping at all.
    ARCHIVE_ROOT = "data/result"
    # Copied out of data/ into the round's folder, in this order. The settings names rather than
    # literal paths, so a relocated data/ directory takes the archive with it. last_upload.json is
    # the important one: it records the task id, which is how the folder is paired back to its
    # task — the folder name alone cannot be, since it is the archive time, not the task's.
    ARCHIVED_PATHS = ("MINER_SUBMISSION_PATH", "LAST_UPLOAD_PATH", "CONTRACT_PATH",
                      "HBB_REFERENCE_PATH")
    # calc.py's summary file. Its presence is what marks a folder as already scored.
    VALIDATION_MARKER = "validation.json"
    # calc.py runs out of process on purpose: it repoints every settings path at its own output
    # directory, which in-process would leave this miner writing its next submission into the
    # archive.
    VALIDATION_SCRIPT = "calc.py"
    # It loads the 135 MB FASTA, builds a k-mer index and fits 15 random forests under KFold —
    # a couple of minutes warm. The cap is loose because it runs long after the upload, and only
    # exists so a hung subprocess cannot hold the worker thread for ever.
    VALIDATION_TIMEOUT_S = 1800

    # Per-cell-type replacements for the constants above. Everything the tuned configuration is
    # built on assumes stage 3's energy clamp: at accessibility 0.77 and above,
    # accessibility * (1.8*gc + 0.6*exp(-d/1500) + offset) exceeds 1.0 near the mutation, so every
    # Cas9 site sits at cut_p 0.99 and every Cas12a site at 0.96. HEK293's 0.35 caps energy near
    # 0.52 — the best Cas9 site reaches 0.9545 and the best Cas12a 0.8745, so the ceilings are
    # unreachable and the cut_p filter has nothing to select on. HEK293 therefore has a dedicated
    # seed-agnostic two-pool builder. These overrides keep the generic config compatible for local
    # reporting and offline tools; the dedicated builder owns its final Cas/mutation distribution.
    #
    # 59 of the backend's 258 tasks are HEK293, so this path is a quarter of all rounds.
    CELL_TYPE_OVERRIDES = {
        "HEK293": {
            "SELECTION": "packed",       # dedicated builder applies its own Cas-specific ranking
            "GC_TOLERANCE": 0.0,         # unread by packed selection
            "CONSTRUCTION": "mh",        # config compatibility only; outcomes are not forced
            "CUT_P_CEILING": False,      # unreachable at HEK293 accessibility
            "CAS_MIX": None,             # dedicated builder balances exact global Cas totals
        },
    }

    # --- seed-agnostic bank builder (non-HEK293) ------------------------------------------------
    # Where the contract carries no seed, the outcome the ordinary construction engineers against
    # seed 0 does not hold when the validator scores under the seed it is later assigned. The bank
    # builder instead picks guides whose *cut* survives every seed in a window: strict Cas9 (cut_p
    # 0.99, so all-window guides exist) plus a min-union Cas12a group (cut_p 0.96, where none do).
    # Measured on a K562 task: mean final ~41-46 across in-window seeds against ~29 for the
    # ordinary seed-0 build, with is_cut's r2 recovering to 1.0 on ~62% of seeds.
    SEED_AGNOSTIC = True
    # HEK293 keeps its own clustered builder: at accessibility 0.35 no strict Cas9 guide exists
    # (best reachable cut_p 0.9545, so 0.9545**900 ~ 6e-19), so this method has nothing to select.
    SEED_AGNOSTIC_SKIP_CELL_TYPES = ("HEK293",)
    # Seconds held back from the upload deadline. This has to cover more than the upload: the
    # builder's time budget bounds only the two scans, while the min-union selection, assembly and
    # coverage check add ~20s after them (measured: a 254s budget produced a 274s build). So the
    # reserve is ~20s post-scan + ~15s for the PUT + margin. Getting this wrong costs the round,
    # not just the hedge.
    SEED_AGNOSTIC_RESERVE_S = 75.0
    # Below this the scan cannot finish, and a truncated scan is *worse* than not trying: the rows
    # it fails to hedge get backfilled from ordinary guides, and a handful of unhedged rows destroy
    # the clean-seed property the whole method rests on.
    SEED_AGNOSTIC_MIN_BUDGET_S = 150.0
    # Unhedged rows tolerated before the ordinary construction is the safer submission. Each fails
    # ~1% (Cas9) to ~4% (Cas12a) of window seeds, so a few cost a little clean fraction while dozens
    # cost all of it.
    SEED_AGNOSTIC_MAX_BACKFILL = 8

    # The constants a contract may override. Anything absent from this tuple is global.
    TUNABLE = ("STRATEGY", "SELECTION", "GC_TOLERANCE", "CONSTRUCTION", "CUT_P_CEILING",
               "CAS_MIX", "BACKFILL_CAS", "VARIANTS", "FLANK", "LENGTHS")

    def __init__(self, config=None):
        super(Miner, self).__init__(config=config)

        self.gen_config = self.base_gen_config()

        # One build at a time, and at most one per task. Every validator broadcasts the same task id
        # with its own presigned URL, and the rows are a deterministic function of the contract, so
        # later broadcasts should reuse the first build and repeat only the upload.
        self._build_lock = asyncio.Lock()
        self._built: tuple[str, list[dict]] | None = None
        # asyncio holds only a weak reference to a running task, so a fire-and-forget
        # ``create_task`` result can be collected mid-flight and the round vanishes with no log
        # line at all. Keeping the handle here until it completes is what makes the submission
        # survive — and what lets the done-callback below report a crash that nothing awaits.
        self._inflight: set[asyncio.Task] = set()
        # Every validator broadcasts the same task, so process_task can run several times per
        # round. One archive/scoring pass at a time, and the losers skip rather than queue.
        self._archive_lock = threading.Lock()

        logger.info(
            f"Generation config: strategy={self.STRATEGY} selection={self.SELECTION} "
            f"construction={self.CONSTRUCTION} variants={self.VARIANTS} flank={self.FLANK} "
            f"lengths={self.LENGTHS} score_locally={self.SCORE_LOCALLY}"
        )
        threading.Thread(target=self._prewarm, name="niome-prewarm", daemon=True).start()

    @classmethod
    def settings_for(cls, contract: dict) -> dict:
        """The generation constants that apply to this contract, after any cell-type override."""
        overrides = cls.CELL_TYPE_OVERRIDES.get(contract.get("cell_type"), {})
        return {name: overrides.get(name, getattr(cls, name)) for name in cls.TUNABLE}

    @classmethod
    def gen_config_for(cls, contract: dict) -> G.GenConfig:
        """The fully resolved config for one contract: constants, overrides, then the knobs that
        can only be keyed once ``rules.cas_systems`` is known.

        The single entry point for anything that needs to reproduce what the miner will build —
        ``_build`` here, and the offline harnesses — so a cell-type override cannot apply on one
        side and not the other.
        """
        chosen = cls.settings_for(contract)
        cfg = G.GenConfig(
            strategy=chosen["STRATEGY"],
            selection=chosen["SELECTION"],
            construction=chosen["CONSTRUCTION"],
            gc_tolerance=chosen["GC_TOLERANCE"],
            backfill_cas=chosen["BACKFILL_CAS"],
            variants=chosen["VARIANTS"],
            flank=chosen["FLANK"],
            lengths=tuple(chosen["LENGTHS"]),
        )
        return G.config_for_contract(cfg, contract, cut_p_ceiling=chosen["CUT_P_CEILING"],
                                     cas_mix=chosen["CAS_MIX"])

    @classmethod
    def base_gen_config(cls) -> G.GenConfig:
        """The static half of the generation config, built from the class constants.

        The contract-dependent knobs — the per-Cas cut_p floors and the cas mix — are keyed by Cas
        system name, and that roster comes from ``rules.cas_systems``, so neither can be resolved
        until a task arrives; ``_build`` fills them in through ``G.config_for_contract``, the same
        call submission.py makes. Everything the prewarm needs (flank, lengths) is already here.

        Exposed as a classmethod so an offline harness can reproduce the miner's exact
        configuration without standing up a wallet and a chain connection.
        """
        return G.GenConfig(
            strategy=cls.STRATEGY,
            selection=cls.SELECTION,
            construction=cls.CONSTRUCTION,
            gc_tolerance=cls.GC_TOLERANCE,
            backfill_cas=cls.BACKFILL_CAS,
            variants=cls.VARIANTS,
            flank=cls.FLANK,
            lengths=tuple(cls.LENGTHS),
        )

    async def forward(self, body: bytes, caller_hotkey: str) -> dict:
        """
        Processes an incoming genomics task request.

        Args:
            body: Raw JSON body bytes from the validator.
            caller_hotkey: Verified hotkey ss58 of the calling validator.

        Returns:
            dict: Response payload (empty acknowledgement).
        """
        try:
            synapse = GenomicsTaskSynapse.model_validate_json(body)
            if synapse.task is None:
                logger.error(f"Task missing from {caller_hotkey}'s request; ignoring")
                return {}

            task = synapse.task
            logger.info(f"Received genomics task {task.id} from {caller_hotkey}")
            logger.info(
                f"[{task.id}] contract_url={task.contract_url} hbb_ref_url={task.hbb_ref_url}"
            )
            if not synapse.presigned_url:
                logger.error(
                    f"[{task.id}] no presigned URL in the request — there is nowhere to upload "
                    "the submission, so this round cannot be scored"
                )
                return {}
            logger.info(f"[{task.id}] upload target {self._url_summary(synapse.presigned_url)}")

            # Fire and forget - run process_task asynchronously without waiting
            handle = asyncio.create_task(self.process_task(task, synapse.presigned_url))
            self._inflight.add(handle)
            handle.add_done_callback(self._on_task_done)
            logger.info(f"[{task.id}] build+upload scheduled in the background; acking validator")

            return {}
        except Exception as e:
            logger.error(f"Forward error: {e}")
            logger.debug(traceback.format_exc())
            return {"error": str(e)}

    def _on_task_done(self, handle: asyncio.Task) -> None:
        """Drop the strong reference and surface anything ``process_task`` failed to catch."""
        self._inflight.discard(handle)
        if handle.cancelled():
            logger.error("Background submission task was cancelled before it finished")
            return
        error = handle.exception()
        if error is not None:
            logger.error(f"Background submission task crashed: {error!r}")
            logger.debug(
                "".join(traceback.format_exception(type(error), error, error.__traceback__))
            )

    async def process_task(self, task: Task, presigned_url: str) -> None:
        """Build this task's dataset and upload it to the validator's presigned URL.

        Called fire-and-forget from ``forward``, so it must neither raise into the request path nor
        block the event loop: the build and the HTTP calls are synchronous and CPU-bound, and go to
        worker threads so ``/forward`` stays answerable while a build is in flight.
        """
        if task is None or not presigned_url:
            logger.error("No task or no presigned URL — nothing to submit")
            return

        started = time.time()
        deadline = self._upload_deadline(presigned_url, started)
        tag = f"[{task.id}]"
        logger.info(
            f"{tag} step 1/5 starting submission; upload deadline in "
            f"{deadline - started:.0f}s"
        )
        # The signed URL is the only way back to this key, and the log deliberately prints it
        # without its signature. Recording it here is what makes scripts/resubmit.py possible at
        # all: a failed upload can be retried by hand while the TTL lasts, which is the only
        # recovery there is — the validator never asks twice within a task id.
        self._record_upload(task.id, presigned_url, deadline, submitted=False)

        try:
            logger.info(f"{tag} step 2/5 fetching contract and HBB reference")
            fetch_started = time.time()
            contract, reference = await asyncio.to_thread(self._fetch_artifacts, task)
            logger.info(
                f"{tag} step 2/5 done in {time.time() - fetch_started:.1f}s "
                f"| seed={contract.get('seed')} cell_type={contract.get('cell_type')} "
                f"mutations={len(contract.get('active_mutations') or [])} "
                f"rules={contract.get('rules')}"
            )

            logger.info(f"{tag} step 3/5 fetching cell-type accessibility table")
            cell_types = await asyncio.to_thread(self._fetch_cell_types)
            logger.info(f"{tag} step 3/5 done | {len(cell_types)} cell types")

            logger.info(f"{tag} step 4/5 building the dataset")
            build_started = time.time()
            async with self._build_lock:
                key = self._build_key(task.id, contract)
                if self._built is not None and self._built[0] == key:
                    rows = self._built[1]
                    logger.info(f"{tag} reusing the {len(rows)}-row build for this task")
                else:
                    rows = await asyncio.to_thread(
                        self._build, contract, reference, cell_types, deadline
                    )
                    # Only the current task's rows are worth keeping: a new task means a new
                    # contract, and the old rows can never be submitted again.
                    self._built = (key, rows)
            logger.info(
                f"{tag} step 4/5 done in {time.time() - build_started:.1f}s | {len(rows)} rows"
            )

            if not rows:
                logger.error(f"Built no rows for task {task.id} — nothing to upload")
                return

            logger.info(
                f"{tag} step 5/5 uploading {len(rows)} rows, "
                f"{deadline - time.time():.0f}s of TTL left"
            )
            await asyncio.to_thread(self._upload, presigned_url, rows, deadline)
            self._record_upload(task.id, presigned_url, deadline, submitted=True, rows=len(rows))
            logger.info(
                f"Submitted {len(rows)} rows for task {task.id} in {time.time() - started:.1f}s "
                f"({deadline - time.time():.0f}s of the URL's TTL to spare)"
            )

            # Housekeeping, not submission: score the previous round now that its seed exists,
            # then archive this one. The PUT has landed and the URL is spent, so nothing below
            # here can cost the round — but it is minutes of CPU, so it goes to a worker thread.
            await asyncio.to_thread(self._post_upload, task.id)
        except Exception as e:
            # Nothing downstream reports a miner failure — a missed upload is indistinguishable from
            # never having been contacted, and there is no retry within a task id. So this log line
            # is the only evidence the round was lost.
            logger.error(f"Failed to submit task {getattr(task, 'id', '?')}: {e}")
            logger.debug(traceback.format_exc())
            logger.error(
                f"{deadline - time.time():.0f}s of the URL's TTL remain — retry by hand with "
                "`python scripts/resubmit.py` while that is positive"
            )

    def _record_upload(
        self,
        task_id: str,
        presigned_url: str,
        deadline: float,
        submitted: bool,
        rows: int | None = None,
    ) -> None:
        """Note the current task's upload target and its outcome, for scripts/resubmit.py.

        Best-effort: a failure to write this must not be what loses the round, so it is logged and
        swallowed. The file holds a signed URL, which is a write capability on one bucket key until
        it expires — ``data/`` is gitignored, and it is worthless a few minutes later.
        """
        try:
            self._persist(
                settings.LAST_UPLOAD_PATH,
                {
                    "task_id": task_id,
                    "presigned_url": presigned_url,
                    "expires_at": deadline,
                    "submitted": submitted,
                    "rows": rows,
                },
            )
        except Exception as e:
            logger.warning(f"Could not record the upload target ({e}); manual retry will need "
                           "the URL from the validator")

    def _post_upload(self, task_id: str) -> None:
        """Score the previous round, then archive this one. Never raises.

        The order is the point. When this runs, ``data/result/<latest>`` is the *previous* round:
        its task has closed, so the backend has stamped its seed and calc.py can score it against
        the stream the validator actually used. This round's folder is written afterwards, which is
        what makes it the next run's "previous round". Scoring a folder before its task closes would
        only measure it at seed 0.
        """
        if not self._archive_lock.acquire(blocking=False):
            logger.info("post-upload: an archive pass is already running; skipping this one")
            return
        try:
            self._score_previous_round(task_id)
            self._archive_round(task_id)
        except Exception as e:
            # A failed archive costs nothing that was not already uploaded, so it is logged and
            # dropped rather than raised into process_task's handler, which reports lost rounds.
            logger.warning(f"post-upload housekeeping failed: {e}")
            logger.debug(traceback.format_exc())
        finally:
            self._archive_lock.release()

    def _latest_archive(self) -> Path | None:
        """The newest round folder. Names are timestamps, so lexicographic order is chronological."""
        root = Path(self.ARCHIVE_ROOT)
        if not root.is_dir():
            return None
        folders = sorted((path for path in root.iterdir() if path.is_dir()), key=lambda p: p.name)
        return folders[-1] if folders else None

    def _score_previous_round(self, task_id: str) -> None:
        """Run calc.py over the newest archive folder, unless it has already been scored.

        calc.py is invoked as a subprocess, pinned to the task id the folder recorded when it was
        archived: it fetches that task from ``/api/v3/tasks`` for the stamped seed, re-runs all five
        validator stages and writes the per-stage detail back into the folder.
        """
        folder = self._latest_archive()
        if folder is None:
            logger.info(f"post-upload: nothing archived under {self.ARCHIVE_ROOT} yet, so there "
                        "is nothing to score")
            return
        if (folder / self.VALIDATION_MARKER).exists():
            logger.info(f"post-upload: {folder} already holds {self.VALIDATION_MARKER}; "
                        "not scoring it again")
            return

        submission = folder / Path(settings.MINER_SUBMISSION_PATH).name
        if not submission.exists():
            logger.warning(f"post-upload: {folder} has no {submission.name} to score")
            return

        recorded = self._read_json(folder / Path(settings.LAST_UPLOAD_PATH).name) or {}
        archived_task = recorded.get("task_id")
        if archived_task and archived_task == task_id:
            # This round's own folder, from an earlier broadcast of the same task or a restart.
            # Its seed is stamped only after the task closes, so scoring it now would measure the
            # submission against a stream no validator will use. It becomes scorable next round.
            logger.info(f"post-upload: {folder} is task {task_id}'s own archive — its seed is not "
                        "stamped yet, so scoring waits for the next round")
            return

        command = [sys.executable, os.path.join(PROJECT_ROOT, self.VALIDATION_SCRIPT),
                   "--folder", str(folder), "--quiet"]
        if archived_task:
            command += ["--task-id", archived_task]
        else:
            # No recorded id: calc.py falls back to the folder name against the task's created_at,
            # then to the contract fingerprint, and refuses outright if what it finds disagrees
            # with the archived contract.
            logger.info(f"post-upload: {folder} records no task id; letting calc.py match it")

        log_path = folder / "validation.log"
        logger.info(f"post-upload: scoring {folder} with {self.VALIDATION_SCRIPT}"
                    + (f" against task {archived_task}" if archived_task else "")
                    + f"; output -> {log_path}")
        started = time.time()
        try:
            with open(log_path, "w") as handle:
                completed = subprocess.run(
                    command,
                    cwd=PROJECT_ROOT,          # every settings path is relative to the repo root
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    timeout=self.VALIDATION_TIMEOUT_S,
                    check=False,
                )
        except subprocess.TimeoutExpired:
            logger.error(f"post-upload: {self.VALIDATION_SCRIPT} exceeded "
                         f"{self.VALIDATION_TIMEOUT_S}s on {folder}; see {log_path}")
            return

        elapsed = time.time() - started
        if completed.returncode != 0:
            logger.error(f"post-upload: {self.VALIDATION_SCRIPT} exited {completed.returncode} "
                         f"after {elapsed:.0f}s — see {log_path}")
            return

        summary = self._read_json(folder / self.VALIDATION_MARKER) or {}
        score = summary.get("score") or {}
        rows = summary.get("rows") or {}
        logger.info(
            f"post-upload: scored {folder} in {elapsed:.0f}s | task {summary.get('task_id')} "
            f"seed {summary.get('scored_seed')} | {rows.get('valid')}/{rows.get('scored')} rows "
            f"valid | final_score {self._number(score.get('final_score'))} = weighted "
            f"{self._number(score.get('total_weighted_score'))} x consistency "
            f"{self._number(score.get('consistency_factor'))} x fidelity "
            f"{self._number(score.get('distribution_fidelity_factor'))}"
        )

    def _archive_round(self, task_id: str) -> Path | None:
        """Copy this round's artifacts into a new ``data/result/<archived at>`` folder.

        Everything needed to re-score the upload later, and nothing else: the rows exactly as sent,
        the broadcast contract and reference they were designed against, and the upload record that
        names the task. Skipped when the newest folder already records this task id, so repeated
        broadcasts of one task archive once.
        """
        latest = self._latest_archive()
        if latest is not None:
            recorded = self._read_json(latest / Path(settings.LAST_UPLOAD_PATH).name) or {}
            if recorded.get("task_id") == task_id:
                logger.info(f"post-upload: task {task_id} is already archived in {latest}")
                return None

        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
        folder = Path(self.ARCHIVE_ROOT) / stamp
        folder.mkdir(parents=True, exist_ok=True)

        copied, missing = [], []
        for name in self.ARCHIVED_PATHS:
            source = Path(getattr(settings, name))
            if not source.exists():
                missing.append(str(source))
                continue
            shutil.copyfile(source, folder / source.name)
            copied.append(source.name)

        if copied:
            logger.info(f"post-upload: archived {', '.join(copied)} for task {task_id} -> {folder}")
        if missing:
            # A missing last_upload.json costs the folder its task pairing; a missing submission
            # makes it unscorable altogether. Either way the round itself is already uploaded.
            logger.warning(f"post-upload: not in data/, so not archived: {', '.join(missing)}")
        return folder

    @staticmethod
    def _read_json(path: Path | str):
        """Read a JSON file, or None if it is absent or unreadable. Used on archive metadata,
        where a half-written or hand-edited file must not take the miner down."""
        try:
            with open(path) as handle:
                return json.load(handle)
        except (OSError, json.JSONDecodeError):
            return None

    @staticmethod
    def _number(value, digits: int = 4) -> str:
        """Format a score that calc.py may not have produced."""
        return f"{value:.{digits}f}" if isinstance(value, (int, float)) else "?"

    def _prewarm(self) -> None:
        """Load the reference — and, if a previous task left its artifacts behind, the k-mer index
        and PAM sites too — before the first task arrives.

        All three caches are task-independent (every task issued so far shares one gene_region and
        one rules block), so this is pure critical-path removal: a warm process only pays for the
        build itself.
        """
        started = time.time()
        try:
            logger.info("Prewarm: loading chr11 reference sequence")
            self._load_sequence()
            logger.info(f"Prewarm: chr11 loaded in {time.time() - started:.1f}s")
            if os.path.exists(settings.CONTRACT_PATH) and os.path.exists(
                settings.HBB_REFERENCE_PATH
            ):
                logger.info("Prewarm: replaying the last task's artifacts to warm the site cache")
                with open(settings.CONTRACT_PATH) as handle:
                    contract = json.load(handle)
                with open(settings.HBB_REFERENCE_PATH) as handle:
                    reference = json.load(handle)
                # cell_types is left empty on purpose: it is carried on the Context but nothing
                # cached here reads it, so warming does not need the backend to answer.
                context = G.build_context(contract, reference, {})
                sites = G.enumerate_sites(context, self.gen_config.flank, self._lengths())
                logger.info(
                    f"Prewarmed chr11 and {len(sites)} PAM sites in {time.time() - started:.1f}s"
                )
            else:
                logger.info(
                    f"Prewarmed chr11 in {time.time() - started:.1f}s; PAM sites will be "
                    "enumerated on the first task"
                )
        except Exception as e:
            logger.error(f"Prewarm failed ({e}); the first task will pay the cold start")
            logger.debug(traceback.format_exc())

    def _seed_agnostic_applies(self, contract: dict, budget_s: float) -> tuple[bool, str]:
        """Whether the bank builder should run for this contract, and why not when it should not."""
        if not self.SEED_AGNOSTIC:
            return False, "disabled by SEED_AGNOSTIC"
        if contract.get("cell_type") in self.SEED_AGNOSTIC_SKIP_CELL_TYPES:
            return False, f"{contract.get('cell_type')} has its own clustered builder"
        if contract.get("seed"):
            # A real seed makes the ordinary construction exact — consistency_factor 1.0 and a
            # score around 290, against ~41 for a hedge against an unknown seed. Never trade that
            # away; the hedge exists only for the case the seed is not known at build time.
            return False, f"contract carries seed {contract['seed']}; the matched build scores far higher"
        if SA.MT._gpu() is None:
            # CPU-only the scan takes ~540s against a 300s TTL, so it would always truncate.
            return False, "no usable GPU; the scan cannot finish inside the upload window"
        if budget_s < self.SEED_AGNOSTIC_MIN_BUDGET_S:
            return False, f"only {budget_s:.0f}s of budget; need {self.SEED_AGNOSTIC_MIN_BUDGET_S:.0f}s"
        return True, ""

    def _build_seed_agnostic(self, contract: dict, reference: dict, cell_types: dict,
                             budget_s: float) -> list[dict] | None:
        """Build from the strict-Cas9 / min-union-Cas12a banks. None means fall back."""
        task = {"content": {"contract": contract, "hbb_reference": reference}}
        cfg = replace(
            SA.SeedAgnosticConfig(),
            time_budget_s=budget_s,
            max_experiments=contract["rules"].get("max_experiments"),
        )
        logger.info(
            f"Build: seed-agnostic bank scan | window {cfg.start_seed}-{cfg.end_seed} "
            f"cas_mix={cfg.cas_mix} budget={budget_s:.0f}s gpu=on"
        )
        started = time.time()
        rows, meta = SA.build_submission(contract, reference, cell_types, task, cfg)
        logger.info(
            f"Build: seed-agnostic done in {time.time() - started:.0f}s | rows={meta['rows']} "
            f"hedged={meta['hedged_rows']} backfilled={meta['backfilled']} "
            f"strict_cas9={meta['strict_cas9_available']} cas12a_pool={meta['cas12a_pool']} "
            f"group={meta['cas12a_group']} union={meta['group_failed_seed_union']} "
            f"clean~{meta['clean_fraction_estimate']:.1%} cas={meta['cas_mix']}"
        )

        expected = cfg.max_experiments or 0
        if expected and meta["rows"] < expected:
            logger.warning(
                f"Build: seed-agnostic produced {meta['rows']}/{expected} rows; "
                "falling back to the ordinary construction"
            )
            return None
        if meta["backfilled"] > self.SEED_AGNOSTIC_MAX_BACKFILL:
            # The hedge is all-or-nothing: measured, 66 unhedged rows took the clean fraction to 0%
            # and the score below the ordinary build. Prefer the build we can reason about.
            logger.warning(
                f"Build: {meta['backfilled']} rows outside the hedge exceeds "
                f"{self.SEED_AGNOSTIC_MAX_BACKFILL}; falling back to the ordinary construction"
            )
            return None
        return rows

    def _build(self, contract: dict, reference: dict, cell_types: dict,
               deadline: float | None = None) -> list[dict]:
        """Generate this task's submission and log what it should be worth.

        Most cell types follow the same sequence as ``submission.build_for_task``. HEK293 instead
        uses the dedicated clustered builder because its provisional seed can be replaced before
        scoring; that builder does not inspect simulated outcomes while choosing designs.
        """
        logger.info("Build: loading chr11 (cached after the first call)")
        self._load_sequence()

        budget = (deadline - time.time() - self.SEED_AGNOSTIC_RESERVE_S
                  if deadline else float(self.SEED_AGNOSTIC_MIN_BUDGET_S))
        applies, why_not = self._seed_agnostic_applies(contract, budget)
        if applies:
            try:
                rows = self._build_seed_agnostic(contract, reference, cell_types, budget)
            except Exception as exc:
                logger.warning(f"Build: seed-agnostic builder failed ({exc}); using the "
                               "ordinary construction")
                logger.debug(traceback.format_exc())
                rows = None
            if rows:
                self._persist(settings.MINER_SUBMISSION_PATH, rows)
                logger.info(f"Build: wrote the submission to {settings.MINER_SUBMISSION_PATH}")
                return rows
        else:
            logger.info(f"Build: seed-agnostic path not used ({why_not})")

        if not contract.get("seed"):
            if contract.get("cell_type") == "HEK293":
                logger.info(
                    "Contract carries no seed — HEK293 row selection is seed-agnostic; only the "
                    "local diagnostic outcome score is provisional"
                )
            else:
                # Stage 3 is seeded from contract.seed, so an unstamped task's construction stops
                # holding the moment the backend assigns a real one. The rows stay valid;
                # consistency_factor does not survive.
                logger.warning(
                    "Contract carries no seed — the outcome construction is provisional and the "
                    "predicted consistency will not hold if the backend restamps the task"
                )

        logger.info("Build: building context (k-mer index over gene_region +/- 50 kb)")
        context = G.build_context(contract, reference, cell_types)
        logger.info(
            f"Build: context ready | mutations={len(context.mutations)} "
            f"cas={context.cas_systems} max_experiments={context.max_experiments}"
        )

        # Resolve the per-contract settings before enumeration so a future cell-type override of
        # FLANK or LENGTHS affects the site pool as well as the downstream generator.
        cfg = self.gen_config_for(contract)
        lengths = tuple(sorted(set(cfg.lengths)))
        logger.info(f"Build: enumerating PAM sites in gene_region +/- {cfg.flank}")
        sites = G.enumerate_sites(context, cfg.flank, lengths)
        logger.info(f"Build: {len(sites)} PAM sites enumerated")

        override = self.CELL_TYPE_OVERRIDES.get(contract.get("cell_type"))
        if override:
            logger.info(
                f"Build: cell type {contract.get('cell_type')} overrides {sorted(override)} — "
                f"its accessibility cannot reach the cut_p ceilings the tuned config assumes"
            )
        logger.info(
            f"Build: config resolved | construction={cfg.construction} selection={cfg.selection} "
            f"gc_tolerance={cfg.gc_tolerance} cut_p_floors={cfg.cut_p_floors} "
            f"cas_share={cfg.cas_share} backfill_cas={cfg.backfill_cas}"
        )
        is_hek293 = contract.get("cell_type") == "HEK293"
        if is_hek293:
            logger.info(
                "Build: generating HEK293 rows from independent seed-agnostic Cas pools "
                "(adaptive fallback plus typed x170 candidate; exact 70/30 Cas and 50/50 strand "
                "margins; Cas-specific length/GC; strict TTTV; global guide uniqueness/entropy; "
                "anonymous seed-free candidate gate)"
            )
            clustered = generate_seed_agnostic_clustered(context, sites, cfg)
            rows, valid, results = clustered.rows, clustered.valid, clustered.results
            cfg = replace(cfg, weight_skew=clustered.weight_skew)
            logger.info(
                f"Build: HEK293 clusters ready | site_clusters={clustered.site_clusters} "
                f"feature_clusters={clustered.feature_clusters} cas={clustered.cas_counts} "
                f"strands={clustered.strand_counts} mutations={clustered.mutation_totals} "
                f"quota_objective={clustered.deterministic_objective:.3f} "
                f"(tws={clustered.deterministic_tws:.3f}, "
                f"fidelity={clustered.deterministic_fidelity:.6f}, "
                f"candidates={clustered.quota_candidates_evaluated}) "
                f"selected={clustered.selected_candidate}"
            )
            gate = clustered.gate_diagnostics
            if gate is not None:
                logger.info(
                    f"Build: HEK293 anonymous gate | estimator={gate.estimator} "
                    f"replicates={gate.replicates} folds={gate.folds} "
                    f"mean_C(fallback={gate.fallback_mean_consistency}, "
                    f"typed={gate.typed_mean_consistency}) "
                    f"mean_final(fallback={gate.fallback_mean_final}, "
                    f"typed={gate.typed_mean_final}) "
                    f"paired_delta_mean={gate.paired_final_delta_mean} "
                    f"se={gate.paired_final_delta_standard_error} "
                    f"lcb={gate.paired_final_delta_lcb} "
                    f"threshold={gate.minimum_lcb_gain} choice={gate.selected_candidate}"
                )
                if gate.typed_candidate_failure:
                    logger.warning(
                        "Build: HEK293 typed candidate unavailable; retained adaptive fallback: "
                        f"{gate.typed_candidate_failure}"
                    )
            for failure in clustered.quota_candidate_failures:
                logger.warning(f"Build: skipped infeasible HEK293 quota candidate: {failure}")
        else:
            if cfg.strategy == "pure":
                # The optimal skew depends on this contract's mutation-weight ratio, which moves
                # task to task; fitting it costs a few ms of surrogate scoring against selection.
                cfg = replace(cfg, weight_skew=G.choose_weight_skew(context, sites, cfg))
                logger.info(f"Build: fitted weight_skew={cfg.weight_skew}")

            logger.info(
                f"Build: generating rows (construction={cfg.construction}, "
                f"{cfg.variants} variants/site)"
            )
            rows, valid, results = G.generate(context, sites, cfg)
            rows = G.order_rows(rows, valid)
        logger.info(f"Build: generated {len(rows)} rows, {len(valid)} of them stage-1 valid")

        strategy_label = "hek293_clustered" if is_hek293 else cfg.strategy
        construction_label = "seed-agnostic" if is_hek293 else cfg.construction
        logger.info(
            f"Built {len(rows)}/{context.max_experiments} rows from {len(sites)} sites "
            f"| strategy={strategy_label} construction={construction_label} "
            f"skew={cfg.weight_skew} cas={dict(Counter(r['cas'] for r in results))} "
            f"outcomes={dict(Counter(r['outcome'] for r in results))}"
        )
        if len(rows) < context.max_experiments:
            # Term 1 is a sum, so a short submission is a straight loss. The backfill should have
            # prevented this; reaching here means the eligible pool is genuinely exhausted.
            logger.warning(
                f"Submission is {context.max_experiments - len(rows)} row(s) under the cap even "
                f"after the backfill — widen FLANK, or relax CUT_P_CEILING / CAS_MIX"
            )
        # The normal generator prints these too, but only to stdout; a miner's evidence is its log.
        # HEK293 deliberately does not force a simulated construction, so only apply the common
        # ID/design checks there. Stage-1 validity is guaranteed by the aligned `valid` row count.
        invariant_results = [] if is_hek293 else results
        for problem in G.check_invariants(rows, invariant_results, cfg, valid):
            logger.warning(f"Submission invariant violated: {problem}")

        if self.SCORE_LOCALLY and len(valid) >= 2:
            logger.info("Build: scoring locally through the validator's own pipeline")
            report = G.score_rows(valid, results, context)
            logger.info(
                f"Predicted score {report['final_score']:.3f} = "
                f"weighted {report['total_weighted_score']:.3f} "
                f"x consistency {report['consistency_factor']:.4f} "
                f"x fidelity {report['distribution_fidelity_factor']:.4f}"
            )

        # Keep a copy of exactly what was sent. This is also the path benchmark_submission reads, so
        # an operator can re-score the upload against the persisted contract and reference.
        self._persist(settings.MINER_SUBMISSION_PATH, rows)
        logger.info(f"Build: wrote the submission to {settings.MINER_SUBMISSION_PATH}")
        return rows

    def _lengths(self) -> tuple[int, ...]:
        """Guide lengths in the form the site cache is keyed on — sorted and deduped, so a prewarm
        and a build hit the same entry."""
        return tuple(sorted(set(self.gen_config.lengths)))

    @staticmethod
    def _load_sequence() -> None:
        """Load chr11, turning a missing FASTA into an ordinary error.

        ``G.load_sequence`` raises SystemExit when ``data/chr11.fa`` is absent — the right reflex for
        a CLI, wrong for a long-lived miner, where it would escape both the prewarm thread's and
        ``process_task``'s ``except Exception`` and take the round down without a log line.
        """
        try:
            G.load_sequence()
        except SystemExit as missing:
            raise RuntimeError(str(missing)) from missing

    def _fetch_artifacts(self, task: Task) -> tuple[dict, dict]:
        """GET the contract and the HBB reference.

        Plain unsigned GETs — the presigning is already in the URL. Both URLs are short-lived, so
        they are fetched per task and never cached; the *contents* are persisted, which is what lets
        a restart prewarm the site cache and an operator re-score a submission locally.
        """
        contract = self._get_json(task.contract_url, "contract")
        reference = self._get_json(task.hbb_ref_url, "hbb_reference")
        self._persist(settings.CONTRACT_PATH, contract)
        self._persist(settings.HBB_REFERENCE_PATH, reference)
        logger.info(
            f"Persisted contract to {settings.CONTRACT_PATH} and reference to "
            f"{settings.HBB_REFERENCE_PATH}"
        )
        return contract, reference

    def _fetch_cell_types(self) -> dict:
        """The accessibility table, read unsigned from the backend.

        Accessibility scales stage 3's energy, which drives the cut and repair draws, so falling back
        to the default 1.0 risks simulating different outcomes than the validator and breaking the
        construction. How much it costs depends on the clamp: energy saturates at 1.0 for
        near-mutation rows at 50% GC, so a real 0.77 was measured to change nothing, while 0.5 drops
        below the clamp and takes consistency_factor from 1.0 to 0.18. Rows stay valid either way —
        it is term 2 of the score that is at risk.
        """
        try:
            logger.info(f"Fetching cell types from {settings.CELL_TYPES_URL}")
            response = requests.get(settings.CELL_TYPES_URL, timeout=settings.TASK_REQUEST_TIMEOUT)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.warning(
                f"Cell-types fetch failed ({e}); falling back to accessibility 1.0, which will "
                "mis-predict stage 3 and likely cost the whole consistency factor"
            )
            return {}

    def _get_json(self, url: str, what: str) -> dict:
        last_error = None
        for attempt in range(self.MAX_RETRIES):
            try:
                logger.info(f"GET {what} ({self._url_summary(url)})")
                response = requests.get(url, timeout=settings.TASK_REQUEST_TIMEOUT)
                response.raise_for_status()
                logger.info(
                    f"GET {what} -> {response.status_code}, {len(response.content) / 1024:.1f} KB"
                )
                return response.json()
            except Exception as e:
                last_error = e
                logger.warning(
                    f"Fetching {what} failed (attempt {attempt + 1}/{self.MAX_RETRIES}): {e}"
                )
                if attempt + 1 < self.MAX_RETRIES:
                    time.sleep(settings.BASE_DELAY_SECONDS * (attempt + 1))
        raise RuntimeError(f"could not fetch {what}: {last_error}")

    def _upload(self, presigned_url: str, rows: list[dict], deadline: float) -> None:
        """PUT the bare JSON array to the validator's bucket.

        The URL must be sent exactly as received — a re-encoded query string is a
        ``SignatureDoesNotMatch`` — and so must the header set, which is *not* free to vary: see
        ``_upload_headers`` for why sending a Content-Type can break the signature outright.
        """
        payload = json.dumps(rows).encode()
        headers = self._upload_headers(presigned_url)
        last_error = None
        logger.info(
            f"PUT {len(payload) / 1024:.1f} KB to {self._url_summary(presigned_url)} "
            f"with headers {headers or '{}'}"
        )

        for attempt in range(self.MAX_RETRIES):
            remaining = deadline - time.time()
            if remaining <= 0:
                raise RuntimeError(
                    "presigned URL expired before the upload went through "
                    f"(last error: {last_error})"
                )
            try:
                logger.info(
                    f"Upload attempt {attempt + 1}/{self.MAX_RETRIES}, {remaining:.0f}s of TTL left"
                )
                response = requests.put(
                    presigned_url,
                    data=payload,
                    headers=headers,
                    timeout=min(remaining, 60),
                )
                response.raise_for_status()
                logger.info(
                    f"Uploaded {len(rows)} rows ({len(payload) / 1024:.1f} KB) -> "
                    f"{response.status_code}, etag {response.headers.get('ETag', '?')}"
                )
                return
            except Exception as e:
                last_error = e
                # S3 answers with XML, and the body is what separates an expired URL from an altered
                # one — the status is 403 either way.
                body = getattr(getattr(e, "response", None), "text", "") or ""
                logger.warning(
                    f"Upload attempt {attempt + 1}/{self.MAX_RETRIES} failed: {e} {body[:300]}"
                )
                if attempt + 1 < self.MAX_RETRIES:
                    time.sleep(settings.BASE_DELAY_SECONDS * (attempt + 1))

        raise RuntimeError(f"upload failed after {self.MAX_RETRIES} attempts: {last_error}")

    @staticmethod
    def _upload_headers(presigned_url: str) -> dict:
        """The headers this PUT may carry, decided by which signing scheme minted the URL.

        The validator presigns ``put_object`` with no ``ContentType``, and botocore answers that
        with a **SigV2** URL (``AWSAccessKeyId`` / ``Signature`` / ``Expires``). SigV2's
        string-to-sign is ``VERB\\n Content-MD5\\n Content-Type\\n Expires\\n resource`` — the
        Content-Type is *in the signature*, signed as the empty string. Sending
        ``Content-Type: application/json`` therefore makes S3 hash a different string than the
        validator did and reject the upload with ``SignatureDoesNotMatch``, which reads like a
        credentials problem and is in fact this header. So V2 URLs get no headers at all.

        SigV4 URLs (``X-Amz-SignedHeaders``) cover only the headers they name, so there a
        Content-Type is required exactly when it was signed — and forbidden otherwise, for the
        same reason. ``requests`` adds no Content-Type of its own for a bytes body, so an empty
        dict really does send none.
        """
        try:
            query = parse_qs(urlparse(presigned_url).query)
        except Exception:
            return {}

        signed = query.get("X-Amz-SignedHeaders")
        if signed:
            names = [name.strip().lower() for name in signed[0].split(";")]
            return {"Content-Type": "application/json"} if "content-type" in names else {}
        return {}

    @staticmethod
    def _upload_deadline(presigned_url: str, now: float) -> float:
        """When S3 stops accepting the PUT.

        The URL carries its own answer, which beats assuming it was minted the instant it arrived:
        ``Expires`` (an absolute epoch) on the SigV2 URLs the validator currently mints, and
        ``X-Amz-Date`` + ``X-Amz-Expires`` on SigV4 ones. Either comes off the *validator's* clock,
        so the result is clamped to our own SUBMISSION_TIMEOUT budget: a slow local clock must not
        talk us into a deadline that has in fact already passed.
        """
        fallback = now + settings.SUBMISSION_TIMEOUT
        try:
            query = parse_qs(urlparse(presigned_url).query)
            if "Expires" in query:
                expires_at = float(query["Expires"][0])
                scheme = "SigV2"
            else:
                signed_at = datetime.strptime(query["X-Amz-Date"][0], "%Y%m%dT%H%M%SZ").replace(
                    tzinfo=timezone.utc
                )
                expires_at = signed_at.timestamp() + int(query["X-Amz-Expires"][0])
                scheme = "SigV4"
            logger.info(
                f"{scheme} presigned URL expires at "
                f"{datetime.fromtimestamp(expires_at, timezone.utc).isoformat()}; "
                f"{expires_at - now:.0f}s of it left on our clock"
            )
            return min(expires_at, fallback)
        except Exception as e:
            logger.warning(
                f"Could not read the URL's expiry ({e}); assuming the local "
                f"SUBMISSION_TIMEOUT of {settings.SUBMISSION_TIMEOUT}s"
            )
            return fallback

    @staticmethod
    def _url_summary(url: str) -> str:
        """host + path of a presigned URL. The query carries the signature, so it stays out of the
        log — the bucket and key are what identify a submission, and the expiry is logged
        separately by ``_upload_deadline``."""
        try:
            parsed = urlparse(url)
            return f"{parsed.scheme}://{parsed.netloc}{parsed.path}"
        except Exception:
            return "<unparseable url>"

    @staticmethod
    def _build_key(task_id: str, contract: dict) -> str:
        """Memo key. The contract is hashed in as well, so a contract that changes under a task id
        rebuilds instead of re-uploading rows designed against the old rules."""
        canonical = json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()
        return f"{task_id}:{hashlib.sha256(canonical).hexdigest()[:16]}"

    @staticmethod
    def _persist(path: str, document) -> None:
        """Write via a temp file and rename: the prewarm on the next restart reads these back, and a
        task arriving mid-write must not leave half a JSON document behind."""
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        temporary = f"{path}.tmp"
        with open(temporary, "w") as handle:
            json.dump(document, handle)
        os.replace(temporary, path)

    async def blacklist(self, caller_hotkey: str) -> bool:
        """
        Determines whether an incoming request should be blacklisted.

        Args:
            caller_hotkey: ss58 hotkey of the caller (already verified by http_auth).

        Returns:
            bool: True if the request should be rejected.
        """
        if caller_hotkey not in self.metagraph.hotkeys:
            if not self.config.blacklist.allow_non_registered:
                logger.debug(f"Blacklisting un-registered hotkey {caller_hotkey}")
                return True

        uid = self.metagraph.hotkeys.index(caller_hotkey)

        if self.config.blacklist.force_validator_permit:
            if not self.metagraph.neurons[uid].validator_permit:
                logger.warning(f"Blacklisting non-validator hotkey {caller_hotkey}")
                return True

        logger.info(f"Allowing recognized hotkey {caller_hotkey} (uid {uid})")
        return False


# This is the main function, which runs the miner.
if __name__ == "__main__":
    with Miner() as miner:
        while True:
            logger.info(f"Miner running... {time.time()}")
            time.sleep(5)
