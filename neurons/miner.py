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
import sys
import threading
import time
import traceback

from datetime import datetime, timezone
from urllib.parse import parse_qs, urlparse

# Add project root to Python path. This must happen before the niome_subnet imports, because running
# this file as a script puts neurons/ on sys.path rather than the repo root.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import niome_subnet.utils.settings as settings  # noqa: E402

from niome_subnet.base.miner import BaseMinerNeuron  # noqa: E402
from niome_subnet.genomics import design  # noqa: E402
from niome_subnet.genomics.model import Task  # noqa: E402
from niome_subnet.protocol import GenomicsTaskSynapse  # noqa: E402

logger = logging.getLogger(__name__)


class Miner(BaseMinerNeuron):
    """
    Miner neuron. Receives genomics tasks from validators via HTTP and processes them.

    The reply to a validator is an empty acknowledgement. The dataset itself travels out of band as
    a PUT to the presigned S3 URL that arrived with the task, and has to land before that URL
    expires — ``SUBMISSION_TIMEOUT`` is 300 s from the moment the validator minted it, and a miner
    that misses the window is indistinguishable from one that was never contacted. There is no
    retry within a task id and no feedback channel, which is why this class scores its own build
    locally before uploading and logs every step.

    What gets built, and why it scores, is :mod:`niome_subnet.genomics.design`.
    """

    MAX_RETRIES = 3

    # Score the build through the validator's own stages 4 and 5 before uploading (~25 s for ten
    # seeds). Validators never report back, so this is the only signal available before the next
    # task. Skipped automatically when the URL's remaining TTL cannot afford it.
    SCORE_LOCALLY = True
    LOCAL_SCORE_SEEDS = 6

    def __init__(self, config=None):
        super(Miner, self).__init__(config=config)

        self.generation_config = design.Config()

        # One build at a time, and at most one per contract. Every validator broadcasts the same
        # task with its own presigned URL, and the rows are a deterministic function of the
        # contract, so later broadcasts reuse the first build and repeat only the upload.
        self._build_lock = asyncio.Lock()
        self._built: tuple[str, list[dict]] | None = None
        # asyncio keeps only a weak reference to a running task, so a fire-and-forget create_task
        # can be garbage-collected mid-flight and the round would vanish without a log line.
        # Holding the handle until it completes is what makes the submission survive.
        self._inflight: set[asyncio.Task] = set()

        logger.info(
            "Generation config: guide_lengths=%s pam_search_flank=%d coordinates_per_cell=%d "
            "guides_per_coordinate=%d score_locally=%s",
            self.generation_config.guide_lengths,
            self.generation_config.pam_search_flank,
            self.generation_config.coordinates_per_cell,
            self.generation_config.guides_per_coordinate,
            self.SCORE_LOCALLY,
        )
        threading.Thread(target=self._prewarm, name="niome-prewarm", daemon=True).start()

    async def forward(self, request_body: bytes, caller_hotkey: str) -> dict:
        """
        Processes an incoming genomics task request.

        Args:
            request_body: Raw JSON body bytes from the validator.
            caller_hotkey: Verified hotkey ss58 of the calling validator.

        Returns:
            dict: Response payload (empty acknowledgement).
        """
        try:
            synapse = GenomicsTaskSynapse.model_validate_json(request_body)
            if synapse.task is None:
                logger.error(f"Task missing from {caller_hotkey}'s request; ignoring")
                return {}

            task = synapse.task
            logger.info(f"Received genomics task {task.id} from {caller_hotkey}")
            if not synapse.presigned_url:
                logger.error(
                    f"[{task.id}] no presigned URL in the request — there is nowhere to upload "
                    "the submission, so this round cannot be scored"
                )
                return {}
            logger.info(f"[{task.id}] upload target {self._url_summary(synapse.presigned_url)}")

            # Fire and forget - run process_task asynchronously without waiting
            background_task = asyncio.create_task(self.process_task(task, synapse.presigned_url))
            self._inflight.add(background_task)
            background_task.add_done_callback(self._on_task_done)
            logger.info(f"[{task.id}] build and upload scheduled; acking the validator now")

            return {}
        except Exception as error:
            logger.error(f"Forward error: {error}")
            logger.debug(traceback.format_exc())
            return {"error": str(error)}

    def _on_task_done(self, background_task: asyncio.Task) -> None:
        """Drop the strong reference and surface anything ``process_task`` failed to catch."""
        self._inflight.discard(background_task)
        if background_task.cancelled():
            logger.error("Background submission task was cancelled before it finished")
            return
        error = background_task.exception()
        if error is not None:
            logger.error(f"Background submission task crashed: {error!r}")
            logger.debug(
                "".join(traceback.format_exception(type(error), error, error.__traceback__))
            )

    async def process_task(self, task: Task, presigned_url: str) -> None:
        """Build this task's dataset and upload it to the validator's presigned URL.

        Called fire-and-forget from ``forward``, so it must neither raise into the request path nor
        block the event loop: the build and the HTTP calls are synchronous and CPU-bound and go to
        worker threads, keeping ``/forward`` answerable while a build is in flight.
        """
        if task is None or not presigned_url:
            logger.error("No task or no presigned URL — nothing to submit")
            return

        started = time.time()
        deadline = self._upload_deadline(presigned_url, started)
        task_tag = f"[{task.id}]"
        logger.info(f"{task_tag} step 1/5 starting; upload deadline in {deadline - started:.0f}s")
        # Recording the target is what makes a manual retry possible: a failed upload can be
        # repeated by hand while the TTL lasts, which is the only recovery there is.
        self._record_upload(task.id, presigned_url, deadline, submitted=False)

        try:
            logger.info(f"{task_tag} step 2/5 fetching the contract and HBB reference")
            contract, reference = await asyncio.to_thread(self._fetch_artifacts, task)
            logger.info(
                f"{task_tag} step 2/5 done | seed={contract.get('seed')} "
                f"cell_type={contract.get('cell_type')} "
                f"mutations={len(contract.get('active_mutations') or [])} "
                f"rules={contract.get('rules')}"
            )

            logger.info(f"{task_tag} step 3/5 fetching the cell-type accessibility table")
            cell_types = await asyncio.to_thread(self._fetch_cell_types)
            logger.info(f"{task_tag} step 3/5 done | {len(cell_types)} cell types")

            logger.info(f"{task_tag} step 4/5 building the dataset")
            build_started = time.time()
            async with self._build_lock:
                memo_key = self._build_key(task.id, contract)
                if self._built is not None and self._built[0] == memo_key:
                    rows = self._built[1]
                    logger.info(f"{task_tag} reusing the {len(rows)}-row build for this contract")
                else:
                    rows = await asyncio.to_thread(self._build, contract, reference, cell_types,
                                                   deadline)
                    self._built = (memo_key, rows)
            logger.info(
                f"{task_tag} step 4/5 done in {time.time() - build_started:.1f}s | {len(rows)} rows"
            )

            if not rows:
                logger.error(f"{task_tag} built no rows — nothing to upload")
                return

            logger.info(
                f"{task_tag} step 5/5 uploading {len(rows)} rows, "
                f"{deadline - time.time():.0f}s of TTL left"
            )
            await asyncio.to_thread(self._upload, presigned_url, rows, deadline)
            self._record_upload(task.id, presigned_url, deadline, submitted=True, rows=len(rows))
            logger.info(
                f"Submitted {len(rows)} rows for task {task.id} in {time.time() - started:.1f}s "
                f"({deadline - time.time():.0f}s of the URL's TTL to spare)"
            )
        except Exception as error:
            # Nothing downstream reports a miner failure, so this log line is the only evidence
            # the round was lost.
            logger.error(f"Failed to submit task {getattr(task, 'id', '?')}: {error}")
            logger.debug(traceback.format_exc())
            logger.error(
                f"{deadline - time.time():.0f}s of the URL's TTL remain — a manual PUT to the "
                f"URL recorded in {settings.MINER_LAST_UPLOAD_PATH} still counts while that is positive"
            )

    # -----------------------------------------------------------------------------------------
    # Build
    # -----------------------------------------------------------------------------------------

    def _build(self, contract: dict, reference: dict, cell_types: dict,
               deadline: float) -> list[dict]:
        """Generate this task's submission and log what it should be worth."""
        context = design.build_context(contract, reference, cell_types)
        logger.info(
            "Build: context ready | mutations=%d cas=%s max_experiments=%d "
            "max_mismatches=%d base_padding=%d",
            len(context.mutations), context.cas_systems, context.max_experiments,
            context.max_mismatches, context.base_padding,
        )
        if not context.seeds():
            # Expected in production: the backend stamps the round seed after the broadcast, so a
            # miner never designs against it. Worth stating once per task, because it is why the
            # design optimises the two seed-independent factors and treats the third statistically.
            logger.info(
                "Build: contract carries no round seed (the backend stamps it before validation), "
                "so outcomes cannot be designed for — see genomics/design.py"
            )

        rows, entries, diagnostics = design.build(context, self.generation_config)
        if not rows:
            logger.error(f"Build produced no rows: {diagnostics.get('error', 'unknown reason')}")
            return []

        logger.info(
            "Built %d/%d rows | term1=%.2f offtarget=%s distinct_feature_vectors=%d",
            diagnostics["rows"], diagnostics["rows_wanted"],
            diagnostics["total_weighted_score"], diagnostics["offtarget_factors"],
            diagnostics["distinct_feature_vectors"],
        )
        logger.info(f"Allocation: {diagnostics['allocation']}")
        if diagnostics["empty_cells"]:
            # Stage 5's geometric mean turns an unoccupied (mutation, cas, strand) cell into a
            # ~1e-9 multiplier on the entire score, so this is close to a lost round.
            logger.error(
                f"{len(diagnostics['empty_cells'])} coverage cell(s) could not be filled: "
                f"{diagnostics['empty_cells']} — widen Config.flank"
            )

        for problem in self._check_invariants(rows, contract):
            logger.warning(f"Submission invariant violated: {problem}")

        # Persist before scoring: the rows are the deliverable, and this file plus the contract and
        # reference beside it are a complete record of what was sent. To re-score it with the
        # validator's own benchmark_submission, copy all three into data/ — deliberately a manual
        # step, because writing them there directly is what would clobber a co-located validator.
        self._persist(settings.MINER_OUTPUT_PATH, rows)
        logger.info(f"Build: wrote the submission to {settings.MINER_OUTPUT_PATH}")

        self._score_locally(rows, context, deadline)
        return rows

    def _score_locally(self, rows: list[dict], context, deadline: float) -> None:
        """Predict the score, if the presigned URL's TTL can afford the forest fits."""
        if not self.SCORE_LOCALLY or len(rows) < 2:
            return
        # Roughly 2.5 s per seed for 15 RandomForest fits, plus headroom for the upload itself.
        remaining_ttl = deadline - time.time()
        if remaining_ttl < self.LOCAL_SCORE_SEEDS * 3 + 45:
            logger.info(
                f"Skipping the local score: {remaining_ttl:.0f}s of TTL left, "
                "which belongs to the upload"
            )
            return
        try:
            seeds = context.seeds() or list(design.SEED_SUPPORT[::150])[:self.LOCAL_SCORE_SEEDS]
            report = design.score_rows(rows, context, seeds=seeds)
            basis = "under the contract's seed" if context.seeds() else \
                   f"averaged over {len(seeds)} sampled seeds"
            logger.info(
                "Predicted score %.3f %s = weighted %.3f x consistency %.4f x fidelity %.4f "
                "(cut rate %.4f)",
                report["final_score"], basis, report["total_weighted_score"],
                report["consistency_factor"], report["distribution_fidelity_factor"],
                report["cut_rate"],
            )
        except Exception as error:
            # A failed prediction must never cost the upload.
            logger.warning(f"Local scoring failed ({error}); uploading anyway")
            logger.debug(traceback.format_exc())

    @staticmethod
    def _check_invariants(rows: list[dict], contract: dict) -> list[str]:
        """The three rules that silently cost rows rather than raising.

        ``truncate_submission`` drops a blank or repeated ``experiment_id`` and anything past
        ``max_experiments``; stage 1 drops a repeated (cas, start, strand, guide). None of that is
        reported anywhere, so a violation shows up only as an unexplained gap between the local
        prediction and what a validator pays.
        """
        problems = []
        experiment_ids = [row.get("experiment_id") for row in rows]
        if any(
            not isinstance(experiment_id, str) or not experiment_id.strip()
            for experiment_id in experiment_ids
        ):
            problems.append("blank or non-string experiment_id")
        if len(set(experiment_ids)) != len(experiment_ids):
            duplicate_count = len(experiment_ids) - len(set(experiment_ids))
            problems.append(f"{duplicate_count} duplicate experiment_id")
        design_keys = [
            (row["cas_system"], row["target_alignment_start"], row["strand"], row["guideRNA"])
            for row in rows
        ]
        if len(set(design_keys)) != len(design_keys):
            duplicate_count = len(design_keys) - len(set(design_keys))
            problems.append(f"{duplicate_count} duplicate (cas, start, strand, guide)")
        row_cap = contract.get("rules", {}).get("max_experiments")
        if row_cap is not None and len(rows) > row_cap:
            problems.append(f"{len(rows)} rows exceeds max_experiments {row_cap}")
        return problems

    def _prewarm(self) -> None:
        """Load the 135 MB reference — and, if a previous task left its artifacts behind, the k-mer
        index and PAM enumeration too — before the first task arrives.

        All three caches are task-independent, so this is pure critical-path removal: a warm
        process pays for nothing but the build itself.
        """
        started = time.time()
        try:
            logger.info("Prewarm: loading the chr11 reference sequence")
            design.load_sequence()
            logger.info(f"Prewarm: chr11 loaded in {time.time() - started:.1f}s")
            if not (os.path.exists(settings.MINER_CONTRACT_PATH)
                    and os.path.exists(settings.MINER_HBB_REFERENCE_PATH)):
                logger.info("Prewarm: no previous task on disk; PAMs will be enumerated on demand")
                return
            with open(settings.MINER_CONTRACT_PATH) as contract_file:
                contract = json.load(contract_file)
            with open(settings.MINER_HBB_REFERENCE_PATH) as reference_file:
                reference = json.load(reference_file)
            # cell_types is deliberately empty: nothing cached here reads it, so warming does not
            # need the backend to answer.
            context = design.build_context(contract, reference, {})
            by_cas_strand = design.enumerate_coordinates(context, self.generation_config)
            logger.info(
                f"Prewarmed chr11, the k-mer index and "
                f"{sum(len(picked) for picked in by_cas_strand.values())} PAM coordinates "
                f"in {time.time() - started:.1f}s"
            )
        except Exception as error:
            logger.error(f"Prewarm failed ({error}); the first task will pay the cold start")
            logger.debug(traceback.format_exc())

    # -----------------------------------------------------------------------------------------
    # Fetching
    # -----------------------------------------------------------------------------------------

    def _fetch_artifacts(self, task: Task) -> tuple[dict, dict]:
        """GET the contract and the HBB reference.

        Plain unsigned GETs — the presigning is already in the URL. The URLs are short-lived so
        they are fetched per task, but the *contents* are persisted, which is what lets a restart
        prewarm the PAM cache and an operator re-score a submission offline.
        """
        contract = self._get_json(task.contract_url, "contract")
        reference = self._get_json(task.hbb_ref_url, "hbb_reference")
        self._persist(settings.MINER_CONTRACT_PATH, contract)
        self._persist(settings.MINER_HBB_REFERENCE_PATH, reference)
        return contract, reference

    def _fetch_cell_types(self) -> dict:
        """The accessibility table, read unsigned from the backend.

        Accessibility is the largest single term in stage 3's energy, which sets the cut
        probability and the repair mix, so a wrong value makes the local prediction describe a
        different simulation than the validator's. The rows stay valid either way — this is the
        prediction's accuracy at risk, and the allocation's, which prices each cell's cut
        probability. A stale table on disk is a better guess than the 1.0 default.
        """
        try:
            response = requests.get(settings.CELL_TYPES_URL,
                                    timeout=settings.TASK_REQUEST_TIMEOUT)
            response.raise_for_status()
            cell_types = response.json()
            self._persist(settings.MINER_CELL_TYPES_PATH, cell_types)
            return cell_types
        except Exception as error:
            logger.warning(
                f"Cell-types fetch failed ({error}); falling back to the last known table"
            )
            try:
                with open(settings.MINER_CELL_TYPES_PATH) as cached_file:
                    return json.load(cached_file)
            except Exception:
                logger.warning(
                    "No cached cell types either; stage 2 will default accessibility to 1.0, so "
                    "the local prediction and the cut-rate pricing will both be optimistic"
                )
                return {}

    def _get_json(self, url: str, artifact_name: str) -> dict:
        last_error = None
        for attempt in range(self.MAX_RETRIES):
            try:
                response = requests.get(url, timeout=settings.TASK_REQUEST_TIMEOUT)
                response.raise_for_status()
                logger.info(
                    f"GET {artifact_name} -> {response.status_code}, "
                    f"{len(response.content) / 1024:.1f} KB"
                )
                return response.json()
            except Exception as error:
                last_error = error
                logger.warning(
                    f"Fetching {artifact_name} failed "
                    f"(attempt {attempt + 1}/{self.MAX_RETRIES}): {error}"
                )
                if attempt + 1 < self.MAX_RETRIES:
                    time.sleep(settings.BASE_DELAY_SECONDS * (attempt + 1))
        raise RuntimeError(f"could not fetch {artifact_name}: {last_error}")

    # -----------------------------------------------------------------------------------------
    # Upload
    # -----------------------------------------------------------------------------------------

    def _upload(self, presigned_url: str, rows: list[dict], deadline: float) -> None:
        """PUT the bare JSON array to the validator's bucket.

        The URL must be sent exactly as received — a re-encoded query string is a
        ``SignatureDoesNotMatch`` — and so must the header set; see ``_upload_headers``.
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
                    f"presigned URL expired before the upload went through "
                    f"(last error: {last_error})"
                )
            try:
                response = requests.put(
                    presigned_url, data=payload, headers=headers, timeout=min(remaining, 60)
                )
                response.raise_for_status()
                logger.info(
                    f"Uploaded {len(rows)} rows ({len(payload) / 1024:.1f} KB) -> "
                    f"{response.status_code}, etag {response.headers.get('ETag', '?')}"
                )
                return
            except Exception as error:
                last_error = error
                # S3 answers with XML, and the body is what separates an expired URL from an
                # altered one — the status is 403 either way.
                error_body = getattr(getattr(error, "response", None), "text", "") or ""
                logger.warning(
                    f"Upload attempt {attempt + 1}/{self.MAX_RETRIES} failed: "
                    f"{error} {error_body[:300]}"
                )
                if attempt + 1 < self.MAX_RETRIES:
                    time.sleep(settings.BASE_DELAY_SECONDS * (attempt + 1))

        raise RuntimeError(f"upload failed after {self.MAX_RETRIES} attempts: {last_error}")

    @staticmethod
    def _upload_headers(presigned_url: str) -> dict:
        """The headers this PUT may carry, decided by which signing scheme minted the URL.

        The validator presigns ``put_object`` with no ``ContentType``, and botocore answers that
        with a SigV2 URL (``AWSAccessKeyId``/``Signature``/``Expires``). SigV2's string-to-sign is
        ``VERB\\n Content-MD5\\n Content-Type\\n Expires\\n resource`` — the Content-Type is *in*
        the signature, signed as the empty string. Sending ``Content-Type: application/json``
        therefore makes S3 hash a different string than the validator did and reject the upload
        with ``SignatureDoesNotMatch``, which reads like a credentials problem and is in fact this
        header. So V2 URLs get no headers at all.

        SigV4 URLs cover only the headers they name in ``X-Amz-SignedHeaders``, so there a
        Content-Type is required exactly when it was signed and forbidden otherwise, for the same
        reason. ``requests`` adds none of its own for a bytes body, so an empty dict sends none.
        """
        try:
            url_query = parse_qs(urlparse(presigned_url).query)
        except Exception:
            return {}
        signed_headers = url_query.get("X-Amz-SignedHeaders")
        if signed_headers:
            signed_names = [name.strip().lower() for name in signed_headers[0].split(";")]
            return {"Content-Type": "application/json"} if "content-type" in signed_names else {}
        return {}

    @staticmethod
    def _upload_deadline(presigned_url: str, now: float) -> float:
        """When S3 stops accepting the PUT.

        The URL carries its own answer, which beats assuming it was minted the instant it arrived:
        ``Expires`` (an absolute epoch) on a SigV2 URL, ``X-Amz-Date`` + ``X-Amz-Expires`` on a
        SigV4 one. Either comes off the *validator's* clock, so the result is clamped to our own
        ``SUBMISSION_TIMEOUT`` budget — a slow local clock must not talk us into a deadline that
        has in fact already passed.
        """
        fallback_deadline = now + settings.SUBMISSION_TIMEOUT
        try:
            url_query = parse_qs(urlparse(presigned_url).query)
            if "Expires" in url_query:
                expires_at = float(url_query["Expires"][0])
                scheme = "SigV2"
            else:
                signed_at = datetime.strptime(
                    url_query["X-Amz-Date"][0], "%Y%m%dT%H%M%SZ"
                ).replace(tzinfo=timezone.utc)
                expires_at = signed_at.timestamp() + int(url_query["X-Amz-Expires"][0])
                scheme = "SigV4"
            logger.info(
                f"{scheme} presigned URL expires at "
                f"{datetime.fromtimestamp(expires_at, timezone.utc).isoformat()}; "
                f"{expires_at - now:.0f}s of it left on our clock"
            )
            return min(expires_at, fallback_deadline)
        except Exception as error:
            logger.warning(
                f"Could not read the URL's expiry ({error}); assuming the local "
                f"SUBMISSION_TIMEOUT of {settings.SUBMISSION_TIMEOUT}s"
            )
            return fallback_deadline

    def _record_upload(self, task_id: str, presigned_url: str, deadline: float,
                       submitted: bool, rows: int | None = None) -> None:
        """Note the current task's upload target and its outcome, so a lost round can be retried.

        Best-effort: a failure to write this must not be what loses the round. The file holds a
        signed URL, which is a write capability on one bucket key until it expires — ``data/`` is
        gitignored and the URL is worthless a few minutes later.
        """
        try:
            self._persist(settings.MINER_LAST_UPLOAD_PATH, {
                "task_id": task_id,
                "presigned_url": presigned_url,
                "expires_at": deadline,
                "submitted": submitted,
                "rows": rows,
            })
        except Exception as error:
            logger.warning(f"Could not record the upload target ({error})")

    # -----------------------------------------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------------------------------------

    @staticmethod
    def _url_summary(url: str) -> str:
        """host + path of a presigned URL. The query carries the signature, so it stays out of the
        log; the bucket and key are what identify a submission."""
        try:
            parsed_url = urlparse(url)
            return f"{parsed_url.scheme}://{parsed_url.netloc}{parsed_url.path}"
        except Exception:
            return "<unparseable url>"

    @staticmethod
    def _build_key(task_id: str, contract: dict) -> str:
        """Memo key. The contract is hashed in as well, so a contract that changes under one task
        id rebuilds instead of re-uploading rows designed against the old rules."""
        canonical = json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()
        return f"{task_id}:{hashlib.sha256(canonical).hexdigest()[:16]}"

    @staticmethod
    def _persist(path: str, document) -> None:
        """Write via a temp file and rename: these are read back on the next restart, and a task
        arriving mid-write must not leave half a JSON document behind."""
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        temporary_path = f"{path}.tmp"
        with open(temporary_path, "w") as output_file:
            json.dump(document, output_file)
        os.replace(temporary_path, path)

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

        logger.debug(f"Allowing recognized hotkey {caller_hotkey}")
        return False


# This is the main function, which runs the miner.
if __name__ == "__main__":
    with Miner() as miner:
        while True:
            logger.info(f"Miner running... {time.time()}")
            time.sleep(5)
