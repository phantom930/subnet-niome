import os
from dotenv import load_dotenv

load_dotenv()

# bittensor 10.x defaults BT_NO_PARSE_CLI_ARGS=true, which skips argument
# parsing and leaves config.neuron as None. Override before importing bittensor.
os.environ.setdefault("BT_NO_PARSE_CLI_ARGS", "false")

# ---- AWS Settings -----
AWS_ACCESS_KEY_ID = os.getenv("AWS_ACCESS_KEY_ID")
AWS_SECRET_ACCESS_KEY = os.getenv("AWS_SECRET_ACCESS_KEY")
AWS_REGION = os.getenv("AWS_REGION")
AWS_S3_BUCKET = os.getenv("AWS_S3_BUCKET")


# ---- General Settings -----
TESTNET_UID = 289
MAINNET_UID = 55

FORWARD_TIMEOUT = 20


# ---- Scoring Settings -----
TOP_MINER_COUNT = 10
SCORE_DISTRIBUTION = [0.3, 0.2, 0.2, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]


# ---- Backend Request -----
BASE_URL = "https://niome-api.genomes.io"
MINER_SCORE_URL = f"{BASE_URL}/api/v3/miners/scores"
MINER_SUBMISSION_URL = f"{BASE_URL}/api/v3/miners/submissions"
TASK_URL = f"{BASE_URL}/api/v3/tasks/current"
CELL_TYPES_URL = f"{BASE_URL}/api/v3/data/cell-types?format=json"


# ---- Data -----
# Per-process instance namespace for the miner's own read/write files. Several hotkeys run as
# separate pm2 apps in one working directory; without this they overwrite each other's submission,
# task artifacts, upload record and local scoring — the *uploads* stay correct (each holds its rows
# in memory and PUTs to its own per-uid S3 key), but the on-disk diagnostics collide. Set
# NIOME_INSTANCE per process (miner.sh sets it to the hotkey name); unset keeps the flat data/ paths
# so standalone tools and single-miner runs are unchanged.
#
# chr11.fa and the bank / k-mer caches deliberately stay shared under data/: they are large,
# read-only or content-keyed (the bank filename already folds in every input, the window included),
# so sharing is correct and avoids duplicating a 130 MB reference seven times.
DATA_DIR = f"data/inst/{os.environ['NIOME_INSTANCE']}" if os.getenv("NIOME_INSTANCE") else "data"

HBB_REFERENCE_PATH = f"{DATA_DIR}/hbb_reference.json"
CONTRACT_PATH = f"{DATA_DIR}/contract.json"
CHR11_PATH = "data/chr11.fa"                                    # shared: large, read-only
MINER_SUBMISSION_PATH = f"{DATA_DIR}/submission.json"
# Where the miner records the task id, presigned URL and outcome of its most recent upload, so
# scripts/resubmit.py can retry a failed one while the URL's TTL lasts.
LAST_UPLOAD_PATH = f"{DATA_DIR}/last_upload.json"
VALID_EXPERIMENTS_PATH = f"{DATA_DIR}/valid_experiments.json"
INVALID_EXPERIMENTS_PATH = f"{DATA_DIR}/invalid_experiments.json"
STAGE3_DATASET = f"{DATA_DIR}/stage3_dataset.json"
STAGE3_SUMMARY_PATH = f"{DATA_DIR}/stage3_summary.json"
FINAL_REWARD_PATH = f"{DATA_DIR}/final_reward.json"
DISTRIBUTION_FIDELITY_PATH = f"{DATA_DIR}/distribution_fidelity_summary.json"
KMER_CACHE_DIR = "data/kmer_cache"                             # shared: content-keyed


# ---- Timeout Values -----
TASK_REQUEST_TIMEOUT = 60  # seconds
BASE_DELAY_SECONDS = 2  # seconds
SUBMISSION_TIMEOUT = 300  # seconds


# ---- Other Settings -----
MAX_TASK_RETRIES = 3
MAX_SUBMIT_RETRIES = 3

WANDB_MAX_LOGS = 60_000

SCORING_SYSTEM = "top"  # "linear", "top"
BURNING_RATE = 0.02
OWNER_HOTKEY = "5DJ5fT174AY8GzbYHnamYQCJd4cTcj2Zf7ogUvBhry1KfYVd"

BASE_BLOCK_NUMBER = 8843300
INTERVAL_BLOCKS = 720
VALIDATION_BLOCK = 600
WEIGHT_SET_BLOCK = 700

FINAL_SUBMISSION_COUNT = 5