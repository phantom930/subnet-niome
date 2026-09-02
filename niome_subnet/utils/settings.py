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
HBB_REFERENCE_PATH = "data/hbb_reference.json"
CONTRACT_PATH = "data/contract.json"
CHR11_PATH = "data/chr11.fa"
MINER_SUBMISSION_PATH = "data/submission.json"
VALID_EXPERIMENTS_PATH = "data/valid_experiments.json"
INVALID_EXPERIMENTS_PATH = "data/invalid_experiments.json"
STAGE3_DATASET = "data/stage3_dataset.json"
STAGE3_SUMMARY_PATH = "data/stage3_summary.json"
FINAL_REWARD_PATH = "data/final_reward.json"
DISTRIBUTION_FIDELITY_PATH = "data/distribution_fidelity_summary.json"
KMER_CACHE_DIR = "data/kmer_cache"


# ---- Miner workspace -----
# The miner keeps its own artifacts out of data/, which belongs to the validation pipeline: every
# stage there communicates through these fixed filenames, and truncate_submission rewrites
# data/submission.json in place. A miner writing data/contract.json would silently replace the
# contract a co-located validator is scoring against — the miner receives an unstamped contract
# (seed 0), so the validator would then score every submission under seed 0.
#
# Only the reference genome is shared, read-only: it is ~130 MB, identical for both roles, and
# downloaded once by scripts/run_validator.sh. Point MINER_GENOME_PATH elsewhere to keep it
# outside data/ as well.
MINER_DIR = os.getenv("NIOME_MINER_DIR", "miner_data")
MINER_CONTRACT_PATH = f"{MINER_DIR}/contract.json"
MINER_HBB_REFERENCE_PATH = f"{MINER_DIR}/hbb_reference.json"
MINER_CELL_TYPES_PATH = f"{MINER_DIR}/cell_types.json"
MINER_LAST_UPLOAD_PATH = f"{MINER_DIR}/last_upload.json"
MINER_OUTPUT_PATH = f"{MINER_DIR}/submission.json"
MINER_GENOME_PATH = os.getenv("NIOME_GENOME_PATH", CHR11_PATH)


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