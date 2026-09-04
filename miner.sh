#!/usr/bin/env bash
# Launch the niome miner hotkeys, each pinned to a disjoint all-HDR clean-band window.
#
# all-HDR's clean band is ~15 seeds and lands wherever NIOME_HDR_WINDOW places it. A coldkey's
# payout is 1-(1-union/900)^3 over its siblings' bands, so the win is making those bands DISJOINT:
# 3 hotkeys on disjoint windows covered 15.2% of rounds against 5.2% for one window three times
# (measured, 2.9x). Assign each hotkey its own window below and keep them from overlapping.
#
# Two ways to run:
#   ./miner.sh                 launch every hotkey as a background child, wait on all (one wrapper)
#   ./miner.sh niome_hotkey1   run just that hotkey in the foreground (exec)
#
# Prefer the second form under pm2 — one pm2 app per hotkey gives independent restart, which the
# whole-wrapper form does not (a single child crash-looping is invisible to pm2 and takes the
# hotkey down silently). Register them with, e.g.:
#   pm2 start ./miner.sh --name miner-h0 -- niome_hotkey
#   pm2 start ./miner.sh --name miner-h1 -- niome_hotkey1
#   pm2 start ./miner.sh --name miner-h2 -- niome_hotkey2

set -euo pipefail

# Absolute venv python: `pm2 restart --update-env` has rewritten PATH before and sent the miner
# into a crash loop when it resolved a different python. Never rely on PATH here.
PY=/root/workspace/subnet-niome/.venv/bin/python
ROOT=/root/workspace/subnet-niome
cd "$ROOT"

EXTERNAL_IP=184.144.255.144

# One row per hotkey:  <wallet-hotkey>  <axon.port>  <axon.external-port>  <NIOME_HDR_WINDOW>
# Windows must be disjoint (no shared seeds) or the siblings' bands overlap and the coverage
# collapses toward a single hotkey's. 200-299 / 500-599 / 800-899 are evenly spread; any three
# non-overlapping 100-seed windows in 100-999 are equivalent, since band position is otherwise free.
HOTKEYS=(
  "niome_hotkey   8091 52760 100-199"
  "niome_hotkey1  8092 52096 200-299"
  "niome_hotkey2  8093 52424 300-399"
  "niome_hotkey3  8094 52069 400-499"
  "niome_hotkey4  8095 52011 500-599"
  "niome_hotkey5  8096 52079 600-699"
  "niome_hotkey6  8097 52799 700-799"
  "niome_hotkey7  8098 52240 800-899"
  "niome_hotkey8  9001 52504 900-999"
)

# EXPERIMENTAL: NIOME_SEED_DEPEND swaps a hotkey's construction for genomics/seed_depend.py, a
# submission pinned to seed 0. It wins the rounds the backend never stamps and scores the ~0.10
# floor on every round that IS stamped, so each listed hotkey gives up its seed window entirely.
#
# The value is NOT read off the task listing: /api/v3/tasks reports `seed: 0` for rounds that were
# in fact stamped (9 list that way, only 4 were scored on it). The reliable test is miners reaching
# cons 1.000 in the score rows. By that test the rate is 2 of 57 rounds since 2026-08-25 = 3.5%.
#
# Siblings compete for the same slots: n hotkeys take ranks 1..n, not n x rank 1. So 1 hotkey is
# worth 0.035 x 0.30 = 0.0105/round against the ~0.0027 a band hotkey contributes, and 3 hotkeys
# 0.035 x 0.70 = 0.0245 against ~0.0081 — each one added is worth less than the last, and the
# fleet's band coverage falls from 9 windows to 6 (117 -> 78 seeds of 900).
#
# The VALUE is the variant index, not a flag. The build is deterministic, so without distinct
# variants siblings would submit byte-identical rows. Measured at variants 0-3: finals spanned
# 338.29-338.39 (0.10 points) while row sets shared only 28-36% of their guides.
# TEMPORARY: all nine on seed-depend. This is -14.7% on the model and taken deliberately —
# siblings 4-9 land in ranks 4..9, worth 0.290 of the curve between them and only on 3.5% of
# rounds (0.0102/round), against the six bands they replace at 0.0161/round on every round:
#   9 band / 0 sd  0.0242    6 band / 3 sd  0.0406    3 band / 6 sd  0.0406    0 band / 9 sd  0.0347
# 3 band / 6 sd is the break-even point if a hedge is wanted back. To revert, shorten this list —
# the hotkeys dropped from it return to their NIOME_HDR_WINDOW bands with no other change.
# Empty = the whole fleet runs the seed-agnostic band construction. Re-enable by listing
# "<hotkey>:<variant>" pairs here; the variant index must differ per hotkey or siblings submit
# byte-identical rows.
SEED_DEPEND_VARIANTS=""

run_one() {
  # <name> <port> <ext_port> <window>. exec so the process replaces this shell — pm2 then
  # supervises the python directly.
  local name=$1 port=$2 ext=$3 win=$4
  echo "starting $name on :$port (ext :$ext) with NIOME_HDR_WINDOW=$win"
  # NIOME_INSTANCE namespaces this hotkey's own read/write files under data/inst/<name>/ so the
  # siblings' submission, task artifacts, upload record and local scoring don't collide (settings.py).
  local sd=""
  for pair in ${SEED_DEPEND_VARIANTS:-}; do
    if [[ "${pair%%:*}" == "$name" ]]; then
      sd="${pair##*:}"
      echo "  ($name is a seed-depend hotkey, variant $sd: pinned to seed 0, no band)"
    fi
  done
  NIOME_INSTANCE="$name" NIOME_HDR_WINDOW="$win" NIOME_SEED_DEPEND="$sd" \
    exec "$PY" neurons/miner.py \
    --netuid 55 \
    --wallet niome_coldkey \
    --wallet-hotkey "$name" \
    --axon.external_ip "$EXTERNAL_IP" \
    --axon.external_port "$ext" \
    --axon.ip 0.0.0.0 \
    --axon.port "$port"
}

# Guard the windows: each must be a valid LO-HI inside 100-999, and no two may overlap. A duplicate
# string is the obvious mistake, but 600-899 vs 800-899 overlap without matching as strings, and
# 900-899 is silently rejected by the miner (falling back to the cell default and re-correlating),
# so the guard parses ranges rather than comparing text. Both cost the whole decorrelation, so they
# are launch-time errors, not warnings.
assert_disjoint_windows() {
  local -a los=() his=() names=()
  for row in "${HOTKEYS[@]}"; do
    read -r name _ _ win <<<"$row"
    if [[ ! "$win" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      echo "ERROR: $name window '$win' is not LO-HH form" >&2; exit 1
    fi
    local lo=${BASH_REMATCH[1]} hi=${BASH_REMATCH[2]}
    if (( lo < 100 || hi > 999 || lo >= hi )); then
      echo "ERROR: $name window $win is outside 100-999 or not increasing" >&2; exit 1
    fi
    local i
    for i in "${!los[@]}"; do
      if (( lo <= his[i] && los[i] <= hi )); then
        echo "ERROR: $name window $win overlaps ${names[i]} window ${los[i]}-${his[i]} — siblings would correlate" >&2
        exit 1
      fi
    done
    los+=("$lo"); his+=("$hi"); names+=("$name")
  done
}
assert_disjoint_windows

if [[ $# -ge 1 ]]; then
  # Run the single named hotkey (foreground / exec) — the per-pm2-app path.
  want=$1
  for row in "${HOTKEYS[@]}"; do
    read -r name port ext win <<<"$row"
    if [[ "$name" == "$want" ]]; then
      run_one "$name" "$port" "$ext" "$win"
    fi
  done
  echo "ERROR: hotkey '$want' is not in the HOTKEYS table" >&2
  exit 1
fi

# No argument: launch every hotkey as a background child and wait on all of them. Simple, but pm2
# sees only this wrapper — a child that dies is not individually restarted. Prefer the per-app form.
pids=()
for row in "${HOTKEYS[@]}"; do
  read -r name port ext win <<<"$row"
  ( run_one "$name" "$port" "$ext" "$win" ) &
  pids+=("$!")
done
# If any child exits, bring the rest down too rather than leaving a partial fleet running.
trap 'kill "${pids[@]}" 2>/dev/null || true' EXIT
wait -n
echo "a hotkey process exited; shutting the rest down" >&2
