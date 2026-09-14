#!/usr/bin/env bash
# HUDEP-2 to 12 contracts. The first 6 re-run identically (seeded rng, same order) which
# doubles as a reproducibility check against the backed-up n=6 result.
set -u
cd /root/workspace/subnet-niome
LOG=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) HUDEP-2 n=12 ==="
CR_CELL=HUDEP-2 CR_N=12 CR_NS=12 .venv/bin/python conj_replicate.py > "$LOG/rep12_HUDEP_2.log" 2>&1
rc=$?; echo "=== $(date -Is) HUDEP-2 n=12 exit $rc ==="
