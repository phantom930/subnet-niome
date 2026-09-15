#!/usr/bin/env bash
# Oracle fleet-coverage run, 40 HEK293 tasks. Do NOT edit while running: bash reads a script
# incrementally by byte offset, so an edit mid-run corrupts the tail.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) oracle 40 tasks ==="
.venv/bin/python conj_oracle.py --tasks 40 > "$L/oracle40.log" 2>&1
rc=$?; echo "=== $(date -Is) exit $rc ==="
