#!/usr/bin/env bash
# 10-hotkey oracle fleet, latest 10 three-seed tasks, scored against each task's own field.
# Do not edit while running: bash reads a script incrementally by byte offset.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) conj_fleet 10 tasks ==="
CF_TASKS=10 .venv/bin/python conj_fleet.py > "$L/fleet10.log" 2>&1
rc=$?; echo "=== $(date -Is) exit $rc ==="
