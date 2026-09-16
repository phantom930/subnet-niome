#!/usr/bin/env bash
# Paired re-run of the same 10 tasks at group 80 / k=8 on every cell type.
# Writes conj_fleet_g80k8.json -- NOT the tracked conj_fleet.json. Do not edit while running.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) conj_fleet g80 k8, 10 tasks ==="
CF_TASKS=10 CF_GROUP=80 CF_K=8 CF_OUT=conj_fleet_g80k8.json \
  .venv/bin/python conj_fleet.py > "$L/fleet_g80k8.log" 2>&1
rc=$?; echo "=== $(date -Is) exit $rc ==="
