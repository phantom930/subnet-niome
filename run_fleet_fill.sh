#!/usr/bin/env bash
# Fill the three tasks that fell off the moving window, at group 80 / k=8, to complete the pairing.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) fill 800d04c4, aa071e5e, 9f14cd8d at g80 k8 ==="
CF_ONLY="800d04c4,aa071e5e,9f14cd8d" CF_GROUP=80 CF_K=8 CF_OUT=conj_fleet_g80k8_fill.json \
  .venv/bin/python conj_fleet.py > "$L/fleet_fill.log" 2>&1
rc=$?; echo "=== $(date -Is) exit $rc ==="
