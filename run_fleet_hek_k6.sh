#!/usr/bin/env bash
# group 80 / k=6 (the SHIPPED HEK293 arm) on the three HEK293 tasks of the pinned 10.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) HEK293 group 80 k=6 ==="
CF_ONLY="920e1837,fe564247,1ab688d3" CF_GROUP=80 CF_K=6 CF_OUT=conj_fleet_hek_g80k6.json \
  .venv/bin/python conj_fleet.py > "$L/fleet_hek_k6.log" 2>&1
rc=$?; echo "=== $(date -Is) exit $rc ==="
