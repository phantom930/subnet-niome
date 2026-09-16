#!/usr/bin/env bash
# HUDEP-2: cut space 900 (whole range), band still the 3 oracle classes at width 150 / stride 30,
# group 80 / k=8 -- paired against the joined-300 cut space already measured.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) HUDEP-2 cut 900, band 300, g80 k8 ==="
CF_ONLY="33c266d5,2bb728d0,aa071e5e,9f14cd8d" CF_CUT_SPAN=900 CF_GROUP=80 CF_K=8 \
  CF_OUT=conj_fleet_hudep_cut900.json \
  .venv/bin/python conj_fleet.py > "$L/fleet_hudep900.log" 2>&1
rc=$?; echo "=== $(date -Is) exit $rc ==="
