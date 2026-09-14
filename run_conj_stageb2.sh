#!/usr/bin/env bash
# Re-run HUDEP-2 (task now pinned) and fill the high-k / group-80 gap on the cells already priced.
set -u
cd /root/workspace/subnet-niome
LOG=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
while pgrep -f "conj_stageb.py" >/dev/null || pgrep -f "conj_grid.py" >/dev/null; do sleep 30; done
echo "=== $(date -Is) HUDEP-2 stage B (pinned task) ==="
CSB_CELL=HUDEP-2 CSB_TOP=12 CSB_NS=20 .venv/bin/python conj_stageb.py > "$LOG/stageb_HUDEP_2.log" 2>&1
rc=$?; echo "=== $(date -Is) HUDEP-2 exit $rc ==="
# k=10 at group 80 was never priced: the w x fid shortlist picked g125 as its k=10 representative.
FORCE="10:80:30:6,10:80:75:6,10:80:150:6,10:80:225:6,9:80:30:6,9:80:150:6"
for cell in K562 CD34+_HSPC HUDEP-2; do
  tag="$(echo "$cell" | tr -d '+' | tr '-' '_')"
  echo "=== $(date -Is) high-k fill $cell ==="
  CSB_CELL="$cell" CSB_NS=20 CSB_FORCE="$FORCE" .venv/bin/python conj_stageb.py > "$LOG/stagebk_$tag.log" 2>&1
  rc=$?; echo "=== $(date -Is) high-k fill $cell exit $rc ==="
done
echo "=== $(date -Is) STAGE B FOLLOW-UP DONE ==="
