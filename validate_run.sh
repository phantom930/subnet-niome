#!/bin/bash
# top-2 configs per cell from the payout ranking, plus the shipped (100:80) baseline
declare -A PAIRS=(
  [K562]="75:80,30:80,100:80"
  [HUDEP-2]="150:80,75:80,100:80"
  [CD34+_HSPC]="150:100,225:80,100:80"
  [HEK293]="75:80,75:100,100:80"
)
for CELL in K562 HUDEP-2 CD34+_HSPC HEK293; do
  for T in $(.venv/bin/python -c "import json;print(' '.join(json.load(open('validate_tasks.json'))['$CELL']))"); do
    echo "########## $CELL  $T ##########"
    SAFE=$(echo "$CELL" | tr -d '+' | tr '-' '_')
    CG_CELL="$CELL" CG_TASK="$T" CG_PAIRS="${PAIRS[$CELL]}" \
      CG_OUT="validate_${SAFE}_${T:0:8}.json" .venv/bin/python -u cell_grid.py
  done
done
