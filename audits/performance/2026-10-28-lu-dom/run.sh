#!/bin/bash
# usage: run.sh WORKTREE BASE_SO NEW_SO OUT.jsonl ; runs each case 5x per binary
set -e
W=$1; BASE=$2; NEW=$3; OUT=$4; PY=/home/keckj/Documents/symox/.venv/bin/python
H=$(dirname "$(readlink -f "$0")")/bench.py
SAVE=$(mktemp -d)
cp "$W/python/symox/_core.abi3.so" "$SAVE/"
: > "$OUT"
CASES=$(cd "$W" && PYTHONPATH=python $PY "$H" list)
for rep in 1 2 3 4 5; do
  for label in base new; do
    if [ $label = base ]; then cp "$BASE" "$W/python/symox/_core.abi3.so"; else cp "$NEW" "$W/python/symox/_core.abi3.so"; fi
    for c in $CASES; do
      for mode in first repeat; do
        (cd "$W" && PYTHONPATH=python timeout 300 $PY "$H" $c $mode) | sed "s/^{/{\"binary\": \"$label\", \"rep\": $rep, /" >> "$OUT"
      done
    done
  done
done
cp "$SAVE/_core.abi3.so" "$W/python/symox/_core.abi3.so"
