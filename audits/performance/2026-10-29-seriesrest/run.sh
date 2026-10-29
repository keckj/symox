#!/bin/bash
# usage: run.sh <python dir> <out file>; 5 fresh processes per call
# (median) and the in-process median of 7 repeated calls.
D=$(cd "$(dirname "$0")" && pwd)
PY=/home/keckj/Documents/symox/.venv/bin/python
: > "$2"
for i in $(seq 0 9); do
  for k in 1 2 3 4 5; do
    timeout 120 env PYTHONPATH="$1" $PY "$D/bench.py" fresh $i
  done | sort -t$'\t' -k2,2g | sed -n 3p | sed 's/^/fresh\t/' >> "$2"
  timeout 300 env PYTHONPATH="$1" $PY "$D/bench.py" repeat $i | sed 's/^/repeat\t/' >> "$2"
done
