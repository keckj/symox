#!/bin/bash
# usage: run.sh <label> <module> <PYTHONPATH or -> <reps> >> raw.txt
# Every line is one fresh process: label index first repeated result-hash integrand.
label=$1; mod=$2; pp=$3; reps=${4:-5}
here=$(cd "$(dirname "$0")" && pwd)
PY=/home/keckj/Documents/symox/.venv/bin/python
cd /tmp
for r in $(seq 1 $reps); do
  for i in $(seq 0 9); do
    if [ "$pp" = "-" ]; then out=$(systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 timeout 120 $PY "$here/bench.py" $mod $i)
    else out=$(systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$pp timeout 120 $PY "$here/bench.py" $mod $i); fi
    echo "$label $out"
  done
done
