#!/bin/bash
# usage: run.sh <label> <module> <PYTHONPATH or -> <reps> > raw-<label>.txt
# Every line is one fresh process: label index seconds result-hash integrand.
label=$1; mod=$2; pp=$3; reps=${4:-7}
here=$(cd "$(dirname "$0")" && pwd)
PY=/home/keckj/Documents/symox/.venv/bin/python
cd /tmp
for r in $(seq 1 $reps); do
  for i in $(seq 0 9); do
    if [ "$pp" = "-" ]; then out=$(timeout 120 $PY "$here/bench.py" $mod $i)
    else out=$(PYTHONPATH=$pp timeout 120 $PY "$here/bench.py" $mod $i); fi
    echo "$label $out"
  done
done
