#!/bin/bash
# run.sh <label> <PYTHONPATH or "sympy"> <outdir>
# first-call: each case in a fresh process, 5 runs; repeat: one process, 3
# runs of median-of-20. All under the 6G memory cap.
set -u
label=$1; pp=$2; out=$3; here=$(cd "$(dirname "$0")" && pwd)
PY=/home/keckj/Documents/symox/.venv/bin/python
lib=symox; [ "$pp" = sympy ] && lib=sympy && pp=/nonexistent
mkdir -p "$out"; cd /tmp
cap="systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0"
for run in 1 2 3 4 5; do
  for i in 0 1 2 3 4 5 6 7 8 9; do
    timeout 300 $cap env PYTHONPATH=$pp $PY $here/bench.py $lib first 0 $i
  done
done > "$out/$label-first.jsonl"
for run in 1 2 3; do
  timeout 900 $cap env PYTHONPATH=$pp $PY $here/bench.py $lib repeat 20
done > "$out/$label-repeat.jsonl"
