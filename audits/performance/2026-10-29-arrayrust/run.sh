#!/bin/bash
# usage: run.sh <worktree> <outdir>; 3 fresh processes per case and mode.
W=$1; O=$2; mkdir -p $O; PY=/home/keckj/Documents/symox/.venv/bin/python
H=$(dirname $(readlink -f $0))/bench.py
CAP="systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0"
for run in 1 2 3; do for i in 0 1 2 3 4 5 6 7 8 9; do
  (cd /tmp && timeout 300 $CAP $PY $H sympy x $i 20) >> $O/sympy.jsonl
  for mode in python native; do
    (cd $W && timeout 300 $CAP env PYTHONPATH=python $PY $H symox $mode $i 200) >> $O/symox-$mode.jsonl
  done
done; done
