#!/bin/bash
# usage: run.sh <label> <python-tree> <out.jsonl>: 3 fresh processes per case
D=$(dirname $0); PY=/home/keckj/Documents/symox/.venv/bin/python
for c in cos3_def cos2_def cos2_indef cos1_def inv2cos_final_false x_pow_ny two_pow_ny three_x_pow_ny x_pow_ny_plus3 x_pow_2y; do
  for i in 1 2 3; do
    systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$2 timeout 300 $PY $D/bench.py $1 $c >> $3
  done
done
