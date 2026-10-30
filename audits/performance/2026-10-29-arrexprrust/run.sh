#!/bin/bash
# usage: run.sh <after-python-dir> <before-python-dir> <outdir>
# 3 fresh processes per case and mode; SymPy 1.14 from /tmp.
AFTER=$1; BEFORE=$2; O=$3; mkdir -p $O; PY=/home/keckj/Documents/symox/.venv/bin/python
H=$(dirname $(readlink -f $0))/bench.py
CAP="systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0"
for run in 1 2 3; do for i in $(seq 0 11); do
  (cd /tmp && timeout 300 $CAP $PY $H sympy $i 20) >> $O/sympy.jsonl
  (cd /tmp && timeout 300 $CAP env PYTHONPATH=$BEFORE $PY $H symox $i 200) >> $O/before.jsonl
  (cd /tmp && timeout 300 $CAP env PYTHONPATH=$AFTER $PY $H symox $i 200) >> $O/after.jsonl
done; done
