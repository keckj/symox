#!/bin/bash
# usage: run.sh <after-python-dir> <before-python-dir> <outdir>
# 3 fresh processes per mode under a 6G address-space cap, from /tmp.
AFTER=$1; BEFORE=$2; O=$3; mkdir -p $O; PY=$HOME/symox-ci/venv/bin/python
H=$(dirname $(readlink -f $0))/bench.py
CAP="prlimit --as=6000000000"
for run in 1 2 3; do
  (cd /tmp && timeout 600 $CAP $PY -I $H sympy 30) >> $O/sympy.jsonl
  (cd /tmp && timeout 600 $CAP env PYTHONPATH=$BEFORE $PY $H symox 200) >> $O/before.jsonl
  (cd /tmp && timeout 600 $CAP env PYTHONPATH=$AFTER $PY $H symox 200) >> $O/after.jsonl
done
