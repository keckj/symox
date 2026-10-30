#!/bin/bash
# run.sh <label> <PYTHONPATH or "sympy"> <outdir>: the ten calls of bench.py
# under the 6G memory cap, three fresh processes per library.
set -u
label=$1; pp=$2; out=$3; here=$(cd "$(dirname "$0")" && pwd)
PY=/home/keckj/Documents/symox/.venv/bin/python
lib=symox; [ "$pp" = sympy ] && lib=sympy && pp=/nonexistent
mkdir -p "$out"; cd /tmp
for run in 1 2 3; do
  timeout 900 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
    env PYTHONPATH=$pp $PY $here/bench.py $lib "$out/$label-$run.json" > /dev/null
done
