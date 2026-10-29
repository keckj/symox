#!/bin/bash
# usage: run.sh <label> <pythonpath|sympy> <outfile>
# Each case runs in 3 fresh processes (first call + median of 5 repeats).
L=$1; PP=$2; OUT=$3; D=$(cd "$(dirname "$0")" && pwd)
for c in shift_a unshift_a unshift_b_build unshift_errors meijer_unshift plan_mixed plan_errors plan_apply reduce_order meijerg_floats hyperexpand_plan; do
  for k in 1 2 3; do
    if [ "$PP" = sympy ]; then
      (cd /tmp && timeout 300 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 /home/keckj/Documents/symox/.venv/bin/python "$D/bench.py" sympy $c 5 | sed "s/^/$L /")
    else
      (cd /tmp && timeout 300 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$PP /home/keckj/Documents/symox/.venv/bin/python "$D/bench.py" symox $c 5 | sed "s/^/$L /")
    fi || echo "$L {\"case\": \"$c\", \"timeout\": true}"
  done
done >> "$OUT"
