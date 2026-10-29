#!/bin/bash
# usage: bench.sh <label> <pythonpath|sympy> <runs> <kind:n>...
label=$1; pp=$2; runs=$3; shift 3
H=$(cd "$(dirname "$0")" && pwd)
for c in "$@"; do k=${c%:*}; n=${c#*:}
  for r in $(seq $runs); do
    if [ "$pp" = sympy ]; then
      (cd /tmp && timeout 900 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env LIB=sympy /home/keckj/Documents/symox/.venv/bin/python $H/nest.py $k $n)
    else
      timeout 900 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env LIB=symox PYTHONPATH=$pp /home/keckj/Documents/symox/.venv/bin/python $H/nest.py $k $n
    fi | sed "s/^/$label /"
  done
done
