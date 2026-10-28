#!/bin/bash
# usage: run.sh <label> <pythonpath-or-"sympy"> <out.jsonl> [timeout] [cases...]
label=$1; pp=$2; out=$3; to=${4:-600}; shift 4
here=$(cd "$(dirname "$0")" && pwd)
PY=/home/keckj/Documents/symox/.venv/bin/python
cases=${@:-$(grep -o "^    '[a-z0-9_]*'" "$here/bench.py" | tr -d "' ")}
for c in $cases; do
  rep=2; case $c in count_5000*|factor_300|factor_expr_300|factor_301) rep=1;; esac
  if [ "$pp" = sympy ]; then
    line=$(cd /tmp && timeout $to systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env BENCH_LIB=sympy $PY "$here/bench.py" $c $rep)
  else
    line=$(cd /tmp && timeout $to systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$pp $PY "$here/bench.py" $c $rep)
  fi
  [ -z "$line" ] && line="{\"case\": \"$c\", \"timeout\": $to}"
  echo "{\"label\": \"$label\", \"data\": $line}" | tee -a "$out"
done
