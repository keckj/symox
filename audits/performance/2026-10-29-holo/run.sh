#!/bin/bash
# run.sh <before-python-dir> <after-python-dir> <out.jsonl>
D=$(dirname "$0"); PY=/home/keckj/Documents/symox/.venv/bin/python
CASES="sin_x0_0 sin_pow3 exp_pow_neg sqrt_x x2_plus_x sinc erfc erf_rewrite_erfc sin_plus_cos float_cos float_param_sin"
: > "$3"
for c in $CASES; do for side in before after; do
  dir=$1; [ $side = after ] && dir=$2
  for r in 1 2 3 4 5; do (cd /tmp && timeout 120 env PYTHONPATH=$dir $PY $D/bench.py $side first $c) >> "$3"; done
  (cd /tmp && timeout 300 env PYTHONPATH=$dir $PY $D/bench.py $side repeat $c) >> "$3"
done; done
