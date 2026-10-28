#!/bin/bash
# usage: run.sh <base_tree> <new_tree> <out.jsonl> [processes]
BASE=$1; NEW=$2; OUT=$3; P=${4:-5}; PY=/home/keckj/Documents/symox/.venv/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
[ -z "$APPEND" ] && : > "$OUT"
CASES=${CASES:-"snf_doc_3x3 invf_doc_3x3 decomp_doc_3x3 decomp_6x6 invf_8x8 decomp_4x6_big invf_12x12 matrix_snf_4x4 round_two_x3m2 round_two_cyclo7 prime_decomp_5 nf_hnf_5x8 nf_hnf_mod_4x4"}
for i in $(seq 1 "$P"); do
  for c in $CASES; do
    for side in base new; do
      T=$BASE; [ $side = new ] && T=$NEW
      (cd /tmp && PYTHONPATH=$T/python $PY "$HERE/bench_case.py" "$c" | sed "s/^{/{\"side\": \"$side\", \"run\": $i, /") >> "$OUT"
    done
  done
done
