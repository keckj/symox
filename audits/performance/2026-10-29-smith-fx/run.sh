#!/bin/bash
# usage: run.sh <base_tree> <new_tree> <out.jsonl> [processes]
BASE=$1; NEW=$2; OUT=$3; P=${4:-5}; PY=/home/keckj/Documents/symox/.venv/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
[ -z "$APPEND" ] && : > "$OUT"
CASES=${CASES:-"zz_snf_doc_3x3 zz_invf_doc_3x3 zz_decomp_doc_3x3 qq_snf_3x3 qq_decomp_6x6 qq_invf_10x10 gf7_invf_10x10 gf101_decomp_8x8 gf65537_snf_20x20 qqx_invf_jordan_4x4 qqx_decomp_jordan_4x4 qqx_snf_char_6x6 qqx_decomp_3x4"}
for i in $(seq 1 "$P"); do
  for c in $CASES; do
    for side in base new; do
      T=$BASE; [ $side = new ] && T=$NEW
      (cd /tmp && PYTHONPATH=$T/python timeout 300 $PY "$HERE/bench_case.py" "$c" | sed "s/^{/{\"side\": \"$side\", \"run\": $i, /") >> "$OUT"
    done
  done
done
