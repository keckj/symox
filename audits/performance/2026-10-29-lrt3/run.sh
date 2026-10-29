#!/bin/bash
# usage: run.sh <base_tree> <new_tree> <out.jsonl> [processes]
BASE=$1; NEW=$2; OUT=$3; P=${4:-5}; PY=/home/keckj/Documents/symox/.venv/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
[ -z "$APPEND" ] && : > "$OUT"
CASES=${CASES:-"lrt_x4p1 lrt_x4p4 lrt_x4pq lrt_x4mx2p1 lrt_mixed lrt_x4_x6p1 lrt_x6p1 lrt_x3p1_x_1 lrt_x3m2 lrt_x5p1 lrt_x6m2 lrt_x8p1 lrt_x2_x4px1 lrt_x5p2 lrt_x7p1"}
for i in $(seq 1 "$P"); do
  for c in $CASES; do
    for side in base new; do
      T=$BASE; [ $side = new ] && T=$NEW
      (cd /tmp && timeout 300 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$T/python $PY "$HERE/bench_case.py" "$c" | sed "s/^{/{\"side\": \"$side\", \"run\": $i, /") >> "$OUT"
    done
  done
done
