#!/bin/bash
# usage: run.sh <base_tree> <new_tree> <out.jsonl> [processes]
BASE=$1; NEW=$2; OUT=$3; P=${4:-5}; PY=/home/keckj/Documents/symox/.venv/bin/python
HERE=$(cd "$(dirname "$0")" && pwd)
[ -z "$APPEND" ] && : > "$OUT"
CASES=${CASES:-"x6p8 x4m2x2p4 x4p2 x8p1 x8p16 x4pxp1 x2_x6p8 x4p1 x4mx2p1 xm1_4p1 x6p27 x6p1 x3pxp1 x3m3xp1 x5p1"}
for i in $(seq 1 "$P"); do
  for c in $CASES; do
    for side in base new; do
      T=$BASE; [ $side = new ] && T=$NEW
      (cd /tmp && timeout 300 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$T/python $PY "$HERE/bench_case.py" "$c" | sed "s/^{/{\"side\": \"$side\", \"run\": $i, /") >> "$OUT"
    done
  done
done
