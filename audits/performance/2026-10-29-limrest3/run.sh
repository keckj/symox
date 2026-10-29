#!/bin/bash
# run.sh <base_pythonpath> <new_pythonpath> <out.jsonl> [n]
# n fresh processes per case and side, base and new interleaved.
set -u
B=$1; N=$2; OUT=$3; R=${4:-5}
PY=/home/keckj/Documents/symox/.venv/bin/python
H=$(dirname "$(readlink -f "$0")")/bench_case.py
CASES="ei erf ugam subfac stirl harm_m1 harm_z1 atan_lim atan_int c_sinx c_e c_gruntz c_ei c_erf c_harm3 c_lgam c_atan"
: > "$OUT"
for i in $(seq 1 "$R"); do
  for c in $CASES; do
    for side in base new; do
      if [ $side = base ]; then P=$B; else P=$N; fi
      line=$(cd /tmp && timeout 300 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH="$P" $PY "$H" "$c" 2>/dev/null | tail -1)
      echo "{\"side\": \"$side\", \"run\": $i, \"data\": ${line:-null}}" >> "$OUT"
    done
  done
done
