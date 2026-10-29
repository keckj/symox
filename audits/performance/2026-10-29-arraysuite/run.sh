#!/bin/bash
# usage: run.sh <python-dir> <runs> <out.jsonl>  - one fresh process per run
D=$(dirname "$(readlink -f "$0")")
for r in $(seq "$2"); do
  timeout 300 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
    env PYTHONPATH="$1" /home/keckj/Documents/symox/.venv/bin/python "$D/bench.py"
done > "$3"
