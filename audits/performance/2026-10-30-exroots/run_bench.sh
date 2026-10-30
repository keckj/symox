#!/bin/bash
# run_bench.sh sympy|symox OUT.jsonl [WORKTREE]
WHO=$1; OUT=$2; WT=${3:-/home/keckj/Documents/symox-worktrees/exroots-20261030i}
H=$(dirname "$0")
: > "$OUT"
for i in $(seq 0 11); do
  if [ "$WHO" = sympy ]; then
    (cd /tmp && timeout 600 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 /home/keckj/Documents/symox/.venv/bin/python "$H/bench.py" sympy $i) >> "$OUT" 2>&1
  else
    (cd "$WT" && timeout 600 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=python /home/keckj/Documents/symox/.venv/bin/python "$H/bench.py" symox $i) >> "$OUT" 2>&1
  fi
done
