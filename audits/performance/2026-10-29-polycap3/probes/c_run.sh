#!/bin/bash
# usage: c_run.sh LIB pythonpath case...
lib=$1; pp=$2; shift 2
for k in "$@"; do
  (cd /tmp && timeout ${CTIMEOUT:-180} systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env LIB=$lib PYTHONPATH=$pp /home/keckj/Documents/symox/.venv/bin/python /home/keckj/Documents/symox-worktrees/.tooling/scratch/polycap3/c_probe.py $k 2>&1 | tail -1) || echo "$k TIMEOUT/exit"
done
