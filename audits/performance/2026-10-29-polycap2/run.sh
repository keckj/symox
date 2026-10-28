#!/bin/bash
# usage: run.sh <python-dir-before> <python-dir-after> <out.jsonl>
# Each case: 3 fresh processes per tree (first call + median of 5 repeats).
here=$(dirname "$(readlink -f "$0")")
py=/home/keckj/Documents/symox/.venv/bin/python
: > "$3"
for case in $($py -c "import sys; sys.path.insert(0,'$here'); import bench; print(' '.join(bench.CASES))"); do
  for tree in before after; do
    dir=$1; [ $tree = after ] && dir=$2
    for run in 1 2 3; do
      line=$(cd /tmp && timeout 120 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$dir $py "$here/bench.py" $case 5 2>/dev/null | tail -1)
      [ -z "$line" ] && line="{\"case\": \"$case\", \"error\": \"timeout or crash\"}"
      echo "{\"tree\": \"$tree\", \"run\": $run, \"data\": $line}" >> "$3"
    done
  done
done
