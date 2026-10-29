#!/bin/bash
# usage: run.sh <python dir of the tree> <label>; prints one line per case
here=$(cd "$(dirname "$0")" && pwd)
for c in prs_200_150:1 prs_150_200:1 prs_120_90:3 prs_60_45:5 prs_40_30:5 subres_60_45:5 poly_prs_60_45:5 prs_small:20 prs_abnormal:20 prs_dense_30:5 factor_2000:0; do
  case=${c%%:*}; rep=${c##*:}
  (cd /tmp && timeout 300 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 env PYTHONPATH=$1 /home/keckj/Documents/symox/.venv/bin/python -I -c "import sys; sys.path.insert(0, '$1'); sys.argv=['b','$case','$rep']; exec(open('$here/bench.py').read())" 2>&1 | tail -1 | sed "s/^/$2 /") || echo "$2 $case TIMEOUT"
done
