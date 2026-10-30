#!/bin/bash
# bench.sh sympy|symox PYTHONPATH OUT [CALLS]
B=/home/keckj/Documents/symox-worktrees/limitsrest-20261030h/audits/performance/2026-10-30-limitsrest/bench.py
: > $3
CALLS=${4:-"lim_pow4 lim_pow4b lim_ei lim_ugamma lim_stirling lim_sinx ser_abs_cdir ser_log_cdir ser_abs_cdirI ser_nonsym ser_logx ser_sin"}
for c in $CALLS; do
  (cd /tmp && env PYTHONPATH=$2 timeout 120 systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 /home/keckj/Documents/symox/.venv/bin/python $B $1 $c 5) >> $3 || echo "{\"call\": \"$c\", \"lib\": \"$1\", \"result\": \"TIMEOUT\"}" >> $3
done
echo done >> $3.done
