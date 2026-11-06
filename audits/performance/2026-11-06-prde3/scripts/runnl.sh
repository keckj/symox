#!/bin/bash
# usage: runnl.sh TAG (runs symox native and stock nonlin corpus)
cd "$(dirname "$0")/../../../.."
R=audits/performance/2026-11-06-prde3/raw
S=audits/performance/2026-11-06-prde3/scripts/random_param.py
PY=~/symox-ci/venv/bin/python
TAG=$1
(PRDE_TOWERS=nonlin PYTHONPATH=python $PY $S $R/nl_${TAG}_symox.json 600 > $R/nl_${TAG}_symox.log 2>&1; echo EXIT $? >> $R/nl_${TAG}_symox.log) &
if [ -n "$2" ]; then
(PRDE_TOWERS=nonlin PRDE_STOCK=1 PYTHONPATH=python $PY $S $R/nl_${TAG}_stock.json 600 > $R/nl_${TAG}_stock.log 2>&1; echo EXIT $? >> $R/nl_${TAG}_stock.log) &
fi
wait
tail -n 3 $R/nl_${TAG}_symox.log
