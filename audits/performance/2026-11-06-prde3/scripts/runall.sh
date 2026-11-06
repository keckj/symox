#!/bin/bash
# usage: runall.sh TAG [stock]  -- default-tower corpus + replays (symox native, optional stock)
cd "$(dirname "$0")/../../../.."
R=audits/performance/2026-11-06-prde3/raw
S=audits/performance/2026-11-06-prde3/scripts/random_param.py
P=audits/performance/2026-11-05-prde
PY=~/symox-ci/venv/bin/python
TAG=$1
export PYTHONPATH=python PRDE_BACKEND_MODULE=symox
($PY $S $R/rp_${TAG}_symox.json 600 > $R/rp_${TAG}_symox.log 2>&1; echo EXIT $? >> $R/rp_${TAG}_symox.log) &
($PY $P/scripts/replay.py $P/raw/calls_tests.json $R/tests_${TAG}_symox.json > $R/tests_${TAG}_symox.log 2>&1;
 $PY $P/scripts/replay.py $P/raw/calls_corpus.json.gz $R/corpus_${TAG}_symox.json > $R/corpus_${TAG}_symox.log 2>&1; echo EXIT $? >> $R/corpus_${TAG}_symox.log) &
if [ -n "$2" ]; then
(PRDE_STOCK=1 $PY $S $R/rp_${TAG}_stock.json 600 > $R/rp_${TAG}_stock.log 2>&1; echo EXIT $? >> $R/rp_${TAG}_stock.log) &
(PRDE_STOCK=1 $PY $P/scripts/replay.py $P/raw/calls_tests.json $R/tests_${TAG}_stock.json > $R/tests_${TAG}_stock.log 2>&1;
 PRDE_STOCK=1 $PY $P/scripts/replay.py $P/raw/calls_corpus.json.gz $R/corpus_${TAG}_stock.json > $R/corpus_${TAG}_stock.log 2>&1) &
fi
wait
tail -n 2 $R/*_${TAG}_*.log
