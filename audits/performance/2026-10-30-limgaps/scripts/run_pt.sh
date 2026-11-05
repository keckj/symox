#!/bin/bash
# run_pt.sh TAG WORKTREE -> raw/pt/TAG_f.{txt,json}: pytest on symox copies of
# the SymPy 1.14 test files, with the native decline trace.
A=$(cd $(dirname $0)/.. && pwd)
W=${WORK:-/tmp/limgaps-work}
cd $W/copies
for f in limits gruntz series nseries order; do
  (TRACE_OUT=$A/raw/pt/$1_$f.json PYTHONPATH=$2/python prlimit --as=6000000000 timeout 1800 ~/symox-ci/venv/bin/python -m pytest -q -p no:cacheprovider test_$f.py > $A/raw/pt/$1_$f.txt 2>&1) &
done
wait
