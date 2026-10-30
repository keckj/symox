#!/bin/bash
# run_tl.sh TAG LIB WORKTREE -> raw/tl/TAG_f.json (+ .trace.json, .stats.json for symox)
# Statement-level evaluation of the SymPy 1.14 series test files.
A=$(cd $(dirname $0)/.. && pwd)
T=$(echo ~/symox-ci/venv/lib/python3*/site-packages/sympy/series/tests)
W=${WORK:-/tmp/limgaps-work}
mkdir -p $W/copies $A/raw/tl
cd $W
[ -f copies/test_limits.py ] || ~/symox-ci/venv/bin/python $A/scripts/mkcopies.py $T
cp $A/scripts/conftest_trace.py copies/conftest.py
for f in limits gruntz series nseries order formal aseries lseries demidovich; do
  if [ $2 = sympy ]; then F=$T/test_$f.py; else F=$W/copies/test_$f.py; fi
  (PYTHONPATH=$3/python prlimit --as=6000000000 timeout 1800 ~/symox-ci/venv/bin/python $A/scripts/tl.py $2 $A/raw/tl/$1_$f.json $F > $A/raw/tl/$1_$f.log 2>&1) &
done
wait
