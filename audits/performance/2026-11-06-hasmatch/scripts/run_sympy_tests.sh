#!/bin/sh
# run_sympy_tests.sh PYTHONDIR OUT.xml : the SymPy 1.14 simplify / meijerint /
# collect / has test files, rewritten for symox, run against PYTHONDIR.
A=$(cd $(dirname $0)/.. && pwd)
W=${WORK:?set WORK}
[ -f $W/copies/test_simplify_trigsimp.py ] || ~/symox-ci/venv/bin/python $A/scripts/mkcopies.py $W/copies
cd $W/copies
PYTHONPATH=$1 prlimit --as=6000000000 ~/symox-ci/venv/bin/python -m pytest -q -p no:cacheprovider --continue-on-collection-errors \
  --junitxml=$2 . > $2.log 2>&1 || true
