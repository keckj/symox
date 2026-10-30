#!/bin/bash
# run.sh <after-python-dir> <before-python-dir> <outdir>: 3 fresh processes
# per mode under a 6G cap, from /tmp. SymPy 50 reps, symox 200.
set -e
H=$(cd "$(dirname "$0")" && pwd)
PYBIN=${PYBIN:-$HOME/symox-ci/venv/bin/python}
mkdir -p "$3"
: > "$3/sympy.jsonl"; : > "$3/before.jsonl"; : > "$3/after.jsonl"
cd /tmp
for r in 1 2 3; do
  prlimit --as=6000000000 $PYBIN "$H/bench_a2m.py" sympy 50 >> "$3/sympy.jsonl" 2>/dev/null
  PYTHONPATH="$2" prlimit --as=6000000000 $PYBIN "$H/bench_a2m.py" symox 200 >> "$3/before.jsonl" 2>/dev/null || true
  PYTHONPATH="$1" prlimit --as=6000000000 $PYBIN "$H/bench_a2m.py" symox 200 >> "$3/after.jsonl" 2>/dev/null
done
