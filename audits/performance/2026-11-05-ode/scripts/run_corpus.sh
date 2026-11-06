#!/bin/sh
# usage: run_corpus.sh WORKTREE OUTDIR
# Replays the 315-case ODE corpus (classify / explicit hint / default dsolve)
# with live SymPy 1.14 and the worktree's symox, using the archived runner
# audits/sympy-1.14.0/corpus-2026-10-26/replay-ode.py (16 jobs, 60 s per call,
# 8 GiB per worker). Score = equal_output + equal_error_observation.
set -e
wt=$(cd "$1" && pwd); out=$2
cases=$(mktemp); zcat "$wt/audits/sympy-1.14.0/corpus-2026-10-26/ode-cases.json.gz" > "$cases"
sha=$(sha256sum "$wt/python/symox/_core.abi3.so" | cut -d' ' -f1)
cd "$wt"
PYTHONPATH="$wt/python" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ~/symox-ci/venv/bin/python audits/sympy-1.14.0/corpus-2026-10-26/replay-ode.py \
  --cases "$cases" --out "$out" --jobs 16 --timeout 60 --native-sha "$sha" > "$out.log" 2>&1
~/symox-ci/venv/bin/python "$wt/audits/performance/2026-11-05-ode/scripts/score.py" "$out"
