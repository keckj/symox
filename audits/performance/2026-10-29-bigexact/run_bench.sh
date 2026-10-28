#!/bin/bash
# usage: run_bench.sh sympy|before|after RUNS OUT
# before: PYTHONPATH of the base tree (58e50206 sources and _core);
# after: this worktree. Each case runs RUNS fresh processes (REP repeats each, default 5).
H=$(cd "$(dirname "$0")" && pwd)
mode=$1; runs=$2; out=$3
PY=/home/keckj/Documents/symox/.venv/bin/python
case $mode in
  sympy) mod=sympy; pp="";;
  before) mod=symox; pp=${BEFORE:-/home/keckj/Documents/symox-worktrees/.tooling/scratch/bigexact/base/python};;
  after) mod=symox; pp=${AFTER:-$H/../../../python};;
esac
CASES=(
  "factorial(100001)" "fibonacci(1000001)" "catalan(100001)" "bernoulli(5002)"
  "legendre(2001, Rational(1, 3))" "lucas(1000001)" "tribonacci(100001)"
  "subfactorial(100001)" "rf(3, 100001)" "factorial2(200001)" "euler(5002)"
  "binomial(300000, 150000)" "gamma(Rational(200001, 2))" "factorial(3000)"
  "bernoulli(1000)" "binomial(3000, 1000)" "legendre(300, Rational(1, 3))"
)
: > "$out"
for c in "${CASES[@]}"; do
  for r in $(seq 1 $runs); do
    (cd /tmp && timeout ${TO:-200} systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
      env PYTHONPATH=$pp $PY $H/bench.py $mod ${REP:-5} "$c") >> "$out" || printf '%s\t%s\tTIMEOUT\n' "$mod" "$c" >> "$out"
  done
done
