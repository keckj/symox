#!/bin/bash
# usage: run_bench.sh sympy|before|after RUNS OUT
# before: PYTHONPATH of the base tree (9e2340a1 sources and its _core);
# after: this worktree. Each case runs RUNS fresh processes (REP repeats each, default 3).
H=$(cd "$(dirname "$0")" && pwd)
mode=$1; runs=$2; out=$3
PY=/home/keckj/Documents/symox/.venv/bin/python
case $mode in
  sympy) mod=sympy; pp="";;
  before) mod=symox; pp=${BEFORE:-/home/keckj/Documents/symox-worktrees/.tooling/scratch/biggcd/base/python};;
  after) mod=symox; pp=${AFTER:-$H/../../../python};;
esac
CASES=(
  "Rational(3**200001, 2**100000)"
  "Rational(3**200001 * 2**7, 2**100000)"
  "gamma(Rational(200001, 2))"
  "igcd(A, B)"
  "gcd(Integer(A), Integer(B))"
  "Rational(A, B)"
  "Integer(A) / Integer(B)"
  "Poly(A*x**2 + B*x + G, x).primitive()"
  "Add(*[Rational(1, 10**40 + k) for k in range(1000)])"
  "functools.reduce(operator.add, [Rational(k, 10**40 + k) for k in range(1, 1001)])"
  "harmonic(20000)"
  "Rational(6, 4)"
  "igcd(12, 18)"
  "Rational(3, 7) * Rational(14, 9) + Rational(5, 11)"
  "Add(*[Rational(1, k) for k in range(1, 101)])"
  "harmonic(100)"
)
: > "$out"
for c in "${CASES[@]}"; do
  for r in $(seq 1 $runs); do
    (cd /tmp && timeout ${TO:-200} systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
      env PYTHONPATH=$pp $PY $H/bench.py $mod ${REP:-3} "$c") >> "$out" || printf '%s\t%s\tTIMEOUT\n' "$mod" "$c" >> "$out"
  done
done
