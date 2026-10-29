# Wave 21 slowcases: nested exponential solve, 1/(x**5+x+1)**3

Release build, Python 3.11 venv, 8-core machine shared with other agents
(load about 7), every run under `timeout` and
`systemd-run -p MemoryMax=6G -p MemorySwapMax=0`. One process per run (first
call, no warm cache). `base` = the 82d674aa extension, `new` = this branch,
`sympy` = SymPy 1.14.0 run from /tmp.

## Harness

- `harness/nest.py {mixed|exp|expp1} n`: `solve(T - 5)` for the tower
  `T = exp(exp(...)+1)+1` (mixed), `solve(E - 2)` and `solve(E + 1)` for
  `E = exp(exp(...(x)))` (exp, expp1).
- `harness/bench.sh label pythonpath|sympy runs kind:n...`: the table below
  (`raw/nest_bench.txt`, 3 runs each).
- `harness/diffsolve.py`: 51 solve inputs, srepr and times
  (`raw/diff_base.json`, `raw/diff_new.json`: 0 srepr differences, total
  0.63 s base, 0.54 s new).
- `harness/tvals.py`: solutions and `is_algebraic` answers printed by SymPy,
  base and new (`raw/tv_*.txt`; base and new identical).
- `harness/ratint.py k`, `harness/split_probe.py k`: the integral and the
  per-factor `log_to_real` prototype (not enabled).

## Nested exponential solve, median of 3 (seconds)

| case | base | new | SymPy 1.14 | base/new |
|---|---|---|---|---|
| mixed 4 | 0.104 | 0.090 | 0.255 | 1.2x |
| mixed 6 | 1.620 | 0.360 | 0.529 | 4.5x |
| mixed 8 | 21.49 | 1.930 | 2.135 | 11x |
| mixed 10 (1 run) | 426 | 58.6 | 77.7 | 7.3x |
| exp 6 | 0.188 | 0.102 | 0.163 | 1.8x |
| exp 8 | 0.948 | 0.367 | 0.193 | 2.6x |
| exp 10 | 7.078 | 2.107 | 1.067 | 3.4x |
| expp1 6 | 0.172 | 0.081 | 0.207 | 2.1x |
| expp1 8 | 1.087 | 0.323 | 0.497 | 3.4x |
| expp1 10 | 6.596 | 1.739 | 6.821 | 3.8x |

Solution sets are identical to base everywhere and to SymPy up to the
pre-existing forms `log(pi) + I*pi/2` (SymPy `log(I*pi)`) and `log(4)`
(SymPy `2*log(2)`).

Profiles (py-spy native): before, 95% of `mixed 6` was the numerical zero
test of `algebraic_proof` while `log` of a numeric complex argument asked
`is_algebraic`; after that fix, the native solve repeated the symbolic
inverse of the same exponential generator at every depth (3586 `solve_rec`
calls at depth 6). Remaining for `exp 10` (cProfile, 3.1 s under the
profiler): 0.64 s native solve, 2.4 s Python simplify of the solution
(1024 nested `simplify` calls through `_simplify_function_arguments`, the
same recursion as SymPy, with symox's Python `has`/`replace` per call). Gap.

## integrate(1/(x**5 + x + 1)**3, x)

- symox base and new: memory exhausted at 6 GB (`raw/symox_k1.log`, exit
  137, already for k = 1).
- SymPy 1.14: k = 3 not finished in 1800 s (`raw/sympy_k3.log`, exit 124);
  k = 1 not finished in 1800 s either (`raw/sympy_k1.log`, exit 124).
- Blowup: Python `log_to_real` general path, `Poly(a.xreplace({u: r_u, v:
  r_v}), x)` with Cardano radicals, for the reducible quintic
  `q(t) = -3381 t**5 + 160 t**3 + 80 t**2 + 15 t + 1` (factors of degree 2
  and 3).
- Prototype `split_probe.py`: per irreducible factor, 0.21 s (k = 1, length
  9694) and 0.84 s (k = 3, length 27728), derivative minus integrand 0 to 30
  digits at 1/3, -5/7, 2. Not enabled (changes the form for reducible `q`
  that SymPy computes; no work estimate yet). Gap.
