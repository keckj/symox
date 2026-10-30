# exroots (wave 36): EX roots, the simplify denominator step, 2026-10-30

Release build, symox before = base `173328b9` (`raw/*_before*`), after = branch
`feat-exroots-20261030i`, reference = SymPy 1.14.0 run from /tmp. Each call
is timed in a fresh process (first call) and as the median of 5 repeats
(`bench.py`, `run_bench.sh sympy|symox OUT [WORKTREE]`). The after run shared
the machine with the two corpus runs; a second run without them
(`raw/bench_after2.jsonl`) is within 10%.

| call | before first / repeat ms | after first / repeat ms | SymPy first / repeat ms | len before / after / SymPy |
|---|---|---|---|---|
| `roots((3125/4 + 3125*sqrt(5)/4)*x**4 + (-125*sqrt(5)/4 - 125/2)*x**2 + 5*sqrt(5)/64 + 75/64)` | 1.13 / 0.83 | 3.34 / 3.09 | 259.57 / 106.44 | 148 / 204 / 204 |
| `roots(sqrt(3)*x**2 - 2*x + sqrt(2))` | 0.56 / 0.27 | 1.38 / 1.07 | 110.73 / 11.17 | 88 / 82 / 82 |
| `roots(x**4 - sqrt(3)*x**2 - 1)` | 0.46 / 0.21 | 0.79 / 0.51 | 115.17 / 13.95 | 120 / 120 / 120 |
| `roots(sqrt(2)*x**4 + 3*x**2 + sqrt(5))` | 0.83 / 0.53 | 1.71 / 1.42 | 150.67 / 30.73 | 226 / 214 / 214 |
| `simplify((2 + sqrt(5))/(50*(1 + sqrt(5))))` | 2.77 / 2.36 | 2.46 / 2.10 | 140.92 / 23.08 | 19 / 19 / 19 |
| `simplify((1 + sqrt(3))/(2 - sqrt(3)))` | 3.15 / 2.71 | 2.77 / 2.46 | 118.76 / 12.38 | 13 / 13 / 13 |
| `simplify((2 + sqrt(5))/(50*(1 + sqrt(5))) + 1/(25*(1 + sqrt(5))))` | 3.62 / 3.17 | 3.33 / 2.99 | 147.62 / 29.35 | 21 / 21 / 21 |
| `simplify(sqrt(-1/(25*(1 + sqrt(5))) + (2 + sqrt(5))/(50*(1 + sqrt(5)))))` | 5.86 / 5.33 | 5.76 / 5.37 | 129.37 / 12.92 | 40 / 40 / 40 |
| `integrate(1/(x**3-2))` | 30.40 / 7.10 | 27.63 / 8.34 | 198.53 / 44.24 | 136 / 136 / 136 |
| `integrate(1/(x**3+x+1))` | 175.43 / 158.26 | 394.51 / 369.70 | 7752.01 / 2084.30 | 7892 / 9543 / 16653 |
| `integrate(1/(x**3+2*x-1))` | 163.79 / 147.80 | 363.12 / 344.07 | 7962.30 / 1961.67 | 8214 / 9987 / 17352 |
| `integrate(1/(x**5+1))` | 34.07 / 16.54 | 51.65 / 32.72 | 1195.35 / 728.36 | 399 / 493 / 493 |

The roots calls are slower than before because they now run SymPy's
`simplify` on the EX linear and quadratic roots (the before outputs were not
SymPy's). The cubic integrates compute both EX constant quotients of
`log_to_atan` and keep the shorter one; they stay 5-20 times faster than
SymPy.

## Outputs

- `rootscmp.py`: 25 `roots` cases (dict and `multiple=True`), raw
  `raw/roots_{before,after,sympy}.txt`. Identical to SymPy: 6/25 before,
  24/25 after. Case 20, `x**4 + 2*x**2 + 3`, keeps the quadratic-in-x**2
  form (`sqrt(-1 + sqrt(2)*I)`), SymPy uses its binomial/trig form; the
  difference predates the wave (`try_decompose` uses `roots_quadratic`).
- `integ.py`: the four integrates with srepr, raw
  `raw/integ_{before,after,sympy}.txt`. Identical to SymPy after:
  `1/(x**3-2)` (already before), `1/(x**5+1)` (new). `1/(x**3+x+1)` and
  `1/(x**3+2*x-1)` differ in form (exact, certified by differentiation in
  tests/test_native_exroots.py, shorter than SymPy's).

## Corpora

- Integration corpus: `{"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1, "SAME": 1171}`, unchanged.
- ODE corpus: `{'cls': 303, 'hint': 279, 'def': 274, 'n': 315}`, unchanged.
