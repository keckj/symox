# solvegaps (wave 31), 2026-10-30

Base 7e8709e3 ("before") against the feat-solvegaps-20261030d build
("after") and SymPy 1.14. Run on the local machine under the 6G memory cap,
while a single-job `cargo test` build was running.

## Ten solve/solveset/unrad calls

`bench.py` runs each call once (cold, in a process that has already
imported the library) and then 5 times (median, warm). `run.sh` runs 3
fresh processes per label, and `summarize.py raw` prints the medians. The
raw JSON is in `raw/`.

| call | SymPy cold / warm (ms) | before cold / warm (ms) | after cold / warm (ms) | after == SymPy |
|---|---|---|---|---|
| solveset 11534 | 181.26 / 31.20 | 14.46 / 13.48 | 15.77 / 14.83 | yes |
| solve x**(2*y)-3 | 82.50 / 50.14 | 6.61 / 5.84 | 8.48 / 7.63 | yes |
| solve x**(2*y)-3*x**y+2 | 47.32 / 34.95 | 2.22 / 2.12 | 5.73 / 5.68 | yes |
| solve x**(3*y)-3 | 126.15 / 67.41 | 4.85 / 4.81 | 10.92 / 10.76 | yes |
| solve sqrt(x)+sqrt(x+1)-3 | 25.41 / 4.82 | 4.34 / 3.98 | 4.59 / 4.04 | yes |
| solveset sqrt(x)+sqrt(x+1)-3 | 11.94 / 11.96 | 0.91 / 0.81 | 0.91 / 0.84 | yes |
| solveset sqrt(x-1)-x+3 | 19.00 / 18.76 | 0.81 / 0.79 | 0.83 / 0.80 | yes |
| solve x**(1/3)+sqrt(x)-1 | 301.36 / 25.73 | 38.24 / 31.47 | 39.04 / 32.08 | yes |
| unrad cov | 3.29 / 2.70 | 0.09 / 0.02 | 0.87 / 0.68 | yes |
| solve x*exp(x)-1 | 40.44 / 28.13 | 5.05 / 4.12 | 4.98 / 4.19 | yes |

- before, solve x**(2*y)-3: `[3**(1/(2*y))]` (SymPy `[3**(1/(2*y)), (-sqrt(3))**(1/y)]`)
- before, solve x**(2*y)-3*x**y+2: `error: NotImplementedError: could not solve (x**y - 2)*(x**y - 1) for x` (SymPy `[1, 2**(1/y)]`)
- before, solve x**(3*y)-3: `[3**(1/(3*y))]` (SymPy `[3**(1/(3*y)), (-3**(1/3)/2 - 3**(5/6)*I/2)**(1/y), (-3**(1/3)/2 + 3**(5/6)*I/2)**(1/y)]`)
- before, unrad cov: `error: NotImplementedError: unrad flags ['cov'] are not supported` (SymPy `(_p**3 + _p**2 - 1, [_p, _p**6 - x])`)

The "after" rows were measured again with the final build, which also
checks parametric candidates with the native checksol (`simplify=False`),
as SymPy's `_solve` does. That check costs about 1.3 ms on
`x**(2*y)-3` and 6.5 ms on `x**(3*y)-3` (three complex candidates); the
earlier "after" medians without it were 6.28 and 4.29 ms warm.

`solve(x**(y/2) + x**y - 6, x)`, measured the same way by hand (3
processes, cold / median of 5 warm): SymPy `[2**(2/y)]` 247.67 / 83.61 ms,
after `[2**(2/y)]` 22.59 / 20.67 ms. The base raised NotImplementedError.

The "before" rows with an error or a missing root are not speed
comparisons. Before this change, `unrad` with `cov` raised at once.
`x**(2*y) - 3*x**y + 2` was not solved.

## checksol differential

`checksol/cslog_plugin.py` is a pytest plugin. It wraps SymPy 1.14's
`checksol` and logs every single-symbol call made by SymPy's own
test_solvers.py and test_solveset.py, and by the ODE corpus
(`.tooling/sspeed/runner.py sympy`), as srepr, flags and result. The logs
are `checksol/pairs_tests.jsonl.gz` (1414 distinct triples) and
`checksol/pairs_ode_new.jsonl.gz` (1363 further triples).
`checksol/diff_checksol.py` replays them through
`_core._checksol_native`. `checksol/summary.json` has the counts:

- old native checksol (`minimal` ignored), SymPy tests: 34 Falses where
  SymPy returns None (all with `minimal=True`);
- new native checksol, SymPy tests: 0 false Falses. The (SymPy, native)
  pairs are:
  - True/True 786, False/False 425, None/None 192;
  - True/None 5, False/None 1;
  - 5 not replayed by the harness (Piecewise pairs and units);
  - native time 2.6 s in all;
- new native checksol, ODE corpus: 0 mismatches. The pairs are:
  - True/True 511, False/False 148, None/None 558;
  - 146 not replayed (applied undefined functions in the srepr).

tests/test_native_solvegaps.py replays 326 of these pairs, with the SymPy
values stored in tests/data/native_checksol_sympy114.json.
