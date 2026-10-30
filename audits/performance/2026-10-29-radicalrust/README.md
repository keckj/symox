# radicalrust (wave 29): native `unrad`, `_solve_radical` and negative-base steps

Ten public calls that go through `unrad`, the `solve` radical fallback, the
`solveset` radical step and the negative-base steps. Builds are release builds.

- before: base 62d4b985, which has the Python resultant `unrad` and the Python
  steps (`base_tree/python` with the base `_core.abi3.so`);
- after: this branch, with native `symox::solve::unrad` and
  `symox::sets::solveset_radical`;
- SymPy 1.14.0, run from /tmp with the project venv.

The runs were made on 2026-10-30 while the two corpora and a single-job
`cargo test` ran on the same 8-core machine (load average about 5). The
absolute numbers are therefore noisy; the ratios are the useful part.

## Harness

- `bench.py <symox|sympy> <first|repeat> [reps] [only]` prints one JSON line
  per case, with the times and the result.
- `run.sh <label> <PYTHONPATH|sympy> <absolute outdir>` measures two things:
  - first call: each case runs in a fresh process, 5 runs, and the median is
    reported;
  - repeated call: one process, 3 runs of a median of 20 calls, and the
    median of the 3 is reported.
  - Every process runs under `timeout` and
    `systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0`.
- Raw numbers are in `raw/{before,after,sympy}-{first,repeat}.jsonl`.
  `python summarize.py before after sympy > summary.md` rebuilds the table.

## Results (ms)

| # | call | before first (ms) | after first (ms) | sympy first (ms) | before repeat (ms) | after repeat (ms) | sympy repeat (ms) |
|---|---|---|---|---|---|---|---|
| 0 | `unrad(sqrt(x) + sqrt(x + 1) + sqrt(1 - x))` | 3.26 | 0.80 | 40.28 | 2.573 | 0.521 | 6.419 |
| 1 | `unrad(root(x, 3) - root(x + 1, 4)/2 + root(x + 2, 3))` | 5.97 | 3.18 | 103.49 | 4.997 | 2.548 | 39.220 |
| 2 | `unrad(sqrt(x) + root(x, 3) + y)` | 4.87 | 1.25 | 27.65 | 4.025 | 0.644 | 3.032 |
| 3 | `solve(sqrt(x) + sqrt(x + 1) - 3, x)` | 6.78 | 4.95 | 135.55 | 5.442 | 3.855 | 4.928 |
| 4 | `solve(sqrt(x) + root(x, 3) - 2, x)` | 9.59 | 9.54 | 148.27 | 7.509 | 7.480 | 14.756 |
| 5 | `solve(root(x + 1, 3) + root(x, 3) - 1, x)` | 12.92 | 13.29 | 189.90 | 10.672 | 11.306 | 47.163 |
| 6 | `solveset(sqrt(2*x - 1) - sqrt(x - 4) - 2, x, S.Reals)` | 5.07 | 1.83 | 138.23 | 3.462 | 1.174 | 25.182 |
| 7 | `solveset(-8*x**2/(9*(x**2 - 1)**(S(4)/3)) + 4/(3*(x**2 - 1)**(S(1)/3)), x, S.Complexes)` | 5.95 | 2.33 | 60.88 | 3.950 | 1.429 | 27.797 |
| 8 | `solveset((-2)**x + 8, x, S.Reals)` | 1.59 | 0.77 | 27.21 | 0.544 | 0.255 | 1.570 |
| 9 | `solveset(3*(-2)**x - 12, x, S.Reals)` | 1.70 | 0.85 | 27.51 | 0.635 | 0.339 | 1.707 |

## Results

- after and SymPy return the same 10 results; the `unrad` cases were compared
  up to the Dummy name.
- before differed on the three `unrad` cases:
  - case 0: the resultant elimination returned `x*(25*x**2 - 40*x + 16)`
    where SymPy returns `x*(5*x - 4)`;
  - cases 1 and 2: it returned a high-degree polynomial in `x` where SymPy
    returns the change of variable `(eq(_p), [_p, _p**k - ...])`.
- The solve and solveset results were the same before and after.
- Speed:
  - `unrad` is 2 to 4 times faster than before and 22 to 50 times faster than
    SymPy on the first call;
  - the solveset radical and negative-base steps are about 2 to 3 times
    faster than before;
  - the `solve` cases are unchanged within noise. Their time is dominated by
    `_unrad_solution` (Python) and the polynomial solve.
