# tsolverust (wave 27): native `_tsolve` steps and `bivariate`

Ten public `solve` calls that go through the `_tsolve` steps, measured with
release builds:

- before: base 31307c81 (Python `_tsolve_steps.py` / `bivariate.py`);
- after-check1: the first version of this branch (native
  `symox::solve::tsolve`, SymPy `checksol` filtering of the candidates);
- after: the final version, where the numerical part of `checksol` is
  cheaper (see below);
- SymPy 1.14.0, run from /tmp with the project venv.

## Harness

- `bench.py <symox|sympy> <first|repeat> [reps] [only]`: prints one JSON line
  per case with the time in ms and the result.
- `run.sh <label> <PYTHONPATH|sympy> <outdir>`:
  - first call: each case runs in a fresh process, 5 runs; the reported value
    is the median;
  - repeated call: one process, 3 runs of a median of 20 calls; the reported
    value is the median of the 3.
  - Every process runs under `timeout` and
    `systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0`.
- `summarize.py raw/` writes `summary.md`.
- Raw numbers: `raw/{before,after-check1,after,sympy}-{first,repeat}.jsonl`;
  `python summarize.py before after-check1 after sympy`.

## Results (ms)

| # | call | before first (ms) | after-check1 first (ms) | after first (ms) | sympy first (ms) | before repeat (ms) | after-check1 repeat (ms) | after repeat (ms) | sympy repeat (ms) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | `solve(x**y - 2, x)` | 6.14 | 6.55 | 5.53 | 165.15 | 4.650 | 5.660 | 4.456 | 31.937 |
| 1 | `solve(-a*x + 2*x*log(x), x)` | 4.01 | 4.52 | 3.94 | 224.21 | 2.629 | 3.750 | 2.979 | 56.576 |
| 2 | `solve(x**(2*y) - 3, x)` | 6.73 | 7.18 | 6.22 | 208.11 | 5.078 | 6.185 | 4.938 | 51.069 |
| 3 | `solve((x**2 - 1)**y - 4, x)` | 11.34 | 12.48 | 10.78 | 235.94 | 9.453 | 11.226 | 9.154 | 82.656 |
| 4 | `solve(x*exp(x**2) - 1, x)` | 53.49 | 215.20 | 32.43 | 541.23 | 33.235 | 143.762 | 16.265 | 321.215 |
| 5 | `solve(sqrt(x)*exp(sqrt(x)) - 1, x)` | 42.37 | 159.31 | 31.12 | 260.90 | 22.183 | 89.278 | 14.973 | 97.915 |
| 6 | `solve(x - exp(-x**2), x)` | 53.66 | 214.61 | 32.30 | 595.72 | 33.203 | 143.136 | 16.360 | 328.366 |
| 7 | `solve(x**a*exp(x) - 1, x)` | 7.27 | 5.70 | 4.93 | 270.10 | 5.661 | 4.930 | 3.871 | 59.993 |
| 8 | `solve(x**a - exp(x), x)` | 13.11 | 9.51 | 8.47 | 291.45 | 10.782 | 8.741 | 7.089 | 81.168 |
| 9 | `solve(sqrt(x) + log(x), x)` | 41.34 | 89.91 | 30.23 | 287.70 | 21.692 | 49.986 | 14.209 | 113.431 |

## Reading

- **Results.** The results are identical in all versions except case 1:
  `solve(-a*x + 2*x*log(x), x)` no longer returns `0`. SymPy also rejects
  that candidate (`log(0)` gives `zoo`), so the output now matches SymPy.
- **Case 2.** SymPy also returns `(-sqrt(3))**(1/y)`; neither the base nor
  this branch does. This is a recorded gap in docs/algorithmic-parity.md.
- **Why after-check1 was slow (cases 4, 5, 6, 9).**
  - It was 2 to 4 times slower than the base. A native profile (py-spy
    `--native`, release build) of case 4 put 1386 of the 1582 samples
    spent in tsolve inside `checksol` -> `evalf` -> the LambertW evaluator.
  - The cause: the numerical test `abs(val.n(18).n(12, chop=True)) < 1e-9`
    first evaluated the residual at 18 digits without chop.
  - A correct candidate's residual is an exact zero, so that evaluation
    raised the working precision up to the 4096-bit cancellation budget,
    re-evaluating LambertW at every step: 50 to 117 ms per candidate.
- **The fix.**
  - The final version evaluates the residual at 18 digits directly under
    the 12-digit chop cutoff. A part that is chopped is proven to be below
    the cutoff, far below 1e-9; a part that is not chopped is resolved to
    18 digits. The absolute test therefore takes the same decision.
  - When the direct evaluation does not give a number, the former two-step
    evaluation runs, so the undecided answers are unchanged.
  - The 40-case candidate corpus gives the same results as after-check1.
  - Cases 4, 5, 6 and 9 are now below the base.
- **Machine load.** The before and after-check1 runs shared the machine
  with other agents and a `cargo test` build; the final after run did not.
  Part of the small-case differences is that load.
