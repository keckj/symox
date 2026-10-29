# tsolverust (wave 27): native `_tsolve` steps and `bivariate`

Ten public `solve` calls that go through the `_tsolve` steps, measured with
release builds:

- before: base 31307c81 (Python `_tsolve_steps.py` / `bivariate.py`);
- after: this branch (native `symox::solve::tsolve`, SymPy `checksol`
  filtering of the candidates);
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
- Raw numbers: `raw/{before,after,sympy}-{first,repeat}.jsonl`.

## Results (ms)

| # | call | before first | after first | SymPy first | before repeat | after repeat | SymPy repeat |
|---|---|---|---|---|---|---|---|
| 0 | `solve(x**y - 2, x)` | 6.14 | 6.55 | 165.15 | 4.650 | 5.660 | 31.937 |
| 1 | `solve(-a*x + 2*x*log(x), x)` | 4.01 | 4.52 | 224.21 | 2.629 | 3.750 | 56.576 |
| 2 | `solve(x**(2*y) - 3, x)` | 6.73 | 7.18 | 208.11 | 5.078 | 6.185 | 51.069 |
| 3 | `solve((x**2 - 1)**y - 4, x)` | 11.34 | 12.48 | 235.94 | 9.453 | 11.226 | 82.656 |
| 4 | `solve(x*exp(x**2) - 1, x)` | 53.49 | 215.20 | 541.23 | 33.235 | 143.762 | 321.215 |
| 5 | `solve(sqrt(x)*exp(sqrt(x)) - 1, x)` | 42.37 | 159.31 | 260.90 | 22.183 | 89.278 | 97.915 |
| 6 | `solve(x - exp(-x**2), x)` | 53.66 | 214.61 | 595.72 | 33.203 | 143.136 | 328.366 |
| 7 | `solve(x**a*exp(x) - 1, x)` | 7.27 | 5.70 | 270.10 | 5.661 | 4.930 | 59.993 |
| 8 | `solve(x**a - exp(x), x)` | 13.11 | 9.51 | 291.45 | 10.782 | 8.741 | 81.168 |
| 9 | `solve(sqrt(x) + log(x), x)` | 41.34 | 89.91 | 287.70 | 21.692 | 49.986 | 113.431 |

## Reading

- **Results.** The results before and after are identical except case 1:
  `solve(-a*x + 2*x*log(x), x)` no longer returns `0`. SymPy also rejects
  that candidate (`log(0)` gives `zoo`), so the output now matches SymPy.
- **Case 2.** SymPy also returns `(-sqrt(3))**(1/y)`; neither the base nor
  this branch does. This is a recorded gap in docs/algorithmic-parity.md.
- **Slower cases (4, 5, 6, 9).** These are slower because the candidates now
  go through the SymPy `checksol` filtering, which the base skipped and SymPy
  performs.
  - A profile of case 4 puts 121 ms of 141 ms inside the native call.
  - The same call with `check=False` takes 17 ms.
  - So the cost is in the native `checksol`; the Python callbacks are not
    the cause. The likely parts are the `simplify` attempt, the 18-digit
    evalf of the LambertW values, and the check that runs at every
    recursion level. They were not isolated further; this is a follow-up.
  - These cases remain faster than SymPy.
- **Machine load.** The runs shared the machine with other agents and a
  `cargo test` build, so differences of about 1 ms in the small cases are
  noise.
