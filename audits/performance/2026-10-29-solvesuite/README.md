# solvesuite: before/after timings

`bench.py` times public `solve` calls on the paths changed by solvesuite,
plus unchanged controls.

- first: median of 5 first calls, each in a fresh process.
- repeat: median of 15 calls with the cache cleared between calls.

- before: base 3d7da725 (`PYTHONPATH=<base tree>/python python bench.py before`)
- after: feat-solvesuite-20261029g (`PYTHONPATH=python python bench.py after`)

Both ran on the same machine under `systemd-run` with MemoryMax=6G. The
load average was about 4 (other worktrees were building), so differences
under about 5% are noise. Raw numbers: before.json, after.json.

| case | first before (ms) | first after (ms) | repeat before (ms) | repeat after (ms) | ratio after/before | result |
|---|---|---|---|---|---|---|
| `solve(Eq(x**2 - 1/(x**2 - 4), 4 - 1/(x**2 - 4)), x)` | 1.266 | 1.589 | 0.827 | 1.073 | 1.30 | `[-2, 2]` -> `[]` |
| `solve(Eq(x + 1/x, 1/x), x)` | 0.881 | 1.159 | 0.530 | 0.732 | 1.38 | `[0]` -> `[]` |
| `solve(Eq(x**2, 4), x)` | 1.026 | 1.019 | 0.586 | 0.589 | 1.00 | same |
| `solve(Eq(1/x + x, 2), x)` | 1.488 | 1.488 | 0.966 | 0.997 | 1.03 | same |
| `solve(x**2 - 4, x)` | 1.003 | 1.003 | 0.552 | 0.549 | 0.99 | same |
| `solve(-a*x + 2*x*log(x), x)` | 3.553 | 3.330 | 2.414 | 2.414 | 1.00 | same |
| `solve(1 - oo*x, x)` | 0.937 | 0.951 | 0.553 | 0.544 | 0.98 | `[nan]` -> `[]` |
| `solve(x - 3, x)` | 0.828 | 0.809 | 0.441 | 0.436 | 0.99 | same |
| `solve([x*y - x - y, x*y - x - y], [x, y])` | 1.078 | 1.086 | 0.736 | 0.701 | 0.95 | same |
| `solve([x**2 + y**2 - 1, x - y], [x, y])` | 1.029 | 1.032 | 0.630 | 0.622 | 0.99 | same |
| `solve((7**x)**x - 2, x)` | 17.862 | 30.178 | 16.227 | 27.860 | 1.72 | `[]` -> `[-sqrt(log(2))/sqrt(log(7)), sqrt(log(2))/sqrt(log(7))]` |
| `solve(x**x - 4, x)` | 5.311 | 5.292 | 4.381 | 4.393 | 1.00 | same |
| `solve(exp(x) - 2, x)` | 4.025 | 3.997 | 3.329 | 3.306 | 0.99 | same |
| `solve(sqrt(A**2 + B**2) - 3, A)` | 13.028 | 8.130 | 12.057 | 7.352 | 0.61 | `[-sqrt(9 - b**2), sqrt(9 - b**2)]` -> `[]` |

Slower cases, all of which now give a different, SymPy-equal result:

- The two `Eq` cases with denominators that cancel in `lhs - rhs` are
  30-38% slower: about 0.2 ms more for `_simple_dens` on both sides plus
  the denominator check. Before, they returned wrong solutions. An `Eq`
  without denominators is unchanged.
- `solve((7**x)**x - 2, x)` takes 1.7x as long: the forced `expand_log`
  retry plus the check against the original equation. Before, it returned
  `[]`, which is wrong.

`solve(sqrt(A**2 + B**2) - 3, A)` with imaginary `A, B` is 1.6x faster. The
assumptions now prove `B**2 < 0`, so the radicand is positive.

Unchanged: the SymPy test_issue_12476 system takes 136.7 s before and
142.0 s after (pytest junit time, same load). A withdrawn variant of the
substitution solver, described in docs/algorithmic-parity.md, made it take
more than 7.5 min, so it was not kept.
