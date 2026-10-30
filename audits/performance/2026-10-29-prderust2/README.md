# prderust2 (wave 28): risch_integrate before/after and against SymPy 1.14

Harness: `bench.py` times one `risch_integrate` call in a fresh interpreter;
`run.sh <label> <module> <PYTHONPATH|-> <reps>` loops over the 10
integrands. The three configurations were interleaved, 5 repetitions each
(`raw.txt`: label, index, seconds, sha1 prefix of `srepr(result)`,
integrand), and `summarize.py raw.txt` gives the medians below.

- before: the tree at 581f7bb1 (`git archive HEAD python`) with the release
  `_core.abi3.so` built from that tree;
- after: this branch, release build from the same tree plus the change;
- SymPy: SymPy 1.14 in /home/keckj/Documents/symox/.venv.

The machine was shared with other builds (8 cores), so differences under
about 3% are noise.

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | n |
|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 39.2 | 40.2 | 1.02 | 67.1 | 0.60 | 5 |
| 1 | `log(x)**2` | 41.9 | 41.8 | 1.00 | 71.9 | 0.58 | 5 |
| 2 | `log(x)**3/x` | 39.0 | 39.4 | 1.01 | 48.9 | 0.81 | 5 |
| 3 | `log(log(x))` | error | 234.4 | - | 137.9 | 1.70 | 5 |
| 4 | `exp(x)*x` | 74.4 | 74.8 | 1.00 | 57.4 | 1.30 | 5 |
| 5 | `x*exp(x)*log(x)` | 178.6 | 180.3 | 1.01 | 107.4 | 1.68 | 5 |
| 6 | `1/(x*log(x))` | 39.2 | 40.5 | 1.03 | 47.9 | 0.84 | 5 |
| 7 | `log(x)**2/x**2` | 43.8 | 43.9 | 1.00 | 89.1 | 0.49 | 5 |
| 8 | `exp(x)/(1+exp(x))**2` | 36.8 | 36.5 | 0.99 | 11.4 | 3.21 | 5 |
| 9 | `log(x)*log(x+1)` | 1739.0 | 1747.2 | 1.00 | 722.3 | 2.42 | 5 |

Reading:

- `log(log(x))` raised `CoercionFailed` before; it now returns
  `x*log(log(x)) + NonElementaryIntegral(-1/log(x), x)`, like SymPy. The
  fix is the `field=True` joint domain in `parallel_poly_from_expr`
  (Python wrapper plumbing; the domain itself is built by the existing
  `construct_domain`). The result hash differs from SymPy only because
  symox `srepr` prints the inner `NonElementaryIntegral` of an `Add` as
  `Integral` (printer gap, the object class is `NonElementaryIntegral`).
- Every other result hash is identical before and after; all but
  `log(x)*log(x+1)` (a different but pre-existing form) are identical to
  SymPy.
- The new native paths (`constant_system` outside the base case and
  `limited_integrate_reduce`) do not change these timings: on these
  integrands only `log(x)*log(x+1)` reaches the primitive-case
  `constant_system` with rational entries (2 calls), and
  `limited_integrate_reduce` is not called by `risch_integrate` in SymPy
  1.14. The time of `log(x)*log(x+1)` is dominated by Python PolyMatrix
  overhead (`_ground_check`, `__getitem__`, `from_list`: about 2.4 s of
  Poly constructions in a profile) and `prde_cancel_liouvillian`, which
  are still Python.

`import symox` (`import-raw.txt`, 3 rounds of 5 fresh processes, medians):
before 207.8 / 206.4 / 208.8 ms, after 202.0 / 205.8 / 206.2 ms; the
import path is unchanged by this branch. With the capture call removed the
median import was 196.8 ms against 198.2 ms with it (15 processes each),
so the wave-27 capture now costs about 1.4 ms, not 6.7 ms. It was kept
eager: a lazy capture at the first prde call could not prove that prde,
rde, risch and PolyMatrix were not patched between `import symox` and that
call without per-module import hooks.

Direct calls of the new native entry points (`kernels.py`, same process,
7 rounds of 200 calls, median per call, `kernels-raw.txt`):

| call | native | Python body | ratio |
|---|---|---|---|
| `limited_integrate_reduce`, base case, 2 terms | 327.8 us | 34805.7 us | 0.009 |
| `constant_system`, primitive case, 3x3 rational | 168.5 us | 1749.1 us | 0.096 |
