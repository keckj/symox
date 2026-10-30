# ratintgen (wave 32): native ratint_ratpart / ratint_logpart outside QQ/ZZ

Harness: `bench.py sympy|symox OUT.json`. SymPy 1.14 runs from /tmp; symox
runs from the worktree root with `PYTHONPATH=python`. For each call, "fresh"
is the first call in a new process (imports excluded) and "repeat" is the
median of 7 calls in one process after a warm-up call. All runs use the
release build and the 6G memory cap.

- `before.json`: the `_core.abi3.so` built from 7e8709e3 (the base of the
  branch, same Python tree).
- `after.json`: the binary of this branch.
- `sympy.json`: the reference SymPy 1.14.

Times are in ms. "b/a" is the repeat time before divided by after. The
results are identical before and after on all 12 calls. Where they differ
from SymPy, the difference was already there in the base binary (see the
gaps in docs/algorithmic-parity.md).

| call | SymPy fresh | SymPy repeat | before fresh | before repeat | after fresh | after repeat | b/a | str equal to SymPy |
|---|---|---|---|---|---|---|---|---|
| `integrate(1/(x**2+I)**2, x)` | 71.5 | 32.02 | 12.0 | 10.19 | 10.2 | 7.90 | 1.29 | yes |
| `integrate(I/(x**2+1)**2, x)` | 57.5 | 21.46 | 0.3 | 0.16 | 0.5 | 0.17 | 0.94 | no (same as before) |
| `integrate((x+I)/(x**3-I), x)` | 55.4 | 18.73 | 2.4 | 1.43 | 2.6 | 1.50 | 0.95 | no (same as before) |
| `integrate(1/(x**2+2*I*x-5), x)` | 51.2 | 16.91 | 0.3 | 0.16 | 0.4 | 0.16 | 0.99 | yes |
| `integrate(1/(x**2+a)**2, x)` | 98.5 | 47.53 | 11.0 | 8.93 | 9.2 | 7.05 | 1.27 | yes |
| `integrate(1/((x-a)*(x-b)), x)` | 158.5 | 82.41 | 14.1 | 11.91 | 13.8 | 11.27 | 1.06 | yes |
| `integrate(x/(x**2+a*x+b)**2, x)` | 291.3 | 178.86 | 24.6 | 22.17 | 21.1 | 19.54 | 1.13 | no (same as before) |
| `integrate(a/(x-a)**3, x)` | 73.5 | 35.95 | 0.3 | 0.13 | 0.3 | 0.13 | 0.99 | no (same as before) |
| `integrate(1/(x-sqrt(2))**2/(x+1), x)` | 844.5 | 567.67 | 165.4 | 142.79 | 167.9 | 125.70 | 1.14 | no (same as before) |
| `integrate(1/(x**2+sqrt(2))**2, x)` | 263.7 | 119.77 | 0.4 | 0.23 | 0.5 | 0.24 | 0.94 | yes |
| `ratint(1/(x**2+a)**3, x)` | 117.7 | 69.83 | 9.6 | 7.91 | 5.6 | 4.69 | 1.69 | yes |
| `ratint((x**2+1)/(x**2+I)**3, x)` | 62.2 | 39.38 | 15.3 | 14.02 | 8.9 | 9.42 | 1.49 | yes |

Calls at about 0.1-0.5 ms never reach `ratint_ratpart`/`ratint_logpart`
over the new domains: for example `I/(x**2+1)**2` pulls the constant out
and integrates over QQ. Their 0.94-0.99 ratios are noise, not slowdowns.
The calls that use the new kernels gain 1.06x to 1.69x. The rest of their
time is spent in the Python `ratint` body, `log_to_real`/`log_to_atan`, the
RootSum and `cancel`.

Kernel-level comparisons (native kernel against the reference Python body
on the same symox process, median of 5 calls):

- ratint_ratpart: QQ_I/ZZ_I 0.5-1.2 ms against 2.8-6.4 ms. QQ(a)/QQ(a,b)
  0.7-1.7 ms against 2.4-3.7 ms. EX 0.4-5.4 ms against 1.4-9.3 ms. Over
  QQ<sqrt(2)> the reference body raises TypeError and the kernel returns
  the SymPy 1.14 value.
- ratint_logpart: 1.0-6.0 ms against 1.4-6.5 ms over ZZ_I, QQ_I, ZZ[a],
  ZZ[a,b], QQ(a) and EX (10-30 percent). The resultant, sqf_list, gcd,
  invert and rem calls were already native, so the gain is the removed
  interpreter and constructor overhead.

## EX `roots` (item 4): attempted, not committed

The code is not in the branch. `roots-ex/roots-ex.patch` is the patch that
was measured (native Rust only, `crates/symox-core/src/poly/roots.rs` and
the `multiple=True` order in `crates/symox-py/src/poly.rs`). Its test file
is kept as `roots-ex/test_native_roots_expression.py.txt`.

Other files in `roots-ex/`:
- `rootscmp.py`: 25 `roots` cases, dict and `multiple=True`.
- Its outputs: SymPy 1.14 (`roots_sympy.txt`), symox at HEAD
  (`roots_before.txt`) and symox with the patch (`roots_symox.txt`).
- `bench_roots.py`: median of 5 calls, in ms. Results in `before.json`
  (HEAD), `after.json` (patch) and `sympy.json`.

On the 25 cases, dict and list both equal SymPy 1.14 in 6 cases at HEAD and
in 22 with the patch. Three cases still differ:
- the two `C(v)` quartics of `1/(x**5+1)`;
- `x**4 + 2*x**2 + 3`, where SymPy gives its trigonometric form (this one
  already differed at HEAD).

| call | HEAD | patch | SymPy 1.14 |
|---|---:|---:|---:|
| `roots(sqrt(3)*x**2 - 2*x + sqrt(2))` | 0.37 | 1.18 | 12.58 |
| `roots(x**4 - sqrt(3)*x**2 - 1)` | 0.23 | 0.55 | 16.14 |
| `roots(sqrt(2)*x**4 + 3*x**2 + sqrt(5))` | 0.62 | 1.51 | 32.57 |
| `roots(x**4 - 5*x**2 + 6)` | 0.12 | 0.11 | 1.84 |
| `integrate(1/(x**3+x+1))` | 185.0 | 336.1 | 7234.6 (one call) |
| `integrate(1/(x**5+1))` | 21.9 | 35.9 | 880.7 |
| `integrate(x/(x**5+1))` | 30.1 | 49.8 | 930.7 |
| `integrate(1/(x**3-2))` | 8.6 | 9.5 | 47.5 |

The full suite with the patch gives 57991 passed and 3 failed:
- `tests/test_native_ratint_dom.py::test_cubic_log_to_real_does_not_grow[1/((x**2+1)*(x**3+2))-9057]`:
  the result has 17199 characters, against a limit of 9057 (SymPy).
- `tests/test_nonlin_hang_parity.py:104`: `nonlinsolve([x**3+y-2, x+y])`
  prints differently.
- `tests/test_solveset_complex_wave4_parity.py:114`:
  `solveset(sin(x)+2*cos(x), x, Complexes)` gives
  `atan(cos(atan(4/3)/2)/sin(atan(4/3)/2))` instead of `atan(2)`. The cause
  is the `(-c)/a` tree on Gaussian quadratics. A second version restricts
  those trees to EX; it was not validated.

The integration corpus with the patch equals the baseline. The patch was
therefore not committed.
