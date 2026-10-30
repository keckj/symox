# prdefrac (wave 30): risch_integrate through the parametric Risch DE

Harness: `bench.py` (one fresh process per call: first call, then the median
of 3 repeated calls; prints the sha1 prefixes of `srepr(result)` and
`str(result)`), `run.sh <label> <module> <PYTHONPATH|-> <reps>` (every
process under `systemd-run -p MemoryMax=6G -p MemorySwapMax=0` and
`timeout 120`, run from `/tmp`), `summarize.py raw.txt` (medians). Raw
numbers: `raw.txt`. Labels: `before` = base tree 28ef4e1b with its release
`_core` build, `after` = this branch (release build), `sympy` = SymPy 1.14 in
the same venv. The three configurations were interleaved, 5 rounds.

## (1) Profiles before the change

cProfile (`prof.py`, one repeated call after a warm-up) and py-spy
(`py-spy record -r 997 --format raw` on a loop of calls), under the caps.
"Native" below is the self time of `symox._core` builtins; py-spy without
`--native` attributes native time to the calling Python frame.

`risch_integrate(log(x)*log(x + 1), x)` (cProfile total 3.06 s, wall
about 1.6 s):

- Python frames 2.17 s self, other builtins (isinstance, hash, getattr, ABC
  checks) 0.72 s, symox native 0.15 s (of which `Poly._from_normalized`
  0.12 s). The computation itself is almost entirely native; the time is
  the Python wrapper around it.
- 37103 `Poly(e, t, domain=K)` constructions (`_poly_constructor.construct`
  2.36 s cumulative, `_options` 1.26 s) made by `PolyMatrix` entry reads:
  `_poly`/`__getitem__` (1.30 s), the ground check of `rref`/`nullspace`
  (`_ground_check` iterates all entries, 0.85 s), `from_list` unify chains
  (1.25 s) and `set_gens` (0.40 s). py-spy: 75 % of the samples inside
  `construct`, 42 % inside `_options`; the self time is `_options` 25 %,
  `construct` 18 %, `sympify` 8 %. Nearly all of these entries are exact
  Integer/Rational constants over QQ or `ZZ(x)`/`QQ(x)`.
- `prde_cancel_liouvillian` is 1.76 s cumulative, but only through those
  matrix reads (its `nullspace` calls are 0.88 s).

`risch_integrate(log(log(x)), x)` (cProfile total 0.38 s, wall about 0.22 s):

- Python frames 0.27 s self, other builtins 0.06 s, symox native 0.016 s.
- 0.114 s self in `_risch_splitfactor_native._plain` (6480 calls): the
  admission scans the 557 keys of `_core.__dict__` in a Python loop on every
  function-metadata read; `_polymatrix_native._unchanged`/`_class_unchanged`
  reach the same scan (0.09 s).
- 1913 `Poly(...)` constructions (0.13 s), as above.

The kernels named by the task were small in both calls: `constant_system`
4 ms (6 and 7 calls), `derivation` 18-36 ms (44 and 88 calls), `is_deriv_k`
6 ms, `is_log_deriv_k_t_radical_in_field` 24 ms.

## (2) Changes (computation location)

- Rust (`crates/symox-py/src/prde_frac.rs`): `_str_keyed(d)`, the "exact
  dict whose keys are exact str" check of the admissions, via
  `PyDict_Next` (no hashing, equality or callbacks); and
  `_function_states_unchanged(records)`, the code/defaults/kwdefaults/
  closure identity check of `_polymatrix_native._unchanged`, through the
  same public CPython getters as `_ddm_function_metadata` (no
  `object.__getattr__` audit event). Both are admission checks; the Python
  loops remain the fallback when the native helpers are absent.
- Admission only (`python/symox/polys/_polymatrix_rational.py`): a
  `PolyMatrix` exact Integer/Rational entry over QQ or a fraction field
  `ZZ(a, ...)`/`QQ(a, ...)` of stock symbols disjoint from its single stock
  generator is built by calling the native constructor
  `Poly._from_normalized(e, gens, domain, None, None)` directly, which is
  exactly what stock `construct` reduces to for that input; the ground check
  of such a matrix is proven without building the entries; `from_list` of
  native `Poly` items over one QQ/ZZ domain object and one generator returns
  that domain without the `unify` chain. No arithmetic moved to Python.

After the change (cProfile): log(x)*log(x+1) 0.82 s total (Python frames
0.51 s, native 0.08 s), the remaining `construct` calls (3242) are
non-constant entries and `parallel_poly_from_expr(field=True)` in
`set_gens`; log(log(x)) 0.21 s total (Python frames 0.10 s), mostly the
splitfactor admission proof.

## (3) Measurements

Medians over 5 fresh processes per configuration and call, release builds,
interleaved (`summarize.py raw.txt`). Calls 0-9 are the test_risch.py /
test_prde.py integrands that go through the parametric Risch DE.

first call in a fresh process

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 36.9 | 26.1 | 0.71 | 62.7 | 0.42 | yes | yes | yes | 5 |
| 1 | `log(x)**2` | 38.2 | 28.4 | 0.74 | 68.3 | 0.42 | yes | yes | yes | 5 |
| 2 | `log(x)**3/x` | 36.6 | 27.5 | 0.75 | 44.2 | 0.62 | yes | yes | yes | 5 |
| 3 | `log(log(x))` | 219.8 | 115.1 | 0.52 | 123.6 | 0.93 | yes | NO | yes | 5 |
| 4 | `exp(x)*x` | 69.2 | 37.2 | 0.54 | 52.6 | 0.71 | yes | yes | yes | 5 |
| 5 | `x*exp(x)*log(x)` | 168.2 | 76.0 | 0.45 | 98.1 | 0.78 | yes | yes | yes | 5 |
| 6 | `1/(x*log(x))` | 36.3 | 26.9 | 0.74 | 44.4 | 0.61 | yes | yes | yes | 5 |
| 7 | `log(x)**2/x**2` | 39.3 | 29.3 | 0.75 | 83.9 | 0.35 | yes | yes | yes | 5 |
| 8 | `exp(x)/(1+exp(x))**2` | 34.2 | 24.6 | 0.72 | 10.3 | 2.38 | yes | yes | yes | 5 |
| 9 | `log(x)*log(x+1)` | 1714.0 | 400.2 | 0.23 | 677.2 | 0.59 | yes | NO | yes | 5 |

repeated call (median of 3 per process)

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 19.1 | 10.4 | 0.54 | 19.5 | 0.53 | yes | yes | yes | 5 |
| 1 | `log(x)**2` | 20.7 | 12.3 | 0.60 | 25.5 | 0.48 | yes | yes | yes | 5 |
| 2 | `log(x)**3/x` | 18.8 | 10.7 | 0.57 | 19.4 | 0.55 | yes | yes | yes | 5 |
| 3 | `log(log(x))` | 194.9 | 94.3 | 0.48 | 112.5 | 0.84 | yes | NO | yes | 5 |
| 4 | `exp(x)*x` | 51.6 | 20.8 | 0.40 | 13.5 | 1.54 | yes | yes | yes | 5 |
| 5 | `x*exp(x)*log(x)` | 150.7 | 59.2 | 0.39 | 53.9 | 1.10 | yes | yes | yes | 5 |
| 6 | `1/(x*log(x))` | 19.0 | 11.2 | 0.59 | 16.4 | 0.68 | yes | yes | yes | 5 |
| 7 | `log(x)**2/x**2` | 22.5 | 13.2 | 0.59 | 38.5 | 0.34 | yes | yes | yes | 5 |
| 8 | `exp(x)/(1+exp(x))**2` | 16.8 | 9.0 | 0.53 | 7.8 | 1.15 | yes | yes | yes | 5 |
| 9 | `log(x)*log(x+1)` | 1648.5 | 377.8 | 0.23 | 655.3 | 0.58 | yes | NO | yes | 5 |

Result: `log(log(x))` 220 -> 115 ms first call (SymPy 124), 195 -> 94 ms
repeated (SymPy 113); `log(x)*log(x + 1)` 1714 -> 400 ms first (SymPy
677), 1649 -> 378 ms repeated (SymPy 655). No call got slower than before.
Still slower than SymPy: `exp(x)/(1+exp(x))**2` first call (25 vs 10 ms)
and the repeated calls of `exp(x)*x` (21 vs 14 ms), `x*exp(x)*log(x)`
(59 vs 54 ms) and `exp(x)/(1+exp(x))**2` (9.0 vs 7.8 ms), all faster than
before. Every result has the same `srepr` as before and prints as
SymPy's (`str`); the two "NO" in "same srepr as SymPy" are the existing
`srepr` gap (a `NonElementaryIntegral` inside an `Add` prints as
`Integral` in symox).

## SymPy test files (run on symox through the SymPy test sources)

test_prde.py: 3 failed, 12 passed before and after the first three commits
(test_prde_special_denom, test_constant_system, test_prde_no_cancel); 1
failed, 14 passed after the follow-up (4) (test_prde_special_denom). test_risch.py: 5 failed, 31
passed before and after (test_gcdex_diophantine, test_laurent_series,
test_recognize_derivative, test_residue_reduce, test_risch_integrate_float).

## Validation

Integration corpus `{"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1,
"SAME": 1171}` and ODE corpus `{'cls': 303, 'hint': 279, 'def': 274, 'n':
315}`, both equal to the 28ef4e1b baselines. `cargo test -p symox
--release`: 793 + 1 unit tests pass, the 5 doctests pass when run with
`--test-threads=1` (2 of them fail to compile in parallel under the 6G cap,
as on the base). clippy correctness: clean. Full suite (24 shards on
demeter plus the local compiler-dependent files): all pass after reruns on demeter of
two load-dependent flakes: `test_dirac_affine_integration.py` (a 5 s
subprocess timeout at load 34) and `test_matrix_inexact_native_inverse.py`
(a 1 ms SIGALRM that fired before the `pytest.raises` block and aborted
the shard; the shard rerun passed 3826 of 3826).

## (4) Follow-up: domain inference of constant_system, remaining gaps

SymPy test_prde.py `test_constant_system` expects `B` and `v` over
`QQ(x)[t]`; symox gave `ZZ(x)[t]`. Causes and fixes:

- Rust, crates/symox-core/src/expand.rs `expand_mul_node`: SymPy's
  `Mul._eval_expand_mul` rebuilds `n/d` from `fraction(expr)` whatever `d`
  is; symox did so only for several denominator factors. In
  `constant_system`, `Rm1 = Ri.applyfunc(D(x)/DAij)` holds
  `2*(x**2/2 + x + 1/2)/(x**2 + 2*x + 1)` unevaluated; SymPy expands it to
  `1` (domain `QQ`), symox left it as a fraction (domain `QQ(x)`), so the
  matrix ring and then `B`, `v` differed. Now `QQ[t]` as in SymPy.
- Python plumbing, polymatrix.py: `PolyMatrix(..., ring=K[gens])` now
  takes `K` and `gens` from the ring (SymPy `from_dm`), and `__eq__` also
  compares the domain (SymPy `DomainMatrix.__eq__`).

Effect: test_constant_system and test_prde_no_cancel pass.

Gaps found, not changed: the `expand_mul` entry point (simp.rs
`expand_mul_deep`) does not rebuild `n/d`; `Mul(x + 1, 1/(x + 1))` held
unevaluated expands to `x/(x + 1) + 1/(x + 1)` (SymPy `1`); two
denominator factors such as `y*(x + 2)` are not expanded.

Not ported (gaps), with the reasons found:

- `constant_system`: its `rref`/`nullspace` over `QQ`, `ZZ(x)`, `QQ(x)`
  already run natively (`_polymatrix_ground`); the Python left is the
  expression-level `derivation` of each entry and the PolyMatrix
  `applyfunc`/subtraction, whose domains are inferred from the expression
  forms (as the bug above shows). A Rust port with identical domains needs
  `expand`/`cancel` parity of every intermediate first.
- `derivation` (exp/primitive), `is_deriv_k`,
  `is_log_deriv_k_t_radical(_in_field)`, `limited_integrate_reduce`
  outside the base case: not attempted in the remaining time (each item
  needs a release build of about 10 minutes plus the full validation).
- `parallel_poly_from_expr(field=True)` / `set_gens` in `log(x)*log(x + 1)`:
  the per-entry `Poly` construction is already native; the cost left is the
  Python option handling around it (`_poly_constructor._options`, 3242
  calls, 0.12 s of 0.82 s profiled) spread over many small calls, not one
  kernel.
- the three `exp` integrands still slower than SymPy: the profile of
  `exp(x)/(1 + exp(x))**2` (66 ms first call) shows 22 ms for the lazy
  first import of mpmath (domains/_methods.py `_inexact_value`, reached
  from `mul_ground`; SymPy imports mpmath at import time) and 26 ms in the
  splitfactor admission (`_risch_splitfactor_native.reduce`, with
  `_risch_native._prove_fields` 14 ms); `exp(x)*x` repeated is also
  dominated by that admission. Neither is on the prde path ported here.
  Skipping the mpmath import in `_inexact_value` for exact elements was
  tried and reverted: the native LLL/SDM cold-start proofs fingerprint that
  function (tests/test_native_sdm_berkowitz.py failed).

Measured after (4) (raw-followup.txt, 3 interleaved fresh processes,
medians, while a one-job cargo build ran, so both sides are slower than in
(3)): `log(log(x))` 140/110 ms first/repeated (SymPy 141/127),
`log(x)*log(x + 1)` 452/425 ms (SymPy 783/742), `exp(x)*x` 44/24 ms (SymPy
62/16), `x*exp(x)*log(x)` 89/71 ms (SymPy 117/65), `exp(x)/(1 + exp(x))**2`
30/10 ms (SymPy 12/9). The result hashes (srepr and str) equal those of (3).

Validation of (4): 48 shards on demeter, all exit 0, 57818 passed, plus 72
local = 57890 passed, no failure; integration corpus `{"COND": 20, "DIFF":
10, "EQUIV": 8, "MISSING": 1, "SAME": 1171}` and ODE corpus `{'cls': 303,
'hint': 279, 'def': 274, 'n': 315}`, both equal to the baselines.
`cargo test -p symox --release`: 793 + 1 unit tests pass, the 5 doctests
pass with `--test-threads=1` (2 fail to compile in parallel under the 6G
cap, as before); clippy correctness clean. SymPy test_risch.py: 5 failed,
31 passed, same tests as before.
