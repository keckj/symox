# Transcendental fallbacks of `solve`, 2026-11-06 (`feat-tsolve-20261106`)

Before this lot, every equation that the native `_core.solve` left with
NotImplementedError went through Python: `_solve_piecewise`, then
`solvers._tsolve_steps.tsolve_step` (the Rust `_tsolve` steps driven by a
PyHost whose nested `solve` calls re-entered the public Python `solve`),
then the `_masked_solution`, `_inverted_solution` and `_unrad_solution`
fallbacks of `_solve_relations` (with `_unrad_resultant`).

The SymPy 1.14 order is now native, inside one native scope
(`crates/symox-core/src/solve/native_tsolve.rs`):

1. the unrad fallback of `_solve` for radicals of the unknown mixed with
   other functions of it (`radical::solve_mixed`), before `_tsolve`;
2. the `_tsolve` steps (`solve/tsolve.rs`): `_invert`, the `Add` factor and
   logcombine branches, the `Pow` branch (`f(x)**g(x) = c` with its
   rational candidates), `sin`/`cos`/`LambertW`, `rebuild(rewrite(exp))`,
   the bivariate block (`_solve_lambert`, then `bivariate_type` with
   `_tsolve(g - u, bivariate=False)` and `_vsolve(p, u)`), and the posify
   retry (flags with `bivariate` popped, `force=False`);
3. the nested `_vsolve` calls of these steps stay native, with the
   reference `_vsolve` semantics: `checksol` only, no filtering by the
   unknown's assumptions (`solve::without_assumption_filter`).

Native `_solve` also gained the multi-generator step of the reference
(`solve.rs::single_function`: `_as_base_q`, then the substitution of the
shallowest function by a Dummy), and its logcombine step now combines only
the terms that depend on the unknown (`log(-pi/2)` keeps its `I*pi`).
Its polynomial step takes the primitive part of a quadratic with symbolic
coefficients before the closed form, as the reference `roots`
(`preprocess_roots`) does (`solve.rs::quadratic_primitive`; the dsolve row
`1st_exact_09` reaches `x**2*y**2 - 2*x*y - x**2*exp(2*C1)` through the
logcombine step).

Anything the native host cannot reproduce exactly declines, and the binding
reruns the previous native solver and leaves the rest to the Python path:
a decline never becomes `[]` or a partial answer. Declines: nested
equations the Python wrapper rewrites first (Floats, hyperbolic functions,
`Abs`, `Piecewise`, held derivatives, integrals, sums, undefined
functions); a nested radical equation left unsolved when `unrad` raises
NotImplementedError (the established `_unrad_resultant` fallback stays
reachable); `nsimplify` of a number that needs `mpmath.identify` and an
undecided `equals`, both taken from Python through hooks (counted below as
`nsimplify`/`equals`). `SYMOX_TSOLVE_TRACE=1` prints the decline reasons
to stderr.

## Corpus

`scripts/cases.py`: 115 equations in five families (`mix_*` radicals mixed
with functions, `sym_*` symbolic exponents, `ts_*` the `_tsolve` tests of
SymPy's `test_solvers.py`, `lw_*` the Lambert tests, `flag_*` the flags
`check=False`, `dict`, `simplify=False`, `_unrad=False`, `Eq`,
`force=False`). `scripts/run_cases.py BACKEND OUT.json --trace` replays
them, one forked child per case, and records which Python fallbacks ran;
`--corpus unrad` replays the 90 radical equations of
`../2026-11-06-unrad`. `scripts/compare.py SYMPY.json SYMOX.json`
classifies the results.

| replay | SAME | EQUIV | SUBSET | SUPERSET | WRONG_EMPTY | DECLINED | EXTRA |
|---|---|---|---|---|---|---|---|
| before (main `e4b1d011`) | 99 | 6 | 3 | 1 | 1 | 5 | 0 |
| after | 106 | 6 | 0 | 2 | 0 | 1 | 0 |

Python computation per case (`python_fallbacks` in the raw files):

| replay | cases through the Python `_tsolve`/fallback chain | cases using the `nsimplify` hook |
|---|---|---|
| before | 66 | - |
| after | 10 | 12 |

The 10 remaining cases:

- `mix_cosh_over_x`, `mix_two_sqrt_exp`, `ts_loop`, `ts_sinh_ni`,
  `lw_sin_ni`: SymPy raises NotImplementedError, symox too. The Python
  chain still runs after the native NotImplementedError and also fails
  (the established fallbacks are kept).
- `sym_x_x4`, `sym_mx`, `sym_x_half`, `sym_x2p1`: the solution is native;
  the Python wrapper's post-processing (`_high_order_roots` /
  `_polynomial_roots`) asks `solve(gen - t, x)`, which raises in native
  code and then in the Python chain, without effect on the result.
- `flag_sym_noforce`: `force=False` is not accepted by the native binding.

Radical corpus of the unrad lot: before 87 SAME, 2 EQUIV, 1 EXTRA; after
87 SAME, 2 EQUIV, 1 EXTRA (unchanged; the 3 cases going through Python are
`posify_z`, `cov_fail`, `cosh_fail`, instead of 4 before).

dsolve corpus (`../2026-11-05-ode`): 315/289/270 before and after, with no
row changed (`raw/ode-before-after.score`, baseline
`~/symox-ci/evidence/integ41-20261106/ode`, after
`~/symox-ci/evidence/tsolve-20261106-ode/ode-after3`).

Files: `raw/sympy.json` (SymPy 1.14 in `~/symox-ci/venv`),
`raw/symox-before.json`, `raw/symox-after.json`,
`raw/unrad-symox-before.json`, `raw/unrad-symox-after.json`,
`raw/compare-before.txt`, `raw/compare-after.txt`, `raw/ode-before-after.score`.

## Remaining differences

- `lw_x3_3x` (`x**3 - 3**x`) and `lw_cos` (`3**cos(x) - cos(x)**3`):
  SUPERSET. The native `lambert::solve` of `solve_rec`, which runs before
  the `_tsolve` steps, also returns the two complex branches
  `-3*LambertW((1 ± sqrt(3)*I)*log(3)/6)/log(3)`. They are valid
  solutions. This predates the lot for `lw_x3_3x`; `lw_cos` now goes
  through the same step after its substitution of `cos(x)`.
- `flag_sym_noforce`: DECLINED (NotImplementedError for `force=False`).
- EQUIV (same values, different forms): `ts_2x10` (`log(10)/log(2)`
  against `1 + log(5)/log(2)`) and five others, all predating the lot.
