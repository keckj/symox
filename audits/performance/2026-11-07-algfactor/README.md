# algfactor: factorization and cancel with algebraic numbers (wave 45)

Lot `feat-algfactor-20261107`, based on main f7888825 (wave 44). Scope: the
gaps left by algpoly. These are natural `factor_list`/`sqf_list` whose
generators include algebraic numbers, `cancel`/`gcd`/`lcm`/`cofactors` with
`greedy=False`, `real_roots`/`all_roots` with `extension=True` on a
multivariate `Poly`, `CRootOf` extensions, and the exception classes of
`factor(..., extension=)` and of polynomial ring conversions.

## Harnesses

- `diffcorpus.py sympy|symox OUT.jsonl [cases...]` (249 cases, in nine groups):
  - natural factor lists (`(x - sqrt(2))*(x + sqrt(2))`, `(x - sqrt(2)*t)**2`,
    `(x - sqrt(2)*pi)*(x + sqrt(2)*pi)`, `frac=True`, explicit generators);
  - constant products (`(1 - sqrt(2))*(1 + sqrt(2))`, `polys=True`);
  - natural square-free lists;
  - Poly methods over the generators `(x, sqrt(2))`;
  - `extension=` / `gaussian=` (univariate and multivariate);
  - `cancel`/`gcd`/`lcm`/`cofactors` with `greedy=False`;
  - roots (`real_roots`, `all_roots`, `nroots`, `ground_roots`, `roots` with
    and without `extension`);
  - `CRootOf` extensions (`QQ.algebraic_field(CRootOf(...))`,
    `minimal_polynomial`, `primitive_element`, `to_number_field`,
    `field_isomorphism`, arithmetic);
  - exception classes.

  Each record keeps the printed result (or the exception) and a `value` that
  `compare.py` uses to tell a representation difference (DIFF) from a wrong
  value (WRONG). A case killed by the 60 s alarm is recorded as `CaseTimeout`.
- `compare.py ref.jsonl X.jsonl` classifies each case:
  - SAME;
  - DIFF;
  - DECLINED (NotImplementedError, or missing);
  - WRONG (the value differs from SymPy's, or symox answers where SymPy raises).
- `upstream.py sympy|symox OUT.jsonl` runs these SymPy 1.14 modules
  statement by statement (a failing statement does not stop its function):
  - `test_polytools.py` (the functions about factor, sqf, cancel, gcd, lcm,
    roots, extensions);
  - `test_factortools.py`;
  - `test_rootoftools.py`;
  - `numberfields/tests/test_minpoly.py`, `test_subfield.py`,
    `test_numbers.py` and `test_utilities.py`.

  For symox, the module text is rewritten to import symox. A SymPy import
  that symox lacks becomes a stub that raises NotImplementedError when used.
  Three statements hang in uninterruptible native code. They are not run
  for symox, and are recorded as failed (`HANG` in the script).

Raw files:
- `raw/sympy.jsonl`: SymPy 1.14 (the shared CI virtual environment);
- `raw/before.jsonl`: the wave-44 `.so` (main f7888825,
  `integ44-20261107` target);
- `raw/after.jsonl`: this branch;
- `raw/compare-*.txt`: the per-case tables;
- `raw/upstream-sympy.jsonl`, `raw/upstream-before.jsonl`,
  `raw/upstream-symox.jsonl`: the statement records.

## Results

| corpus | before | after |
|---|---|---|
| algfactor corpus (249) | SAME 179, DIFF 64, DECLINED 3, WRONG 3 | SAME 248, DIFF 1, DECLINED 0, WRONG 0 |
| upstream statements (1362; SymPy passes 1360) | 948 | 955 |
| `tests/test_algfactor_parity.py` (live comparisons) | — | 73 passed |

By group (algfactor corpus):

| group | cases | before SAME / DIFF / DECLINED / WRONG | after SAME / DIFF / DECLINED / WRONG |
|---|---:|---|---|
| natural factor lists | 42 | 21 / 19 / 2 / 0 | 42 / 0 / 0 / 0 |
| constant products | 11 | 3 / 8 / 0 / 0 | 11 / 0 / 0 / 0 |
| natural square-free lists | 19 | 10 / 9 / 0 / 0 | 19 / 0 / 0 / 0 |
| Poly methods | 5 | 5 / 0 / 0 / 0 | 5 / 0 / 0 / 0 |
| `extension=` / `gaussian=` | 58 | 56 / 2 / 0 / 0 | 58 / 0 / 0 / 0 |
| `greedy=False` cancel/gcd/lcm | 37 | 20 / 14 / 0 / 3 | 37 / 0 / 0 / 0 |
| roots | 29 | 25 / 3 / 1 / 0 | 29 / 0 / 0 / 0 |
| `CRootOf` extensions | 35 | 34 / 1 / 0 / 0 | 34 / 1 / 0 / 0 |
| exception classes | 13 | 5 / 8 / 0 / 0 | 13 / 0 / 0 / 0 |

The three WRONG cases before the lot were wrong values from
`gcd(x**2 - 2, x + sqrt(2), greedy=False)`, `lcm(...)` and
`gcd(x**2 - 2*t**2, x - sqrt(2)*t, greedy=False)`. The pre-lot code treated
`sqrt(2)` as a generator and returned 1, or the product. SymPy computes over
`EX`.

By module (upstream statements; SymPy-recorded statements):

| module | statements | SymPy | before | after |
|---|---:|---:|---:|---:|
| test_polytools | 518 | 517 | 457 | 464 |
| test_factortools | 216 | 216 | 6 | 6 |
| test_rootoftools | 205 | 204 | 160 | 160 |
| test_minpoly | 159 | 159 | 147 | 147 |
| test_subfield | 149 | 149 | 91 | 91 |
| test_numbers | 93 | 93 | 67 | 67 |
| test_utilities | 22 | 22 | 20 | 20 |
| total | 1362 | 1360 | 948 | 955 |

The run fixes these statements: `test_factor#113` (constant product
`factor_list(s2p*s2n)`), `#115`, `#121` (`sqrt(2)*pi` and natural
generators), `test_real_roots#19`, `test_all_roots#17` (`extension=True`
on a multivariate Poly), and `test_cancel#22`, `#25` (`greedy=False`).
There are no regressions.

Validation:

- suite: SUITE CLEAN on the first run, with 60 893 passed, 14 xfailed and
  45 skipped, out of 60 952 collected in 701 files. One more xfail than the
  baseline comes from the hash-seed-dependent, `strict=False` case
  `test_checkinfsol_matches_sympy[x*(f(x).diff(x)) + 1 - f(x)**2]`;
- `cargo test --release`: 876 + 1 + 0 + 5;
- clippy correctness: clean;
- ODE corpus: 315/289/270;
- the related test files: 9054 passed.

## Fixes (where the computation lives)

| fix | language | files |
|---|---|---|
| reconstruction certificate by expansion when the generator certificate cannot decide (algebraic numbers as generators) | Rust | `poly/checked_factors.rs` |
| natural `factor_list`/`sqf_list` componentwise, in public `Mul.args` order, square-free groups multiplied by multiplicity as `_symbolic_factor_list` does (global and fraction sides) | Rust | `poly/api.rs`, `symox-py/src/poly_checked_factors.rs` |
| `factor_list`/`sqf_list` of a constant product (`(1 - sqrt(2))*(1 + sqrt(2))`), with `polys=True` | Rust | `poly/api.rs` (`checked_constant_factor_list`), `symox-py/src/poly_checked_factors.rs` |
| `frac=True` lists with numeric generators | Rust | `symox-py/src/poly_checked_factors.rs` |
| `PolynomialError("a polynomial expected, got ...")` for a non-polynomial without `frac` | Rust | `symox-py/src/poly_checked_factors.rs` |
| `greedy=False` for `gcd`/`lcm`/`cofactors`/`cancel` (the algebraic and transcendental coefficients stay in the ground domain, over `EX`) | Rust | `poly/convert.rs` (`nongreedy_generators`), `symox-py/src/poly.rs` |
| `CoercionFailed` for a coefficient outside `QQ<a>` (`factor(x**2 - 2*t**2, x, extension=sqrt(2))`) | Rust | `poly/algebraic.rs` |
| `ZZ[t](1/t)`, `QQ_I[t](1/t)`, `from_sympy`: ValueError "expected an expression convertible to a polynomial in Polynomial ring ..."; `convert`: CoercionFailed | Rust | `poly/domains.rs`, `symox-py/src/domains.rs` |
| `real_roots`/`all_roots` with `extension=True` on a non-algebraic Poly: the Poly is rebuilt from its expression with `extension=True`, as the reference does | Python (routing only) | `python/symox/polys/polytools.py` (`_univariate_for_roots`) |

`CRootOf` extensions already matched SymPy on main. The algpoly gap came
from `test_polytools` statements that SymPy itself fails. The 35 corpus
cases are kept as a guard.

## Deliberate deviations (to arbitrate)

- `sqf_list` of a product now follows SymPy's `_symbolic_factor_list`: the
  components of equal multiplicity are multiplied together. So
  `sqf_list((x**2 - 1)*(x + 1))` is `(1, [(x**3 + x**2 - x - 1, 1)])`, which
  is not square-free, like SymPy (`sqf` already did this on main). One
  assertion of `tests/test_fraction_factor_lists.py` that required global
  square-freeness now compares live with SymPy.
- `minpoly r*rc` (`r`, `rc` two conjugate `CRootOf` of `x**3 + x + 3`): SymPy
  hits the 60 s alarm, and symox answers `y**3 - y**2 - 9`. This is checked
  by hand: the product of the three roots is `-3`, so `r*rc = -3/r2` with
  `r2` the third root. Substituting `r2 = -3/y` into `x**3 + x + 3` and
  clearing denominators gives `3*y**3 - 3*y**2 - 27`, which is irreducible.

## Remaining gaps

- algpoly item 5 (representatives): `monic` over `QQ<a>[t]`, and
  `gcd`/`lcm`/`sqf_list` over `QQ<a>(t)`. They give equal values with
  another numerator/denominator pair. Every Poly output still goes through
  `gaussian_ground::normalize_poly`. Fixing this needs Poly coefficients
  that keep their pair, and a port of the unnormalized subresultant PRS. Not
  done in this lot.
- `real_roots` over `QQ_I`/gaussian with `extension`: the sorting goes
  through the Python `which_real_roots`. The corpus matches SymPy, but the
  computation is not native.
- `EX[t](1/t)` and `GF(p)[t](1/t)` raise NotImplementedError (SymPy:
  ValueError). `RR[t](1/t)` keeps its pre-lot class.
- The `CoercionFailed` message of a gaussian ground differs from SymPy's
  wording (same class).
- `factor(x/t + 1, x, domain=ZZ[t])`: the exception differs.
- `nroots` has no `extension` parameter in SymPy 1.14, so there is nothing
  to align.
- Three native computations cannot be interrupted, so the upstream runner
  skips them:
  - `factor(g, gaussian=True)` with exponent 200000 (`test_factor_large#9`);
  - `cancel(q, _signsimp=False)` (`test_cancel#80`);
  - `minimal_polynomial` of `sum(sqrt(1..11))` (`test_minimal_polynomial_sq#9`).
- Upstream failures outside the algebraic-extension scope:
  - private APIs: ring `dup_*`/`dmp_*` (`test_factortools`),
    `_torational_factor_list`, `field_isomorphism` internals,
    `is_isomorphism_possible`, `_minpoly_*`;
  - `factor_large` resource limits;
  - `gcd_list`/`lcm_list` edge cases;
  - `cofactors(8, 6)` and `gcd(3.0, 9.0)` ("can't initialize without
    generators");
  - `gcdex(..., auto=False)` does not raise;
  - `factor` with a symbolic exponent `i`, `deep=True`, `expand=False`,
    lists;
  - `Piecewise` inside `cancel`;
  - `same_root` domain errors;
  - `cancel_modulus`.
