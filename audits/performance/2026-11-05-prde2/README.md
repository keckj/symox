# prde: parametric Risch DE above the base level in Rust (2026-11-06)

Base: branch `feat-prde2-20261105` on cab152e3 (main 15e75d3a + 88545232
"cancel(include) keeps unified QQ" + cab152e3). Parity notes:
`docs/algorithmic-parity.md`, section "prde (towers, parametric)".
Reference: SymPy 1.14 (`~/symox-ci/venv`). Machine: demeter, shared with
other agents, so timings are noisy.

## What is computed where

| Item | Computation | Python role |
| --- | --- | --- |
| `param_rischDE(fa, fd, G, DE)` on an exp/primitive level: `prde_normal_denom` (`weak_normalizer`, `gcdex_diophantine`), `prde_special_denom` (exp and primitive arms, `parametric_log_deriv` of the level below), `prde_linear_constraints`, `constant_system`, `poly_linear_constraints`, the block matrices, `nullspace`, the final `h` and `M` | Rust (`symox-py/src/prde_param.rs`), driving the native Poly/Expr objects and the Rust `DDM` nullspace | dispatch (`_prde_native._tower_matrix`), Dummy creation from the native trace, PolyMatrix construction from `(rows, cols, items, gens, domain)`, raising the stock exceptions |
| `param_poly_rischDE`: `prde_spde`, `bound_degree` (base, exp and primitive arms, including the `db == da` primitive branch that calls `limited_integrate` one level down), `prde_no_cancel_b_large`, `prde_no_cancel_b_small` at the base level, `prde_cancel_liouvillian`, the recursion | Rust (`prde_param.rs`) | same |
| `prde_cancel_liouvillian(b, Q, n, DE)`: the `i` loop, `DecrementLevel`, the lower-level `param_rischDE` (native base kernel or the tower kernel), `set_gens`, the block matrix | Rust (`prde_param.rs`) | hook at the top of the stock body (`prde.py`) |
| `limited_integrate(fa, fd, G, DE)`: `limited_integrate_reduce`, the `param_rischDE` call, the `constant_system`/nullspace selection of the solution | Rust (`prde_param.rs`) | dispatch |
| base-level calls made by the tower kernels (`param_rischDE` base case, `is_log_deriv_k_t_radical_in_field`, `parametric_log_deriv_heu`) | Rust (existing `prde.rs`, `prde_tower.rs`, `prde_log.rs` kernels) | none |
| `real_roots` in `bound_degree` (clear_denoms, primitive part, retract to ZZ) | Rust (`prde_param.rs`) | none |
| tan / nonlinear levels, `prde_no_cancel_b_small` above the base level, `prde_no_cancel_b_equal`, the declines listed in the parity notes | Python (stock body), recorded as gaps | |

All computation of the four functions above the base level is in Rust; no
Python algorithm was added or extended. The Rust code drives the native
`symox._core.Poly`/`Expr` methods through PyO3 (the same objects the stock
body manipulates), keeps the PolyMatrix as `(rows, cols, entries, gens,
domain)` in Rust and computes the nullspace with the Rust `DDM`.

## Scripts (`scripts/`)

- `random_param.py OUT [COUNT] [SEED] [TIMEOUT]`: seeded random calls
  (`random.Random(20261106)`) of the four functions on ten towers (exp,
  log, exp(2x), log(x+1), exp(x**2/2), log then exp, exp then log,
  log(log(x)), exp then a primitive, tan), with random coefficients in the
  lower variables (some rational). Each row is
  `[index, kind, spec, result, seconds, dummies, native]`, `native` being
  the number of native tower answers of the four functions during the
  call; a call over TIMEOUT seconds (default 20) is recorded as
  `['timeout']`. `$PRDE_ONLY` restricts the run to some indices.
- `compare_rp.py SYMPY NATIVE STOCK [SHOW]`: three-way classification of
  each call (same as SymPy, pre-existing difference shared with the stock
  symox body, native-only difference, timeouts), Dummy counts and times.
- `cases/gen.py BACKEND`, `cases/more2.py`: the handcrafted regression
  cases of `tests/test_native_prde_param.py`, printed with their srepr,
  domains and Dummy counts.
- The replays of the 142 calls recorded from the SymPy prde/Risch test
  copies and the 485 calls recorded from `risch_integrate` on the
  integration corpus use `../2026-11-05-prde/scripts/replay.py` and
  `compare.py` (inputs in `../2026-11-05-prde/raw/calls_*`).

`PRDE_STOCK=1` sets `_prde_native._ROOT = None`, which disables every native
prde admission (stock symox bodies, the comparison baseline).

## Recorded calls against SymPy 1.14 (`raw/r11_*`)

| set | calls | SymPy | native | stock symox | main 2e623e66 |
| --- | ---: | ---: | ---: | ---: | ---: |
| SymPy test copies | 142 | 3.095 s | 0.551 s | 2.735 s | 1.212 s |
| integration corpus | 485 | 37.113 s | 16.882 s | 34.302 s | 44.866 s |

Native against SymPy: 0 result differences, 0 Dummy-count differences on
both sets. The stock symox body on this base has the same results (1 extra
Dummy on the first call of a process, lazy initialisation). Main 2e623e66
differed from SymPy on 29 test calls and 254 corpus calls by the ZZ/QQ
domain of returned Polys; 88545232, on which this lot is based, removed
those differences (they are not part of this lot).

Native coverage (`STATS`, nested calls included): corpus 61
`param_rischDE`, 61 `param_poly_rischDE`, 61 `prde_cancel_liouvillian` and
61 `limited_integrate` tower answers; tests 13, 9, 10 and 5. With
`SYMOX_PRDE_DEBUG=1` the tower kernel reports no decline on either set. In the test
calls, records 119 to 122 run the stock body: their DE has `log(5)` in the
ground domain (`ZZ[x,_t0,log(5)]`), which the admission (plain symbol
generators only) rejects.

## Random calls against SymPy 1.14 (`raw/rp_*.json.gz`)

600 draws, 584 calls (16 draws had a zero denominator): 144
`param_rischDE`, 140 `param_poly_rischDE`, 154 `prde_cancel_liouvillian`,
146 `limited_integrate`. The native kernels answered in 461 calls (92, 110,
132 and 127 per kind); the others are tan levels, stock raises and the
declines of the parity notes.

| class (`compare_rp.py`) | calls |
| --- | ---: |
| same result as SymPy | 536 |
| differs from SymPy, identical in the stock symox body | 31 |
| differs from SymPy only natively | **0** |
| SymPy over 20 s, native and stock answer | 4 |
| SymPy and stock over 20 s, native answers | 6 |
| all three over 20 s | 7 |

No Dummy-count difference on calls whose results agree. Time on the 567
calls without timeout: SymPy 186.41 s, native 19.19 s, stock 146.76 s.

The 31 pre-existing differences (same in the stock body, not touched by
this lot):
- 12: generator order of the fraction-field domain of the returned
  PolyMatrix (`QQ(t0,x)` against SymPy's `QQ(x,t0)`), entries equal;
- 8: `TypeError` message names the class `symox._core.Poly` where SymPy
  says `Poly` (`'>' not supported between instances of ...`);
- 5: `TypeError('cannot convert object of type NoneType ...')` where SymPy
  raises `SympifyError: None` (the stock `derivation` adds `None` to a Poly;
  also the `risch_prim` case of the tests);
- 2: a differently scaled `h` basis (480, 514), same matrix;
- 4: the constraint matrix `M` differs (47, 241, 407, 573), `h` identical:
  the stock symox PolyMatrix nullspace/rref returns another row form.

TIMEOUT_RERUN

## Handcrafted cases (`raw/cases_*.txt.gz`)

`cases/gen.py` (the cases of `tests/test_native_prde_param.py` and more):
outputs identical to SymPy except `risch_prim` (`TypeError` against
`SympifyError`, pre-existing, kept as a strict xfail).

## SymPy test copies

test_prde.py, test_rde.py, test_risch.py copied from SymPy 1.14 with
`sympy` replaced by `symox`: 5 failed, 56 passed, the same list with the
stock bodies (test_prde_special_denom, test_gcdex_diophantine,
test_laurent_series, test_recognize_derivative, test_residue_reduce).
