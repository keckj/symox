# prde3: the remaining parametric Risch DE paths in Rust (2026-11-06)

Base: branch `feat-prde3-20261106` on main e4b1d011. Parity notes:
`docs/algorithmic-parity.md`, section "prde3 (tan and non-linear levels,
no-cancel cases above the base)". Reference: SymPy 1.14 (`~/symox-ci/venv`).
Machine: demeter, shared with other agents, so timings are noisy. "before"
is main e4b1d011 (same `.so` as the base of the branch), "after" is the
final `.so` of this lot.

## What is computed where

| Item | Computation | Python role |
| --- | --- | --- |
| `param_poly_rischDE` on tan and non-linear levels: `prde_no_cancel_b_large`, `prde_no_cancel_b_small` above the base level (the `DE.d.degree() >= 2` arm), `prde_no_cancel_b_equal` (the stock `NotImplementedError`), the "non-linear and hypertangent" `NotImplementedError`, `bound_degree` on these levels | Rust (`symox-py/src/prde_param.rs`) | dispatch, Dummy creation, PolyMatrix construction, raising the stock exceptions |
| `prde_special_denom` tan arm: `real_imag`, `recognize_log_derivative` (with the `real_roots` integrality test), the two `parametric_log_deriv` calls one level down, `n = min(n, s/2)`; `order_at` and the `ValueError` of unknown cases | Rust (`prde_param.rs`) | none |
| `param_rischDE` with no polynomial solution (`v == 0`): SymPy's `wl[:m] + wl[-v:]` (the whole `wl`) and `Matrix(len(W), m + v, items)` keeping the first `len(W)*(m + v)` items, the nullspace of the empty matrix (identity) | Rust (`prde_param.rs`; the kernel declined there before) | none |
| `parametric_log_deriv_heu` on ground arguments (tan arm with ground real/imaginary parts): the stock `None` after its `Dummy('c1')` | Rust (`prde_param.rs`, `heu`) | none |
| public `prde_no_cancel_b_large` / `prde_no_cancel_b_small` above the base level | Rust (`prde_param.rs`, kinds `*_tower` in `prde.rs`) | `_prde_native._no_cancel` routes above-base calls to the tower kernel and writes back `Q` |
| kkinv `PolynomialError` of SymPy `residue_reduce` (case 44 of `tests/test_native_prde_tower.py`) | Rust (`symox-core/src/integrate/prde_tower.rs` `Stop::Polynomial`, `symox-py/src/prde_tower.rs` `kkinv_error`, status 3) | `_prde_native._tower` raises `symox.polys.polyerrors.PolynomialError(payload)` |
| generator order of `construct_domain` / PolyMatrix domains (`QQ(x,t0)`, SymPy `_sort_gens`) | Rust (`symox::poly::convert::gen_key` in `prde_const.rs`; `_core._composite_generators`) | `polys/constructor.py._infer_domain` asks the native generators for the order (routing, no algorithm) |
| `Poly + None` raises `SympifyError` (SymPy `_polifyit`) | Rust (`symox-py/src/poly.rs`, binop) | none |

No Python algorithm was added or extended: the Python changes are the
dispatch of the public `prde_no_cancel_*` calls above the base level, the
status-3 exception and the generator order taken from the native
`_construct_composite`.

## Scripts (`scripts/`)

- `random_param.py OUT [COUNT] [SEED] [TIMEOUT]`: the prde2 harness
  (`../2026-11-05-prde2/scripts/random_param.py`) with three changes:
  `Timeout` is a `BaseException` (an `Exception` raised by SIGALRM inside a
  native kernel callback is turned into a decline, so the prde2 harness let
  some calls run past the limit, e.g. rp 196 at 125.8 s "before");
  `PRDE_TOWERS=nonlin` draws on 12 tan / non-linear towers (tan(x),
  tan(2x), tan(x**2/2), `D t0 = t0**2`, `1 + x*t0**2`, `t0**3 + 1`,
  `t0 + 1`, log then tan, exp then tan, tan then a primitive, tan then
  exp(tan), log then a non-linear level); `PRDE_PROGRESS=1` prints each
  index. Rows are `[index, kind, spec, result, seconds, dummies, native]`.
- `compare_rp.py SYMPY NATIVE STOCK [SHOW]`: three-way classification
  (gzip aware).
- `runnl.sh TAG [stock]`, `runall.sh TAG [stock]`: the nonlin corpus, the
  default corpus and the two replays of `../2026-11-05-prde`.
- `gen3.py MODULE`: 19 handcrafted tan / non-linear calls (output in
  `raw/gen3_*.json`).
- `mform.py MODULE INDICES OUT [SOURCE]`: the public forms of the returned
  `M` (`to_Matrix` srepr, domain, `M[i, j]`).

`PRDE_STOCK=1` disables every native prde admission (stock symox bodies).

## Handcrafted tan / non-linear calls (`raw/gen3_*.json`)

19 calls (`param_poly_rischDE` 10, `param_rischDE` 5, `limited_integrate`
1, public `prde_no_cancel_b_small` 2 and `prde_no_cancel_b_large` 1): 19/19
identical to SymPy (results, exception types and messages, Dummy counts).
Answers: 13; `NotImplementedError("non-linear and hypertangent cases have
not yet been implemented")`: 4; `NotImplementedError("prde_no_cancel_b_equal()
is not yet implemented.")`: 1; `ValueError("case must be one of ...")`: 1.
Every call is a native answer (`STATS`: `param_poly_rischDE_tower` 20,
`param_rischDE_tower` 10, `limited_integrate_tower` 2,
`prde_no_cancel_b_small_tower` 4, `prde_no_cancel_b_large_tower` 2, two
runs per call). They are the `CASES` of `tests/test_native_prde3.py`,
with two random calls (nonlin 190 and 454) on the `v == 0` path.

## Random calls on tan / non-linear towers (`raw/nl_*.json.gz`)

600 draws (`PRDE_TOWERS=nonlin`, seed 20261106, 20 s), 582 calls: 126
`param_rischDE`, 147 `param_poly_rischDE`, 165 `prde_cancel_liouvillian`,
144 `limited_integrate`.

| | before | after |
| --- | ---: | ---: |
| calls with a native answer | 64 | **376** |
| (param_rischDE / param_poly / cancel_liouvillian / limited) | 7 / 20 / 11 / 26 | 84 / 121 / 35 / 136 |
| same result as SymPy | 555 | **573** |
| differs from SymPy, identical in the stock body | 24 | 6 |
| differs from SymPy only natively | 0 | **0** |
| over 20 s in SymPy | 3 | 3 (2 answered natively) |
| time without timeouts: SymPy / native / stock | 138.89 / 112.76 / 126.73 s | 138.89 / **12.27** / 111.17 s |

No Dummy-count difference on calls whose results agree.

Over 20 s (`raw/nl_long_*.json`, reruns at 600 s): 75
(`param_poly_rischDE`) native 22.6 s against SymPy 495.8 s, 151
(`limited_integrate`) native 3.0 s against SymPy 113.9 s, both identical to
SymPy; 350 over the limit in both (native rerun at 300 s).

The 6 remaining differences (all in the stock body too):
- 124: `CoercionFailed` message (stock path);
- 490, 554: `PolynomialError` message prints `t0**(-2)` where SymPy prints
  `1/t0` (stock path);
- 545: `TypeError` names `symox._core.Poly` (stock path);
- 212, 598: the raw entries of the returned PolyMatrix hold uncancelled
  zeros (`_items`); every public form (`to_Matrix`, domain, `M[i, j]`,
  equality, repr) is identical to SymPy (`mform.py`).

## Random calls on the default towers (`raw/rp_*.json.gz`)

The prde2 draws (600 draws, 584 calls) against `../2026-11-05-prde2/raw/rp_sympy.json.gz`.

| | before | after |
| --- | ---: | ---: |
| calls with a native answer | 461 | **506** |
| (param_rischDE / param_poly / cancel_liouvillian / limited) | 92 / 110 / 132 / 127 | 103 / 129 / 132 / 142 |
| same result as SymPy | 536 | **554** |
| differs from SymPy, identical in the stock body | 31 | 13 |
| differs from SymPy only natively | 0 | **0** |
| over 20 s in SymPy | 17 | 17 |
| time without timeouts: SymPy / native / stock | 186.41 / 23.48 / 159.22 s | 186.41 / **17.72** / 139.03 s |

The 18 differences removed: 12 PolyMatrix domain generator orders
(`QQ(t0,x)` → `QQ(x,t0)`), 5 `SympifyError` (stock `derivation` adding
`None` to a Poly), 1 constraint-matrix form (241). The 13 left:
- 8 (48, 60, 74, 126, 203, 346, 380, 568): `TypeError` message names
  `symox._core.Poly` where SymPy says `Poly`;
- 3 (47, 407, 573): uncancelled zeros in the raw entries of `M`; public
  forms identical to SymPy (`raw/mform_rp_*.json`);
- 2 (480, 514): a differently scaled `h` basis, same `M` (prde2 gap).

Timeouts: 196 ran 125.8 s natively "before" because the prde2 harness let
SIGALRM be swallowed in a kernel callback; it is now reported as a timeout
(see the rerun below). The other calls keep their class.

Rerun of 196 (`raw/rp_long_symox_196.json`, 600 s, build before the
`v == 0` change): native 117.0 s,
result identical to the "before" run (`None`); SymPy is over 20 s there
(not rerun).

## Recorded calls against SymPy 1.14 (`raw/tests_*`, `raw/corpus_*`)

Replays of `../2026-11-05-prde/raw/calls_*` compared with `r3_*_sympy`:

| set | calls | SymPy | native before | native after | stock before | stock after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| SymPy test copies | 142 | 3.095 s | 0.905 s | 0.530 s | 3.638 s | 2.826 s |
| integration corpus | 485 | 37.113 s | 27.437 s | 18.631 s | 60.765 s | 34.220 s |

Native before and after: 0 result differences, 0 Dummy-count differences
on both sets (stock body: 1 extra Dummy on the first call of a process, lazy
initialisation). Native coverage is unchanged (corpus: 61 tower answers
per function; tests: 13 / 9 / 10 / 5): these calls hold no tan or
non-linear level, and records 119 to 122 of the tests still run the stock
body because `log(5)` is a ground generator (`ZZ[x,_t0,log(5)]`). The time
differences are machine noise.

## SymPy test copies

test_prde.py, test_rde.py, test_risch.py copied from SymPy 1.14 with
`sympy` replaced by `symox`: 5 failed, 56 passed, natively and with the
stock bodies, unchanged (test_prde_special_denom, test_gcdex_diophantine,
test_laurent_series, test_recognize_derivative, test_residue_reduce).

## Native declines on the tan / non-linear corpus (`raw/nl_after_debug.log`)

`SYMOX_PRDE_DEBUG=1` on the final build lists each kernel decline (nested
declines included). Of the 582 calls, 206 run the stock body; SymPy
answers only 1 of them (417, below), the others raise or time out in
SymPy too, and the stock body reproduces them. Decline events (a call
can decline at several nested levels):
- 130 `prde_cancel_liouvillian` called directly on a tan or non-linear
  level with `n >= 0`: SymPy raises `UnboundLocalError` (`ba` is only bound
  in the primitive and exp arms); its message depends on the Python
  version, so the kernel declines (`prde_param.rs` `cancel_liouvillian`);
- 43 `real_roots` whose preprocessed polynomial is not over ZZ (stock
  `NotImplementedError: sorted roots not supported over ...`);
- 61 stock-body exceptions met while transliterating
  (`PolynomialError: 1/t0 contains an element of the set of generators`,
  `TypeError: cannot determine truth value of Relational`, `CoercionFailed`,
  `ExactQuotientFailed`): the kernel declines on an unmodelled Python
  exception and the stock body raises the same one;
- 7 `parametric_log_deriv_heu` base-kernel declines in the tan arm of
  `prde_special_denom`: SymPy's `real_imag` leaves the current monomial
  in the coefficient domain (`Poly(-t0**3 + 2*t0, x, domain='ZZ[t0]')`)
  and the base kernel only admits ZZ/QQ ground coefficients; one of them
  (417) is a SymPy answer computed by the stock body;
- 4 `derivation` terms that are not polynomials, 2 step caps.

## Remaining gaps

- the tan arm of `prde_special_denom` when `real_imag` keeps `t` in the
  coefficient domain and the base `parametric_log_deriv_heu` kernel
  declines (1 SymPy answer in 582 calls, stock body);
- `prde_cancel_liouvillian` on tan/non-linear levels (`UnboundLocalError`
  of SymPy, stock body);
- the `v == 0`, `W` non-empty truncation of `Matrix(len(W), m, items)` is
  ported but met by no recorded or random call (scan of the SymPy runs,
  both corpora); the stock symox `PolyMatrix` raises `ValueError` on a
  list longer than `rows*cols` where SymPy truncates;
- non-symbol ground generators (`log(5)`, `exact_ground` declines): stock
  body; 4 recorded test calls, none in the corpora;
- `TypeError` messages naming `symox._core.Poly` (the PyO3 class name;
  changing `tp_name` touches pickling and the module layout);
- raw PolyMatrix entries with uncancelled zeros (not publicly visible);
- differently scaled `h` basis (480, 514);
- `CoercionFailed` / `PolynomialError` message spellings on stock paths;
- the stock symox body still computes the kkinv gcd of case 44 (only the
  native kernel raises SymPy's `PolynomialError`).

