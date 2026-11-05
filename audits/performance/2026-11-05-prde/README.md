# prde: is_log_deriv_k_t_radical_in_field and parametric_log_deriv_heu on exp/primitive towers in Rust (2026-11-05)

Base: main 2e623e66. Parity notes: `docs/algorithmic-parity.md`, section
"prde (towers)". Reference: SymPy 1.14 (`~/symox-ci/venv`). Machine:
demeter, shared with other agents, so timings are noisy.

## What is computed where

| Item | Computation | Python role |
| --- | --- | --- |
| `is_log_deriv_k_t_radical_in_field` on exp/primitive levels (any `case` string): `cancel`, `residue_reduce` over K(t) (subresultant PRS, `splitfactor_sqf`, the `H` loop with `gcd`/`invert`), `real_roots` checks, `p in k`, the exp branch (`parametric_log_deriv` one level down), the primitive branch (recursion), the base branch, the final `cancel` | Rust (`symox-core/src/integrate/prde_tower.rs`) | dispatch, Dummy creation from the native trace, raising the stock exceptions |
| `parametric_log_deriv_heu` (and `parametric_log_deriv`, which only wraps it) on exp/primitive levels, `splitfactor` included | Rust (`prde_tower.rs`) | same |
| Dummy trace of the base-case kernel (`z`, one per `real_roots`, `i`, `c1`) | Rust (`prde_log.rs`) | creation from the trace |
| `Poly.free_symbols` with the symbols of a typed composite domain | Rust (`symox-py/src/poly.rs`) | none |
| splitfactor admission loops (scopes, function states, lookups) | Rust (`symox-py/src/risch_sf_admission.rs`) | call order in `_risch_splitfactor_native.reduce` |
| `limited_integrate`, `param_rischDE`, `prde_cancel_liouvillian` above the base level | Python (stock body), recorded as a gap | |
| `z`/`c1` given, explicit exp case on a primitive level, kkinv `PolynomialError`, stock `real_roots` raises | Python (stock body), recorded | |

## Scripts (`scripts/`)

- `codec.py`: encoding of values (srepr of every expression, Poly domains,
  exceptions as `['raise', type, message]`) and decoding of recorded calls.
- `integrands.py`, `record_calls.py`: record every call of the prde entry
  points made by the SymPy prde/Risch test copies (`raw/calls_tests.json`,
  142 calls, 5 of them with a `ZZ_I` domain the codec does not decode)
  and by `risch_integrate` on the integration corpus
  (`raw/calls_corpus.json.gz`, 485 calls).
- `replay.py IN OUT [KIND...]`: replays them in `$PRDE_BACKEND_MODULE`
  (`sympy` or `symox`; `PRDE_STOCK=1` disables the native prde admission),
  storing the encoded result, the time and the number of Dummies created.
- `random_towers.py OUT N`: seeded random calls on exp/primitive/tan towers
  (`random.Random(20261105)`): `fa/fd` built from logarithmic derivatives
  of random radicals plus noise, random `case` strings, random heu inputs;
  533 distinct calls out of 600 draws. Each row also records whether the
  native tower kernel answered.
- `compare.py REF OTHER [SHOW]`: call-by-call comparison (results, Dummy
  counts, times).

## Random towers against SymPy 1.14 (`raw/rt_*.json`)

533 calls (429 `is_log_deriv_k_t_radical_in_field`, 104
`parametric_log_deriv_heu`). The native tower kernel answers 529 of them
(419 + 110 counted, nested calls included); it declines 4, which run the
stock body:

| rt | reason | SymPy |
| ---: | --- | --- |
| 17, 438 | explicit `case='exp'` on a primitive level | None |
| 71 | stock `real_roots` raise over QQ(...) | `NotImplementedError` |
| 146 | kkinv gcd of `residue_reduce` | `PolynomialError` |

Differences against SymPy: 1 result (rt 146: SymPy raises
`PolynomialError` in the kkinv gcd, the symox stock body computes the gcd
and goes on to the tan `NotImplementedError`; the same difference exists
before this lot), 2 Dummy counts (rt 136: lazy first-call Dummy of the
symox process, see below; rt 427: SymPy's `solve` creates 4 temporary
`subs` Dummies in `is_constant`, the symox `solve` does not; same in the
stock body). Stock symox has the same 3 differences (rt 7 instead of 136
for the lazy Dummy).

Times (sum over calls): SymPy 43.42 s, stock 4.43 s, native 0.45 s.

## Recorded calls: SymPy tests and integration corpus (`raw/r3_*`, `raw/r4_*`)

`r3_*_sympy`: SymPy 1.14; `r4_*_symox`: this branch; `r4_*_main`: main
2e623e66 (python tree of main with the wave-38 build); `r4_*_stock`: native
prde admission disabled (`PRDE_STOCK=1`). The corpus files are gzipped;
`compare.py` and `replay.py` read `.json.gz` directly.

| part | calls | result diffs vs SymPy | same vs main 2e623e66 | same vs stock | SymPy | native | main |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: |
| tests | 142 | 29 | all results | all results | 3.10 s | 1.12 s | 1.21 s |
| corpus | 485 | 254 | all results | all results | 37.11 s | 42.67 s | 44.87 s |

The result differences are all pre-existing (identical in main and in the
stock bodies) and lie in `param_rischDE`, `param_poly_rischDE` and
`prde_cancel_liouvillian`: 272 are a ZZ/QQ ground domain of returned Poly
objects, 11 (tests) a differently scaled basis of the solution space. The
Dummy-count differences against SymPy (31 tests, 112 corpus) are also
pre-existing (the base-case kernels of earlier waves create fewer internal
Dummies in `param_rischDE`/`limited_integrate`); against main, this lot
creates more Dummies (closer to SymPy) in 31 test and 20 corpus calls.

In the tests, every tower call of the two functions is answered natively
(7 `is_log_deriv`, 3 heu); the 6 non-native calls are base-level calls with
`ZZ[t]` coefficients (a prdelog admission gap). The corpus only reaches the
two functions at the base level. The corpus time is dominated by the
non-base `limited_integrate`/`param_rischDE`, still in Python.

## Splitfactor admission (`raw/splitfactor_*.json`)

`risch.splitfactor` per call (minimum over 5 x 100 calls, `exp(x)*log(x)`
extension):

| input | SymPy | main | branch |
| --- | ---: | ---: | ---: |
| `QQ` `2t+1` | 1023 µs | 3874 µs | 2837 µs |
| `QQ` `t^2 - 3t/2` | 2254 µs | 7798 µs | 5852 µs |
| `ZZ(x)` `t+x` | 1857 µs | 339 µs | 339 µs |

The native loops remove the Python iteration over the 50 namespaces, 400
function states and 3 000 lookups. What remains on QQ inputs is the shared
`_prove_fields` protocol, which calls `_metadata` about 395 times per
`splitfactor` call; each call rescans `_core.__dict__` (612 entries) with
`_str_keyed`, plus the `_polymatrix_native._class_unchanged` checks.
Scanning once per call would change the authentication contract (Python
code runs between two reads); it is recorded as a gap for user validation.

## Validation

- `cargo test --release`: 858 passed (1 ignored) + 1 + 0 + 5 (baseline
  852 + 1 + 5; the 6 extra tests are new `prde_tower` unit tests).
- `cargo clippy -- -D clippy::correctness`: clean.
- `tests/test_native_prde_tower.py`: 120 passed, 1 xfailed (strict, case 44:
  the kkinv `PolynomialError` gap).
- SymPy prde/Risch test copies: 6 failed / 55 passed, unchanged from main.
- Full suite (`audits/performance/2026-10-30-smallgaps/run_suite.py`, run
  with `XLA_PYTHON_CLIENT_PREALLOCATE=false`, see the jax note in the
  build memo): SUITE CLEAN, 60 155 passed, 6 xfailed, 45 skipped,
  60 206 collected, 683 files, 0 failures (baseline 60 020 / 60 085).
