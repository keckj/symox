# prdelog: base-case is_log_deriv_k_t_radical_in_field, parametric_log_deriv_heu, limited_integrate in Rust (2026-10-30)

Base: main 21068f42. Parity notes: `docs/algorithmic-parity.md`, section
"prdelog". Machine: demeter (112 cores), shared with other agents; load
average 13 to 67 during the measurements, so timings are noisy and are
given as minima over repetitions.

## What is computed where

| Item | Computation | Python role |
| --- | --- | --- |
| `is_log_deriv_k_t_radical_in_field`, base case | Rust (`symox-core/src/integrate/prde_log.rs`) | admission, argument passing |
| `parametric_log_deriv_heu` (and `parametric_log_deriv`), base case | Rust (`prde_log.rs`) | admission, argument passing |
| `limited_integrate`, base case | Rust (`prde_log.rs`, native `param_rischDE`) | admission, argument passing |
| `param_rischDE` with `h == []` and `W == []` (`C = eye(m)`) | Rust (`symox-core/src/integrate/prde.rs`) | none (was a decline) |
| splitfactor admission precheck | Rust (`symox-py/src/risch_splitfactor.rs`) | call order in `_risch_splitfactor_native.reduce` |
| non-base cases, non-ZZ/QQ coefficients, explicit `z`/`c1`, declines | Python (stock body), recorded as gaps | |

## Scripts (`scripts/`)

- `integrands.py`: the 114 integrands of the profile corpus.
- `record_calls.py OUT.json`: records the base-case calls of the three
  functions made by `risch_integrate` on the corpus (38 distinct calls) and
  adds 900 seeded random base-case calls (`random.Random(20261030)`).
- `replay_calls.py <sympy|symox|symox-stock> CALLS OUT`: replays every
  call, 3 repetitions, 20 s alarm, stores the srepr-level encoding of the
  result (or the exception class) and the minimum time. `symox-stock`
  disables the three native entries (Python bodies).
- `compare_calls.py CALLS SYMPY SYMOX STOCK`: call-by-call comparison and
  per-kind times.
- `splitfactor_overhead.py OUT`: per-call cost of `risch.splitfactor` and of
  the stages of the native admission on the splitfactor calls of the first
  60 integrands.
- `prof_targets.py <backend> nowrap`: `risch_integrate` on the 114
  integrands, 30 s alarm each, time per integrand, exceptions, and the native
  counters. `symox-stock` = the state before this lot (three entries
  disabled, precheck removed); `symox-noprecheck` = new entries, precheck
  removed. The wrapped mode (without `nowrap`) defeats the native protocol
  check and is not used for symox.
- `sum_min.py RUN...`: sum over integrands of the minimum over runs,
  excluding integrands that timed out.
- `make_tests.py`: generates the `CASES` table of
  `tests/test_native_prdelog.py` from the replay files.

## Differential replay against SymPy 1.14 (`raw/replay_*.json`, `raw/compare.txt`)

938 calls; "eq" counts identical encodings (srepr of every value, Poly
domains, None, exception class). Times are sums of per-call minima, in s.

| kind | part | n | eq SymPy | eq stock | SymPy | native | stock |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| heu | random | 300 | 300 | 300 | 5.673 | 1.333 | 5.561 |
| is_log_deriv | corpus | 10 | 10 | 10 | 0.037 | 0.021 | 0.028 |
| is_log_deriv | random | 300 | 300 | 300 | 7.948 | 0.889 | 7.300 |
| limited | corpus | 28 | 28 | 28 | 0.647 | 0.087 | 0.093 |
| limited | random | 300 | 300 | 300 | 3.643 | 0.766 | 1.145 |

0 differences. Native counters: is_log_deriv 310/310, heu 300/300,
limited_integrate 327/328; the declined call has `ZZ[log(2)]`
coefficients and runs the stock body (same result). Totals: SymPy 17.95 s,
stock 14.13 s, native 3.10 s.

## Splitfactor admission overhead (`raw/splitfactor_overhead_*.json`)

63 splitfactor calls from `risch_integrate` (first 60 integrands); none has
a QQ Poly (domains `ZZ`, `ZZ[x]`, `ZZ(x)`, `ZZ[_t0]`, `ZZ(_t0, x)`, ...), so
the kernel admits none of them, before and after. Mean per call, in µs:

| | `risch.splitfactor` | `_risch_splitfactor_native.reduce` |
| --- | ---: | ---: |
| before (authentication, then kernel conditions) | 3906 | 3437 |
| after (Rust precheck first) | 454 | 3.7 |

Before, `reduce` spent 905 µs in `_states_match` (416 function states) and
2497 µs in `_prove_fields`. After, these stages are not reached when the
precheck fails. QQ inputs still go through the full authentication.

## risch_integrate on the 114 integrands (`raw/prof/`)

Four backends run concurrently, 3 runs each, `nowrap`. Sum of the
per-integrand minima over the 3 runs, excluding the one integrand that hits
the 30 s alarm on every backend (`exp(x)*log(exp(x) + 1)`):

| backend | s |
| --- | ---: |
| SymPy 1.14 | 6.494 |
| symox before (stock) | 3.551 |
| symox, new entries, no precheck | 3.377 |
| symox after | 2.333 |

The exception lines (8, including the timeout) are identical on all twelve
runs. Most of the gain comes from the splitfactor precheck (about 1 s, a
roughly constant 40 ms on each exp/log integrand); the native base-case
entries contribute about 0.17 s, because most `limited_integrate` time in
the corpus is in non-base calls (55 of 89 calls), which stay in Python.
Native counters for the symox run: is_log_deriv 31, limited_integrate 36,
heu 0 (the corpus reaches `parametric_log_deriv_heu` only in non-base
cases).

## SymPy test files

Copies of SymPy 1.14 `test_prde.py` and `test_risch.py` with `sympy`
replaced by `symox`: test_prde 1 failed (`test_prde_special_denom`), 14
passed; test_risch 5 failed (`test_gcdex_diophantine`,
`test_laurent_series`, `test_recognize_derivative`, `test_residue_reduce`,
`test_risch_integrate_float`), 31 passed. The same tests fail on the base.

## Validation

Raw summaries in `raw/validation/`.

- Full suite, `scripts/run_suite24.py` (24 shards; identity: base
  21068f42 plus the working diff, `.so` sha256 4ee1188b...):
  59 617 tests collected, 59 177 passed, 44 skipped, 4 xfailed, 36
  failures in the first pass. The machine was heavily loaded (the corpora,
  cargo and other agents ran at the same time):
  - 28 failures in 21 files were subprocess timeouts or wall-clock
    assertions (`TimeoutExpired`, `'timeout' == 'completed'`, an elliptic_pi
    timing at 5.1 s against 5 s). Rerun file by file: all pass.
    `test_sympy_parity_audit.py` needs the suite environment
    (`SYMOX_CI_PY`); it passes 13/13 there and fails the same 3 tests on
    the base outside it.
  - `tests/extra` (shard 23) stopped on a `MemoryError` inside pytest's
    junit writer under the 6 GB address-space cap (7 failures, 358 tests not
    run). Rerun under `systemd-run MemoryMax=6G`: 367 passed, 16 skipped.
  - 4 failures of
    `test_native_lll_classifier_callbacks.py::test_rank_exception_constructor_observes_original_reduction`
    are order-dependent: the file passes 8/8 alone, and the same shard file
    list fails the same 4 tests on a clean 21068f42 worktree (without
    `test_native_prdelog.py`). Pre-existing, not caused by this lot.
- `cargo test -p symox --release`: 836 passed, 1 ignored; 1 + 5 doctests
  pass. Run under `systemd-run MemoryMax=16G`: under `prlimit
  --as=6000000000` rustc aborts with an LLVM out-of-memory error while
  building the lib test (lto, codegen-units=1).
- `cargo clippy -- -D clippy::correctness`: one error, pre-existing on the
  base (`approx_constant` in `crates/symox-core/src/cephes_gamma.rs:137`,
  not touched here). With `-A clippy::approx_constant` added, both crates
  pass.
- Integration corpus (`tools/integration_corpus.py`, MemoryMax=4G, fresh
  SymPy reference): `{"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1,
  "SAME": 1171}`, the baseline.
- ODE corpus, `scripts/ode_corpus.py` (symox only, scored against the
  stored SymPy 1.14 observations `corpus-2026-10-26/ode-570f`, 16 jobs,
  60 s): branch `{cls: 315, hint: 283, def: 267}`, base 21068f42
  `{cls: 315, hint: 285, def: 267}`, same script at the same time. The 5
  rows that differ between the two runs are all a 60 s DEADLINE on one
  side: `lie_group_13`/`lie_group_14` (hint) rerun on both builds give the
  same text in 12 to 20 s; `separable_reduced_02`/`_11` take 43 to 78 s on
  both builds. No output differs where both runs finished. This scorer does not
  reproduce the coordinator's 303/279/274 baseline numbers (they come from
  another scorer, not found), so the comparison here is against the base build.

## Remaining gaps

See the parity section: non-base cases of the three functions (Python),
coefficients outside ZZ/QQ, explicit `z`/`c1`, kernel declines, the
`PolyMatrix` truncation difference, and the unchanged QQ splitfactor
authentication.
