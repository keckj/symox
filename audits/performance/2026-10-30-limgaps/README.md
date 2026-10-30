# limgaps: limit and Order gaps, Gruntz declines (2026-10-30)

Branch `feat-limgaps-20261030k`, base main 21068f42. Reference SymPy 1.14.0
(`~/symox-ci/venv`). Machine demeter (112 cores, shared with other agents),
so the timings are indicative only. The parity notes are in
`docs/algorithmic-parity.md`, section "limgaps (wave 36)".

"before" is the base tree with a native module that has only the decline
trace (`_core._limit_trace`) added. "after" is this branch.

## Fixes and where the computation runs

| Fix | Computation |
| --- | --- |
| Gruntz: merge exponential factors whose exponents cancel (`exp(E)`) | Rust `limits/gruntz.rs` |
| `loggamma(a)`, `a -> -oo`: SymPy's heuristics value `zoo`/`oo`/0 | Rust `limits/gamma_heuristic.rs` |
| Bounded `sin`/`cos` integer powers in `bounded_affine` | Rust `limits/mod.rs` |
| Bounded factor times a part with a limit, next to `floor`/`ceiling` | Rust `limits/mod.rs` (`bounded_times_limit`) |
| `frac`/`floor`/`ceiling` series at an infinite argument (bounds + rest, `nan`) | Rust `series/accum_value.rs` |
| Order of sums: `factor` + `extract_leading_order` | Rust `series/order_dominant.rs` |
| Order: bounded trig factors, `PoleError` branch, containment by limit | Rust `series/order.rs` |
| Order at different points: `NotImplementedError` | Rust `series/order.rs`, `crates/symox-py/src/series.rs` |
| `O(e).subs(x, cos(x))` at `+-oo*I`, multivariable error | Rust `series/order_point.rs`, `subs.rs` |
| Decline trace (audit only) | Rust `limits/trace.rs` |
| Build `AccumBounds(lo, hi) + rest` from the native result | Python (host node construction only) |
| `nseries`/`series` call the native accum value after `PoleError` | Python (routing) |
| `oscillating_limit` tries native `bounded_affine` first | Python (routing) |
| `_bounded_affine` adds the native terms with the `AccumBounds` operators | Python (host node construction only) |

No Python algorithm was added or extended.

## SymPy test files

### Statement-level replay (`scripts/tl.py`, `scripts/run_tl.sh`)

Each test assertion that mentions limit/series/Order is evaluated, and its
left-hand side is compared with SymPy's output text. The table counts
statements equal to SymPy (`scripts/compare.py before after`).

| file | statements | before | after |
| --- | ---: | ---: | ---: |
| test_limits | 554 | 519 | 520 |
| test_gruntz | 107 | 95 | 95 |
| test_series | 83 | 44 | 44 |
| test_nseries | 186 | 164 | 165 |
| test_order | 210 | 157 | 172 |
| test_formal | 85 | 59 | 59 |
| test_demidovich | 54 | 53 | 53 |
| test_aseries | 9 | 5 | 5 |
| test_lseries | 5 | 5 | 5 |
| total | 1293 | 1101 | 1118 |

17 statements fixed and 0 broken. The full list is in the output of
`compare.py`.

### pytest on symox copies (`scripts/mkcopies.py`, `scripts/run_pt.sh`)

The copies replace `sympy` imports with `symox`; names that are missing
raise `NotImplementedError`.

| file | before passed | after passed | failed after |
| --- | ---: | ---: | ---: |
| test_limits | 140 | 141 | 32 |
| test_gruntz | 19 | 19 | 19 |
| test_series | 26 | 26 | 14 |
| test_nseries | 51 | 52 | 17 |
| test_order | 30 | 37 | 20 |

## Gruntz declines

A decline is the native `_limit` returning undecided. The counts are
distinct (reason, limit) pairs (`scripts/declines.py`).

| source | before | after |
| --- | ---: | ---: |
| statement replay (all 9 files) | 92 | 90 |
| pytest copies (5 files) | 79 | 77 |

On the statement replay, 68 declines are "series in mrv variable", both
before and after (before: `sin` 24, `cos` 23, `tan` 10, `floor` 10, `cot` 1).
SymPy's `gruntz` fails there too. The value is an `AccumBounds`, which is
a host node, so these limits go on to the oscillating front end. Its route
on the statement replay (`raw/tl/*_*.stats.json`):

| route | before | after |
| --- | ---: | ---: |
| native `bounded_affine` (Rust) | 2 | 19 |
| Python walker | 57 | 38 |
| held | 31 | 29 |

## Gap cases (`scripts/cases.py`, `raw/cases_*.txt`)

37 cases: 11 equal to SymPy before, 33 after. The 4 that still differ:

- `factorial(x)/x` and `gamma(x)/x` at `-oo`: held (intentional, en
  attente de validation utilisateur);
- `limit(Abs(x + y)**2, x, 0)`: the `Abs`/`sign(c)` rule (intentional,
  kept);
- `floor(1/sin(p)).series(p)`: a gap.

Single-call times in ms (after process start; SymPy's times include its
caches warming up):

| case | before | after | SymPy |
| --- | ---: | ---: | ---: |
| `limit(((x**(x+1)+(x+1)**x)/x**(x+1))**x, x, oo)` | 22.9 (held) | 8.5 | 900.3 |
| `frac(x).series(x, oo)` | 1.8 (PoleError) | 1.9 | 47.8 |
| `limit(loggamma(x), x, -oo)` | 1.6 (held) | 0.5 | 31.2 |
| `O(exp(1/x) + 1)` | 0.5 | 1.7 | 62.3 |
| `O(log(x) + log(log(x)))` | 3.1 | 0.8 | 69.4 |
| `O(x**2, (x, oo)).subs(x, cos(x))` | 6.7 | 0.7 | 172.2 |
| `O(x + factorial(x), (x, oo))` | 0.5 | 0.7 | 531.0 |
| `limit(cos(x)*floor(x)/x, x, oo)` | 4.9 | 1.9 | 4.5 |
| `limit(sin(x)**2, x, oo)` | 4.6 | 1.2 | 2.3 |
| `limit(atan(x) + sin(x), x, oo)` | 6.9 | 2.5 | 74.0 |
| `floor(1/p + p).series(p)` | 0.3 (PoleError) | 1.0 | 12.2 |

A "before" value marked held or with an error was not equal to SymPy.

The other comparison files are `raw/fr_*` (frac/floor series), `raw/ord*_*`
(Order forms) and `raw/lg_*` (loggamma). Each symox output equals SymPy's,
except `floor(1/sin(p)).series(p)` in `fr`. The `ord`/`ord2` symox files
add a time column.

## Validation

- Full `tests/` suite, 24 shards (`scripts/run_suite.py`): 59 586
  collected and run, 59 399 passed, 83 skipped, 4 xfailed, 100 failed.
  Coverage is complete: every collected test ran.
  - 96 failures are in `tests/extra` (jax/torch codegen, CUDA interop).
    They are `CUDA_ERROR_OUT_OF_MEMORY` and `CUDA error: unknown error`
    while other agents held the GPU. A rerun of that shard alone gave 366
    passed, 16 skipped, 1 failed; that test
    (`test_jax_codegen.py::test_custom_functions_and_standalone_source`)
    passes alone.
  - 4 failures are
    `test_native_lll_classifier_callbacks.py::test_rank_exception_constructor_observes_original_reduction[*]`.
    They are order dependent: the file passes alone (8 passed), and the
    same shard file list fails the same way with the base module.
- `cargo test -p symox --release`: 836 passed, 1 ignored; 1 + 5 (doctests)
  passed.
- `cargo clippy -- -D clippy::correctness`: 1 error. It is
  `approx_constant` (`EULER_GAMMA`) in `crates/symox-core/src/cephes_gamma.rs:137`,
  a file this lot does not touch (last change e0a95719 on main), reported
  by the current clippy. With `-A clippy::approx_constant` added: clean
  (exit 0, both crates).
- Integration corpus (`tools/integration_corpus.py`):
  `{"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1, "SAME": 1171}`, the
  baseline.
- ODE corpus (`corpus-r11` `replay-ode.py`, 315 cases):
  - classify 315 equal to SymPy;
  - hint 282 equal;
  - default 263 equal.
  - Against the smallgaps run, which is at the baseline, the symox outputs
    show 0 per-case differences (`scripts/ode_compare.py`). The only change
    is `fact_06`, a resource-limited observation (inconclusive) that moved
    between hint and default.
- `tests/test_native_limgaps.py`: 85 passed.

## Reproduce

From the worktree root, with `S` a scratch directory:

```
WORK=$S/work bash audits/performance/2026-10-30-limgaps/scripts/run_tl.sh after symox $PWD
WORK=$S/work bash audits/performance/2026-10-30-limgaps/scripts/run_tl.sh sympy sympy $PWD
WORK=$S/work bash audits/performance/2026-10-30-limgaps/scripts/run_pt.sh after $PWD
cd audits/performance/2026-10-30-limgaps
python3 scripts/compare.py before after
python3 scripts/declines.py after raw
cd ../../..
PYTHONPATH=$PWD/python ~/symox-ci/venv/bin/python audits/performance/2026-10-30-limgaps/scripts/cases.py symox
~/symox-ci/venv/bin/python audits/performance/2026-10-30-limgaps/scripts/cases.py sympy
```

`scripts/run_suite.py` (with `validation_plugin.py`) is the 24-shard runner
of the smallgaps lot. `scripts/ode_compare.py` compares two `replay-ode.py`
output directories case by case.
