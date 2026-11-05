# smallgaps audit, 2026-10-30

Corrections: Float Bareiss determinant, checksol issue 17799, second
"base as E" branch of `Pow.__new__`, units in checksol. lambdarr
gamma/poly2 was measured, with no code change. Reference: SymPy 1.14.0 in
`~/symox-ci/venv`. Before = main 21068f42 (scratch copy of
`python/symox` with the main `.so`). After = branch
`feat-smallgaps-20261030k`.

Where each computation runs:

| Correction | Computation |
| --- | --- |
| Float Bareiss `det` (default and `method='bareiss'`) | Rust (`matrix/linalg.rs`, `matrix/ratfunc.rs`) |
| Issue 17799, evalf ValueError during argument conversion | Rust (`evalf_mp.rs`; `symox-py/src/special.rs` only marks the error class) |
| Issue 17799 in checksol | Rust (`solve/checksol.rs`) |
| Second E branch of `Pow.__new__` | Rust (`ops.rs`) |
| Units in checksol (`as_independent(Unit)`) | Rust. Python collects the Quantity atoms and passes them as an argument (`_solver_checks.py`, `symox-py/src/tsolve.rs`) |
| lambdarr gamma/poly2 | no change: no Python frame on the call path |
| `cephes_gamma.rs` `#[allow(clippy::approx_constant)]` | lint attribute only. The error predates this change; the literal is correctly rounded Euler gamma |

## Scripts

- `det_cases.py <module> <out.json>`: 77 Float/symbolic matrices (sizes 1
  to 7, precisions 15 and 50, including SymPy's `detdata.py` matrices). It
  records `srepr(M.det(method=m))` for None, bareiss, berkowitz, lu, bird,
  laplace and domain-ge.
- `cmp_det.py ref.json[.xz] other.json[.xz]`: counts identical and different results
  per method. The raw outputs are stored xz-compressed
  (`raw/det_*.json.xz`).
- `bench.py <module> det_cases.py <out.json>`: median timings of the
  smallgaps cases.
- `markov6.py <module>`: the Float Markov chain of the markovrust gap.
- `ode_compare.py <base dir> <after dir>`: per-case comparison of two
  `replay-ode.py` outputs (corpus-r11 script).
- `run_suite.py <worktree> <evidence dir>` and `validation_plugin.py`: the
  full `tests/` suite in 24 file shards with collection/execution coverage
  proof. Change from the native-sdm-irref version: the collection and the
  `tests/extra/` shard run under `systemd-run --user --scope -p
  MemoryMax=6G` instead of `prlimit --as=6000000000`. Those files import
  torch, jax and cuda, which reserve more address space than `--as` allows,
  so pytest raised MemoryError while collecting. The other shards keep
  `prlimit`.
- checksol replays use the scripts of `../2026-10-30-checksol/`
  (`run_cases.py`, `cmp.py`), with the logged calls of SymPy's
  `test_solvers.py` and `test_solveset.py` and the 5304 generated cases.
- lambdarr timings use `../2026-10-29-lambdarr/bench.py` (cases
  `gamma_f64` and `_float`).

## Results

Float determinant, exact `srepr` against SymPy (`raw/det_cmp_*.txt`):

| method | before identical | after identical |
| --- | --- | --- |
| default | 43/77 (1 ran out of memory) | 77/77 |
| bareiss | 43/77 (1 ran out of memory) | 77/77 |
| berkowitz | 76/77 | 76/77 |
| lu | 75/77 | 75/77 |
| bird | 31/77 | 31/77 |
| laplace | 35/77 | 35/77 |
| domain-ge | 33/77 | 38/77 |

Markov chain, seed 7 (`markov6.py`): `(eye(6) - P).det()` was exactly 0
before. It is now `Float('-2.073958280700135e-17', precision=53)`, and
`fundamental_matrix()[0, 0]` is `2.04792404838000e+15`, as in SymPy.

checksol:

- 1897 distinct logged calls (`raw/logged_{sympy,before,after}.json`).
  Before: 1888 agree, 1 decision differs (issue 17799). After: 1889 agree,
  0 decisions differ. In both runs, 7 differences remain in the mapping
  values written back (print order of radicals, public `simplify`), and 1
  call cannot be rebuilt (setup).
- 5304 generated cases: 5304 agree, before and after.

Probes against SymPy:

- `checksol(-erf(x**(1/3))**pi + I, x, erfinv(I**(1/pi))**3)` is False and
  `solve(-erf(x**(1/3))**pi + I, x)` is `[]`.
- `gamma(erfinv(2)).n()` is unevaluated, and `sin(erfinv(2)).n()` raises
  ValueError, as in SymPy.
- `(2*I)**(x/log(2*I))` is `exp(x)`, and `(-2)**(x/(log(2) + I*pi))` stays a
  power.
- `checksol(cm*(x - 2), x, Float('2.0000000001'))` is True and
  `checksol(2*cm*x, x, 1)` is False.
- Remaining differences, all documented as gaps in
  `docs/algorithmic-parity.md`:
  - `(I**(1/pi))**pi` ;
  - `8*kilo*newton + x + y` (Prefix representation) ;
  - `solve(1/x - 1/(2*cm))` (Quantity is a free symbol) ;
  - issue 14860.

SymPy's own tests: `test_determinant.py` and `test_solvers.py` copied, with
private imports replaced by public calls and `test_issue_12114`
deselected because it aborts on both builds (`raw/sympy_copied_tests_*.txt`):

| build | failed | passed |
| --- | --- | --- |
| before | 55 | 200 |
| after | 54 | 201 |

The only change is `test_issue_17799`, which now passes. The other
failures are the same in both builds. They include `test_units` and
`test_issue_14860` (units representation), plus `test_issue_14517` and
`test_det`, which fail for other reasons (`det(iszerofunc="test")` raises
TypeError instead of ValueError).

Corpora:

- Integration corpus (`tools/integration_corpus.py`): `{"COND": 20,
  "DIFF": 10, "EQUIV": 8, "MISSING": 1, "SAME": 1171}`, the same as the
  baseline.
- ODE corpus (`replay-ode.py`, 16 jobs), base and after run on the same day
  (`raw/ode_compare.txt`):
  - no case has a different symox observation ;
  - 3 cases hit a resource limit in the base run and finished in the after
    run (`fact_06` default, `lie_group_13` and `lie_group_14` hint). These
    are inconclusive ;
  - counts: classify 315 equal in both runs ; hint 279 equal before and
    281 after ; default 263 equal before and 264 after.

## Timings

Median seconds per call, machine under load (load average 45 to 64 of 112
cores). Files: `raw/bench_*`.

| case | before | after | SymPy |
| --- | --- | --- | --- |
| det third4 bareiss | 0.000024 | 0.000023 | 0.001012 |
| det n4_k3 bareiss | 0.000149 | 0.000169 | 0.001299 |
| det n5_k3 bareiss | 0.003708 | 0.019214 | 0.185691 |
| det n6_k3 bareiss | 0.016291 | 0.080965 | 1.022628 |
| det n7_k0 bareiss | 0.000066 | 0.000073 | 0.003537 |
| det prec50 bareiss | 0.000029 | 0.000031 | 0.000861 |
| checksol 17799 | 0.051392 | 0.022520 | 0.001299 |
| solve 17799 | 0.031383 | 0.024748 | 0.083787 |
| checksol units 2cmx | 0.000853 | 0.000062 | 0.000063 |
| checksol units float | 0.000939 | 0.000066 | 0.000142 |
| pow E branch 2I | 0.000019 | 0.000102 | 0.000007 |
| pow E branch none | 0.000005 | 0.000006 | 0.000006 |
| pow plain | 0.000005 | 0.000006 | 0.000007 |

- Symbolic Float determinants (`n*_k3`) are about 5 times slower than
  before. The reason is the RR `cancel` that SymPy's `dotprodsimp`
  applies. Before, those entries stayed unsimplified and the results
  differed from SymPy. They are still 10 times faster than SymPy.
- `pow E branch 2I` now evaluates `factor_terms`, `log`, `sign` and `im`
  for the second branch. This happens only after the structural
  pre-check, so `pow plain` is unchanged. It is slower than SymPy on this
  case.
- `checksol 17799` is faster than before but still 17 times slower than
  SymPy. The time goes to the failing mpmath evaluations.
- lambdarr (main build, which this branch does not change; median per
  call, symox against SymPy):

  | case | symox | SymPy |
  | --- | --- | --- |
  | gamma float64, 10 elements | 0.87 us | 0.65 us |
  | gamma float64, 1000 elements | 11.4 us | 11.3 us |
  | gamma float64, 100000 elements | 1.12 ms | 1.06 ms |
  | `gamma(0.5)` scalar | 0.39 us | 0.74 us |
  | poly2 scalar | 0.37 us | 0.29 us |

## Validation

- `tests/test_smallgaps_parity.py`: 30 tests, all pass.
- `cargo test -p symox --release`: 836 passed, 1 ignored ; 1 + 5 doctests
  pass.
- `cargo clippy -- -D clippy::correctness`: clean, after the
  `cephes_gamma.rs` allow. Before, it failed on that function with
  `approx_constant`, an error that predates this change.
- Full suite (`run_suite.py`, 24 shards), on commit d73e7dc4: SUITE CLEAN.
  59531 tests collected and finished in 677 files ; 59467 passed, 60
  skipped, 4 xfailed, 0 failed ; 356 s at load average 12 to 20.
  Three earlier runs of the same build, under higher load, were not clean,
  each for a different test unrelated to this change:
  - run 1 (load about 60, 888 s) had 3 wall-clock timeouts:
    `test_fuzz_harness::test_worker_recovers_from_death_and_timeout`,
    `test_fuzz_harness::test_finite_run_and_deduplication` and
    `test_dirac_affine_integration::test_large_numeric_powers_are_preflighted_but_cheap_large_powers_remain_supported`.
    The fuzz deduplication test is flaky on the base build too (1 failure
    in 4 runs) ;
  - run 2: a 1 ms SIGALRM in
    `test_matrix_inexact_native_inverse::test_pending_interruption_propagates_from_native_elimination`
    fired before `pytest.raises` was entered and stopped that shard ;
  - run 3: the jax subprocess of
    `tests/extra/test_jax_codegen.py::test_custom_functions_and_standalone_source`
    exited with status 1 (the shard peaked at 4.5 GB, under its 6 GB cap).

  Each of these tests passes when run alone, and each passed in the other
  runs.
