# odepart: native particular solutions of the linear ODE hints (2026-10-30)

Branch `feat-odepart-20261030k`, based on main `21068f42`. Reference:
SymPy 1.14.0. Machine: demeter (112 cores, shared with other agents).

"before" is main's Python tree with main's `_core.abi3.so`
(sha256 `438c4ebb...`); "after" is this branch (`5b3834d2...`);
"stock" is the after build with `ode_uc_trialset`, `ode_uc_solve` and
`ode_vop_solve` deleted from `_core`, so the reference bodies of
`nonhomogeneous.py` run.

## What is native

| Step (SymPy 1.14 `nonhomogeneous.py`) | Where |
| --- | --- |
| `_undetermined_coefficients_match`: `test`, trial set, `x**m` resonance shift | Rust (`symox::ode::undetermined::trial_set`) |
| `_solve_undetermined_coefficients`: substitution, coefficient system, exact solve, `gsol.rhs + psol` | Rust (`symox::ode::undetermined::solve_undetermined`) |
| `_solve_variation_of_parameters`: Wronskians, integrands, held integrals, `psol` sum, `homogen_sol.rhs + psol` | Rust (`symox::ode::variation::variation_of_parameters`, `matrix::wronskian`, `binders::held`) |
| `simplify` + `trigsimp` of the Wronskian and of `psol` (`simplify=True`) | Python (public symox functions, passed to the kernel as callbacks) |
| evaluation of the held integrals (`_handle_Integral`, `integrate`) | existing path (unchanged) |
| Euler wrappers in `single.py` (indicial equation, `posify`, `subs`, `expand`, envelope) | Python (unchanged) |
| admission, decline, host-class check, `Eq(f(x), rhs)` | Python (`_odesolve_native.py`, hooks in `nonhomogeneous.py`) |

Decline = `None` from the kernel; the reference body that follows each
hook then runs unchanged.

## Parity fix in `matrix::ratfunc::dotprodsimp` (Rust)

The 2x2 and 3x3 determinants, hence `wronskian`, go through
`dotprodsimp`. SymPy's version replaces every subexpression that is not an
`Add`, `Mul` or `Pow` by a Dummy (`nonalg_subs_dummies`), cancels, and
keeps the expanded form when `cancel` changes nothing. symox cancelled the
expression itself, so `exp(-2*x)` was turned into a denominator:

    wronskian([exp(-x)*sin(3*x), exp(-x)*cos(3*x)], x)
    SymPy 1.14:   -3*exp(-2*x)*sin(3*x)**2 - 3*exp(-2*x)*cos(3*x)**2
    main symox:   (-3*sin(3*x)**2 - 3*cos(3*x)**2)*exp(-2*x)
    this branch:  -3*exp(-2*x)*sin(3*x)**2 - 3*exp(-2*x)*cos(3*x)**2

The native version now does the Dummy substitution (Dummies numbered in
SymPy's first-encounter order, which is their generator order in
`cancel`), cancels, and substitutes back. Tests: Rust
`dotprodsimp_exp_negative_argument_is_opaque`, Python
`test_wronskian_keeps_exp_negative_argument_as_factor` (Wronskian and
`dsolve(..., hint='nth_linear_constant_coeff_variation_of_parameters_Integral',
simplify=False)` of `f'' + 2*f' + 10*f = exp(x)/x`, both equal to SymPy's
strings). The kernel is shared, so the stock runs below include the fix
too; the before runs (main) do not. All numbers below, except the before
columns and the timing tables, are from the build with the fix
(`5b3834d2...`).

## Native against the reference bodies

- `trialset_exhaustive.py`: every monomial `c*x**k*exp(a*x)*T(b*x)`,
  `k <= 4`, with no homogeneous equation and with resonant and
  non-resonant operators, plus sums of two monomials: 1340 cases,
  0 mismatches, 0 declines (`raw/trialset-exhaustive.json`).
- `native_vs_reference.py 300 20261030` (rerun on `5b3834d2...`): random characteristic polynomials
  and forcing terms; `_undetermined_coefficients_match` and `dsolve` with
  the UC and VOP hints, `simplify` True/False, kernels on and removed:
  300/300 identical (match and all four hint variants). Kernel calls:
  `ode_uc_trialset` 900 native / 0 declined, `ode_uc_solve` 600 / 0,
  `ode_vop_solve` 375 / 0 (`raw/native-vs-reference.jsonl`; seed 1, 60
  cases: 60/60, `raw/native-vs-reference-seed1.jsonl`).

## Differential test against SymPy 1.14 (runner.py, compare.py)

The `undet` and `vop` families of `gen_eqs.py` (seed 20261030, the cases
of the odesolve audit), hinted, `simplify=False` and default hint. SymPy
rows: `raw/diff-ref-sympy-*.jsonl` (copied from the odesolve audit).
`vop_016` skipped (`ODE_SKIP`): symox exceeds 6 GB inside integration,
before and after; SymPy times out at 30 s. `vop_013` skipped in the after
and stock runs: after the `dotprodsimp` fix (see below) its
`simplify=False` call exceeds 6 GB inside integration (before the fix
it returned held integrals); the hinted and default calls still complete
with the before result (65 s and 48 s, run alone,
`f'' + 6*f' + 16*f = log(x)`), and the `_Integral` hint with
`simplify=False` gives SymPy's held integrals. SymPy times out at 30 s on
all three fields. The
before run still contains it (29 rows); compare.py joins on the case
name.

| Family | field | after == SymPy | before == SymPy | after == before | after == stock |
| --- | --- | ---: | ---: | ---: | ---: |
| undet | hint | 55/60 | 55/60 | 60/60 | 60/60 |
| undet | nosimp | 59/60 | 59/60 | 60/60 | 60/60 |
| undet | default | 54/60 | 54/60 | 60/60 | 60/60 |
| vop | hint | 25/28 | 25/28 | 28/28 | 28/28 |
| vop | nosimp | 19/28 | 19/28 | 28/28 | 28/28 |
| vop | default | 26/28 | 26/28 | 28/28 | 28/28 |

The remaining differences are all identical before the lot: cubic radical
forms of the homogeneous part (`undet_004`, `undet_023`, `undet_040`),
the `collect` partial exponent match (`undet_007`), SymPy's
`NotImplementedError` on casus irreducibilis roots (`undet_014`,
`undet_017`), and, for `vop`, 14 differences:

- SymPy times out at 30 s and symox answers: `vop_000`, `vop_003`,
  `vop_019`, `vop_028` (`simplify=False`);
- `Ei(-3*x)` where SymPy gives `Ei(3*x*exp_polar(I*pi))`: `vop_001`,
  `vop_026` (all three fields), from `integrate`;
- forms of evaluated integrals with `sqrt` roots: `vop_000` (hint),
  `vop_014`, `vop_025` (`simplify=False`);
- `vop_023` (`simplify=False`): the held integrals built by the kernel are
  SymPy's, character for character (checked with the `_Integral` hint);
  `Integral.doit` of the unevaluable integral then returns the integrand
  denominator `x*(sin(3*x)**2 + cos(3*x)**2)` where SymPy returns
  `x*sin(3*x)**2 + x*cos(3*x)**2`. This normalisation belongs to the
  integration path.

## SymPy test_single examples (examples.py)

The four hints (`nth_linear_undetermined_coefficients`,
`nth_linear_var_of_parameters`, `euler_undetermined_coeff`,
`euler_var_para`), 67 examples (`var_of_parameters_15` skipped: symox
exceeds 6 GB, before and after), `hint` and `simplify=False`:

| Field | after == SymPy | before == SymPy | after == before | after == stock |
| --- | ---: | ---: | ---: | ---: |
| hint | 62/67 | 62/67 | 67/67 | 67/67 |
| simplify=False | 59/67 | 59/67 | 67/67 | 67/67 |

Differences (all identical before): `undet_38`/`undet_39` (symbolic
quadratic roots; symbolic coefficients, declined), `euler_undet_03`,
`euler_undet_06`, `var_of_parameters_13` (`simplify` forms),
`var_of_parameters_01/02/06/08/13/14`, `euler_var_05` (unsimplified forms
of evaluated integrals; the held integrals are SymPy's).

## ODE corpus

`audits/sympy-1.14.0/corpus-2026-10-26/replay-ode.py`, live SymPy and
symox, 16 jobs, 60 s per call; score = `equal_output` +
`equal_error_observation` (`raw/ode-replay/`):

| | cls | hint | def |
| --- | ---: | ---: | ---: |
| before | 315 | 284 | 267 |
| after | 315 | 283 | 265 |

945/945 rows exactly once in both runs (after: native sha256
`5b3834d2...`, the build with the `dotprodsimp` fix). The three rows that
are equal before and not after are all `inconclusive_resource_limit`
after, and in each of them the SymPy side hit the 60 s deadline:
`fact_06` (default; SymPy deadline), `lie_group_13` (default) and
`lie_group_14` (hint) (deadline on both sides). SymPy's time does not
depend on the symox build, so these are machine-load effects, not
regressions. On symox's side, `lie_group_13` and `lie_group_14` were
rerun alone, before and after interleaved, twice (`scripts/lg13.py`,
`scripts/lg14.py`, `raw/ode-replay/lie-group-rerun.txt`): same results,
8.6/16.1 s before vs 8.3/10.0 s after (`lie_group_13`, default) and
14.0/13.5 s vs 16.3/13.5 s (`lie_group_14`, hint), load average about 20. No row changed from
equal to different. The baseline given for this lot was 303/279/274; it
was measured on an earlier state of main, and the before and after runs
here use the same script on the same machine.

## Timings (bench.py, ms)

After and stock rerun on `5b3834d2...` (with the `dotprodsimp` fix);
SymPy and before from the same day. Median of 3 fresh processes; "first"
is the first call after import,
"repeat" the median of 5 further calls in the same process. Results are
identical before, stock and after on every case; after equals SymPy
except `vop_sec/False` (integration form, see above).

| Case/simplify | SymPy first | before first | stock first | after first | SymPy repeat | before repeat | after repeat |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| uc_2nd/True | 318.29 | 74.48 | 73.87 | 62.89 | 170.39 | 68.30 | 59.53 |
| uc_2nd/False | 202.84 | 15.42 | 15.11 | 5.40 | 109.00 | 12.60 | 3.61 |
| uc_3rd/True | 267.41 | 70.71 | 71.23 | 64.89 | 150.66 | 66.96 | 60.94 |
| uc_3rd/False | 146.63 | 12.73 | 12.90 | 6.34 | 75.93 | 10.57 | 3.06 |
| uc_res4/True | 458.30 | 98.50 | 90.95 | 79.71 | 364.66 | 93.86 | 73.70 |
| uc_res4/False | 302.35 | 22.83 | 21.03 | 6.10 | 214.66 | 17.63 | 3.37 |
| uc_poly/True | 125.56 | 28.15 | 38.94 | 27.20 | 76.32 | 25.27 | 21.90 |
| uc_poly/False | 87.26 | 7.84 | 8.65 | 4.22 | 47.86 | 6.01 | 2.64 |
| vop_sec/True | 323.97 | 32.61 | 32.84 | 31.93 | 126.39 | 28.93 | 28.22 |
| vop_sec/False | 858.38 | 6.02 | 5.97 | 5.45 | 437.75 | 4.14 | 3.87 |
| vop_expx/True | 291.56 | 30.24 | 30.17 | 29.55 | 86.45 | 26.62 | 26.27 |
| vop_expx/False | 106.06 | 5.13 | 5.03 | 4.39 | 44.67 | 3.17 | 2.85 |
| vop_3rd/True | 379.02 | 49.34 | 43.46 | 43.29 | 122.24 | 42.08 | 39.66 |
| vop_3rd/False | 98.27 | 4.75 | 4.98 | 4.52 | 41.67 | 3.07 | 2.76 |
| euler_uc/True | 509.47 | 35.88 | 35.46 | 34.16 | 312.29 | 32.48 | 30.86 |
| euler_uc/False | 123.48 | 12.95 | 9.87 | 7.95 | 66.02 | 7.45 | 5.73 |
| euler_vop/True | 532.98 | 34.87 | 31.57 | 30.72 | 287.40 | 29.15 | 27.61 |
| euler_vop/False | 122.09 | 16.88 | 5.63 | 5.34 | 57.60 | 4.16 | 3.46 |

Helpers alone (median of 5 calls, `simplify=False` cases):

| Case, helper | SymPy | before | stock | after |
| --- | ---: | ---: | ---: | ---: |
| uc_2nd, `_undetermined_coefficients_match` | 5.95 | 5.62 | 5.63 | 0.04 |
| uc_2nd, `_solve_undetermined_coefficients` | 44.83 | 3.20 | 3.17 | 0.17 |
| uc_3rd, `_undetermined_coefficients_match` | 5.99 | 4.32 | 4.36 | 0.04 |
| uc_3rd, `_solve_undetermined_coefficients` | 62.28 | 3.28 | 3.33 | 0.21 |
| uc_res4, `_undetermined_coefficients_match` | 20.11 | 8.56 | 8.54 | 0.04 |
| uc_res4, `_solve_undetermined_coefficients` | 161.19 | 5.29 | 5.14 | 0.25 |
| uc_poly, `_undetermined_coefficients_match` | 1.36 | 1.54 | 1.58 | 0.03 |
| uc_poly, `_solve_undetermined_coefficients` | 13.07 | 1.86 | 2.01 | 0.09 |

The undetermined-coefficients hints gain 2-5x (repeat) with
`simplify=False`; with `simplify=True` the time is dominated by `odesimp`/`constantsimp`
(Python). Variation of parameters is unchanged: its time is the
`simplify`/`trigsimp` callbacks and the integration, both outside this
lot.

## Validation

Build `5b3834d2...` (with the `dotprodsimp` fix).

- `tests/` in 24 file shards (`scripts/run_shards.py`, no memory cap on
  the shards): 59594 collected, 59521 passed, 60 skipped, 5 xfailed,
  8 failed. The 8 failures are environmental and pass on rerun:
  `tests/test_downstream_audit.py` (3; they need symox installed in the
  interpreter, which this worktree must not do: 14 passed in a private
  venv with the worktree installed), `tests/extra/test_matrix_log_jax.py`
  (4, "No BLAS support for stream" under load) and
  `tests/test_fuzz_harness.py::test_finite_run_and_deduplication`
  (timeout under load); rerun: 32 passed.
- `tests/test_native_odepart.py`: 92 passed, 1 xfailed (strict xfail,
  the `f'' + f = 1/cos(x)` integral form).
- `cargo test -p symox --release`: 838 passed, 1 ignored; 1 passed;
  5 passed.
- `cargo clippy -- -D clippy::correctness`: one error,
  `approx_constant` (`EULER_GAMMA`) in `cephes_gamma.rs`, present on main
  and outside this lot; nothing in the files of this lot.
- Integration corpus (`audits/sympy-1.14.0/corpus-2026-10-26`):
  `{"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1, "SAME": 1171}`,
  the baseline.
- ODE corpus, differential test, examples and native-vs-reference: above.

## Files

- scripts/trialset_exhaustive.py, scripts/native_vs_reference.py: native
  against the reference bodies.
- scripts/gen_eqs.py, scripts/runner.py (`ODE_SKIP`), scripts/compare.py:
  differential test (from the odesolve audit; stock mode deletes the new
  kernels).
- scripts/examples.py: SymPy test_single examples (symox needs a copy of
  `sympy/solvers/ode/tests/test_single.py` with `sympy` replaced by
  `symox`).
- scripts/bench.py: timings.
- scripts/lg13.py, scripts/lg14.py: reruns of two ODE corpus rows.
- scripts/run_shards.py: full suite in file shards.
- raw/: outputs of the above.
