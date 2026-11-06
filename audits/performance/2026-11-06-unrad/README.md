# Radical equations in `solve`, 2026-11-06 (`feat-unrad-20261106`)

The branches of SymPy 1.14's `_solve` that an algebraic equation with
radicals of the unknown takes are now native
(`crates/symox-core/src/solve/radical.rs`). They run in the reference order,
after the polynomial/product/trig/symbolic-power steps and before the
`_tsolve` steps (`solve_by_generator`, `invert`, Lambert, ...):

1. one generator `g` of `Poly(f_num)` contains the unknown
   (`sqrt(y)` for `y**(3/2) - x`): the roots of the polynomial in `g`, then
   the inversion `_vsolve(g - t, symbol)` of every root, through `unrad`
   and then the `_tsolve` inversion, as the reference `tsolve` pass does;
2. several generators (`y` and `sqrt(y**3)` for the ODE reproducer
   `-2*y/sqrt(y**3) - C1 - x`): the unrad fallback, with the native
   `unrad` (`solve/unrad.rs`) and a native host, then the nested
   `_vsolve` of the radical-free equation (with its change of variable when
   there is one). The nested solutions are simplified per factor as the
   default `_solve` does, unless a factor has more than two roots.

Before the lot, the native solver raised NotImplementedError on the second
kind, and the Python wrapper reached the `_tsolve` step first, which returned
`[]` for `-2*y/sqrt(y**3) = C1 + x`: a wrong claim that there is no solution.
The first kind kept a single root.

Where the native branch cannot decide (`unrad` raises or returns None, the
radical-free equation is not solved), it returns "not taken", never `[]`.
The later steps then run as before. The candidates are checked by the
existing `solve_checked`: denominators, the unknown's assumptions, and
`checksol` (`check=False` disables both, as in the reference). The values the
reference leaves unsimplified (unrad results, more than two generator roots)
are recorded natively, and the Python wrapper reads them through
`_core._solve_unsimplified(eq, x)`. This is dispatch only. The `_unrad=False`
flag is accepted.

The native `solveset` keeps its own `_solve_radical` step (`radical_step`):
its generic call to `solve` runs without these branches
(`solve::without_radical_branch`). The reference `_solveset` does not use
`_solve`'s unrad branch, and an undecided candidate must stay in a
ConditionSet (`test_issue_11534`).

## Corpus

`scripts/cases.py` holds 90 radical equations:

- the ODE reproducer and its variants (`check=False`, `dict`, `Eq`, positive
  unknown);
- the lot's completeness examples;
- the `solve` calls of `test_unrad1`, `test_unrad_slow` and
  `test_unrad_fail` in SymPy's `test_solvers.py` (Lamar, purplemath,
  posify, covariance cases);
- issues touching radicals;
- parameters, assumptions, `check=False`, `set`, `dict` and `force`.

`scripts/run_cases.py` replays them with one backend, one forked child per
case and a 60 s timeout. `scripts/compare.py` classifies the results against
SymPy 1.14.

| replay | SAME | EQUIV | SUBSET | SUPERSET | WRONG_EMPTY | DECLINED | EXTRA |
|---|---|---|---|---|---|---|---|
| before (main `5ebd60c0`) | 63 | 12 | 7 | 1 | 4 | 2 | 1 |
| after | 87 | 2 | 0 | 0 | 0 | 0 | 1 |

Changed rows (before → after):

- WRONG_EMPTY → SAME: `ode_fact`, `ode_fact_dict`, `ode_fact_eq`,
  `posify_z6i`.
- SUBSET → SAME:
  - `ode_fact_nocheck`, `sqrt_y3_minus`: the kept root 0;
  - `pow32`, `pow32_nocheck`, `pow32_pos`, `pow23`, `pow52`: all the
    generator roots.
- SUPERSET → SAME: `ode_fact_pos`. The reference check rejects the
  candidate for a positive unknown.
- DECLINED → SAME:
  - `composite_y`: the `_unrad=False` flag;
  - `mixed_set`: `set=True`.
- EQUIV → SAME (now in SymPy's form): `ode_fact_cube`, `cbrt_lin`,
  `cbrt_lin_nocheck`, `sqrt_quad`, `sqrt_quad_nocheck`, `sqrt_quad_pos`,
  `sqrt_quad_simpl`, `unrad_some`, `param_cbrt`, `mixed_force`.

Remaining:

- EQUIV `cov_not_sym`, `ratio_den`: equal values, different printed forms
  (simplification of the nested radicals).
- EXTRA `cov_fail`: SymPy raises NotImplementedError ("multiple
  generators"); the native branch declines too (`_core.solve` raises). The
  established Python resultant fallback (`_unrad_resultant`) then returns
  the checked `CRootOf` root, as before the lot (an extension, unchanged).

ODE corpus (`audits/performance/2026-11-05-ode/scripts/run_corpus.sh`):
315/287/268 before, 315/289/270 after (`raw/ode-before.score`,
`raw/ode-after.score`). fact_13 and fact_16 (hint and default) are now equal
to SymPy. The default dsolve of `2nd_nonlinear_autonomous_conserved_03/04`
goes from `inconclusive_resource_limit` to `worker_failure_review_required`
(memory exhaustion). `solve` now solves the `_is_special_case_of` system as
SymPy does, and the `checksol` that follows integrates a large expression.
SymPy exceeds its deadline on both rows, so neither verdict is scored.

Time over the 90 cases: SymPy 67.8 s, symox 6.7 s before, 6.8 s after.

## Reproduce

```sh
PYTHONPATH=<wt>/python ~/symox-ci/venv/bin/python scripts/run_cases.py symox raw/symox-after.json
~/symox-ci/venv/bin/python scripts/run_cases.py sympy raw/sympy.json
~/symox-ci/venv/bin/python scripts/compare.py raw/sympy.json raw/symox-after.json
```

`raw/` holds the three replays and both comparisons.
