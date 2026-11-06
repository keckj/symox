# ode3 (wave 42): native nth_algebraic, factorable and 2nd_linear_bessel

Branch `feat-ode3-20261106` from main `e4b1d011`. Reference: SymPy 1.14.0
(`~/symox-ci/venv`).

## What was ported

| Reference (SymPy 1.14 `solvers/ode/single.py`) | Rust | Binding |
|---|---|---|
| `NthAlgebraic._matches` | `symox::ode::algebraic::nth_algebraic` | `_core.ode_nth_algebraic` |
| `Factorable._matches` | `symox::ode::algebraic::factorable` | `_core.ode_factorable_match` |
| `Factorable._get_general_solution` | `symox::ode::algebraic::factorable_solutions` | `_core.ode_factorable_solutions` |
| `SecondLinearBessel._matches` | `symox::ode::bessel::bessel_match` | `_core.ode_bessel_match` |
| `SecondLinearBessel._get_general_solution` | `symox::ode::bessel::bessel_solution` | `_core.ode_bessel_solution` |

Python keeps the admission and the call (`_classify_native.py`, a few
lines in each solver of `single.py`); the reference bodies stay as the
fallback when a kernel declines (`None`).

## Scripts

- `scripts/gen_ode3.py`: SymPy 1.14 values (`classify_ode`, `dsolve` with
  the hint) of the 17 equations of `tests/test_native_ode3.py`.
- `scripts/per_example.py <dir> <out.json>`: per-example pass/fail of the
  `nth_algebraic`, `factorable`, `2nd_linear_bessel` and `lie_group`
  example lists of SymPy's `test_single.py` (copy with `sympy` replaced by
  `symox`), each example in a forked child with a 300 s limit, the same
  check as `_ode_solver_test` (classification, expected solution,
  `checkodesol`).
- `scripts/native_share.py <dir> <out.json>`: on every example of
  `test_single.py`, counts the wrapper calls that the native kernels decide
  and those they decline (`classify_ode`, then `dsolve` with the example
  hint for the three ported hints).
- ODE corpus: `audits/performance/2026-11-05-ode/scripts/run_corpus.sh`
  and `score.py` (unchanged).

## Results

All runs on the same machine, one at a time, SymPy 1.14.0 from
`~/symox-ci/venv`. "Before" is main `e4b1d011` (its Python tree and its
`_core` build).

SymPy `test_single.py` examples, per example (`raw/sympy-examples-*.json`):

| hint | before | after |
|---|---|---|
| nth_algebraic | 23/23 | 23/23 |
| factorable | 20/20 | 20/20 |
| 2nd_linear_bessel | 12/12 | 12/12 |
| lie_group (not ported) | 18/20 | 18/20 |

The two lie_group failures (lie_group_07 `AttributeError: not a
relational`, lie_group_09 different solution) are the same before and
after.

Copied SymPy `test_single.py` (39 tests) and `test_ode.py` (36 tests)
(`raw/sympy-tests-*.log`): 53 passed, 20 failed, 1 skipped, 1 xfailed
before and after, with the same failing tests (14 in `test_single.py`, 6
in `test_ode.py`, none of them for the three ported hints).
`test_nth_algebraic`, `test_factorable`, `test_slow_examples_factorable`
and `test_2nd_linear_bessel_equation` pass.

Native share (`raw/native-share.json`): over the 315 examples of
`test_single.py` (`classify_ode` on each, then `dsolve` with the example
hint for the three hints), 791 kernel calls, all decided natively and none
declined: `nth_algebraic` 340, `factorable_match` 377,
`factorable_solutions` 23, `bessel_match` 37, `bessel_solution` 14. No
timeout or crash.

ODE corpus, 315 rows (`raw/corpus-*.tar.gz`; `run_corpus.sh` then
`score.py`): classify/hint/default 315/289/270 before and after. Per
prefix, before = after: `algeb` 23/22/23, `fact` 20/19/19, `bessel`
12/12/12, `lie_group` 20/16/14. The only row whose verdict changed is
`2nd_nonlinear_autonomous_conserved_03` default (worker failure ->
resource limit, both unscored: the known `checksol` memory blowup of
`_is_special_case_of`). algeb_04 (hint) is the order of `solve` on a
multi-generator equation; fact_11 hits the resource limit in both runs.

`tests/test_native_ode3.py`: 73 tests. The 17 equations against SymPy 1.14;
native match objects and solutions equal to the reference bodies' (native
kernels hidden) for each solver on each equation; proof that the native
path is taken (the reference `collect`/`fraction` replaced by failures);
kernel values, declines and error propagation.

Validation of the branch (evidence in `~/symox-ci/evidence/ode3-20261106`):
full suite 60 425 passed, 0 failed ("SUITE CLEAN", 8 xfailed, 45 skipped,
60 478 collected, 691 files; a first pass, kept in `suite-run1`, failed
one 5 s subprocess timeout of `test_dirac_affine_integration.py` under
machine load, which passes alone and in the second pass); `cargo test
--release` 871 + 1 + 0 + 5; `cargo clippy --release -- -A clippy::all -D
clippy::correctness` clean.
