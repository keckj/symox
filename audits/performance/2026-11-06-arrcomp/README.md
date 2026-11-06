# arrcomp (wave 41): native Kronecker, matrix-expression, trace and comprehension models

Branch feat-arrcomp-20261106, base main 5ebd60c0. The port is described in
docs/algorithmic-parity.md, section "arrcomp (wave 41)". SymPy reference:
1.14.0 in `~/symox-ci/venv`.

"Before" was measured on the worktree at the base, with the release build
of 5ebd60c0. "After" uses this branch's release build
(`python/symox/_core.abi3.so` of the worktree).

## Measurements

| measure | before | after |
|---|---|---|
| SymPy array/matrix-expression test functions (21 files, 201) passing | 165 | 182 |
| failing (out of 201; 2 more are SymPy xfails) | 34 | 17 |
| Kronecker cases (`kronecker_cases.py`), identical to SymPy | - | 7/8 |
| `convert_array_to_matrix` stock fallback (527 cases) | 6 | 6 (unchanged) |

Per-file summaries: `raw/sympy_tests_before.txt`, `raw/sympy_tests_after.txt`.
Failing tests: `raw/failed_before.txt`, `raw/failed_after.txt`. No test
that passed before fails after.

`raw/kronecker_sympy.txt` and `raw/kronecker_symox.txt` are the outputs of
`kronecker_cases.py` under SymPy and symox. The only difference (case 4) is
the deliberate one: `combine_kronecker(KroneckerProduct(X, Y)**x)` stays
unevaluated in symox. SymPy distributes a symbolic power over the factors,
which is false for non-invertible factors (see the parity document).

## Reproduction

From the worktree, with `PY=~/symox-ci/venv/bin/python`:

- `PYTHONPATH=python $PY sympy_arr_tests.py $PWD outdir` runs the SymPy
  test files, rewritten to `symox`. Copied unchanged from
  `2026-11-05-arrexpr`.
- `PYTHONPATH=python $PY kronecker_cases.py symox` runs the Kronecker cases
  under symox; `$PY kronecker_cases.py sympy`, run from outside the
  worktree, runs them under SymPy.
