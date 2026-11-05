# arrexpr (wave 39): stock fallback of arr2mat, native `_chain_array_to_matrix`, array tests

Branch feat-arrexpr-20261105, base 2e623e66. The port is described in
docs/algorithmic-parity.md, section "arrexpr (wave 39)". SymPy reference:
1.14.0 in `~/symox-ci/venv`.

The "before" figures were measured on the worktree at the base commit,
with a release build of 2e623e66, before any change. The "after" figures
use this branch's release build (`python/symox/_core.abi3.so` of the
worktree).

## Measurements

| measure | before | after |
|---|---|---|
| `convert_array_to_matrix` corpus (527 cases), native path | 441 | 521 |
| same, stock (vendored SymPy) fallback | 86 | 6 |
| same, agreeing with SymPy | 517 | 520 |
| SymPy array test functions (21 files, 201) | 159 | 165 |
| arrindex corpus (168 cases) | 162 | 167 |
| `_chain_array_to_matrix` | Python | Rust (Python fallback kept) |

Corpus passes after (agree/total): rtd 526, hadamard 524, removable,
remove_id and multi_id 527, diag2contr 524, editor 312/312.

The 6 cases still on the stock path (`raw/decline.json`) each need a
matrix expression that holds an array operand:

- Hadamard product of an array expression: arr34, arr308, arr336;
- `DiagMatrix` of a non-matrix: arr15, arr145;
- transpose of a non-matrix: arr249.

At the base, 86 cases fell back (`raw/decline_before.json`). Most of them
were products of an array expression, `MatAdd` of non-matrices, or traces
of such products.

SymPy tests newly passing (no test that passed before fails; see
`raw/sympy_tests_before.txt` and `raw/sympy_tests_after.txt`):

- test_array_derivatives.py::test_array_derivative_construction;
- test_array_expressions.py::test_arrayexpr_contraction_construction;
- test_deprecated_conv_modules.py::test_deprecated_conv_module_results;
- matrices test_derivatives.py::test_derivatives_matrix_norms,
  test_derivatives_of_complicated_matrix_expr and
  test_derivatives_of_hadamard_expressions.

`_chain_array_to_matrix`: native and Python forms are identical on 3000
fuzzed array expressions and on the 78 chains the matrix-derivative tests
produce.

## Cost of routing `matrix_derivative` through the array path

SymPy 1.14's `_matrix_derivative` always goes through
`convert_matrix_to_array`, then `array_derive`, then
`convert_array_to_matrix`. symox now does the same (all three native), and
keeps its previous chain as the fallback on errors. It also declines
(an admission check) when a scalar leaf of the converted expression depends
on `x`. SymPy 1.14 would differentiate that leaf as zero: `exp(Trace(X))`,
`X[0, 1]` or `DotProduct(a, a)`. The previous chain gives the correct
explicit derivative there. Per derivative:

- native array path: 0.5 to 1.7 ms;
- previous lines kernel: 0.2 to 36 µs;
- SymPy 1.14 on a cold cache: 12 to 30 ms.

The profile is dominated by the Python `matrix_sort_key` callback of the
codec's ordering. Moving it to Rust is a performance follow-up.

## Validation

- Full suite (`audits/performance/2026-10-30-smallgaps/run_suite.py`): SUITE
  CLEAN, with 60 054 passed, 45 skipped and 5 xfailed out of 60 104. The
  base had 60 020 passed out of 60 085. The summary is in
  `raw/suite_summary.json`.
  - The run sets `XLA_PYTHON_CLIENT_PREALLOCATE=false` in the environment;
    the code and tests are unchanged.
  - Without it, `tests/extra/test_jax_codegen.py::test_custom_functions_and_standalone_source`
    fails with `CUDA_ERROR_OUT_OF_MEMORY` in its subprocess, because the
    shard's own JAX pre-reserves 75% of the GPU. The ratintex branch fails
    the same way, so the failure does not come from this change.
- `cargo test --release`: 852 passed (1 ignored), plus 1 and 5.
- Clippy correctness (`-A clippy::all -D clippy::correctness`): no error.

## Reproduction

From the worktree, with `PY=~/symox-ci/venv/bin/python` and
`A=audits/performance/2026-10-30-arr2mat`:

- `PYTHONPATH=python $PY $A/diff_a2m.py symox out.jsonl`, then
  `$PY $A/compare.py <sympy.jsonl> out.jsonl 30`, for the corpus. The SymPy
  side is `$PY $A/diff_a2m.py sympy`.
- `PYTHONPATH=python $PY decline.py $PWD decline.json`, for the reasons of
  the fallback.
- `PYTHONPATH=python $PY sympy_arr_tests.py $PWD outdir`, for the SymPy test
  files.
- `PYTHONPATH=python $PY audits/performance/2026-10-30-arrindex/corpus.py symox`,
  compared against `sympy_expected.json` of that folder, for arrindex.
