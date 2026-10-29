# arrayrust (wave 26): native matrix-lines derivative

Harness: `bench.py <sympy|symox> <native|python> <case|all> [reps]` and
`run.sh <worktree> <outdir>` (3 fresh processes per case and mode, each
reporting the first call and the median of repeated calls; 20 repetitions for
SymPy, 200 for symox). Release build of this tree. "before" is the same
binary with the native kernel hidden from `matrix_derivative` (the Python
matrix-lines/array algorithm that ran before this change); "after" is the
native kernel (`crates/symox-core/src/matrixexpr/lines.rs`). SymPy 1.14 from
/tmp. Values below are medians over the 3 processes, in ms. Raw numbers in
`raw/*.jsonl`. All 10 results are identical in the three modes (strings
compared, SymPy included). Measured on 2026-10-29 with cargo test (and other
agents) running on the same machine (absolute times are noisy; ratios are
interleaved runs).

| # | case | SymPy first / median | before first / median | after first / median | after/before (first / median) | SymPy/after (median) |
|---|---|---|---|---|---|---|
| 0 | `(x.T*x).diff(x)` | 49.839 / 8.404 | 0.538 / 0.143 | 0.323 / 0.021 | 0.60 / 0.15 | 402x |
| 1 | `(x.T*A*x).diff(x)` | 53.879 / 12.472 | 0.610 / 0.195 | 0.386 / 0.041 | 0.63 / 0.21 | 306x |
| 2 | `Trace(A*X*B).diff(X)` | 47.981 / 7.006 | 0.638 / 0.399 | 0.172 / 0.027 | 0.27 / 0.07 | 255x |
| 3 | `Trace(X**3).diff(X)` | 60.589 / 17.896 | 1.249 / 0.984 | 0.171 / 0.038 | 0.14 / 0.04 | 465x |
| 4 | `Trace(X.T*A*X).diff(X)` | 55.347 / 13.734 | 1.032 / 0.730 | 0.221 / 0.049 | 0.21 / 0.07 | 278x |
| 5 | `(x.T*X**-1*y).diff(X)` | 66.656 / 22.014 | 0.591 / 0.184 | 0.346 / 0.033 | 0.58 / 0.18 | 659x |
| 6 | `(Trace(X)*Trace(A*X)).diff(X)` | 56.610 / 13.755 | 1.239 / 0.918 | 0.280 / 0.056 | 0.23 / 0.06 | 244x |
| 7 | `Trace(A*X*B*X).diff(X)` | 59.737 / 17.441 | 1.048 / 0.743 | 0.250 / 0.062 | 0.24 / 0.08 | 283x |
| 8 | `(A*X).diff(X)` (rank 4, declines) | 45.559 / 4.437 | 0.715 / 0.289 | 0.689 / 0.284 | 0.96 / 0.98 | 16x |
| 9 | `convert_array_to_matrix(convert_matrix_to_array(A*X*B))` (not ported) | 47.047 / 5.813 | 0.563 / 0.401 | 0.565 / 0.402 | 1.00 / 1.00 | 14x |

No case got slower. Case 8 (rank-4 result) declines in a structural check
(no unit dimension, about 3 microseconds) and runs the Python array
algorithm; case 9 is unchanged because the conversions are not ported (gap,
docs/algorithmic-parity.md "arrayrust"). `symox._matrix_derivatives` is
imported before timing in both modes, so first calls compare the
algorithms, not a lazy import.

Run: `bash run.sh <worktree> <outdir>`.

## SymPy 1.14 array/matrix-expression test files

Passed tests per file (arraysuite r6 baseline -> this tree; a collection
error counts the file as 0 passed, as in the baseline).

| file | before | after |
|---|---|---|
| test_applyfunc.py | 3 | 3 |
| test_array_comprehension.py | 0 | 0 |
| test_array_derivatives.py | 0 | 0 |
| test_arrayexpr_derivatives.py | 0 | 0 (collection error) |
| test_array_expressions.py | 0 | 0 (collection error) |
| test_arrayop.py | 9 | 9 |
| test_as_explicit.py | 0 | 0 (collection error) |
| test_convert_array_to_indexed.py | 0 | 0 |
| test_convert_array_to_matrix.py | 0 | 0 (collection error) |
| test_convert_indexed_to_array.py | 0 | 0 (collection error) |
| test_convert_matrix_to_array.py | 0 | 0 (collection error) |
| test_deprecated_conv_modules.py | 0 | 0 |
| test_derivatives.py | 0 | 0 (collection error) |
| test_hadamard.py | 7 | 8 |
| test_immutable_ndim_array.py | 18 | 18 |
| test_kronecker.py | 11 | 11 |
| test_matexpr.py | 23 | 25 |
| test_mutable_ndim_array.py | 15 | 15 |
| test_ndim_array_conversions.py | 1 | 1 |
| test_ndim_array.py | 4 | 5 |
| test_trace.py | 6 | 6 |
| total | 97/201 | 101/201 |
