# arr2mat (wave 28, third part): native convert_array_to_matrix and _EditArrayContraction

Branch feat-arr2mat-20261030k, base 21068f42. The port is described in
docs/algorithmic-parity.md, section "arr2mat (wave 28, third part)".

## Binaries

- "after": this branch, built with
  `cargo build --release -p symox-py`. The binary is
  `python/symox/_core.abi3.so` of the worktree.
- "before": the base 21068f42. Its `python/` tree is a clean checkout of
  21068f42, and its `_core.abi3.so` is a release build of the same commit
  taken from another worktree (`arrindex`, clean at 21068f42). Before use,
  that binary was checked to export the same `_core` symbols as the
  21068f42 Python tree expects. The base was not rebuilt in its own
  worktree.
- SymPy: 1.14.0 in `~/symox-ci/venv`.

## Timings

Harness: `bench_a2m.py <sympy|symox> <reps>` prints, for each case, the
first call and the median of `reps` calls in ms, with the printed result.
`run.sh <after-python> <before-python> <outdir>` runs 3 fresh processes per
mode under a 6 GB cap, from /tmp, with 50 reps for SymPy and 200 for
symox. The table gives the median of the 3 per-process medians. Raw data:
`raw/sympy.jsonl`, `raw/before.jsonl`, `raw/after.jsonl`.

M, N and P are k x k MatrixSymbols, a and b are k x 1, I = Identity(k).
`cm` is convert_matrix_to_array, `ca` is convert_array_to_matrix, and
`tp`/`con`/`dia` are the canonicalizing constructors.

| case | SymPy 1.14 (ms) | before (ms) | after (ms) | SymPy/after | before = SymPy? |
|---|---|---|---|---|---|
| `ca(cm(M*N))` | 3.8779 | 0.0192 | 0.0953 | 40.7x | yes |
| `ca(cm((M*N).T))` | 3.8923 | 0.0196 | 0.0982 | 39.6x | yes |
| `ca(cm(M*N*P))` | 5.8191 | 0.0244 | 0.1408 | 41.3x | yes |
| `ca(cm(a.T*M*b))` | 5.9030 | 0.0250 | 0.1390 | 42.5x | yes |
| `ca(cm(Trace(M*N)))` | 3.4184 | 0.0937 | 0.2013 | 17.0x | yes |
| `ca(cm(hadamard(M,N)))` | 0.3035 | 0.1457 | 0.0157 | 19.3x | no: `FunctionMatrix(k, k, Lambda(...))` |
| `ca(cm(M*N + N.T*P))` | 7.9110 | 0.0434 | 0.2016 | 39.2x | yes |
| `ca(tp(a,b))` | 0.1495 | 0.0435 | 0.0120 | 12.5x | no: unconverted |
| `ca(dia(tp(M,a),(1,2)))` | 6.1746 | 0.0010 | 0.1422 | 43.4x | no: unconverted |
| `ca(con(tp(M,I,N),(1,2),(3,4)))` | 4.5186 | 0.0239 | 0.1080 | 41.8x | yes |
| `ca(con(tp(M,N,P),(0,2),(1,5)))` | 5.7985 | 0.0247 | 0.1399 | 41.4x | yes |
| `_EditArrayContraction` round trip | 0.4958 | NotImplementedError | 0.0657 | 7.5x | no: stub |

After this change, all 12 results are identical to SymPy 1.14.

Before, the conversion was a short Python pattern match. It was not
SymPy's algorithm: it returned a wrong form or left the input unconverted
on 3 of the 11 conversions, and the editor was a stub. On the 8 cases it
handled like SymPy, it was 2x to 7x faster than the port, which runs
SymPy's full pass pipeline.

The port is 7x to 43x faster than SymPy 1.14. A profile of
`ca(cm(M*N))` puts 85% of the time in the Rust call: canonical forms,
`MatMul(...).doit()` and Identity recognition. The Python encoding and
rebuild take the rest.

## Differential corpus

`diff_a2m.py <sympy|symox> <out.jsonl>` generates 600 inputs from a fixed
seed:

- 400 random contractions, diagonals, permutations and sums over
  MatrixSymbols, Identity, ArraySymbol, OneArray and scalars;
- 200 `convert_matrix_to_array` images of random matrix expressions.

For each input it records str(result) or the exception type for:

- `convert_array_to_matrix` and the passes `_array2matrix`,
  `_remove_trivial_dims`, `identify_hadamard_products`,
  `identify_removable_identity_matrices`, `remove_identity_matrices`,
  `_array_contraction_to_diagonal_multiple_identity` and
  `_array_diag2contr_diagmatrix`;
- the `_EditArrayContraction` state (contraction indices, `_ArgE`
  elements and indices, relative positions, diagonal count), followed by
  `merge_scalars` and `to_array_contraction`.

In symox mode, a result computed by the stock fallback (symox._arrayexpr)
is marked ` [stock]`. `compare.py` tabulates the agreement. Raw data:
`raw/diff_sympy.jsonl`, `raw/diff_symox.jsonl`, `raw/compare.txt`.

The 73 inputs that cannot be constructed raise the same exception in both
libraries (no input differs). Over the 527 others:

| pass | agree / total | native (agree) | stock (agree) |
|---|---|---|---|
| convert | 517/527 | 441 (440) | 86 (77) |
| a2m | 517/527 | 441 (440) | 86 (77) |
| rtd | 526/527 | 500 (500) | 27 (26) |
| hadamard | 521/527 | 481 (481) | 46 (40) |
| removable | 527/527 | 508 (508) | 19 (19) |
| remove_id | 527/527 | 508 (508) | 19 (19) |
| multi_id | 527/527 | 508 (508) | 19 (19) |
| diag2contr | 524/527 | 524 (524) | 3 (0) |
| editor | 312/312 | 299 (299) | 13 (13) |

The one native disagreement is arr47, a scalar placement in a held MatMul:

- SymPy: `Trace(I*Trace(a*c))*Trace(I)`
- symox: `Trace(Trace(a*c)*I)*Trace(I)`

The stock-path inputs are those where SymPy puts an array object
(OneArray, ArraySymbol) inside a matrix operation. The random generator
draws these often: 3 of its 14 factors are array objects. The
disagreements on that path (arr12, arr15, arr34, arr127, arr145, arr158,
arr249, arr300, arr308, arr336) are listed in `raw/compare.txt` and in the
parity notes.

`pick_cases.py` selects the corpus cases pinned by
tests/test_native_arr2mat.py. It keeps native-path agreements whose input
string rebuilds the same expression with the plain constructors.

## SymPy array-expression tests

The 21 SymPy 1.14 test files of `sympy/tensor/array/tests`,
`sympy/tensor/array/expressions/tests` and the matrix-expression files
that use array expressions were run with imports rewritten to symox.
Results per file: `raw/sympy_tests_before.txt` and
`raw/sympy_tests_after.txt`. In total, 120 of 201 tests passed before and
151 after:

| file | before | after |
|---|---|---|
| expressions/test_array_expressions.py | 13/25 | 24/25 |
| expressions/test_as_explicit.py | 0/2 | 2/2 |
| expressions/test_convert_array_to_matrix.py | collection error | 16/16 |
| matrices/expressions/test_applyfunc.py | 3/5 | 4/5 |
| matrices/expressions/test_derivatives.py | 5/12 | 6/12 |
| the 16 other files | unchanged | unchanged |

With the stock fallback disabled (a pytest plugin that makes
`_stock.to_stock` raise), test_convert_array_to_matrix.py still passes
16/16, test_array_expressions.py 24/25 and test_as_explicit.py 2/2.

The remaining failure of test_array_expressions.py
(`test_arrayexpr_contraction_construction`) needs
`_get_contraction_tuples` (convert_indexed_to_array), which is outside this
lot.

## Full suite

`run_suite.py <repo> <outdir> 24` runs every file under tests/ (677 files)
in 24 disjoint shards and sums the junit results. Shards 0 to 22 run under
a 6 GB address-space cap. Shard 23 is tests/extra, which runs without the
cap: under it, jax aborted the process (`Fatal Python error: Aborted`, no
junit), so it was rerun separately without the cap.

| part | passed | failed | skipped |
|---|---|---|---|
| shards 0-22 (first run) | 59153 | 4 | 48 |
| tests/extra (rerun without the cap) | 367 | 0 | 16 |

The 4 failures of the first run:

- `tests/test_matrix_calculus.py::test_trace_cycles_rectangular_products_and_explicit_lowering`
  asserted `trace(A*B) == trace(B*A)` and `trace(A*B) == trace((A*B).T)`.
  SymPy 1.14 gives False for both, and this lot aligned `Trace(MatMul)`
  with `MatMul._eval_trace`. The test now pins SymPy's result. After the
  change the file passes 21/21. The alignment awaits user validation; see
  the parity notes.
- Three tests of `tests/test_downstream_audit.py` start a subprocess with
  `PYTHONPATH` set to `tools/downstream` only, so `symox` must be
  installed in the interpreter. It is not installed in the CI venv, and
  the rules forbid a shared install. They fail the same way on the base
  21068f42. With a private venv overlay in the scratchpad (a `.pth` file
  pointing to this worktree's `python/`), the file passes 14/14.

After these two points, the suite has 0 failures. Shards 0-22 hold
59157 tests (passed and failed). Without the 87 new tests of
tests/test_native_arr2mat.py, that leaves 59070, against 59068 passed for
the main reference without tests/extra. The difference of 2 was not
investigated; the base suite was not rerun here. No second full run was
made after the test change. Raw data: `raw/suite_after_summary.json`,
`raw/suite_shards.json`, `raw/suite_after_extra.txt`.
