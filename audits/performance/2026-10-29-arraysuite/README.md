# arraysuite timings (2026-10-29)

Before: base 3d7da725 (`git archive 3d7da725 python` with the same release
`_core.abi3.so`; no Rust change). After: feat-arraysuite-20261029g.

- `bench.py`: per case, first call and median of 10 repeated calls in one
  fresh process (symox import excluded). The `Flatten` and
  `determinant.det` cases exist only after (new names).
- `run.sh <python-dir> <runs> <out.jsonl>`: one process per run, 6G cap.
- `summarize.py before.jsonl after.jsonl`: median over the 5 runs
  (`summary.md`). Raw: `before.jsonl`, `after.jsonl`.

The committed change is API surface only (`Flatten`, determinant-module
names); the matrix-derivative and conversion rows are controls and are
unchanged within noise (the 1.50 vs 0.63 ms first call of `(A*A).diff(A)` is
one-off import order; repeated 0.554 vs 0.543 ms).

Not committed, for the record: running SymPy 1.14's Python array-expression
algorithm for `MatrixExpr.diff` (an attempt made during this task, dropped
because computation must be ported to Rust) cost 0.8x-1.3x of SymPy and
10x-70x the current path, e.g. `(x.T*x).diff(x)` 0.14 ms -> 9.8 ms
repeated (SymPy 8.3 ms). A native port of the canonicalization of
ArrayContraction/PermuteDims is the precondition for that path.
