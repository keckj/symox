# arrexprrust (wave 28): native rank-four matrix derivatives

Harness: `bench.py <sympy|symox> <case> <reps>` times one derivative call
(median over `reps` calls after a first call; `first_ms` is that first
call). `run.sh <after-python> <before-python> <outdir>` runs 3 fresh
processes per case and mode under the 6G memory cap, from /tmp. SymPy 1.14
does 20 reps per process and symox 200. "before" is the base dba30517 tree
with its binary, "after" is this branch. Raw output is in `raw/`
(`sympy.jsonl`, `before.jsonl`, `after.jsonl`, with the printed results).
The table gives the median of the 3 per-process medians, in ms. A, B, X
are 3x3 MatrixSymbols, N and M are n x n MatrixSymbols with a
symbolic n, x is a 3x1 MatrixSymbol.

| # | call | SymPy 1.14 | before | after | before/after | SymPy/after |
|---|---|---|---|---|---|---|
| 0 | (A*X).diff(X) | 4.5368 | 0.2679 | 0.0945 | 2.84x | 48.0x |
| 1 | (X*A).diff(X) | 4.4836 | 0.2651 | 0.0932 | 2.84x | 48.1x |
| 2 | (A*X*B).diff(X) | 6.5320 | 0.2779 | 0.1003 | 2.77x | 65.1x |
| 3 | X.T.diff(X) | 0.5378 | 0.2127 | 0.0846 | 2.51x | 6.4x |
| 4 | (X*X).diff(X) | 8.7343 | 0.5187 | 0.2008 | 2.58x | 43.5x |
| 5 | (X*A*X).diff(X) | 12.4177 | 0.5214 | 0.2121 | 2.46x | 58.6x |
| 6 | (Trace(X)*A).diff(X) | 4.4722 | 0.5946 | 0.0737 | 8.07x | 60.7x |
| 7 | (A*X**-1*B).diff(X) | 22.2212 | 0.3210 | 0.1147 | 2.80x | 193.8x |
| 8 | (N.T*M*N).diff(N) | 12.8989 | 0.5629 | 0.2317 | 2.43x | 55.7x |
| 9 | (X**3).diff(X) | 17.5151 | 0.7308 | 0.2781 | 2.63x | 63.0x |
| 10 | (x.T*x).diff(x) (rank 0, control) | 8.7443 | 0.0209 | 0.0213 | 0.98x | 410.6x |
| 11 | Trace(A*X*B).diff(X) (rank 2, control) | 7.3697 | 0.0284 | 0.0283 | 1.01x | 260.8x |

Cases 0-9 go through the new native rank-four kernel
(`_core._matrix_lines_rank_four`); the controls 10-11 keep the existing
native rank <= 2 kernel and are unchanged (within noise). The results
are identical before and after for every case, and equal to SymPy 1.14
for every case except 9, where the terms come in a different order
(a pre-existing gap from ArrayContraction canonicalization, see
docs/algorithmic-parity.md, "arrexprrust (wave 28)"). The remaining
per-call cost after the change is mostly the Python construction of the
`PermuteDims`/`ArrayTensorProduct`/`ArrayAdd` objects.

Correctness probes (scratch, not committed): 60 rank-four cases
give identical output on the base and on this branch; 51 equal SymPy
1.14 exactly and are pinned in tests/test_native_arrexprrust.py.
