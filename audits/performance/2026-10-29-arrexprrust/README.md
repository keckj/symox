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

## Second part (2026-10-30): native canonical forms and convert_matrix_to_array

Harness: `bench_canon.py <sympy|symox> <reps>` (median ms per call after a
first call), 3 fresh processes per mode under the 6G cap, from /tmp;
SymPy 1.14 does 50 reps, symox 200. "before" is the base tree, where
`_array_tensor_product`/`_array_diagonal` do not exist (TypeError on the
call) and `convert_matrix_to_array` was the old Python conversion, which
returned different results (ArraySymbols, unconverted Trace and Hadamard
products). Raw output: `raw/canon/`. The table gives the median of the 3
per-process medians, in ms; M, N, P, Q are k x k MatrixSymbols, a is k x 1.

| # | call | SymPy 1.14 | before | after | SymPy/after | result vs SymPy |
|---|---|---|---|---|---|---|
| 0 | `con(con(tp(M,N,P),(1,2)),(0,1))` | 0.2708 | TypeError | 0.1143 | 2.4x | same |
| 1 | `pd(con(tp(M,N,P),(1,2)),[3,2,1,0])` | 0.3249 | TypeError | 0.0968 | 3.4x | same |
| 2 | `tp(Q,con(tp(M,N),(1,2)))` | 0.1997 | TypeError | 0.0720 | 2.8x | same |
| 3 | `dia(pd(tp(M,N),[1,2,3,0]),(0,1))` | 0.1772 | TypeError | 0.0782 | 2.3x | same |
| 4 | `cm(M*N*P)` | 0.3888 | 0.2585 | 0.0252 | 15.4x | same |
| 5 | `cm((M*N).T)` | 0.7788 | 0.2145 | 0.0255 | 30.5x | same |
| 6 | `cm(Trace(M*N))` | 0.3355 | 0.0027 | 0.0542 | 6.2x | same |
| 7 | `cm(a.T*M*a)` | 0.9408 | 0.2779 | 0.0329 | 28.6x | same |
| 8 | `cm(hadamard(M,N,P))` | 0.1741 | 0.0023 | 0.0266 | 6.5x | same |

Before, cases 6 and 8 returned the input unchanged (`Trace(M*N)`,
`HadamardProduct(M, N, P)`), and 4, 5, 7 returned ArraySymbol-based or
unsorted results, so the "before" times there are not comparable. In the
canonical-form calls the Rust call is a small part of the time; most of
it is encoding and rebuilding the Python objects. A first profile showed
`_shape` re-sympifying every dimension (11 us per dimension); `_shape` now
keeps native expressions as they are (sympify returned them unchanged),
which halved these calls and also speeds up the rank-four matrix
derivatives (same harness as above, after: case 0 0.0488 ms, case 6
0.0476 ms, case 9 0.0846 ms). The controls are unchanged:
`(x.T*x).diff(x)` 0.0215 ms, `Trace(A*X*B).diff(X)` 0.0281 ms.

Final validation on cbe35bbf: rshards suite 57819 passed, 0 failed
(57752 before); `cargo test -p symox --release --lib` 797 passed, clippy
correctness clean; integration corpus {"COND": 20, "DIFF": 10,
"EQUIV": 8, "MISSING": 1, "SAME": 1171} and ODE corpus {'cls': 303,
'hint': 279, 'def': 274, 'n': 315}, both equal to the baselines.
