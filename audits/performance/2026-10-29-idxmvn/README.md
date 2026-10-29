# idxmvn wave: conjugate, re, im, Abs of abstract matrices (2026-10-29)

Harness: `bench.py` (12 public calls, 3x3 MatrixSymbols A, B and the Symbol x).
First call: median over 5 fresh processes of the first call after the
imports and symbol setup. Repeated: median of 7 runs of 200 calls, per call.
Release build, CPython 3 from the symox venv, run from /tmp under
`systemd-run --scope -p MemoryMax=6G`. The machine was shared with a
single-job `cargo test` build during the run.

- before: `git archive` of base 65a3f2af (python/) with the base release
  `_core.abi3.so` (`raw-before.jsonl`).
- after: branch feat-idxmvn-20261029e, release build (`raw-after.jsonl`).

Commands:

    PYTHONPATH=<before tree>/python python bench.py before > raw-before.jsonl
    PYTHONPATH=<worktree>/python python bench.py after > raw-after.jsonl

| Call | First before (us) | First after (us) | Repeated before (us) | Repeated after (us) | Before result | After result |
|---|---:|---:|---:|---:|---|---|
| `conjugate(A)` | 35.8 | 31.5 | 3.1 | 3.2 | `MatrixConjugate(A)` | `Adjoint(A.T)` |
| `conjugate(x*A)` | 95.6 | 62.2 | 6.4 | 6.1 | `conjugate(x)*MatrixConjugate(A)` | `Adjoint((x*A).T)` |
| `conjugate(A**2)` | 52.6 | 43.8 | 4.7 | 4.8 | `MatrixConjugate(A)**2` | `Adjoint(A.T)**2` |
| `str(sin(A - A.T))` | 381.4 | 399.9 | 20.8 | 22.7 | `sin(-A.T + A)` | `sin(A - A.T)` |
| `srepr(sin(2*A - B))` | 435.1 | 443.0 | 42.8 | 45.0 | `sin(MatAdd(MatMul(Integer(-1), MatrixSymbol("B", Integer(3),` | `sin(MatAdd(MatMul(Integer(2), MatrixSymbol("A", Integer(3), ` |
| `re(A)` | 246.4 | 390.9 | 3.2 | 54.1 | `!! NotImplementedError: re of abstract matrices is not imple` | `(1/2)*Adjoint(A.T) + (1/2)*A` |
| `im(x*A)` | 264.3 | 458.2 | 6.1 | 74.9 | `!! NotImplementedError: im of abstract matrices is not imple` | `-(I/2)*(x*A + (-conjugate(x))*Adjoint(A.T))` |
| `Abs(A)` | 219.0 | 284.5 | 3.1 | 11.0 | `!! NotImplementedError: Abs of abstract matrices is not impl` | `(A*Adjoint(A.T))**(1/2)` |
| `Abs(x*A*B)` | 252.0 | 401.9 | 6.9 | 49.6 | `!! NotImplementedError: Abs of abstract matrices is not impl` | `Abs(x)*(A*Adjoint(A.T))**(1/2)*(B*Adjoint(B.T))**(1/2)` |
| `re(log(A))` | 86.3 | 364.1 | 5.7 | 24.3 | `re(log(A))` | `log((A*Adjoint(A.T))**(1/2))` |
| `Heaviside(A) error` | 309.0 | 1100.3 | 7.8 | 421.2 | `!! NotImplementedError: Heaviside of abstract matrices is no` | `!! ValueError: Function defined only for Real Values. Comple` |
| `x + A error` | 32.5 | 31.7 | 2.9 | 3.0 | `!! TypeError: Expected a matrix expression or explicit matri` | `!! TypeError: Mix of Matrix and Scalar symbols` |

Reading:

- Before, `re`, `im`, `Abs` and `Heaviside` of a MatrixSymbol raised
  NotImplementedError, and `re(log(A))` stayed unevaluated. Their "before"
  times are the time to raise, so they are not comparable speeds: the "after"
  rows compute the SymPy 1.14 result.
- `conjugate` now builds SymPy's `Adjoint(Transpose(M))` natively at the same
  cost as the old MatrixConjugate node.
- Ordering MatAdd terms in the printer (`str`, `srepr` of a function of a
  MatAdd) costs about 2 us on a 20-45 us call.
- `Heaviside(A)` builds `im(A)` and its repr for the error message, as SymPy
  does. That is the 0.4 ms.
- `x + A` raises SymPy's TypeError at the same cost as the old message.
