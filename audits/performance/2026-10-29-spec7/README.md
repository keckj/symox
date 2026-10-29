# spec7: scalar functions of matrix expressions, Str Markov states

Harness: `bench.py` (one case, `fresh` = one call in a new process,
`repeat` = median of repeated calls), `drive.py` (7 fresh processes, median,
plus one repeat run per case). Release builds. Before is the base tree
(82d674aa, extension `_core.base.so`), after is this branch (build5).

    systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
      .venv/bin/python drive.py before <base>/python before.json
    systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0 \
      .venv/bin/python drive.py after <worktree>/python after.json

Raw numbers: `before.json`, `after.json` (seconds). The machine was shared
with other builds (load average 5-7), so differences under about 15% on
fresh calls are noise.

| case | before fresh | after fresh | before repeat | after repeat |
|---|---|---|---|---|
| `sin(A)` | refused NotImplementedError | 179.7 us | refused NotImplementedError | 5.2 us |
| `cos(-A)` | refused NotImplementedError | 213.9 us | refused NotImplementedError | 9.2 us |
| `exp(log(A))` | 85.0 us | 79.4 us | 6.3 us | 6.9 us |
| `sqrt(A)` | 36.9 us | 29.0 us | 1.6 us | 7.0 us |
| `sin(x).subs(x,A)` | 136.4 us | 148.1 us | 20.0 us | 20.2 us |
| `Function(f)(A)` | 100.1 us | 111.9 us | 4.9 us | 5.1 us |
| `2*x*A` | 54.0 us | 65.9 us | 3.7 us | 3.8 us |
| `A**2*B` | 25.8 us | 30.2 us | 2.0 us | 2.1 us |
| `ZeroMatrix**(1/2)` | 22.2 us | 26.5 us | 7.1 us | 6.5 us |
| `sin(x)*cos(x)+x**2` | 111.2 us | 106.4 us | 13.3 us | 13.3 us |
| `srepr(Tuple)` | 108.8 us | 100.8 us | 19.8 us | 18.8 us |
| `Markov Str P` | refused NotImplementedError | 1583.7 us | refused NotImplementedError | 945.1 us |
| `Markov int P` | 1396.2 us | 1408.5 us | 767.3 us | 782.7 us |

Notes:

- `sin(A)`, `cos(-A)` and the `Str` Markov query were refused before and are
  now built as in SymPy 1.14.
- `sqrt(A)` and `exp(log(A))` built wrong objects before (`sqrt(A)` was a scalar
  `Pow`, `exp(log(A))` a scalar `Expr`); they are now the reference `MatPow`
  and `A`. The repeat cost of `sqrt(A)` (1.6 -> 7.0 us) is the cost of
  building the `MatPow`.
- `ZeroMatrix**(1/2)` is now the ZeroMatrix (it was `MatPow(0, 1/2)`).
- The unchanged paths (`subs`, `Function(f)(A)`, MatMul, scalar arithmetic,
  srepr, integer Markov chains) are within noise.
