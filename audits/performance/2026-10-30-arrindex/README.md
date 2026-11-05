# arrindex (arrexprrust, third part): native array_derive, convert_indexed_to_array, convert_array_to_indexed

Files:

- `corpus.py`: differential corpus (57 `array_derive`, 86
  `convert_indexed_to_array`/`_convert_indexed_to_array`, 25
  `convert_array_to_indexed` cases), serialized structurally.
  `python corpus.py <sympy|symox>` prints the JSON results.
- `sympy_expected.json`: SymPy 1.14 results; `symox_before.json`: the base
  tree (21068f42 with its binary); `symox_after.json`: this branch.
- `bench.py <sympy|symox> <reps>`: one JSON line per timed call (median
  over `reps` calls after a first call). `run.sh <after-python>
  <before-python> <outdir>` runs 3 fresh processes per mode under the 6G
  address-space cap, from /tmp; SymPy 1.14 does 30 reps per process, symox
  200. `table.py raw` prints the table below (median of the 3 per-process
  medians, in ms). Raw output: `raw/`.

## Results against SymPy 1.14

| section | cases | before = SymPy | after = SymPy |
|---|---|---|---|
| array_derive | 57 | 0 | 57 |
| convert_indexed_to_array | 86 | 26 | 82 |
| convert_array_to_indexed | 25 | 4 | 23 |

Before, `array_derive` returned an unevaluated `_ArrayDerivative` for
every case and the Python `convert_indexed_to_array` raised
NotImplementedError outside a few patterns. The 6 remaining differences
come only from the order of commutative factors in symox's core Mul/Add
(docs/algorithmic-parity.md, section arrindex): `M[i, 0]*N[0, j]`,
`M[i, j] + v[i, 0]`, `IA[i, 0]*IB[i]` and `C[2]*B[0, 1]`-style products.

## Timings

| # | call | SymPy 1.14 | before | after | SymPy/after | before/after | result vs SymPy |
|---|---|---|---|---|---|---|---|
| 0 | `array_derive(X, X)` | 0.0942 | 0.0044 | 0.0940 | 1.0x | 0.0x | same |
| 1 | `array_derive(X.T, X)` | 0.1740 | 0.0045 | 0.0950 | 1.8x | 0.0x | same |
| 2 | `array_derive(X*Y*X, X)` | 2.6876 | 0.0061 | 0.3530 | 7.6x | 0.0x | same |
| 3 | `array_derive(X**-1, X)` | 1.3775 | 0.0064 | 0.1714 | 8.0x | 0.0x | same |
| 4 | `array_derive(tp(A3, B), B)` | 0.5084 | 0.0052 | 0.1637 | 3.1x | 0.0x | same |
| 5 | `array_derive(con(tp(X, Y, X), (1, 2), (3, 4)), X)` | 2.1921 | 0.0190 | 0.3211 | 6.8x | 0.1x | same |
| 6 | `c(Sum(M[i,j]*N[j,k], j))` | 1.2342 | 0.2417 | 0.3172 | 3.9x | 0.8x | same |
| 7 | `c(Sum(M[i,j]*N[j,k]*P[k,l], j, k))` | 1.8204 | 0.3293 | 0.4178 | 4.4x | 0.8x | same |
| 8 | `c(M[i,j]*N[k,l] + N[i,j]*M[k,l], [i,j,k,l])` | 2.1898 | NotImplementedError | 0.3076 | 7.1x | - | same |
| 9 | `c(E[i,j] + F[j,i])` | 0.2906 | 0.2063 | 0.3418 | 0.9x | 0.6x | same |
| 10 | `a2i(con(tp(M, N), (1, 2)), [i, j])` | 0.3016 | 0.2812 | 0.2654 | 1.1x | 1.1x | same |
| 11 | `a2i(pd(tp(M, N), [1, 0, 3, 2]), [i, j, k, l])` | 0.1405 | 0.1577 | 0.2126 | 0.7x | 0.7x | same |
| 12 | `a2i(dia(tp(M, N), (0, 2)), [i, j, k])` | TypeError | 0.1472 | TypeError | - | - | same |

Median of 3 process medians, ms (2026-10-30, shared 112-core host, load
between 13 and 40). "result vs SymPy" compares the after result with SymPy's.

- "before" for `array_derive` only builds the unevaluated `_ArrayDerivative`
  object, so its time is not a derivative; the before/after ratio is not
  meaningful for rows 0-5.
- Row 8: the base tree raised NotImplementedError. Row 9: the base tree
  returned the unconverted expression.
- Row 12: SymPy 1.14 raises `TypeError: list indices must be integers or
  slices, not NoneType` (a bug in its ArrayDiagonal index plan); symox now
  raises the same error, the base tree returned an expression.
- Rows 9 and 11 are slower than SymPy (0.9x, 0.7x): the encoding of the
  objects and the rebuilding of the result in Python dominate these small
  inputs.

M, N, P are
d x d MatrixSymbols, X, Y k x k, A a (k, k, k) and B a (k, k) ArraySymbol,
E (3, 2) and F (2, 3) ArraySymbols.
