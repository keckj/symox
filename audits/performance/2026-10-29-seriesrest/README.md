# seriesrest (wave 22) timings, 2026-10-29

Release builds of the base (c7289185) and new `_core.abi3.so`, same Python
sources otherwise (`git archive c7289185 python` for the base). Harness:
`bench.py` (10 calls) and `run.sh <python dir> <out>`:

- fresh: 5 fresh processes per call, one call each, the median wall time of
  the call (import excluded);
- repeat: 7 calls in one process, the median.

Raw results: `base.tsv`, `new.tsv` (columns: mode, call, seconds, result).

| call | base fresh | new fresh | base repeat | new repeat | base result | new result |
|---|---|---|---|---|---|---|
| floor(x+1/2)-floor(x) at oo | 1.64 ms | 2.60 ms | 0.57 ms | 1.13 ms | Limit | AccumBounds(-1/2, 3/2) |
| ceiling(x+1/2)-ceiling(x) at oo | 1.72 ms | 2.64 ms | 0.58 ms | 1.11 ms | Limit | AccumBounds(-1/2, 3/2) |
| floor(x)-x at oo | 1.56 ms | 2.14 ms | 0.50 ms | 0.72 ms | Limit | AccumBounds(-1, 0) |
| Mod(x, 3) at oo | 1.69 ms | 2.32 ms | 0.47 ms | 0.78 ms | Limit | AccumBounds(0, 3) |
| floor(x)/x at oo | 1.73 ms | 1.18 ms | 0.57 ms | 0.24 ms | Limit | 1 |
| x*floor(3/x)/2 at 0+ | 1.89 ms | 1.28 ms | 0.68 ms | 0.32 ms | Limit | 3/2 |
| Mod(x, b) at n*b+ | 1.47 ms | 0.83 ms | 0.56 ms | 0.17 ms | Limit | 0 |
| (-1)**x*2**(-x) at oo | 1.89 ms | 1.26 ms | 0.74 ms | 0.38 ms | Limit | 0 |
| (-1/4)**x at oo | 2.55 ms | 1.92 ms | 1.11 ms | 0.87 ms | Limit | 0 |
| Piecewise((x, y<1), (x**5, True)).series(x, oo, 3) | 0.62 ms | 2.68 ms | 0.13 ms | 1.21 ms | NotImplementedError | the Piecewise |

The base did not compute any of these calls. Its times are the time it took
to give up: an unevaluated Limit, or NotImplementedError for the series. The
AccumBounds results now take about 1 ms longer than giving up did, because
they need the bounded-affine pass. Every call now returns SymPy 1.14's value
in a few milliseconds.
