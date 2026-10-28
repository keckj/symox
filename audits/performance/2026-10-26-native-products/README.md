# Exact dense matrix kernels, 26 October 2026

These paired release measurements compare `d9b14020` (Python DDM arithmetic) with `576a9f76` (native DDM multiplication and QQ Gauss–Jordan). Both immutable packages ran on Demeter, pinned to CPU 24, with one BLAS/OMP thread. Three fresh processes per release alternated order. The table uses the median of the three process medians. Each process uses distinct preconstructed inputs; public output allocation remains inside each timed call. These are selected DDM operations, not a whole-library performance claim.

| Case | Before | After | Before / after |
| --- | ---: | ---: | ---: |
| Native Python `int`, 4×4 accumulation | 7.26 µs | 6.81 µs | 1.07× |
| Native Python `int`, 32×32 accumulation | 1.103 ms | 0.437 ms | 2.52× |
| `ZZ`, 8×8 multiplication | 2.385 ms | 0.0412 ms | 58.0× |
| `ZZ`, 32×32 multiplication | 150.326 ms | 0.731 ms | 205.6× |
| `QQ`, 8×8 multiplication | 2.438 ms | 0.0790 ms | 30.8× |
| `QQ`, 24×24 multiplication | 64.713 ms | 1.703 ms | 38.0× |
| `QQ`, 16×16 nonzero accumulation | 19.364 ms | 0.503 ms | 38.5× |
| `QQ`, 8×8 RREF | 0.666 ms | 0.0434 ms | 15.3× |
| `QQ`, 20×20 RREF | 7.814 ms | 0.287 ms | 27.3× |
| `QQ`, 12×24 RREF, rank 10 | 7.326 ms | 0.385 ms | 19.0× |

The small integer case shows little gain; no broad speedup is inferred from it. Small operations include binding, conversion and publication costs. No debug-build timings are used.

Every one of the **756 products** was checked against an independent exact SymPy matrix product, including the original nonzero accumulation target and row identities. Every one of the **216 RREFs** was checked against an independent exact SymPy RREF, including pivots, rank and the QQ output domain. All value digests match across all six processes for each case. Native workers block SymPy imports during timed calculations; reference verification starts after all timers stop.

Raw per-process observations, exact-check counts and digests are in `matmul-*.json` and `rref-*.json`. `comparison.json` includes full binary hashes and unrounded measurements. The actual measurement and orchestration scripts are archived alongside them. Source/build/test identity is recorded in [the full validation](../../sympy-1.14.0/validation-2026-10-26-products.json).

The QQ RREF port covers DDM and explicit dense Gauss–Jordan. DomainMatrix's automatic ZZ/QQ route can still choose a sparse fraction-free Python algorithm; its port is a separate task. Custom domains, aliases and changed protocols retain complete fallbacks. Unsupported computations are not counted as faster answers.
