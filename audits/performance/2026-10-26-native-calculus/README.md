# Sparse calculus and fraction-free reduction, 26 October 2026

Paired release measurements compare `3ef4ccbc` with `146f6917` on Demeter CPU24, one BLAS/OMP thread. Three fresh processes per release alternate order. Every case uses nine distinct preconstructed inputs; output construction, public method selection and conversions remain inside each timed call. The table uses medians of process medians.

| Case | Before | After | Before / after |
| --- | ---: | ---: | ---: |
| EX three-term derivative, degree 1,000, order 1 | 379 µs | 6.79 µs | 55.8× |
| EX three-term derivative, degree 4,000, order 1 | 1,577 µs | 6.87 µs | 229.5× |
| EX three-term primitive, degree 4,000, order 1 | 3,046 µs | 10.05 µs | 303.1× |
| EX public univariate eval(0), degree 4,000, control | 184.2 µs | 185.8 µs | 0.99× |
| ZZ fraction-free dense, 24×36 | 108.822 ms | 1.770 ms | 61.5× |
| QQ fraction-free dense, 16×24 | 32.371 ms | 0.858 ms | 37.7× |
| ZZ DomainMatrix AUTO, 32×48, selected FF/sparse | 29.973 ms | 2.624 ms | 11.4× |
| QQ DomainMatrix AUTO, 24×36, selected CD/sparse | 18.602 ms | 6.348 ms | 2.93× |
| ZZ AUTO, 8×12, selected GJ/sparse, control | 1.024 ms | 1.030 ms | 0.99× |

The controls show no gain. The public univariate EX eval adapter still reconstructs an expression before substitution, bypassing the new backend zero-evaluation path in this release; its separate constructor/payload adapter fix is pending. AUTO GJ/sparse remains an existing unported algorithm. CD includes denominator clearing and other Python orchestration, so its whole-call gain is smaller than the new reduction kernel's selected dense gains.

All **216 independent sparse-calculus checks** and **270 independent fraction-free checks** pass. Calculus inputs use exact rational coefficients, so the timings exclude the separate float-grouping/precision corrections where the old result was wrong. Both releases finish and return identical output digests for every timed case. Matrix numerators divided by their nonzero denominator equal independent exact SymPy RREFs, with matching pivots, full rank and domains. Selected AUTO methods are recorded per input. Native calculations forbid SymPy imports; verification starts after all timers.

Raw samples, full binary hashes, exact-check counts, digests and actual runners are included. No resource refusal or incorrect old result is counted as a speedup. Huge-coordinate capability tests are separate from these paired timings. [Complete Demeter validation](../../sympy-1.14.0/validation-2026-10-26-calculus.json) covers all 46,647 Python tests, 661 Rust unit tests, 5 doctests and Clippy. No whole-library performance or global parity claim is made.
