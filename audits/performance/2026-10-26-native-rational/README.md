# Native rational reduction, 26 October 2026

Paired measurements compare immutable releases `576a9f76` (previous general reduction) and `3ef4ccbc` (native QQ Horowitz–Ostrogradsky selection). Three fresh processes per release alternate order on Demeter CPU24, with one BLAS/OMP thread. Each case contains nine distinct preconstructed proper fractions; output construction and cancellation remain inside the public `ratint_ratpart` call. The table uses medians of process medians.

| Case | Before | After | Before / after |
| --- | ---: | ---: | ---: |
| Repeated linear pole, order2, plus quadratic | 2.683 ms | 0.205 ms | 13.1× |
| Repeated linear pole, order6, plus quadratic | 6.980 ms | 0.535 ms | 13.1× |
| Repeated linear pole, order12, plus quadratic | 15.475 ms | 3.447 ms | 4.49× |
| Square-free denominator, degree12 | 3.699 ms | 0.180 ms | 20.6× |

All **216 exact independent checks** pass: the rational part's derivative plus the residual equals the input fraction; the residual denominator is square-free; both components equal the independent SymPy reduction. Every case's output digest agrees across all six processes. No reference library is imported during the native calculations; verification occurs after the timers. No debug timings or result caches are used.

The native kernel selects stock proper QQ/ZZ fractions with denominator degree at most128 and preserves complete existing routes outside that selection. These cases do not establish a speedup for the remaining logarithmic-construction algorithm, other domains or whole integration.

Raw process samples, binary hashes and digests are in `before-*.json`, `after-*.json` and `comparison.json`; actual runners are included. The additional seed-based 96-case differential was run on a private debug candidate before the release build, and is labelled accordingly in `fuzz-exact.json`. Its script is included for reproduction, not as a performance measurement. [Full Demeter validation](../../sympy-1.14.0/validation-2026-10-26-rational.json) covers all 46,421 collected Python tests, Rust tests, Clippy and the corrected compiler PATH. No global parity is certified.
