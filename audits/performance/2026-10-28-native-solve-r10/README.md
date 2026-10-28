# Integrated SOLVE release observations, 28 October 2026

Accepted corrected R10 native **96d342cf5fd2349f01ba38ec33a9aeec5669baada3de36e66ed75588d5447a04**, source **82feae2c9d065fd338e9439693f9c5b2bd1663b683ac2a9e411937b923e38754**, is compared with accepted R9; Risch additionally includes accepted R8. All arms pass complete identified Demeter Python/Rust/core and binding-clippy gates. The failed initial R10 snapshot is excluded.

Public QQ triangular solve, including DDM and dense/sparse DomainMatrix, multiple right-hand sides and consistent overdetermined systems. Independent Fraction A*X=B and rank=n establish unique complete solutions. Other-domain and error controls remain outside timings.

**1,152 exact certificates agree over 18 cases and 144 fresh processes.** All exact values, domains, aliases and input state are verified after all timers, with runtime SymPy imports blocked. Timed processes perform a real first public operation after construction, followed by 9 calls on the same object; no preceding route query, certificate or artificial cache priming occurs. Separate instrumented processes run after every timing process. Three interleaved rounds pin CPU 24; arm/case order alternates (three-arm Latin square for Risch). Imports and construction are outside timers. All source maps and native hashes are unchanged before/after. Agents acknowledge quiet compute before the sequential campaigns; ordinary host load is retained and the host is not dedicated.

R9/R10 ratios span **0.886–3.946× first**, **1.035–4.889× repeated**. A ratio above 1 means R10 is faster. Every slower sample remains in the table and [explicit regression list](regressions.json); no whole-engine gain or global parity is inferred.

| Case | First R9 → R10 (µs) | Repeated R9 → R10 (µs) | First/repeated R9÷R10 |
| --- | ---: | ---: | ---: |
| ddm-2x2-small-rhs1 | 145.8 → 118.9 | 58.1 → 53.6 | 1.227× / 1.085× |
| domainmatrix_dense-2x2-small-rhs1 | 130.5 → 119.0 | 56.6 → 53.9 | 1.096× / 1.052× |
| domainmatrix_sparse-2x2-small-rhs1 | 145.2 → 152.3 | 78.0 → 75.4 | 0.953× / 1.035× |
| ddm-2x2-rhs3 | 167.6 → 120.8 | 91.4 → 58.1 | 1.387× / 1.574× |
| domainmatrix_dense-2x2-rhs3 | 165.5 → 121.9 | 92.6 → 56.4 | 1.357× / 1.641× |
| domainmatrix_sparse-2x2-rhs3 | 199.6 → 146.5 | 118.5 → 78.9 | 1.362× / 1.501× |
| ddm-4x4-rhs1 | 208.1 → 134.0 | 125.5 → 69.9 | 1.553× / 1.794× |
| domainmatrix_dense-4x4-rhs1 | 209.2 → 142.3 | 126.7 → 74.7 | 1.470× / 1.696× |
| domainmatrix_sparse-4x4-rhs1 | 236.2 → 174.5 | 165.5 → 101.2 | 1.354× / 1.636× |
| ddm-4x2-tall-rhs3 | 276.2 → 129.4 | 170.4 → 64.3 | 2.134× / 2.651× |
| domainmatrix_dense-4x2-tall-rhs3 | 249.5 → 139.7 | 174.1 → 64.9 | 1.786× / 2.681× |
| domainmatrix_sparse-4x2-tall-rhs3 | 278.1 → 164.4 | 193.3 → 94.1 | 1.691× / 2.055× |
| ddm-8x8-rhs3 | 1210.7 → 306.9 | 1060.8 → 226.2 | 3.946× / 4.689× |
| domainmatrix_dense-8x8-rhs3 | 1223.0 → 747.5 | 1089.5 → 222.8 | 1.636× / 4.889× |
| domainmatrix_sparse-8x8-rhs3 | 1233.3 → 788.4 | 1127.8 → 274.1 | 1.564× / 4.115× |
| ddm-8x4-tall-rhs2 | 550.0 → 182.6 | 435.8 → 106.9 | 3.012× / 4.077× |
| domainmatrix_dense-8x4-tall-rhs2 | 536.1 → 193.7 | 441.2 → 112.7 | 2.768× / 3.914× |
| domainmatrix_sparse-8x4-tall-rhs2 | 580.3 → 654.7 | 476.8 → 141.9 | 0.886× / 3.361× |

The safe raw archive contains all scripts/fixtures, process samples/stderr, identities, commands and selection/certificate observations. Wrong, incomplete, held, DEBUG and refused baselines are never counted as gains. This bounded operation comparison does not establish import, end-to-end program or whole-integrator acceleration.
