# Integrated HNF release observations, 28 October 2026

Accepted corrected R11 native **e7b2975992bcffd2fd992ca11ab93f033231efe11cb679269c6ccbe906aa4bd4**, source **5f89ef4eb3c40b63fb53afb216cbb53cecb468c352a06a12d4dca53200326c35**, is compared with accepted corrected R10 native **96d342cf5fd2349f01ba38ec33a9aeec5669baada3de36e66ed75588d5447a04**, source **82feae2c9d065fd338e9439693f9c5b2bd1663b683ac2a9e411937b923e38754**. Both arms pass complete identified Demeter Python/Rust/core and binding-clippy gates. Initial failed snapshots are excluded.

Public ordinary ZZ column HNF over 18 exact fixtures, including empty, deficient, tall/wide and large-integer cases. Independently reconstructed unimodular V certifies A*V=[0|H], zero discarded columns and the defining HNF conditions. Corrected R10 uses its complete original HNF route; R11 native admission is positively checked separately.

**1,152 exact certificates agree over 18 cases and 144 fresh processes.** Exact values, domains, aliases and input state are checked after all timers, with runtime SymPy imports blocked. Each timed process performs its real first public operation after construction, followed by nine calls on the same object, without preceding route queries, certificates or artificial cache priming. Separate instrumented route processes run after all timing processes. Three interleaved rounds pin CPU 24 with alternating arm/case order; imports and construction are outside timers. Both complete source maps and native hashes remain unchanged before/after. Agents acknowledged quiet compute before these sequential campaigns; ordinary host load remains recorded and the host is not dedicated.

R10/R11 ratios span **0.642–15.471× first**, **0.553–18.482× repeated**. A ratio above 1 means R11 is faster. Every slower median remains in the table and [explicit regression list](regressions.json); no whole-engine gain or global parity follows.

| Case | First R10 → R11 (µs) | Repeated R10 → R11 (µs) | First/repeated R10÷R11 |
| --- | ---: | ---: | ---: |
| matrix-2x2-small | 16277.8 → 16911.0 | 428.1 → 520.0 | 0.963× / 0.823× |
| domainmatrix_dense-2x2-small | 182.3 → 284.1 | 121.5 → 219.7 | 0.642× / 0.553× |
| domainmatrix_sparse-2x2-small | 193.2 → 299.4 | 129.3 → 226.6 | 0.645× / 0.571× |
| matrix-2x4-wide | 16658.4 → 16924.6 | 736.7 → 704.0 | 0.984× / 1.046× |
| domainmatrix_dense-2x4-wide | 296.7 → 293.2 | 249.7 → 222.4 | 1.012× / 1.123× |
| domainmatrix_sparse-2x4-wide | 309.0 → 305.9 | 242.0 → 229.4 | 1.010× / 1.055× |
| matrix-4x4-signed | 17607.5 → 17443.7 | 1677.0 → 1079.1 | 1.009× / 1.554× |
| domainmatrix_dense-4x4-signed | 870.2 → 297.1 | 787.3 → 246.4 | 2.930× / 3.196× |
| domainmatrix_sparse-4x4-signed | 862.3 → 313.4 | 795.6 → 240.6 | 2.752× / 3.307× |
| matrix-4x2-tall | 16468.8 → 16954.1 | 597.0 → 669.9 | 0.971× / 0.891× |
| domainmatrix_dense-4x2-tall | 209.7 → 291.3 | 154.1 → 223.5 | 0.720× / 0.689× |
| domainmatrix_sparse-4x2-tall | 217.4 → 293.6 | 158.2 → 229.6 | 0.741× / 0.689× |
| matrix-8x8-dense | 26463.6 → 19794.2 | 10109.8 → 3502.7 | 1.337× / 2.886× |
| domainmatrix_dense-8x8-dense | 6958.3 → 459.5 | 6929.2 → 374.9 | 15.144× / 18.482× |
| domainmatrix_sparse-8x8-dense | 7381.1 → 477.1 | 6838.1 → 390.4 | 15.471× / 17.514× |
| matrix-8x8-rank6 | 24376.6 → 19951.2 | 7976.6 → 3258.0 | 1.222× / 2.448× |
| domainmatrix_dense-8x8-rank6 | 5448.9 → 391.2 | 5377.8 → 312.2 | 13.929× / 17.227× |
| domainmatrix_sparse-8x8-rank6 | 5244.0 → 396.4 | 5234.6 → 324.0 | 13.228× / 16.155× |

The safe raw archive includes scripts/fixtures, every process sample/stderr, identities, commands and route/certificate observations. Wrong, incomplete, held, DEBUG and refused baselines are never counted as gains. These bounded operation comparisons do not establish import, end-to-end program or whole-engine acceleration.
