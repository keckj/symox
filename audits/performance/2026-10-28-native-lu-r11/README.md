# Integrated LU release observations, 28 October 2026

Accepted corrected R11 native **e7b2975992bcffd2fd992ca11ab93f033231efe11cb679269c6ccbe906aa4bd4**, source **5f89ef4eb3c40b63fb53afb216cbb53cecb468c352a06a12d4dca53200326c35**, is compared with accepted corrected R10 native **96d342cf5fd2349f01ba38ec33a9aeec5669baada3de36e66ed75588d5447a04**, source **82feae2c9d065fd338e9439693f9c5b2bd1663b683ac2a9e411937b923e38754**. Both arms pass complete identified Demeter Python/Rust/core and binding-clippy gates. Initial failed snapshots are excluded.

Public QQ LU over the identical 18 fixtures of the earlier LU harness. Both arms positively enter the native path in separate processes. Exact factorization, permutation, aliases, publication and error controls remain certified. This measures the current metadata proof and its costs.

**1,152 exact certificates agree over 18 cases and 144 fresh processes.** Exact values, domains, aliases and input state are checked after all timers, with runtime SymPy imports blocked. Each timed process performs its real first public operation after construction, followed by nine calls on the same object, without preceding route queries, certificates or artificial cache priming. Separate instrumented route processes run after all timing processes. Three interleaved rounds pin CPU 24 with alternating arm/case order; imports and construction are outside timers. Both complete source maps and native hashes remain unchanged before/after. Agents acknowledged quiet compute before these sequential campaigns; ordinary host load remains recorded and the host is not dedicated.

R10/R11 ratios span **0.576–0.819× first**, **0.347–0.831× repeated**. A ratio above 1 means R11 is faster. Every slower median remains in the table and [explicit regression list](regressions.json); no whole-engine gain or global parity follows.

| Case | First R10 → R11 (µs) | Repeated R10 → R11 (µs) | First/repeated R10÷R11 |
| --- | ---: | ---: | ---: |
| ddm-2x2-pivot | 68.2 → 117.4 | 24.2 → 69.7 | 0.581× / 0.347× |
| domainmatrix_dense-2x2-pivot | 72.0 → 123.8 | 25.2 → 69.1 | 0.582× / 0.364× |
| domainmatrix_sparse-2x2-pivot | 106.0 → 156.9 | 48.1 → 92.5 | 0.675× / 0.520× |
| ddm-2x4-wide | 68.9 → 119.7 | 24.8 → 67.5 | 0.576× / 0.368× |
| domainmatrix_dense-2x4-wide | 74.9 → 124.4 | 26.9 → 69.8 | 0.602× / 0.386× |
| domainmatrix_sparse-2x4-wide | 112.9 → 158.7 | 51.8 → 93.4 | 0.711× / 0.555× |
| ddm-4x4-pivot | 81.2 → 136.2 | 34.9 → 78.6 | 0.596× / 0.444× |
| domainmatrix_dense-4x4-pivot | 84.9 → 142.4 | 37.3 → 79.4 | 0.596× / 0.470× |
| domainmatrix_sparse-4x4-pivot | 121.3 → 185.2 | 71.1 → 112.7 | 0.655× / 0.630× |
| ddm-4x2-tall | 78.7 → 126.7 | 30.0 → 74.0 | 0.621× / 0.406× |
| domainmatrix_dense-4x2-tall | 89.2 → 131.7 | 33.7 → 76.4 | 0.677× / 0.440× |
| domainmatrix_sparse-4x2-tall | 121.5 → 169.9 | 61.9 → 105.1 | 0.715× / 0.589× |
| ddm-8x8-full | 161.0 → 221.5 | 106.4 → 146.0 | 0.727× / 0.729× |
| domainmatrix_dense-8x8-full | 167.4 → 222.3 | 107.0 → 146.4 | 0.753× / 0.731× |
| domainmatrix_sparse-8x8-full | 233.9 → 289.4 | 174.0 → 217.1 | 0.808× / 0.802× |
| ddm-8x8-rank6 | 130.9 → 186.4 | 85.9 → 122.4 | 0.702× / 0.702× |
| domainmatrix_dense-8x8-rank6 | 137.8 → 192.3 | 83.0 → 125.4 | 0.716× / 0.662× |
| domainmatrix_sparse-8x8-rank6 | 209.2 → 255.5 | 154.6 → 186.2 | 0.819× / 0.831× |

The safe raw archive includes scripts/fixtures, every process sample/stderr, identities, commands and route/certificate observations. Wrong, incomplete, held, DEBUG and refused baselines are never counted as gains. These bounded operation comparisons do not establish import, end-to-end program or whole-engine acceleration.
