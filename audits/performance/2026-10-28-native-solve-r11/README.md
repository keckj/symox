# Integrated SOLVE release observations, 28 October 2026

Accepted corrected R11 native **e7b2975992bcffd2fd992ca11ab93f033231efe11cb679269c6ccbe906aa4bd4**, source **5f89ef4eb3c40b63fb53afb216cbb53cecb468c352a06a12d4dca53200326c35**, is compared with accepted corrected R10 native **96d342cf5fd2349f01ba38ec33a9aeec5669baada3de36e66ed75588d5447a04**, source **82feae2c9d065fd338e9439693f9c5b2bd1663b683ac2a9e411937b923e38754**. Both arms pass complete identified Demeter Python/Rust/core and binding-clippy gates. Initial failed snapshots are excluded.

Public QQ triangular solve over the identical 18 fixtures of the earlier solve harness, including multiple right-hand sides and consistent overdetermined systems. Independent Fraction A*X=B and rank=n certify complete unique solutions. Both arms positively enter the native path in separate processes.

**1,152 exact certificates agree over 18 cases and 144 fresh processes.** Exact values, domains, aliases and input state are checked after all timers, with runtime SymPy imports blocked. Each timed process performs its real first public operation after construction, followed by nine calls on the same object, without preceding route queries, certificates or artificial cache priming. Separate instrumented route processes run after all timing processes. Three interleaved rounds pin CPU 24 with alternating arm/case order; imports and construction are outside timers. Both complete source maps and native hashes remain unchanged before/after. Agents acknowledged quiet compute before these sequential campaigns; ordinary host load remains recorded and the host is not dedicated.

R10/R11 ratios span **0.593–2.012× first**, **0.376–0.752× repeated**. A ratio above 1 means R11 is faster. Every slower median remains in the table and [explicit regression list](regressions.json); no whole-engine gain or global parity follows.

| Case | First R10 → R11 (µs) | Repeated R10 → R11 (µs) | First/repeated R10÷R11 |
| --- | ---: | ---: | ---: |
| ddm-2x2-small-rhs1 | 133.1 → 224.5 | 55.0 → 141.1 | 0.593× / 0.390× |
| domainmatrix_dense-2x2-small-rhs1 | 140.7 → 224.6 | 56.4 → 141.7 | 0.626× / 0.398× |
| domainmatrix_sparse-2x2-small-rhs1 | 171.2 → 252.4 | 80.3 → 159.4 | 0.678× / 0.504× |
| ddm-2x2-rhs3 | 139.1 → 232.2 | 55.2 → 146.8 | 0.599× / 0.376× |
| domainmatrix_dense-2x2-rhs3 | 144.3 → 229.1 | 57.2 → 141.7 | 0.630× / 0.404× |
| domainmatrix_sparse-2x2-rhs3 | 171.9 → 250.4 | 78.7 → 169.8 | 0.687× / 0.463× |
| ddm-4x4-rhs1 | 152.3 → 242.3 | 69.5 → 156.2 | 0.629× / 0.445× |
| domainmatrix_dense-4x4-rhs1 | 160.7 → 239.9 | 70.0 → 154.1 | 0.670× / 0.455× |
| domainmatrix_sparse-4x4-rhs1 | 181.8 → 278.6 | 103.7 → 181.4 | 0.653× / 0.571× |
| ddm-4x2-tall-rhs3 | 143.1 → 225.9 | 65.5 → 148.8 | 0.634× / 0.441× |
| domainmatrix_dense-4x2-tall-rhs3 | 160.7 → 236.2 | 65.8 → 158.0 | 0.680× / 0.416× |
| domainmatrix_sparse-4x2-tall-rhs3 | 185.2 → 270.0 | 102.1 → 175.9 | 0.686× / 0.581× |
| ddm-8x8-rhs3 | 322.0 → 414.8 | 225.5 → 307.5 | 0.776× / 0.733× |
| domainmatrix_dense-8x8-rhs3 | 785.7 → 412.7 | 223.7 → 305.9 | 1.904× / 0.731× |
| domainmatrix_sparse-8x8-rhs3 | 798.4 → 467.5 | 272.9 → 362.9 | 1.708× / 0.752× |
| ddm-8x4-tall-rhs2 | 205.6 → 281.3 | 108.9 → 190.7 | 0.731× / 0.571× |
| domainmatrix_dense-8x4-tall-rhs2 | 198.6 → 293.0 | 111.9 → 195.9 | 0.678× / 0.571× |
| domainmatrix_sparse-8x4-tall-rhs2 | 673.8 → 334.9 | 142.1 → 229.0 | 2.012× / 0.621× |

The safe raw archive includes scripts/fixtures, every process sample/stderr, identities, commands and route/certificate observations. Wrong, incomplete, held, DEBUG and refused baselines are never counted as gains. These bounded operation comparisons do not establish import, end-to-end program or whole-engine acceleration.
