# Integrated LLL release observations, 28 October 2026

Accepted corrected R10 native **96d342cf5fd2349f01ba38ec33a9aeec5669baada3de36e66ed75588d5447a04**, source **82feae2c9d065fd338e9439693f9c5b2bd1663b683ac2a9e411937b923e38754**, is compared with accepted R9; Risch additionally includes accepted R8. All arms pass complete identified Demeter Python/Rust/core and binding-clippy gates. The failed initial R10 snapshot is excluded.

The same 18 exact LLL fixtures as the accepted R9 replay. First-call mathematics remains the established Python path; separate fresh probes confirm native selection only after actual cache state allows it. This measures a metadata-proof repair, not a completed cold LLL port.

**1,152 exact certificates agree over 18 cases and 144 fresh processes.** All exact values, domains, aliases and input state are verified after all timers, with runtime SymPy imports blocked. Timed processes perform a real first public operation after construction, followed by 9 calls on the same object; no preceding route query, certificate or artificial cache priming occurs. Separate instrumented processes run after every timing process. Three interleaved rounds pin CPU 24; arm/case order alternates (three-arm Latin square for Risch). Imports and construction are outside timers. All source maps and native hashes are unchanged before/after. Agents acknowledge quiet compute before the sequential campaigns; ordinary host load is retained and the host is not dedicated.

R9/R10 ratios span **0.977–1.120× first**, **0.964–1.070× repeated**. A ratio above 1 means R10 is faster. Every slower sample remains in the table and [explicit regression list](regressions.json); no whole-engine gain or global parity is inferred.

| Case | First R9 → R10 (µs) | Repeated R9 → R10 (µs) | First/repeated R9÷R10 |
| --- | ---: | ---: | ---: |
| ddm-2x2-delta3_4 | 532.2 → 479.9 | 196.3 → 196.9 | 1.109× / 0.997× |
| domainmatrix_dense-2x2-delta3_4 | 540.2 → 485.2 | 205.6 → 201.6 | 1.113× / 1.020× |
| domainmatrix_sparse-2x2-delta3_4 | 564.3 → 529.0 | 238.7 → 244.5 | 1.067× / 0.976× |
| ddm-2x4-delta3_4 | 519.2 → 463.5 | 206.4 → 207.7 | 1.120× / 0.994× |
| domainmatrix_dense-2x4-delta3_4 | 523.3 → 477.2 | 220.2 → 213.5 | 1.096× / 1.032× |
| domainmatrix_sparse-2x4-delta3_4 | 561.2 → 514.3 | 257.2 → 240.4 | 1.091× / 1.070× |
| ddm-4x4-delta3_4 | 1540.7 → 1464.5 | 277.5 → 270.3 | 1.052× / 1.027× |
| domainmatrix_dense-4x4-delta3_4 | 1537.8 → 1470.1 | 280.4 → 271.4 | 1.046× / 1.033× |
| domainmatrix_sparse-4x4-delta3_4 | 1576.1 → 1532.5 | 325.2 → 328.4 | 1.028× / 0.990× |
| ddm-4x6-delta3_4 | 1886.4 → 1815.6 | 313.5 → 317.2 | 1.039× / 0.989× |
| domainmatrix_dense-4x6-delta3_4 | 1889.6 → 1829.9 | 318.7 → 314.4 | 1.033× / 1.014× |
| domainmatrix_sparse-4x6-delta3_4 | 1979.5 → 1906.1 | 367.0 → 365.2 | 1.039× / 1.005× |
| ddm-6x6-delta3_4 | 3854.9 → 3712.5 | 641.8 → 631.2 | 1.038× / 1.017× |
| domainmatrix_dense-6x6-delta3_4 | 3866.2 → 3655.4 | 626.8 → 629.1 | 1.058× / 0.996× |
| domainmatrix_sparse-6x6-delta3_4 | 3871.9 → 3963.8 | 693.0 → 693.2 | 0.977× / 1.000× |
| ddm-6x8-delta3_4 | 3858.4 → 3847.2 | 771.0 → 765.2 | 1.003× / 1.008× |
| domainmatrix_dense-6x8-delta3_4 | 3822.6 → 3745.5 | 782.5 → 812.0 | 1.021× / 0.964× |
| domainmatrix_sparse-6x8-delta3_4 | 3963.0 → 3885.2 | 834.4 → 836.7 | 1.020× / 0.997× |

The safe raw archive contains all scripts/fixtures, process samples/stderr, identities, commands and selection/certificate observations. Wrong, incomplete, held, DEBUG and refused baselines are never counted as gains. This bounded operation comparison does not establish import, end-to-end program or whole-integrator acceleration.
