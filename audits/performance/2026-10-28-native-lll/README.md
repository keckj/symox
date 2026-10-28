# Public LLL release comparison, 28 October 2026

Accepted R6r1 native cc6d0210/source 80c98ada is compared with accepted R7 native 40cb0e7a/source 34632bf7. Both releases pass their identified full Demeter suites. This is a bounded public LLL comparison, not a whole-engine performance claim.

**1,152 exact certificates agree**, across 18 cases from 6 moderate matrices (2/4/6 rows, square/rectangular) and 3 public entries: DDM, dense DomainMatrix and sparse DomainMatrix. Original unmodified reference and exact-floor reference produce identical B/T on every fixture; known lossy-rounding/assertion cases are excluded from speed comparisons. Independent Fraction checks certify T*A=B, determinant +/-1, size reduction and Lovasz conditions; input matrices remain unchanged. SymPy imports are blocked in both measured arms.

Three interleaved rounds use 108 fresh timing processes pinned to CPU 24. Every first lll_transform call follows actual import/construction with no earlier LLL, QQ cache priming or result certificate; nine repeated calls follow on the same object. All output conversion/certificates occur after every timer in that process. Medians of process medians are reported below. Separate 36 instrumented fresh processes run after all timings: every R7 first call makes no native LLL call, and every second call selects native LLL. The tracing process is never a timing sample.

**Repeated 4/6-row calls gain 4.20–5.60×; tiny 2-row calls gain 1.04–1.18×. Every first-call median regresses**, by about 1.9–21.1% (about 66–287 microseconds): the old complete Python algorithm still runs with added admission checks. This cost is retained as open architecture/performance debt; no first-call acceleration is claimed.

| Public case | First call R6r1 → R7 (µs) | First cost | Repeated R6r1 → R7 (µs) | Repeated gain |
| --- | ---: | ---: | ---: | ---: |
| ddm-2x2-delta3_4 | 460.7 → 552.7 | +20.0% | 210.2 → 178.5 | 1.18× |
| domainmatrix_dense-2x2-delta3_4 | 476.5 → 576.8 | +21.1% | 219.1 → 190.6 | 1.15× |
| domainmatrix_sparse-2x2-delta3_4 | 481.2 → 581.2 | +20.8% | 244.6 → 217.5 | 1.12× |
| ddm-2x4-delta3_4 | 443.8 → 526.7 | +18.7% | 201.6 → 193.6 | 1.04× |
| domainmatrix_dense-2x4-delta3_4 | 446.5 → 536.1 | +20.1% | 205.6 → 192.4 | 1.07× |
| domainmatrix_sparse-2x4-delta3_4 | 498.4 → 564.7 | +13.3% | 248.1 → 231.6 | 1.07× |
| ddm-4x4-delta3_4 | 1438.6 → 1516.9 | +5.4% | 1178.0 → 250.5 | 4.70× |
| domainmatrix_dense-4x4-delta3_4 | 1434.7 → 1595.1 | +11.2% | 1169.2 → 261.9 | 4.46× |
| domainmatrix_sparse-4x4-delta3_4 | 1454.0 → 1641.0 | +12.9% | 1256.1 → 298.7 | 4.20× |
| ddm-4x6-delta3_4 | 1768.0 → 1874.5 | +6.0% | 1485.1 → 291.8 | 5.09× |
| domainmatrix_dense-4x6-delta3_4 | 1774.2 → 1892.3 | +6.7% | 1499.5 → 294.3 | 5.09× |
| domainmatrix_sparse-4x6-delta3_4 | 1841.1 → 1993.3 | +8.3% | 1562.8 → 352.0 | 4.44× |
| ddm-6x6-delta3_4 | 3684.5 → 3776.8 | +2.5% | 3380.8 → 607.2 | 5.57× |
| domainmatrix_dense-6x6-delta3_4 | 3702.9 → 3839.3 | +3.7% | 3409.7 → 608.7 | 5.60× |
| domainmatrix_sparse-6x6-delta3_4 | 3784.3 → 3930.8 | +3.9% | 3473.9 → 668.8 | 5.19× |
| ddm-6x8-delta3_4 | 3801.2 → 3874.0 | +1.9% | 3455.4 → 755.5 | 4.57× |
| domainmatrix_dense-6x8-delta3_4 | 3703.7 → 3990.5 | +7.7% | 3445.8 → 756.8 | 4.55× |
| domainmatrix_sparse-6x8-delta3_4 | 3770.0 → 3900.4 | +3.5% | 3497.3 → 817.3 | 4.28× |

All source manifests and binary hashes are checked before/after for both immutable packages, and all output digests agree between arms. The 62-second campaign starts after Python/Rust/clippy gates and private agent activity end. Ordinary Demeter host load is retained in raw identity/summary (1-minute load 4.20 before, 2.22 after); the CPU is pinned, not an isolated machine. No new source edit, incomplete baseline, DEBUG timing or wrong result is counted as a gain.

The safe raw archive includes scripts, fixtures, all 144 process records/stderr, commands, both source identities and all first/repeated samples. summary.json exposes every comparison and instrumented route outcome.
