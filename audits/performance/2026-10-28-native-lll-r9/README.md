# Integrated LLL release observations, 28 October 2026

Accepted R8 native 7932ea23/source 40be6712 is paired with accepted R9 native c3f413bd/source d81c42c9. Both identified packages pass complete Demeter Python/Rust/clippy gates. These are selected public operations; no global parity or whole-engine acceleration is claimed.

**1,152 exact certificates agree** over **18 cases and 144 fresh processes**. Independent Fraction checks, output digests/domains, mutation/alias conventions and input/extension state are verified outside every timer; runtime SymPy imports are blocked. Unsupported/error controls are preserved, never counted as a gain. See exact raw records for the full case protocol.

Three interleaved rounds pin each process to CPU 24. A real first call follows actual package import and construction without an earlier operation, selection query or certificate. Nine repeated calls follow on the same object. All result conversion/certificates occur after all timers. Independent fresh instrumented processes follow all timing processes; their traces are not timing samples. Both entire source maps and binaries match before/after. Source identities, per-process samples and ordinary host load remain recorded; the host is not dedicated. The three campaigns run sequentially after all agents acknowledge compute stopped and complete gates have ended.

The first-call ratios span 0.977–1.075: several cases improve modestly, two 6-row DDM controls are about 2.3% slower. Every repeated-call median regresses 2.8–15.3% (roughly 10–43µs) versus R8, so no warm improvement is claimed. Every separate probe observes no native LLL on the first call and native LLL on the second; genuinely cold mathematics remains Python. The added read-only proof avoids the heavy cold rejection but has a measurable warm cost, retained as open work.

| Case | First R8 → R9 (µs) | First ratio | Repeated R8 → R9 (µs) | Repeated ratio |
| --- | ---: | ---: | ---: | ---: |
| ddm-2x2-delta3_4 | 556.1 → 539.4 | 1.031× | 183.0 → 211.0 | 0.867× |
| domainmatrix_dense-2x2-delta3_4 | 562.7 → 531.0 | 1.060× | 183.9 → 210.1 | 0.875× |
| domainmatrix_sparse-2x2-delta3_4 | 596.2 → 569.9 | 1.046× | 218.9 → 237.3 | 0.923× |
| ddm-2x4-delta3_4 | 528.6 → 501.9 | 1.053× | 187.9 → 205.9 | 0.913× |
| domainmatrix_dense-2x4-delta3_4 | 537.7 → 516.7 | 1.041× | 191.7 → 220.1 | 0.871× |
| domainmatrix_sparse-2x4-delta3_4 | 570.7 → 532.3 | 1.072× | 224.7 → 242.2 | 0.928× |
| ddm-4x4-delta3_4 | 1615.2 → 1519.4 | 1.063× | 253.9 → 276.3 | 0.919× |
| domainmatrix_dense-4x4-delta3_4 | 1613.8 → 1531.4 | 1.054× | 258.5 → 280.2 | 0.923× |
| domainmatrix_sparse-4x4-delta3_4 | 1639.3 → 1589.1 | 1.032× | 301.5 → 330.1 | 0.913× |
| ddm-4x6-delta3_4 | 1907.6 → 1867.3 | 1.022× | 303.7 → 313.2 | 0.970× |
| domainmatrix_dense-4x6-delta3_4 | 1912.0 → 1874.7 | 1.020× | 295.6 → 320.0 | 0.924× |
| domainmatrix_sparse-4x6-delta3_4 | 1988.3 → 1983.4 | 1.002× | 345.4 → 367.5 | 0.940× |
| ddm-6x6-delta3_4 | 3831.6 → 3920.7 | 0.977× | 608.7 → 638.6 | 0.953× |
| domainmatrix_dense-6x6-delta3_4 | 4131.3 → 3842.0 | 1.075× | 621.6 → 657.2 | 0.946× |
| domainmatrix_sparse-6x6-delta3_4 | 4043.2 → 3818.3 | 1.059× | 672.8 → 693.4 | 0.970× |
| ddm-6x8-delta3_4 | 3875.2 → 3952.9 | 0.980× | 747.9 → 772.0 | 0.969× |
| domainmatrix_dense-6x8-delta3_4 | 3934.4 → 3839.7 | 1.025× | 754.4 → 797.6 | 0.946× |
| domainmatrix_sparse-6x8-delta3_4 | 4021.6 → 3977.2 | 1.011× | 815.3 → 838.1 | 0.973× |

Safe raw archive contains scripts/fixtures, all process results/stderr, commands, identity and route observations. No DEBUG timing, wrong/held/incomplete baseline, missing solution or resource refusal is a speedup. Ratios use median-of-process-medians; construction/imports and certificates are outside timers, so import or complete-program gains are not claimed.
