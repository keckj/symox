# Integrated LU release observations, 28 October 2026

Accepted R8 native 7932ea23/source 40be6712 is paired with accepted R9 native c3f413bd/source d81c42c9. Both identified packages pass complete Demeter Python/Rust/clippy gates. These are selected public operations; no global parity or whole-engine acceleration is claimed.

**1,152 exact certificates agree** over **18 cases and 144 fresh processes**. Independent Fraction checks, output digests/domains, mutation/alias conventions and input/extension state are verified outside every timer; runtime SymPy imports are blocked. Unsupported/error controls are preserved, never counted as a gain. See exact raw records for the full case protocol.

Three interleaved rounds pin each process to CPU 24. A real first call follows actual package import and construction without an earlier operation, selection query or certificate. Nine repeated calls follow on the same object. All result conversion/certificates occur after all timers. Independent fresh instrumented processes follow all timing processes; their traces are not timing samples. Both entire source maps and binaries match before/after. Source identities, per-process samples and ordinary host load remain recorded; the host is not dedicated. The three campaigns run sequentially after all agents acknowledge compute stopped and complete gates have ended.

Selected 8-row public LU calls gain 4.03–6.42× first and 5.07–9.63× repeated; selected 4×4 calls gain 1.42–1.82× first and 1.96–3.22× repeated. Tiny 2×2 first calls retain 9–17µs extra cost (about 15–26%); sparse 2×2 repeated also costs about 6.7%. These costs remain open. This measures LU decomposition, not the still-Python triangular solver.

| Case | First R8 → R9 (µs) | First ratio | Repeated R8 → R9 (µs) | Repeated ratio |
| --- | ---: | ---: | ---: | ---: |
| ddm-2x2-pivot | 55.6 → 70.1 | 0.793× | 24.6 → 23.8 | 1.033× |
| domainmatrix_dense-2x2-pivot | 63.4 → 72.8 | 0.872× | 27.3 → 25.8 | 1.059× |
| domainmatrix_sparse-2x2-pivot | 93.7 → 110.5 | 0.848× | 47.3 → 50.5 | 0.936× |
| ddm-2x4-wide | 67.5 → 68.7 | 0.983× | 34.0 → 24.9 | 1.367× |
| domainmatrix_dense-2x4-wide | 73.1 → 71.3 | 1.025× | 37.5 → 26.6 | 1.407× |
| domainmatrix_sparse-2x4-wide | 104.1 → 111.4 | 0.935× | 61.2 → 52.5 | 1.165× |
| ddm-4x4-pivot | 151.0 → 82.8 | 1.823× | 111.0 → 34.4 | 3.223× |
| domainmatrix_dense-4x4-pivot | 153.9 → 86.9 | 1.770× | 118.3 → 36.9 | 3.205× |
| domainmatrix_sparse-4x4-pivot | 199.0 → 139.8 | 1.423× | 146.3 → 74.8 | 1.958× |
| ddm-4x2-tall | 87.3 → 79.9 | 1.093× | 51.7 → 30.1 | 1.717× |
| domainmatrix_dense-4x2-tall | 91.2 → 82.9 | 1.100× | 54.7 → 31.7 | 1.726× |
| domainmatrix_sparse-4x2-tall | 136.2 → 120.9 | 1.126× | 84.8 → 62.4 | 1.359× |
| ddm-8x8-full | 870.0 → 162.4 | 5.355× | 823.0 → 102.6 | 8.019× |
| domainmatrix_dense-8x8-full | 881.4 → 168.5 | 5.231× | 819.7 → 109.5 | 7.486× |
| domainmatrix_sparse-8x8-full | 958.9 → 237.9 | 4.031× | 889.9 → 175.5 | 5.070× |
| ddm-8x8-rank6 | 859.4 → 142.5 | 6.033× | 808.8 → 83.9 | 9.634× |
| domainmatrix_dense-8x8-rank6 | 874.8 → 136.4 | 6.416× | 800.9 → 84.3 | 9.504× |
| domainmatrix_sparse-8x8-rank6 | 910.0 → 205.9 | 4.420× | 855.0 → 149.1 | 5.733× |

Safe raw archive contains scripts/fixtures, all process results/stderr, commands, identity and route observations. No DEBUG timing, wrong/held/incomplete baseline, missing solution or resource refusal is a speedup. Ratios use median-of-process-medians; construction/imports and certificates are outside timers, so import or complete-program gains are not claimed.
