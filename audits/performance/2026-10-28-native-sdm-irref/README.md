# Sparse QQ public field GJ release comparison, 28 October 2026

Accepted R7 native40cb0e7a/source34632bf7 is compared with accepted R8 native7932ea23/source40be6712. Both releases pass their complete identified Demeter suites. These are selected public sparse GJ observations, not a whole-engine claim.

**1,152 exact certificates agree** across18 public cases (six matrices × SDM.rref, DomainMatrix GJ and DomainMatrix AUTO), with108 timing processes and36 separate untimed route probes. Independent Fraction elimination certifies complete rows, pivots, rank and selected nullspace witnesses; output order, key/scalar identity relationships and unchanged input containers are checked outside all timers. All output digests/aliases/choices agree between releases. SymPy imports are blocked. Tiny, sparse8/24-row, rank-deficient, moderately dense and million-column coordinate cases retain their full shapes without dimension-sized densification.

Three interleaved rounds use fresh processes pinned to CPU24 for each case/arm. The actual first rref follows import/construction without any earlier rref, choice query or certificate; nine repeated calls follow. Every certificate/conversion is after all timings in that process. The separate probes verify native GJ on BOTH first and repeated R8 calls. Two AUTO controls choose CD and never enter the field GJ kernel; their near-unchanged times are controls. Process medians are retained individually, with median-of-process-medians below.

Selected24-row sparse/rank-deficient public calls gain **2.26–2.70× on first calls and2.67–3.34× repeated**. Million-column coordinate cases gain1.64–1.84× first and2.15–2.89× repeated, without densification. The tiny2×3 first SDM/GJ calls regress from61.2/76.2 to78.5/96.7µs (**+17.3/+20.5µs, +28.3/+26.9%**), although their repeated calls gain1.35–1.53×. This first-call regression remains explicit; no acceleration is claimed for those tiny first calls or AUTO/CD controls.

| Public case | First R7 → R8 (µs) | First ratio | Repeated R7 → R8 (µs) | Repeated ratio |
| --- | ---: | ---: | ---: | ---: |
| sdm-tiny-2x3 | 61.2 → 78.5 | 0.78× | 26.8 → 17.5 | 1.53× |
| domainmatrix_gj-tiny-2x3 | 76.2 → 96.7 | 0.79× | 39.2 → 29.0 | 1.35× |
| domainmatrix_auto-tiny-2x3 | 16596.2 → 16522.6 | 1.00× | 262.6 → 264.9 | 0.99× |
| sdm-sparse-8x12 | 156.0 → 137.3 | 1.14× | 115.3 → 57.6 | 2.00× |
| domainmatrix_gj-sparse-8x12 | 176.3 → 156.6 | 1.13× | 133.7 → 76.7 | 1.74× |
| domainmatrix_auto-sparse-8x12 | 265.7 → 220.9 | 1.20× | 188.5 → 127.7 | 1.48× |
| sdm-sparse-24x40 | 1140.3 → 422.5 | 2.70× | 1076.1 → 322.7 | 3.34× |
| domainmatrix_gj-sparse-24x40 | 1165.8 → 483.1 | 2.41× | 1097.6 → 364.4 | 3.01× |
| domainmatrix_auto-sparse-24x40 | 1229.0 → 543.7 | 2.26× | 1156.7 → 433.0 | 2.67× |
| sdm-rankdef-24x40 | 986.8 → 376.0 | 2.62× | 941.9 → 295.3 | 3.19× |
| domainmatrix_gj-rankdef-24x40 | 1038.8 → 423.9 | 2.45× | 993.6 → 324.3 | 3.06× |
| domainmatrix_auto-rankdef-24x40 | 1108.6 → 486.4 | 2.28× | 1032.4 → 385.6 | 2.68× |
| sdm-moderate-12x16 | 4813.8 → 1987.1 | 2.42× | 4642.0 → 1885.9 | 2.46× |
| domainmatrix_gj-moderate-12x16 | 4880.1 → 2084.2 | 2.34× | 4783.2 → 1923.8 | 2.49× |
| domainmatrix_auto-moderate-12x16 | 17898.1 → 17700.3 | 1.01× | 1469.5 → 1477.7 | 0.99× |
| sdm-large-coordinates-8x1048577 | 318.7 → 173.3 | 1.84× | 279.3 → 96.8 | 2.89× |
| domainmatrix_gj-large-coordinates-8x1048577 | 362.0 → 197.1 | 1.84× | 310.0 → 115.1 | 2.69× |
| domainmatrix_auto-large-coordinates-8x1048577 | 439.4 → 268.6 | 1.64× | 368.8 → 171.8 | 2.15× |

All complete source maps and native binaries are verified before/after, with both snapshot identities and normal host load in raw records. The 60.5-second campaign runs after all Python/Rust/clippy gates and acknowledgements that all three agents have stopped compute; ordinary services remain active, so this is not a dedicated isolated host. No DEBUG timing, wrong/incomplete result, missing solution or resource refusal counts as a gain.

The safe raw archive includes scripts, exact fixtures, functional non-timed R7 validation, all 144 process results/stderr, every sample and route/identity checks.
