# Rational selector release controls and remaining overhead, 27 October 2026

Two immutable campaigns compare accepted R2 and R3b, then R3b with only the selector optimization overlaid. The optimized package uses the same **8a8a** native binary as R3b and a separately identified adapter (753ad09a); it is **not the later R5 release**. Source identities, adapter source, runners, fixtures, raw outputs and portable archives are included in each directory.

The [initial campaign](initial/RESULTS.txt) uses twelve fresh processes: three per backend/plain-or-vector variant, alternating releases on CPU24 during a confirmed quiet window. The [confirmation](optimized/RESULTS.txt) uses eighteen fresh processes, adding the guard-only overlay. Each process takes five rounds of twenty batches, with inputs preconstructed and exact checks/spies outside timers. Native workers block runtime SymPy imports. The campaigns perform **684 and 1,026 exact checks**, with nineteen independent LRT trace/Atan Euclid-derivative certificates and identical digests, including plain versus vector. Native selection and the complete R2/vector fallback are explicitly recorded.

Medians below are per pair, median of three process medians in the confirmation window. Ratio means R2 divided by optimized latency.

| Variant | Family | R2, us | R3b, us | Optimized overlay, us | R2 / optimized |
| --- | --- | ---: | ---: | ---: | ---: |
| plain | lrt/quadratic | 367.16 | 377.30 | 356.49 | 1.030x |
| plain | lrt/quartic | 442.97 | 453.85 | 428.72 | 1.033x |
| plain | lrt/nonconstant_numerator | 492.55 | 503.05 | 479.36 | 1.028x |
| plain | lrt/multiplicity_groups | 791.10 | 809.40 | 790.23 | 1.001x |
| plain | lrt/distinct_arguments | 756.41 | 766.98 | 745.38 | 1.015x |
| plain | atan/linear_branch | 92.35 | 92.14 | 69.19 | 1.335x |
| plain | atan/quadratic | 104.42 | 103.11 | 80.22 | 1.302x |
| plain | atan/quartic | 137.48 | 137.75 | 113.71 | 1.209x |
| plain | atan/distinct_arguments | 124.13 | 124.63 | 99.85 | 1.243x |
| vector | lrt/quadratic | 1172.26 | 476.71 | 409.80 | 2.861x |
| vector | lrt/quartic | 1307.05 | 558.23 | 492.70 | 2.653x |
| vector | lrt/nonconstant_numerator | 1395.99 | 608.17 | 539.66 | 2.587x |
| vector | lrt/multiplicity_groups | 2598.94 | 907.78 | 839.16 | 3.097x |
| vector | lrt/distinct_arguments | 1733.85 | 871.22 | 802.80 | 2.160x |
| vector | atan/linear_branch | 89.19 | 157.32 | 100.39 | 0.888x |
| vector | atan/quadratic | 119.35 | 167.06 | 109.77 | 1.087x |
| vector | atan/quartic | 241.50 | 198.27 | 144.00 | 1.677x |
| vector | atan/distinct_arguments | 224.47 | 181.30 | 132.19 | 1.698x |

The first campaign exposed real small-case regressions: vector Atan linear increased from about86 to152us and quadratic from115 to167us. These results remain unchanged in the initial raw archive. The optimization improves all eighteen R3b family/variant medians in the confirmation. It distinguishes genuine Ore layers from vector forwarding, shares metadata inspections within a call and retains all callback checks and native selections, without a mutable global cache.

One measured regression remains against R2: **vector linear Atan89.187 ->100.389us**, about **+11.2us (+12.6%)**, versus157.324us for R3b in that same window. Vector LRT gains2.16x–3.10x against R2; the other vector Atan cases gain1.09x–1.70x. Attributing the remaining overhead to authentication/binding is an inference from selected routes, not an isolated component profile. No case was excluded or redirected to Python to manufacture a gain.

[Initial manifest](initial/manifest.json), [confirmation manifest](optimized/manifest.json), both complete archives and all member hashes were verified before copying. These are observations on these inputs, with construction/imports outside timers; they establish neither whole-engine speed nor mathematical completeness. Correctness validation and the eventual integrated source identity are separate gates.
