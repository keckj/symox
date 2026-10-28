# Paired LRT release measurements, 27 October 2026

Before: release `cb5c0c98`, source `db859db2`. After: accepted release `daef4243`, source `8f254811`. Full identities are in [manifest.json](manifest.json). Both are immutable packages validated independently; [the final full-suite gate](../../sympy-1.14.0/validation-2026-10-27-lrt-limits.json) is separate from performance.

Three fresh processes per release alternate CB5/R2, R2/CB5, CB5/R2 on CPU24 of Demeter during a quiet window. Inputs are constructed before timing; each process takes five rounds of 50 batches. The table gives the median of the three process medians, per input pair. The eight-input batch contains distinct arguments. Imports, construction, selection spies and certification are outside timers. The other single-pair families measure repeated calls after warmup; there is no cold-start or whole-engine claim.

| LRT family | Before, ms/pair | After, ms/pair | Observed ratio |
| --- | ---: | ---: | ---: |
| Quadratic | 1.096 | 0.358 | 3.06x |
| Quartic | 1.232 | 0.441 | 2.79x |
| Nonconstant numerator | 1.314 | 0.487 | 2.70x |
| Multiple groups | 2.542 | 0.788 | 3.23x |
| distinct arguments | 1.667 | 0.739 | 2.26x |

All workers block runtime SymPy imports. Every exact polynomial coefficient, domain and pair representative is compared before and after timing with independently prepared SymPy 1.14 fixtures: **144 full pair checks**, identical digests across all six workers. **36 native selections** are verified outside timers. Twelve additional exact rational trace certificates prove that the derivative of the represented logarithmic sum is f/g; these run locally after measurements. Matching representatives retain the established branch conventions. No incorrect-old result is a performance baseline.

[Six raw reports and summary](summary.json), [fixtures](fixtures.json), [native worker](worker.py), [execution runner](run.py), [exact trace certificates](exact-trace-certificates.json) and [certificate reproducer](certify-fixtures.py) are preserved unchanged. [The small portable archive](lrt-cb5-r2-release-20261027.tar.gz) and checksums support replay. These five families gain 2.26x to 3.23x; other domains, complete rational integration and global mathematical parity are outside that claim.
