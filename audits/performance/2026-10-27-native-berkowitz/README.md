# Public dense Berkowitz release measurements, 27 October 2026

Before: immutable `cb5c0c98` / `db859db2`. After: accepted `daef4243` / `8f254811`. [Full identities, protocol and checksums](manifest.json) and [raw medians](summary.json) are preserved. This benchmark calls **DDM.charpoly()**, including its public adapter, on three distinct matrices in each family. All input domains, matrices and expected coefficients are constructed before timers.

Three fresh processes per release alternate CB5/R2, R2/CB5, CB5/R2 on CPU24 of Demeter during a quiet window. Each process takes five rounds of three batches (nine matrix calls per timed round); the table gives the median of the three process medians, per matrix. These are repeated calls after warmup, with the exact small sample sizes recorded. Imports, construction, fixture certification and spies are outside timers. No cold-start or full matrix-engine speed claim is inferred.

| Domain and size | Before, microseconds/matrix | After, microseconds/matrix | Observed ratio |
| --- | ---: | ---: | ---: |
| ZZ, 4x4 | 89.12 | 17.21 | 5.18x |
| ZZ, 8x8 | 469.06 | 40.88 | 11.47x |
| ZZ, 12x12 | 1356.57 | 127.21 | 10.66x |
| QQ, 4x4 | 100.51 | 17.95 | 5.60x |
| QQ, 8x8 | 714.96 | 52.16 | 13.71x |

All 15 matrices have independent exact SymPy 1.14 characteristic-polynomial and Cayley-Hamilton certificates. The six native workers block runtime SymPy imports, check the complete coefficients before and after timing (**180 complete coefficient-vector comparisons**) and preserve all digests. The final release confirms 45 native selections outside timers; original row and element identities remain unchanged. The fixtures use both signs and rational denominators; no incorrect or uncomputed old result supplies a baseline.

[Fixtures](fixtures.json), [fixture reproducer](prepare-fixtures.py), [worker](worker.py), [execution runner](run.py) and all six raw process reports are included. The [full accepted release gate](../../sympy-1.14.0/validation-2026-10-27-lrt-limits.json) separately tests domains, callbacks, aliases and error behavior. Observed gains of 5.2x–13.7x apply to these five dense ZZ/QQ families; sparse, other domains and global mathematical parity are outside that scope.
