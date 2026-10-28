# Integration replay of accepted solve/Risch/LLL release, 28 October 2026

Native **96d342cf5fd2349f01ba38ec33a9aeec5669baada3de36e66ed75588d5447a04**, sources **82feae2c9d065fd338e9439693f9c5b2bd1663b683ac2a9e411937b923e38754** (3,299 files), the [accepted complete R10 gate](../validation-2026-10-28-solve-risch-lll.json).

All **1,210 inputs per backend finish exactly once** (2,420 calls), with 16 identities each and no errors/timeouts, under the unchanged 60-second/8-GiB policy. Native workers block runtime SymPy imports; reference workers use SymPy 1.14. The exact replay runner/helper are unchanged from R2/R3. **All 1,210 native output trees are identical to R2**. The separate 21 parameter and 10 E1 controls pass on this release.

Raw counts remain 1,164 exact matches, 36 differences and 10 held on both sides. Exact/domain review retains 20 conditional values with full original Integral fallback, 10 principal-E1 reference corrections, 5 Ci equivalences and 1 affine-nonzero capability gap. Point substitutions are regression controls; the analytical/domain proofs are separate. This broad corpus does not replace component domain/precision/completeness regressions.

[Manifest and command identities](replay-manifest.json), [coverage](coverage-check.json), [R2 delta](r2-delta.json), [exact review](integration-triage-r10.json), and [complete portable evidence](integration-r10-96d342cf-20261028.tar.gz) are retained. Every manifest entry was verified in the safe archive. Earlier R2/R3 archives remain unchanged. No global parity or speedup follows from these counts.
