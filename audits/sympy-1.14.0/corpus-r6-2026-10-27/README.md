# Integration replay of accepted rational/AUTO/sparse release, 27 October 2026

Native **cc6d021087f8c3d8404011ecf59947d071dfa7ef9bb865c2ae7bd44431cc09ab**, sources **80c98adae719e1cdfc2e7621a234d0ade781cdb3bc26c4285360263b0ca719ae** (3,018 files), the [accepted complete R6r1 gate](../validation-2026-10-27-matrix-rational.json).

All **1,210 inputs per backend finish exactly once** (2,420 calls), with 16 identities each and no errors/timeouts, under the unchanged 60-second/8-GiB policy. Native workers block runtime SymPy imports; reference workers use SymPy 1.14. The exact replay runner/helper are unchanged from R2/R3. **All 1,210 native output trees are identical to R2**. The separate 21 parameter and 10 E1 controls pass on this release.

Raw counts remain 1,164 exact matches, 36 differences and 10 held on both sides. Exact/domain review retains 20 conditional values with full original Integral fallback, 10 principal-E1 reference corrections, 5 Ci equivalences and 1 affine-nonzero capability gap. Point substitutions are regression controls; the analytical/domain proofs are separate. This broad corpus does not replace component domain/precision/completeness regressions.

[Manifest and command identities](replay-manifest.json), [coverage](coverage-check.json), [R2 delta](r2-delta.json), [exact review](integration-triage-r6.json), and [complete portable evidence](integration-r6-cc6d0210-20261027.tar.gz) are retained. Every manifest entry was verified in the safe archive. Earlier R2/R3 archives remain unchanged. No global parity or speedup follows from these counts.
