# Integration replay of accepted LU/Risch/LLL release, 28 October 2026

Native **c3f413bd8a3bbc4b04aa9d7ebc28c69d9aba2170f0b702b09c2a6136e25f0a0c**, sources **d81c42c90dd3bcc21e0a3bad4acf6a8fffcba60194f564e9764da8e1443cd54a** (3,225 files), the [accepted complete R9 gate](../validation-2026-10-28-lu-risch-lll.json).

All **1,210 inputs per backend finish exactly once** (2,420 calls), with 16 identities each and no errors/timeouts, under the unchanged 60-second/8-GiB policy. Native workers block runtime SymPy imports; reference workers use SymPy1.14. The exact replay runner/helper are unchanged from R2/R3. **All1,210 native output trees are identical to R2**. The separate21 parameter and10 E1 controls pass on this release.

Raw counts remain1,164 exact matches,36 differences and10 held on both sides. Exact/domain review retains20 conditional values with full original Integral fallback,10 principal-E1 reference corrections,5 Ci equivalences and1 affine-nonzero capability gap. Point substitutions are regression controls; the analytical/domain proofs are separate. This broad corpus does not replace component domain/precision/completeness regressions.

[Manifest and command identities](replay-manifest.json), [coverage](coverage-check.json), [R2 delta](r2-delta.json), [exact review](integration-triage-r9.json), and [complete portable evidence](integration-r9-c3f413bd-20261028.tar.gz) are retained. Every manifest entry was verified in the safe archive. Earlier R2/R3 archives remain unchanged. No global parity or speedup follows from these counts.
