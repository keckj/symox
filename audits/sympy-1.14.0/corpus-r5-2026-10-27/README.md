# Integration replay of accepted PolyMatrix/log_to_real release, 27 October 2026

Native **532386bda46046d3d21c99f49c1b863fac62bdee5c357d03cf5123bde6913d4b**, sources **2955abf8b02c1e520ddb921d17e56614b8f4c5d391b1fb87730118881665dbe0** (2,879 files), the [accepted complete R5r1 gate](../validation-2026-10-27-polymatrix-real-log.json).

All **1,210 inputs per backend finish exactly once** (2,420 calls), with 16 identities each and no errors/timeouts, under the unchanged 60-second/8-GiB policy. Native workers block runtime SymPy imports; reference workers use SymPy 1.14. The exact replay runner/helper are unchanged from R2/R3. **All 1,210 native output trees are identical to R2**. The separate 21 parameter and 10 E1 controls pass on this release.

Raw counts remain 1,164 exact matches, 36 differences and 10 held on both sides. Exact/domain review retains 20 conditional values with full original Integral fallback, 10 principal-E1 reference corrections, 5 Ci equivalences and 1 affine-nonzero capability gap. Point substitutions are regression controls; the analytical/domain proofs are separate. This broad corpus does not replace component domain/precision/completeness regressions.

[Manifest and command identities](replay-manifest.json), [coverage](coverage-check.json), [R2 delta](r2-delta.json), [exact review](integration-triage-r5.json), and [complete portable evidence](integration-r5-532386bd-20261027.tar.gz) are retained. Every manifest entry was verified in the safe archive. Earlier R2/R3 archives remain unchanged. No global parity or speedup follows from these counts.
