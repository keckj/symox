# Integration replay of accepted HNF/metadata/diffgeom release, 28 October 2026

Native **e7b2975992bcffd2fd992ca11ab93f033231efe11cb679269c6ccbe906aa4bd4**, sources **5f89ef4eb3c40b63fb53afb216cbb53cecb468c352a06a12d4dca53200326c35** (3,436 files), the [accepted complete R11 gate](../validation-2026-10-28-normalforms-risch-lll.json).

All **1,210 inputs per backend finish exactly once** (2,420 calls), with 16 identities each and no errors/timeouts, under the unchanged 60-second/8-GiB policy. Native workers block runtime SymPy imports; reference workers use SymPy 1.14. The exact replay runner/helper are unchanged from R2/R3. **All 1,210 native output trees are identical to R2**. The separate 21 parameter and 10 E1 controls pass on this release.

Raw counts remain 1,164 exact matches, 36 differences and 10 held on both sides. Exact/domain review retains 20 conditional values with full original Integral fallback, 10 principal-E1 reference corrections, 5 Ci equivalences and 1 affine-nonzero capability gap. Point substitutions are regression controls; the analytical/domain proofs are separate. This broad corpus does not replace component domain/precision/completeness regressions.

[Manifest and command identities](replay-manifest.json), [coverage](coverage-check.json), [R2 delta](r2-delta.json), [exact review](integration-triage-r11.json), and [complete portable evidence](integration-r11-20261028.tar.gz) are retained. Every manifest entry was verified in the safe archive. Earlier R2/R3 archives remain unchanged. No global parity or speedup follows from these counts.
