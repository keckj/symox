# Homogeneous symbolic-power Gruntz leading terms — frozen delivery

Private source `/home/keckj/.cache/symox-validation/native-limit-symbolic-20261026`.
No source in main, the parent's exploratory copy, or the frozen QQ/DDM stages
was changed. This copy started from the parent's native-limit-20261026 input;
the final patch is against the real immutable02fa baseline, excluding every
exploratory GUARD/TRACE edit. Python is validation-only immutable02fa overlay;
no runtime Python implementation changed.

## Provenance and integration

- Every one of the328 original native source-input.json entries matches
  native-derivative-20261026/after before this patch. Only the two owned limit
  files differ in the final input cohort.
- Native baseline cohort manifest: `0b5e227c9e9f9d477bb1f94103b9b378ad9d0e4b34e03bd03cd8cc210fff5d1f`.
- Native final cohort manifest: `16b9fbbe42a035c32d2619b823b62a4768a423e82088ee62c8fe096787cff270`.
- All784 validation Python files match immutable02fa (binary/cache excluded):
  manifest `6525c6c0bb1256a4974ee5ac477916e2d1c37d806f84ddaf1a0ba8dee6b4a164`.
- New test hash: `65563e5272ab5b4765105e2875b58be0be4c2e0aefc2de02bf96bcf3f72e076d`.
- Patch `symbolic-power-limit.patch` SHA256:
  `12086902e8c7e18dd8fc801e43aa1a7d581692d154911b76af2c6e75745c92e6`.
- apply--check succeeds on immutable02fa. Per-file hashes in file-hashes.json.
- Final private DEBUG extension:
  `58b25015c4ba16905c361fe8ed536ff72d500ec63cfe9329abb08b2cb076ed11`.
- Demeter private tree `/home/keckj/symox-ci/wt/native-limit-symbolic-20261026`;
  target `/home/keckj/symox-ci/targets/native-limit-symbolic`, Python3.10 venv.
- Final native manifest and new test hash were independently reproduced there.
- Rust source unchanged since build-v2; only Python tests were extended later.
- No local Cargo, maturin, shared installation, commit or performance claim.

Files: limits/gruntz.rs, limits/tests.rs, and new
`tests/test_gruntz_homogeneous_symbolic_powers.py`.

**Integration also needs the owner to remove the existing strict xfail decorator
from test_limwrong_parity.py::test_issue_22334_value.** This file was deliberately
not edited because ownership covered newly created Python tests only. Its body
passes with --runxfail; leaving the decorator in the integrated full suite will
produce the expected XPASS(strict).

## General algorithm and domain proof

The existing symbolic MRV transformation replaces each symbolic w**s by
P*w**(s-s0), P fresh and positive. p_degree already proves exact homogeneity.
The original calculation on g(P,w) remains first, and every old success returns
immediately with the exact old coefficient/exponent grouping.

Previously, failure to eliminate P in an otherwise successful series result
returned None before the error fallback. The new code records both that None
and the original Err, then can try g(1,w). A rational homogeneous degree m
certifies the identity g(P,w)=P**m*g(1,w). The final exponent is e+m*s0. If the
new proof is unavailable, the exact previous None/Err is returned, retaining
the established leadterm/fallback paths. No expression-specific pattern or
residual expansion was added; p_degree and the rest of Gruntz are unchanged.

A finite final coefficient alone is insufficient. Series coefficient arithmetic
can cancel (k+1)/(k+1), hiding the source denominator's zero at k=-1. Before
normalization's series answer is accepted, each reciprocal/nonintegral-power
base gets a finite, nonzero leading-coefficient certificate. Nested powers are
checked before any enclosing cancellation. The final coefficient must also be
finite, nonzero and independent of w/P. Unknown parameter poles retain the old
result rather than receiving an unconditional finite continuation.

For noninteger powers, the positive factor P itself preserves principal
arguments. The series answer additionally needs a safe branch: positive leading
coefficient, certified nonreal leading coefficient, or a base proven real/with
nonnegative imaginary part near positive w. On a negative real leading base,
the latter preserves principal argument +pi; an unproved/lower-side approach
cannot silently use the upper-side coefficient. All tests use general domain
properties, not named values, symbols or expressions.

## Validation on Demeter

- Final **43 new Python tests pass** (`native-limit-symbolic-new-final.log`).
- **8 Rust limit tests pass** (`native-limit-symbolic-rust-v2.log`), including
  four new tests: shifted symbolic family, parameter poles even after leading
  coefficient cancellation, reciprocal coefficient proof, and branch sides.
- Existing targeted limits/sequence/Heaviside/Gamma tests: **253 unique cases
  proved passing**, including the old22334 body run with --runxfail. This is
  across the archived phases below, not a claim of one full-suite run.
- Original targeted v2 run:287 passed,1 environment failure,1 deselected. Its
 36 new cases passed; the old subprocess test failed because relative
 PYTHONPATH=python vanished under cwd=/tmp. Rerun with absolute PYTHONPATH:
 42 passed (41 new at that phase plus that existing subprocess test). The
 old22334 test then passed separately with --runxfail. Two final shifted-family
 cases were added; final new-file run43 passes. No Rust changed during these
 Python phases. No full suite was repeated.
- Inputs cover positive/nonnegative k, arbitrary positive symbolic shift gaps,
  positive/negative numeric shift gaps, k=-p-1 with p>0, exact real/complex
  scalar coefficients, concrete rational/complex exponents, singular k=-1,
  unknown real/negative/complex k, a=b degeneracy, and source poles hidden by
  numerator leading-factor cancellation.
- Public principal-branch results: sqrt(-1-I/n) gives -I/(k+1),
  sqrt(-1+I/n) gives I/(k+1), real-negative approach gives I/(k+1),
  off-cut coefficient sqrt(1+I) stays sqrt(1+I)/(k+1).
- Float coefficients at1/53/139/140/211bits retain their original exact _mpf_
  payload and precision metadata; no comparison tolerance is used.
- Two initial test failures were direct SymPy1.14 TypeErrors for concrete k=I
  and k=1+I, not native result failures. These are now separate strict native
  results with an analytic certificate: (1+w)**k tends1 and the derivative of
  (1+w)**(k+1)-1 at0 equals k+1. Native gives (1-I)/2 and (2-I)/5 exactly.

No test assertion for an already computed native result was weakened. All
exploratory eprintln/TRACE/GUARD statements are absent from the final patch.

## Remaining boundaries

This closes22334 for certified parameter domains, not all symbolic limits.
Generic k, generic real/negative/complex k and unproved shift differences stay
held where the new coefficient/domain proof is undecidable; those domains can
include k=-1 or equal shifts. Unknown branch-side certificates also retain old
paths. The new normalization admits rational homogeneous degrees only; other
degrees keep their prior implementation. No old success is placed behind these
new restrictions. No claim is made about whole Gruntz completeness, the separate
Derivative02fa failures, or a full integrated suite; the parent owns that proof.
