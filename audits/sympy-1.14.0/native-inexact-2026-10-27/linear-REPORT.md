# Native solve_lin_sys inexact reduction — frozen delivery

Stage: `/home/keckj/.cache/symox-validation/native-inexact-linear-20261027/work`.
Proofs and patch: sibling `proof/`. No main source, received package, shared environment or installed extension was modified. No benchmark was run.

## Provenance

The base matches all 328 build inputs of immutable r1b `/home/keckj/.cache/symox-validation/native-derivative-r1b-20261027/after`; all production Python/Rust and tests also match that snapshot. The absent five `.validation` runner metadata files do not affect code. See `r1b-base-provenance.json`. Root snapshot identity is `db859db2e5949ea295c362ae9b06b1b7819f0f28da0753790edda29b45445d2f` (release native `cb5c0c98572e4b05eadca9f234fb6df3bbabea93a202e4a2ac292f5541b5748e`). We did not claim the transitory main extension as a baseline.

Our reproducible file cohort, including unchanged auxiliary files, has these manifest hashes:

- Base: `a13553bd4efa8b6079a7e62e29b08a66c7df425483bc7b1983c6c72e53eb9deb`.
- Final (2891 files): `62d8dbd2086febd4708dca0fd8d18baf1c4809045feff11561ef4616621541c5`.
- Patch: `1f434b41dadb4fefccb3849eb33adf7ed16e5fe7192aa3605e44dfcf345f5fef`.
- Private DEBUG native: `acd1ed8ff7d63c38c88bb2a56bd9957adf06ee014d0dc7f04bad252d3d7e3801`.

`identity.json` contains every changed file's base/final hash. `source-final.sha256` was checked successfully in full on Demeter after the final tests. `git apply --check inexact-linear.patch` passes on base.

## Six changed files

1. `crates/symox-core/src/poly/ddm.rs`: typed inexact Gauss-Jordan kernel plus two tests.
2. `crates/symox-py/src/poly_ddm.rs`: strict native lists/DomainElements binding and registration.
3. `crates/symox-py/src/domains.rs`: only expose the existing `element_error` to sibling bindings, preserving existing exception classes/messages.
4. `python/symox/polys/solvers.py`: dispatch the elimination loop, preserving the original complete loop as fallback and leaving input conversions and answer construction in place.
5. `python/symox/polys/_inexact_linear.py`: protocol eligibility only, no mathematical algorithm.
6. `tests/test_native_inexact_linear_system.py`: 203 tests.

**Integration:** merge `poly_ddm.rs` by hunks with the newly received Berkowitz binding, preserving `_ddm_berk` and its registration. This delivery's file is intentionally based on immutable r1b, which precedes Berkowitz. Do not replace the whole current binding file. No other frozen stage was changed.

## Mathematical and protocol contract

The existing matrix inverse kernel has magnitude partial pivoting and forced pivot zeros/ones, so it cannot reproduce this solver. The new kernel instead reuses `DomainElement.binary`, the same scalar operation as the established Python arithmetic. It scans the first nonzero pivot, calculates one/pivot, multiplies every column by that reciprocal, and eliminates each other row as `a-factor*b` in column order. It does not rationalize Floats, fuse operations, omit earlier columns, force residuals to zero, change p to presentation precision, change the selected domain, or substitute RR/CC for a composite.

Plain RR/CC and parameter polynomial/fraction domains retain their complete selected domain. A non-unit polynomial-ring pivot has the same error as before; it is not promoted to a fraction field. The existing augmented-column consistency test and public result assembly remain unchanged. A binding `None` denotes only an unhandled protocol/input and invokes the old full calculation. Arithmetic errors and interruptions propagate separately; they never become an inconsistent system or empty answer.

Original `K.convert` calls remain before selection; original `K.to_sympy`, ring factories, gens, and result conversions remain after reduction. Native selection proves stock scalar arithmetic, domain methods, installed Gaussian division wrapper/delegate, relevant builtins and attribute probes. Customized operators, conversions, getters, module builtin callbacks and scalar probes retain the old calculation with the same callback order. Foreign domains/storage, host-dependent or nonfinite coefficient payloads also preserve the full fallback. Nonfinite inputs rejected by the original conversion still raise that original error before dispatch.

## Demeter evidence

Private tree: `/home/keckj/symox-ci/wt/inexact-linear-20261027`; private target `/home/keckj/symox-ci/targets/inexact-linear`. Build used Cargo directly with `PYO3_PYTHON=/home/keckj/symox-ci/venv/bin/python`, then copied only that private `.so` into that private package. No maturin or shared installation.

- `cargo build -p symox-py`: pass (`inexact-linear-build-v1.log`).
- `cargo test -p symox inexact_linear_tests`: **2 passed** (`inexact-linear-rust-v1.log`). Includes cancellation/error distinction, domain refusal, first-nonzero pivot and interruption propagation.
- Final new file: **203 passed** (`inexact-linear-final-new.log`).
- Existing linear-system file: **21 passed**, together with the then-current 197 new cases (`inexact-linear-final.log`: 218 passed).
- Existing compact-conjugates/cartesian Heurisch and inexact-domain-followup files: **120 passed**, in the earlier focused run (`inexact-linear-focused-v3.log`). That run also recorded two new certificate assertions incorrectly parsing the displayed p1 string; these were fixed to inspect the exact stored rational value. No production code changed for that fix, and all final 203 new tests pass.

Thus **344 distinct Python tests** and 2 Rust tests passed; this is a union of bounded runs, not a claim of one combined or full-suite run. The parent owns full integration/release validation.

The new tests compare every RREF entry's Float mantissa/exponent/precision, pivot list, domain, and unchanged input, against the unchanged old reduction loop using identical converted inputs. Cases span p=1/53/139/140/211, real/complex coefficients, first tiny pivot versus a larger candidate, row swaps, rank deficiency, inconsistent/underdetermined systems, exact dyadic reconstruction, mixed rationals, tiny values, composite rings/fields, raw result factories, error classes/text and callback traces. Eight complete RREF comparisons against SymPy use the same exact stored binary inputs (RR at all five precisions; CC at 1/53/139). Independent dyadic solution certificates cover all ten real/complex precision combinations. Real Heurisch calls prove native dispatch and preserve primitive bits plus exact differentiated identities at 53/140 bits; two pre-existing unknown Heurisch outcomes remain unknown.

## Explicit remaining gaps and scope

`reference-complex-audit.json` and its reproducer record a **pre-existing** discrepancy of CC scalar arithmetic versus SymPy for a two-row rational-coefficient witness at p140 and p211. Every native entry is bit-identical to the pre-port loop, but some components differ from SymPy. Inputs and all differing output components are recorded as exact binary rationals; no presentation precision is used as arithmetic precision. This port preserves the old arithmetic rather than quietly changing it. It does not establish universal CC scalar parity.

The existing Heurisch unknown for `exp(Float('0.5')*(1+I)*x)` remains unknown at p53 and p140. It is not counted as a solved integral. Plain imaginary exponential and real exponential cases have independently certified successful primitives.

PolyMatrix ground multiplication/RREF/nullspace, exact-domain solve_lin_sys orchestration, sparse Berkowitz and other audit debts are outside this patch. Protocol fallback and unsupported scalar operations remain the established complete paths, not reduced-capability substitutes. There is no speed claim without release measurement.
