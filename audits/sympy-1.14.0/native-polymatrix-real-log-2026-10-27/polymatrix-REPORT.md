# PolyMatrix exact ground RREF/nullspace — frozen R3-based delivery

PolyMatrix's `_rref_rows` and nullspace construction were real Python mathematical loops used by PRDE, despite ordinary Matrix.rref and QQ DDM elimination already having native implementations. This patch routes stock exact ground PolyMatrix elimination and basis construction into Rust, sharing the existing DDM first-pivot engine. No main or previously frozen stage was modified. No performance measurement or speed claim was made.

Final stage: `/home/keckj/.cache/symox-validation/native-polymatrix-ground-r3-20261027/work`.
Patch and all evidence: sibling `proof/`.
The first private port remains frozen at `native-polymatrix-ground-20261027`; use this newer R3-based delivery for integration.

## Identity and integration

The base is immutable Root R3 `/home/keckj/.cache/symox-validation/native-inexact-20261027/after`. Every one of its 2649 source entries and its original release binary were independently checked against its manifest. The inherited validation metadata was removed from the working copy so it cannot be mistaken for final patched validation metadata.

- Root R3 manifest identity: `1ad26aa2fbe980b159a50f38583dce6fc0a301045b92194803207c059ef07875`.
- Root R3 release native: `8a8a8e126be486e31b83b7cae2c306e54ab21e9d569f01bbe46adeb8bb661f3a`.
- Our base code-only manifest: `3d2001dabd72b11532d964a5ad44275fec9161d6668677ed7a09b02bfb1525aa`.
- Final code-only manifest (2651 files): `7cbe05e88d576b83fef4d7fac903f9d000da4b851bc3148490ff68008cad59d0`.
- Patch: `7c4f1582707e7cf4938dd94a9a445ba4578b80d65641ac902ca34a9e21278d33`.
- Final private DEBUG native: `a3b0fe87684dfdafded74c5b37de0563cdf4610a131a37a4ce542961df9bd5da`.

`identity.json` lists every changed file's base/final hash. `r3-provenance.json` records Root's identity and complete base verification. `source-final.sha256` was checked in full on Demeter after the final tests. `git apply --check` passes on base.

Six changed files:

1. `crates/symox-core/src/poly/ddm.rs`: shared generic exact Field reduction and basis construction, existing QQ first-pivot delegation, four new Rust tests.
2. `crates/symox-py/src/poly_ddm.rs`: strict PolyMatrix representation/field checks, native coefficient preparation, result construction and error/fallback boundary.
3. `python/symox/polys/polymatrix.py`: dispatch around the original complete RREF/nullspace paths; result reshaping/factories only.
4. New `python/symox/polys/_polymatrix_native.py`: protocol authentication only, no mathematical algorithm.
5. `python/symox/__init__.py`: short end-of-package registration of the final cancellation/sympification protocols.
6. New `tests/test_polymatrix_native_ground.py`: 86 tests.

**Integration instructions:** preserve current changes by applying hunks. Removing the two insertions from the final `poly_ddm.rs` reproduces the R3 binding byte for byte, hash `58ef0495f1414258ff66cc2b807d4c1ed2a3f8309ba195fcea6681188aff121b`: Berkowitz and prior inexact bindings remain intact. The registration hook belongs at the **very end** of the current `__init__.py`, after any newer export fixes or installations. Do not overwrite those later changes with the older complete file.

## Actual call paths and supported native scope

`integrals/prde.py` imports `PolyMatrix as Matrix`; `constant_system` calls its rref, and parametric Risch routines call its nullspace. Before this patch those calls stayed in PolyMatrix's Python loops rather than using the existing native DDM or Matrix routines.

The native scope is QQ and rational function fields with ZZ/QQ base and formal independent commutative native Symbol parameters, including several parameters. The selected domain is retained. In particular, a real PRDE `constant_system` fixture infers **ZZ(y)** rather than QQ(y); the same integer-polynomial fraction backend handles it exactly without replacing its parent. The test proves this actual PRDE call enters the native binding, then checks both the independent reference and the original complete Python route.

RR/CC, algebraic grounds, EX, other field representations, non-symbol parameters, host-dependent inputs, modified protocols and Float payloads retain the complete old calculation. Float detection traverses the stored expressions, so a Float hidden inside an Add/Mul is not rationalized merely because the declared ground is exact. No refusal is counted as a solved computation or an empty nullspace.

Public `_ground_check` stays in its original location before dispatch, preserving Poly construction, domain validation and errors. Input storage and its entry objects are not mutated. Output class, ground, generators, shape and ignored rref options remain unchanged. Python retains adaptation and factories; elimination, pivot arithmetic and basis construction execute in Rust.

## Shared algorithm and preserved operation/publication order

The exact engine has two explicit modes:

- Existing DDM first-pivot mode retains inverse-then-multiply normalization from the pivot column onward, explicit pivot-zero publication, omitted preceding columns, row swaps and per-cell publication order. Its wrapper maps the generic events back to the existing event type. Existing partial pivoting remains on its previous code path unchanged.
- PolyMatrix mode retains its first nonzero pivot, **e/p on every column**, followed by **a-c*b on every column**, rather than introducing a different fraction-free/modular sequence. Only exact native field arithmetic is admitted.

The nullspace engine enumerates free columns in ascending order, sets that free entry to one, and uses negated reduced pivot entries. Its empty result denotes proven full column rank. Unhandled representations have a separate None route. Exact rational function arithmetic uses the existing RF sparse polynomial algorithms; its ambient parameter count is scoped with RAII and restored on success, refusal or interruption.

Checks occur at entry, pivot search, rows and cells, and basis construction. Signal errors propagate as errors, never as successful refusals or empty bases. Existing individual coefficient operations retain their own granularity. The separate previously reported CC exact-search checkpoint concern is not modified by this patch.

A typed exponent-storage cancellation encountered in the new rational-function computation returns to the original full expression calculation, rather than reducing capacity. Unexpected errors/panics are not blanket-caught. Shape values that cannot fit usize also request the old route: public zero-row RREF with `N=2**80` columns still succeeds. Nullspace checks only actual Rust vector representation limits and preserves the complete fallback when its actual dense output cannot be represented. No arbitrary work quota was introduced.

## Protocol authentication

At the end of package initialization the adapter records the canonical cancellation wrapper chain, its function code/closure delegates and referenced helper identities, relevant builtins and native expression classes. It also authenticates the sympify delegate used for basis zero/one values and its converter registry. Registry checks compare raw stored entries by identity without invoking user key equality. Class metadata is read through captured native type descriptors, avoiding speculative user getters.

Custom `_entry`, `_cancel`, `_like`, instance reduction overrides, subclasses, list-subclass storage, Expr operators, public cancel, wrapped-cancel closure delegates, module builtins, sympify delegates and registered integer converters preserve their complete previous routes and callback order. Merely copying a wrapper name/metadata cannot authenticate a user callback. The tests assert both no speculative callback during selection and identical actual callback traces against the old route.

## Evidence on Demeter

Private tree: `/home/keckj/symox-ci/wt/polymatrix-ground-r3-20261027`.
Private target: `/home/keckj/symox-ci/targets/polymatrix-ground-r3`.
Build uses Cargo directly with PYO3_PYTHON `/home/keckj/symox-ci/venv/bin/python`. Only the private package's .so was copied; no maturin or shared installation occurred.

- `cargo build -p symox-py`: passes (`polymatrix-ground-r3-build.log`).
- Final focused Python command: **423 passed**, comprising **86 new** and **337 existing** (`polymatrix-ground-r3-focused-final.log`). It includes the existing DDM publication/alias tests, DomainMatrix algorithms, Risch parity and the received native Berkowitz tests.
- `cargo test -p symox poly::ddm`: **12 passed** (four new plus eight existing; `polymatrix-ground-r3-rust-final.log`). These include exact RF reconstruction, distinct publication sequences, empty/invalid shape distinction and interruption propagation.
- `cargo clippy -p symox-py`: exit success (`polymatrix-ground-r3-clippy-final.log`). Repository-wide pre-existing warnings remain; there are no diagnostics in the newly added functions. The two diagnostics in ddm.rs concern its unchanged legacy partial-pivot loops.
- All 2651 final source entries verified (`polymatrix-ground-r3-source-check.log`).

The exact tests compare full RREFs, pivots, canonical ordered bases, shapes, classes and parents with SymPy1.14. Every basis satisfies **A*v=0** in independent exact rational/function arithmetic, with basis count equal to columns minus rank. The old complete Python path remains a second differential oracle. Inputs include zeros, rectangular/empty matrices, rank deficiency, row swaps, several formal parameters, and coefficients beyond 521 bits. The PRDE constant-system and parametric polynomial-Risch fixtures prove the intended upstream calls use the kernel and preserve their exact solutions.

No benchmark was performed during the concurrent full suites/corpora. Pytest/build durations are validation logs, not speed claims.

## Replay

Stage a clean copy of `work/` without any inherited .so, copy `proof/source-final.sha256` to its root, and run `proof/replay-demeter.sh` on Demeter. The script accepts SYMOX_STAGE, SYMOX_TARGET and SYMOX_PYTHON overrides, verifies the source manifest, builds/copies a private ABI module, confirms the imported binary matches that build, and runs the exact focused/Rust/clippy commands. `finalize.py` regenerates the patch/manifests from base/work and asserts that removing the binding insertions restores base and that the protocol hook is last.

## Explicit remaining gaps

This is a first significant PolyMatrix ground port, not a claim that PolyMatrix or Risch/PRDE is wholly native. PolyMatrix multiplication, construction/validation and other orchestration remain in Python. Nonselected domains keep their existing complete paths; no general RR/CC rounding parity is claimed for PolyMatrix's old expression route.

Nested QQ(y)(z) is already rejected by the public Poly constructor in the original `_ground_check` for nonempty PolyMatrix input. The generic RF engine can represent flattened functions, but this earlier public gap is not bypassed or counted as a native gain. The independent immutable-baseline probe in `baseline-reference-gaps.log` records it.

SymPy PolyMatrix with zero columns raises ZeroDivisionError while iterating its ground check in this environment (`baseline-reference-gaps.log`). Symox's prior mathematically correct empty-shape behavior is retained and independently checked; the reference failure is not copied or advertised as a new capability. A full integrated release/suite remains the parent's gate after applying the patch.
