# Native DDM Bareiss determinant — frozen delivery

## Scope and provenance

Private stage `../ddm_determinant`, immutable base `../ddm_determinant_base`.
Base is the preceding frozen QQ clearing lot, not a copy of the subsequent
main Derivative/conditional-integral work. Integrate only this four-file patch
after the QQ patch; never replace main with this stage. Main, QQ and all other
frozen stages were left unchanged.

- Base source: 2778 files, SHA256 `8e7a4c81ee29ad5eb8904f5b61c8e4bfa44d37548a4f42eb0a645247a92b95f4`.
- Final source: 2779 files, SHA256 `a278a51048e0a0c520708cc85c57e4e1c1fa0852e68fa2c02679863e15582de7`.
- Final source manifest independently reproduced on Demeter, same count/hash.
- Patch SHA256 `23eacf1ca6a8013acdcf47de0155be2e0064a79fecc3c2d3bcea412644017932`; apply--check passes on the immutable base.
- DEBUG extension SHA256 `c00fadc53d6ad7b8522b53c2fea5b3e2884e1fedb3a2b823e19b8492d1d74a3d`.
- Private Demeter tree `/home/keckj/symox-ci/wt/ddm-determinant-20261026`.
- Rust compiled once from source-v1; subsequent source changes added the new Python test only. Both Rust files and the dense adapter are unchanged since that successful build.
- Native build: cargo build -p symox-py, PYO3_PYTHON=/home/keckj/symox-ci/venv/bin/python, CARGO_TARGET_DIR=/home/keckj/symox-ci/targets/expression-sparse. Copied resulting library only into the private Python package.
- No local Cargo, maturin, shared installation, commit, benchmark or main mutation.

## Diagnostic and implementation

Public Matrix.det was already Rust: exact rational matrices use RatMat/Bareiss
and other determinant algorithms are native as well. The actual Python gap was:

    DomainMatrix.det -> DDM.det -> dense.ddm_idet
    SDM.det -> existing dense conversion -> DDM.det -> dense.ddm_idet

The existing matrix/ff.rs Bareiss kernel removes the pivot row and compacts the
matrix at each step. It cannot reproduce the public in-place ddm_idet state,
row identities, swaps and untouched lower triangle. This port therefore adds
poly/ddm::idet using existing exact Coeff/Int/Rat arithmetic and the same
fraction_free::Failure typed errors. It reproduces the original DDM algorithm
with an event for each existing-row swap and each individual cell assignment.
No alternative pivot rule, modular shortcut, denominator clearing, different
arithmetic grouping or class-level factory replacement was introduced.

The binding parses exact native ZZ/QQ cells only, rejects foreign/aliased rows
before any mutation, then publishes each completed write immediately to the
original PyList row. Original untouched coefficient objects and existing row
objects survive. The returned determinant is a fresh typed element, matching
the original final sign multiplication. Exact division errors and signal errors
propagate directly, preserving the published prefix.

The Python change is a small stock-protocol dispatch before the complete old
body. Existing FF domain/operator and DDM list protocol guards are reused.
Every unsupported native representation/domain falls back to that original
body; no newly unsupported result is reported as a speedup.

## Validation on Demeter

- New Python file: **82 passed** on first execution.
- Combined focused matrix run: **588 passed** in19.86s:
  new determinant82, preceding QQ clear_denoms97, native matmul/RREF/FF and
  existing DomainMatrix algorithm tests. Log `ddm-determinant-combined-v1.log`.
- Targeted core Rust: **2 passed**, exact mutation/write ordering and pivot swap;
  empty/non-square shapes and direct interruption propagation.
- 48 seeded exact ZZ/QQ direct-kernel comparisons against SymPy1.14 ddm_idet
  using the same native scalar objects, plus independent SymPy rational
  determinant certificates. Tests compare full mutated cell state, original
  row permutation, coefficient/return alias partition and domain parents.
- Ten cases with151–219bit coefficient scale, rational denominators and
  multiple pivot swaps. Empty, singular and vanishing leading pivot cases.
- Stock dense/sparse DomainMatrix.det route spies prove actual native use and
  prove that original input values/stored coefficient identities are unchanged.
- Complete old-path comparison for aliased rows, foreign row/outer subclasses,
  raw Python integer/mixed cells, ragged/wide/tall/zero-column representations.
- Callback trace parity for exquo, one, truth, multiplication, subtraction and
  DDM.__setitem__. An exquo callback exception preserves the exact same exception
  object and the already-published mutation prefix.
- RR, polynomial-ring and finite-field cases retain their original algorithms;
  public non-square errors stay unchanged.

All new tests pass. No intermediate Python failure occurred. No full suite or
speedup claim is made; the parent owns integrated release/full Demeter proof.

## Retained boundaries

Native admission is square matrices with distinct builtin-list rows and stock
native ZZ/QQ cells, either in a builtin outer list or the registered stock DDM.
Aliased rows, foreign callbacks/classes, raw/mixed Python cells, unsupported
shapes and every other domain keep the original complete implementation.
This is a native execution boundary, not a new mathematical refusal.

Dense and sparse Berkowitz/characteristic-polynomial routines remain Python,
including ddm_berk and sdm_berk. The public Matrix determinant was already native
and is not counted as a gain. SDM.det still uses its established dense conversion;
no new sparse determinant algorithm or memory-capacity claim is made. Direct
low-level nonsquare ddm_idet's historical behavior is retained by falling back;
normal public DDM/DomainMatrix callers continue to reject non-square matrices.
