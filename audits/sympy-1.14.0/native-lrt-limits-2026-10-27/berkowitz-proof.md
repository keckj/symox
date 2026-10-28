# Native dense DDM Berkowitz / exact characteristic polynomial — frozen delivery

Private source `../ddm_berkowitz`; immutable base `../ddm_berkowitz_base`, copied
from the previously delivered DDM determinant stage. Main and all received
packages were left unchanged. Apply only this five-file patch after the DDM
stage; do not replace main with this tree (later Derivative, conditions, limits
and other work are not part of this older base).

## Provenance

- Base:2779 source files, manifest `a278a51048e0a0c520708cc85c57e4e1c1fa0852e68fa2c02679863e15582de7`.
- Final:2780 source files, manifest `9877ed82d9f268049bb0c7c6591bef5fbdf129901393d64313e548de23ba358a`.
- Final manifest independently reproduced on Demeter with the same count/hash.
- Patch `ddm-berkowitz.patch` SHA256 `72f46d64ceff2410a5288bda4872991a47aa8ae33e1033195350a7d90cacfe04`.
- apply--check passes on the immutable base. Per-file base/final hashes are in file-hashes.json.
- Private DEBUG extension SHA256 `687d44aacc544369db3b681193ab1611c985ae4d093a595cd12c12f8d43b7f3b`.
- Demeter tree `/home/keckj/symox-ci/wt/ddm-berkowitz-20261026`.
- Built once with cargo build -p symox-py, PYO3_PYTHON pointing to the private
  Python3.10 venv, target `.../targets/expression-sparse`. Both core files and
  the Rust binding have not changed since that build; subsequent changes were
  Python protocol validation and tests. The library was copied only into this
  private Python package.
- No local Cargo, maturin, shared installation, commit or timing campaign.

Files: matrix/ff.rs, matrix/linalg.rs, py/poly_ddm.rs,
Python matrices/dense.py and new test_native_ddm_berkowitz.py.

## Actual Python gap and reused native implementation

The uncovered path was DDM.charpoly -> dense.ddm_berk, also reached by dense
DomainMatrix.charpoly_berk and suitable charpoly blocks. Sparse SDM.charpoly
calls the separate sdm_berk, which this port does not change.

The native RatMat::charpoly already implements exact rational characteristic
polynomials using one common denominator D, integer Berkowitz on B=D*A, and
coefficient rescaling c_k(B)/D**k. This port reuses that algorithm intact.
The existing integer Berkowitz summation, zero filtering, operand order and
coefficient order are unchanged. Checked variants add fallible interruption
checkpoints to the existing kernels; the old infallible entry points delegate
with Infallible and preserve their existing values and operation order.

The typed binding admits exact native ZZ/QQ scalar cells in existing dense
storage. It converts the two native exact integer representations without
passing through Expr, Python arithmetic or Float, runs RatMat's checked kernel,
and returns the same [[c0],...,[cn]] shape with typed ZZ/QQ elements. Output rows
are distinct; coefficient objects and their alias partition match the old
algorithm. Input rows/cells, including aliased input rows, are read-only.

The Python change is dispatch/protocol adaptation only. Native eligibility
protects the domain/scalar operators, recursive ddm_berk protocol, Python
imatmul function/code, its native multiplication hook and the old builtin
lookups. Custom hooks keep the complete recursive old implementation.
DDM metadata is read through a validated builtin instance-dict descriptor;
poisoned shape/domain/rows/cols and custom metadata descriptors retain the old
path without introducing getter calls. Foreign/mixed/ragged storage similarly
keeps its old result or precise error. No Unsupported default or dense resource
quota was added.

## Validation on Demeter

- Final new file: **92 passed**, ddm-berkowitz-new-final-v2.log.
- Combined focused DDM/FF and public charpoly run: **803 passed** (87 new at
  that phase plus716 existing), ddm-berkowitz-combined-v1.log. Five metadata
  regression tests and the final metadata guard were added afterward; all92
  new tests then passed. No Rust changed after the successful build, and no
  full suite was repeated.
- Core ff/linalg targeted tests: **13 passed**, including two new tests and
 11 existing. Checked kernel results, rational scaling, empty matrices and
 direct interruption errors are verified.
-44 seeded ZZ/QQ matrices (sizes0–7), complete coefficient-vector equality
 against SymPy1.14 and an independent Matrix.charpoly reference. Exact
 Cayley-Hamilton Horner identities certify P(A)=0. Cases include zero matrices,
 singularity, valid aliased rows, coefficient/return domains, output shapes,
 fresh rows, input immutability and full alias partitions.
- Large exact coefficients through211bits and pairwise coprime denominators.
- Native route spies prove real use by stock DDM.charpoly,
 DomainMatrix.charpoly_berk and DomainMatrix.charpoly.
- Trace parity for one/zero/negation/multiplication/addition/sum,
 Python/native multiplication hooks, recursive Berkowitz hooks and custom DDM
 indexing. Callback exception object identity is retained.
- Complete fallback comparisons for foreign rows/outer containers, tuples,
 raw Python integers, mixed coefficients, inconsistent scalar domains,
 ragged/non-square storage and poisoned metadata.
- RR/CC at1/53/139/140/211bits keep the previous recursive grouping and exact
 Float component _mpf_/precision payloads; CC uses nonzero imaginary input.
- GF(4), GF(5), polynomial-ring and EX fallback behavior matches the original
 division-free algorithm, including the ring with zero divisors.
- SDM.charpoly still succeeds with every new dense/native conversion forbidden
 by a spy: no new sparse densification was introduced.

All tests passed. No native mathematical result assertion was relaxed. The
final alias helper was strengthened after the combined run to use contiguous
identity classes even for repeated source objects; the final92 include that
stronger certificate. Test log runtimes are verification durations, not a
benchmark or a speedup claim.

## Retained boundaries

Only stock exact dense ZZ/QQ goes through the new route. Other domains and
custom or invalid representations retain the complete old algorithm. This is
an execution boundary, not a new mathematical refusal. RR/CC precision and
rounding behavior remain on their established recursive path.

Sparse sdm_berk remains Python; SDM.charpoly is not counted as a native gain.
Existing DomainMatrix.charpoly SCC/block selection, preparation and result
assembly remain their old code, including any densification it already chose.
The entire high-level charpoly call is not claimed to be native. No timings,
large-sparse capability claim, whole-library parity or full integrated suite
claim is made; the parent owns release integration and full Demeter validation.
