# Symbolic matrices

[README](../README.md)

`MatrixSymbol`, `MatAdd`, `MatMul`, `MatPow`, `Transpose`, `Adjoint`,
`Inverse`, `Identity`, `ZeroMatrix`, `OneMatrix`, `HadamardProduct`,
`HadamardPower` and `KroneckerProduct` share a native typed matrix expression
tree. Products preserve their order and validate dimensions, including symbolic
dimensions. Matrix elements participate in scalar arithmetic, substitutions,
cached free symbols, pickling and SymPy conversion. `as_explicit()` materializes
matrices with concrete dimensions for the existing exact algorithms.
`A.I` aliases the inverse. Public coefficient views distinguish scalar
`as_coeff_Mul()` from matrix `as_coeff_mmul()` and `MatMul.as_coeff_matrices()`;
the latter preserve matrix factor order, and `as_coeff_mmul()` retains an
unevaluated matrix product even for a single factor. `as_coeff_mul()` supports
rational coefficient selection and dependency partitioning. Instance-based
`has()` queries cross scalar/matrix parents without expanding shared graphs.
Dependency partitions must retain dimension-compatible matrix products.
Dimensionless `GenericIdentity` and `GenericZeroMatrix` results are explicitly
unsupported in empty dependency partitions and singleton `MatAdd` tails.
`DotProduct(a, b)` holds an ordered scalar contraction of row or column vectors
with equal lengths. It is bilinear over complex entries; request conjugation
explicitly for a Hermitian pairing. `.doit()` evaluates finite contractions or
retains a symbolic Sum for a symbolic length. Scalar differentiation, rewrites,
pickling, SymPy interchange and finite code generation preserve both matrix
parents. Empty vectors contract to zero.
Dense matrices support exact `exp()` and `analytic_func(f, x)`. Polynomial
functions evaluate directly on the matrix, independently of spectral solving;
other functions use exact Hermite interpolation, including repeated eigenvalues
and defective matrices. Poles raise `ValueError`, and unresolved exact spectra
raise `NotImplementedError`. Results are ordinary scalar expressions usable by
the existing code printers, including PyTorch/JAX JIT compilation and gradients.
Dense `log()` applies the existing scalar logarithm branch to each eigenvalue
through the same exact Hermite calculus. For example,
`Matrix([[2, 1], [0, 2]]).log()` gives `[[log(2), 1/2], [0, log(2)]]`;
negative real eigenvalues use `log(abs(value)) + I*pi`. This is a primary
matrix logarithm with those scalar determinations, including on the negative
real axis, rather than a promise of analyticity across the scalar branch cut.
An explicit zero eigenvalue or nonfinite input raises `ValueError`. Generic
parameter formulas apply at finite values with defined denominators; recompute after
specializing to singular or coincident-eigenvalue cases. `exp(log(A)) = A` is
the reconstruction identity; `log(exp(A)) = A` and logarithms of products are
not unconditional identities. No numerical eigensolver is used.
`log(simplify=False)` and `log(simplify=None)` disable optional final
cancellation; necessary internal simplifications and exact certificates remain.
The default final simplifier is native `cancel`. A supplied callable instead
postprocesses the resulting matrix once, with the original mutability preserved;
`True` and other truthy noncallables raise `TypeError`.
Generated scalar logs retain the backend's real/complex dtype semantics. If a
parameter can make a log argument negative, pass a complex input or request
`dtype=torch.complex128` / `dtype=jax.numpy.complex64` (as appropriate) when
calling the generated function to evaluate the complex branch.
Matrix-function operations check conservative work, storage, coefficient and
expansion bounds before large intermediates, including held numerical powers
and spectral substitutions. Identity and nilpotent shortcuts retain supported
large powers.
`DiagMatrix` turns a row or column vector into a square diagonal matrix;
`DiagonalOf` extracts the main diagonal as a column vector, and `DiagonalMatrix`
retains a matrix's shape while zeroing its off-diagonal entries. These nodes
preserve empty shapes, support symbolic dimensions, and lower to finite
PyTorch/JAX outputs alongside the other matrix expressions.
`FunctionMatrix(rows, cols, Lambda((i, j), body))` defines entries through a
native two-argument symbolic function. It supports scalar derivatives,
substitutions, explicit materialization and finite PyTorch/JAX output generation,
including JIT compilation and gradients. Its `.lamda` property follows SymPy's
spelling. Existing native symbolic function heads are accepted; strings and
arbitrary Python callables are rejected.
`A.applyfunc(Lambda(z, body))` and `ElementwiseApplyFunction` apply a scalar
function to every matrix entry while retaining the parent matrix graph.
They support identity and composition, transformations, expansion, pickling
and finite PyTorch/JAX outputs. Scalar differentiation includes the chain rule
and explicit parameters in the Lambda body. General non-holomorphic complex
chain derivatives remain unsupported.
Hadamard operations act element by element and require matching matrix shapes;
Kronecker products multiply dimensions and preserve factor order. These nodes
also support finite matrix output generation for PyTorch and JAX.
`BlockMatrix` validates block row heights and column widths, and exposes a typed
immutable grid through `.blocks`. `block_collapse` combines compatible block
sums, products and nonnegative integer powers while retaining unsupported
operations. `BlockDiagMatrix` also supports rectangular diagonal blocks.
Two-by-two blocks expose `schur(mat="A")` and `schur(mat="D")`, plus
`LDUdecomposition()` and `UDLdecomposition()`. They preserve rectangular outer
blocks when the selected diagonal pivot is invertible and square. Singular
Schur complements and empty blocks are supported. Generalized inverses and
B/C pivots are explicitly unsupported.
`MatrixSlice` and matrix slice indexing follow Python clipping and stride rules,
including negative steps, omitted bounds, empty ranges and nested slices.
Scalar indices preserve a one-row or one-column matrix when combined with a
slice. Symbolic strides require a known nonzero sign; normalized affine ranges
retain their exact lengths through reconstruction, transformations and pickling.
Matrix hashes, symbol discovery, rewrites, sort keys and pickles preserve shared
DAG nodes, including derivatives of large matrix powers. The versioned pickle
codec records scalar and matrix nodes directly, preserving held arithmetic,
function metadata and arbitrary-precision numbers without printing expressions.
Scalar differentiation preserves product order, supports integer powers and
inverses, and uses shared matrix subexpressions for large powers. `Trace` and
`Determinant` retain abstract parents while applying supported exact identities;
determinant derivatives with concrete dimensions use the adjugate, including
at singular matrices. Complex non-holomorphic matrix derivatives that cannot
be represented correctly remain explicit errors.

```python
import symox as sm

A, B = sm.MatrixSymbol("A", 2, 3), sm.MatrixSymbol("B", 3, 2)
(A*B).shape                      # (2, 2)
(A*B)[0, 0]                      # finite scalar contraction
f = sm.extra.torch_lambdify((A, B), A*B, compile=True)
g = sm.extra.jax_lambdify((A, B), A*B, jit=True)
```

PyTorch and JAX functions accept nonempty `MatrixSymbol` inputs with concrete
dimensions. Input matrix dimensions occupy the last two axes; leading batch
axes broadcast with scalar inputs. Generated matrix outputs follow the same
convention and support the backend's autodiff and compilation. These operations
currently expand to scalar entries at generation time. Matrix inputs with
symbolic dimensions and matrix-valued inputs to
the Slang/SYCL/CUDA function exporters are not yet supported.
Trace and determinant outputs with concrete dimensions are lowered to scalar
entries for native and Torch/JAX generation. Ordinary `.doit()` preserves
abstract matrix symbols; explicit materialization is confined to generation.

`Matrix` provides mutable dense matrices; `ImmutableMatrix` and
`ImmutableDenseMatrix` provide immutable, hashable matrices. Constructors
accept nested sequences, NumPy arrays, SymPy matrices and nested matrix
blocks. Arithmetic, slicing, entry transformations, symbolic derivatives,
Jacobian/Hessian construction, exact determinant/inverse, row reduction,
linear solves, characteristic polynomials and supported eigenproblems run
in Rust. Cholesky, LDL, QR and diagonalization are also available.

`rank_decomposition`, `pinv(method="RD")`, `pinv_solve`, triangular solves,
`QRsolve`, `cholesky_solve`, `LDLsolve` and `solve_least_squares` use exact
native algorithms. Complex matrix solves use Hermitian conjugation.
Row reduction, rank, subspace, LU, determinant and inverse methods accept custom
zero tests. Pivot and characteristic coefficient transformations are also
available; see [matrix zero tests and callbacks](matrix-zero-tests.md).
Entrywise substitution supports simultaneous rules and one-shot pair iterators;
see [expression substitution](expression-methods.md#substitution).
`inv` also accepts `CH`, `LDL` and `QR` methods. `pinv(method="ED")` uses
exact spectral projectors for supported Hermitian eigenproblems.
`singular_values`, spectral norms (`norm(2)` / `norm(-2)`),
`condition_number`, `jordan_form` and `jordan_cells` are available for supported
exact eigenproblems. QR solves require full column rank. Methods requiring a
property that cannot be established symbolically raise an explicit error.
`singular_value_decomposition()` returns a condensed exact `(U, S, V)` with
`A = U*S*V.H` and orthonormal columns. The native implementation certifies the
spectral equations and reconstruction, handles rectangular, complex and
rank-deficient inputs, and preserves matrix mutability. It bounds intermediate
work and rejects floating-point inputs or undecidable rank and spectra rather
than using a numerical fallback. The resulting expressions feed all scalar
code-generation targets, including batched PyTorch and JAX execution.
`inv_mod(m)` computes exact inverses over positive composite or prime moduli,
using Bézout row operations and certifying both matrix products modulo `m`.
It preserves mutability, supports exactly integral floats without materializing
large dyadic integers, and shares resource checks with Sparse delegation before
densification. Modulo one returns the zero matrix of the trivial ring, avoiding
SymPy 1.14's error on nonempty inputs. Nonpositive moduli and noninteger entries
are rejected.
`permute`, `permute_rows`, `permute_cols`, `permuteFwd` and `permuteBkwd`
return new matrices with SymPy's permutation conventions: array form lists the
source indices, while `direction` reverses the order of a sequence of cycles.
`elementary_row_op` and `elementary_col_op` implement scaling, swaps and adding
left scalar multiples. These operations preserve shape and mutability. Sparse
transformations operate on stored entries and finite permutation support,
including axes too large for dense storage. If a coefficient such as infinity
turns implicit zeros into `nan`, the affected line is materialized only after
a resource check. Mutable `row_op`, `col_op` and `zip_row_op` apply callbacks
sequentially, allow callbacks to read the matrix, and keep earlier updates if a
later callback fails. Transformation indices must be integers; negative indices
are accepted by the mutable callbacks, but not elementary operations or
permutation specifications. Conservative native allocation and arithmetic
limits raise `NotImplementedError` before a large transformation is attempted.
Dense, immutable and sparse matrices also provide SymPy's auxiliary
`MatrixBase` methods. These include `berkowitz*`, `LUdecomposition_Simple`,
`LUdecompositionFF`, `per`, `inverse_BLOCK`, `cramer_solve`, `rref_rhs`,
`upper_hessenberg_decomposition`, `bidiagonal_decomposition`,
`lower/upper_triangular`, `rot90`, `vech`, `hat`, `vee`, `dual`, `diagonal(k)`,
`equals`, `pow(method=...)` and the container protocol (`key2ij`, `copyin_*`,
`fill`, `row_add`, `todok`, `table`, ...). The constructors `jordan_block`,
`companion`, `wilkinson`, `irregular` and `orthogonalize` are available, as are
`liupc` and `row_structure_symbolic_cholesky` on sparse matrices. The top-level
helpers are `GramSchmidt`, `casoratian`, `jordan_cell`, `randMatrix` (same
seeded output as SymPy), `banded`, `rot_axis*`, `rot_ccw_axis*`, `rot_givens`,
`symarray`, `list2numpy`, `matrix2numpy`, `matrix_multiply_elementwise`,
`MutableMatrix` and `MatrixBase`, which is an isinstance facade.
Householder results are exact; bases and signs can differ from SymPy's.
The bidiagonal tuple is `(U, B, R)` with `A = U*B*R`, matching SymPy's
returned right factor (`V.H`). Zero pivots and complex reflectors are supported.
Classic dense constructors are bounded to 100,000 entries. `irregular(0, ...)`
rejects nonempty input instead of looping indefinitely. Symbolic powers of
unipotent matrices use the native analytic-function engine. Cramer's rule
rejects mismatched right-hand-side row counts with `ShapeError`, including
cases where SymPy 1.14 instead raises an incidental `IndexError`. `equals` follows SymPy's heuristics conservatively and may return
`None` where SymPy decides. `connected_components_decomposition()` and
`strongly_connected_components_decomposition(lower=True)` return factors
`(P, B)` satisfying `P.T * B * P == A` after explicit evaluation. `B` is a
native `BlockDiagMatrix` or `BlockMatrix`; `P` is a native symbolic
`PermutationMatrix`. Its storage and composition do not require a dense
permutation factor. [Permutation matrices](matrix-permutations.md) documents
this node and held row/column `MatrixPermute` operations.
`to_DM(domain=None, **options)` returns a `DomainMatrix` from
`symox.polys.matrices`, storing native domain elements. Inference covers
integer/rational and Gaussian coefficients, polynomial rings and rational-function fields;
`field=True` chooses a field, and `extension=True` enables supported algebraic
extensions. Numerical coefficients select `RealField` or `ComplexField` at the
maximum source floating-point precision. General expressions use `EXRAW`;
function expressions such as `sin(x)` are not polynomial-ring generators.
An explicit native domain also supports finite fields. Conversion back with
`to_Matrix()`, domain changes, dense/sparse list or dictionary views, native
coefficient arithmetic and determinants are supported. See
[domain matrix algorithms](domain-matrices.md) for elimination, solving,
nullspaces and characteristic polynomials.
Not supported: `DeferredVector`, `MatrixKind`, and the `Basic`/`Expr` methods of immutable
matrices.

`lll(delta=Rational(3, 4))` reduces an integer row-lattice basis;
`lll_transform()` returns `(B, T)` with `T*A=B` and unimodular integer `T`.
The native deterministic algorithm follows SymPy's reduction order and
`floor(mu+1/2)` tie convention. It certifies size reduction and Lovasz
inequalities by fresh exact Gram-Schmidt, and checks an integer two-sided
inverse of `T` to certify that the whole lattice is preserved. Inputs must
have independent rows, at most as many rows as columns, and integer entries;
exactly integral floats are converted to integers, with bit limits checked
before compact dyadic values are materialized. Other entry domains and held
arithmetic are rejected without evaluation. The rational parameter must satisfy
`1/4 < delta < 1`; float parameters follow SymPy's rational-field coercion.
Zero-row-count shapes preserve their width and return a `0x0` transform.
A single all-zero row raises `DMRankError` because it is not a basis, fixing
SymPy 1.14's accidental acceptance of that case. Shape/range/rank errors use
`symox.polys.matrices.exceptions` (`DMShapeError`, `DMValueError`, `DMRankError`),
and noninteger domains raise `CoercionFailed`. Resource bounds cover integer
and rational products, normalization, storage and the certificates, with
`NotImplementedError` when a conservative limit is exceeded. Dense and sparse
results preserve their mutability; Sparse checks the native limits before
densifying. No SymPy runtime or floating-point reduction heuristic is used.
The SymPy algorithm attribution and notices are retained in `LICENSE-SYMPY`.

```python
from symox import Matrix, ImmutableMatrix, symbols

x = symbols("x")
a = Matrix([[x, 1], [1, x]])
a.det()                          # x**2 - 1
a.inv()
a.gauss_jordan_solve(Matrix([1, 2]), freevar=True)
ImmutableMatrix(a)               # rejects mutation; usable as a dictionary key
```

`dot` supports Hermitian conjugation conventions; `diagonalize` supports
`reals_only`, `sort` and `normalize`. Matrix pickles preserve entry
assumptions and unevaluated arithmetic. Immutable results preserve their
immutability through arithmetic, slicing and decompositions; `as_mutable()`
returns an independent editable copy.

`SparseMatrix` / `MutableSparseMatrix` and `ImmutableSparseMatrix` use actual
dictionary-based sparse storage. Construction, indexing, arithmetic,
matrix products, nonnegative integer powers, joins and transposes operate on
stored entries.
`nnz`, `todok`, NumPy/SymPy conversion and pickling are supported. Advanced
factorizations and negative powers explicitly convert to the native dense
algorithms; they are not sparse numerical factorization routines.
The dense delegation also covers `exp`, `log`, `analytic_func`,
`singular_value_decomposition`, `inv_mod`, `lll` and `lll_transform`, converting
matrix results back to sparse storage with their mutability preserved. These methods reject oversized
sparse shapes before allocating a dense input.

Eigenproblems beyond the supported exact algebraic cases and unsupported
constructor/algorithm options are rejected explicitly. Custom zero and pivot
simplification callbacks are supported by the operations described in
[matrix zero tests](matrix-zero-tests.md). `log` accepts its separate final
matrix postprocessor.
Matrix `evalf` supports the same arbitrary precision and substitutions as
scalar expressions.

```bash
python tests/compare_matrices.py
python -m pytest tests/test_matrix_compat.py tests/compare_matrices.py
```
