# Matrices over coefficient domains

[README](../README.md) · [Matrices](matrices.md) · [Domains](domains.md)

`symox.polys.matrices` ships the reference `domainmatrix`, `ddm`, `sdm`,
`dense`, `rref`, `lll`, `normalforms`, `linsolve`, `dfm` and `domainscalar`
modules, ported over native typed coefficients. `DomainMatrix` wraps a dense
(`DDM`) or sparse (`SDM`) representation, as in SymPy; arithmetic does not
convert coefficients into scalar expression matrices or import SymPy. The hot
kernels for the stock domains run in Rust (see below); other domains, and any
modified protocol, run the ported reference algorithm.

```python
from symox import ZZ, QQ, GF
from symox.polys.matrices import DM, DomainMatrix

A = DM([[2, -1], [-1, 2]], ZZ)
numerator, denominator = A.inv_den()
assert A * numerator == DomainMatrix.eye(2, ZZ).to_dense() * denominator
assert A.det() == ZZ(3)
assert A.charpoly() == [ZZ(1), ZZ(-4), ZZ(3)]

reduced, pivots = A.rref()
assert reduced.domain == QQ
assert pivots == (0, 1)

modular = A.convert_to(GF(3))
assert modular.rank() == 1
kernel = modular.nullspace()
assert (modular * kernel.transpose()).is_zero_matrix
```

## Construction and arithmetic

Construction accepts rectangular rows or a dictionary of row dictionaries,
an explicit shape and a native domain. `DM(rows, domain)` infers the shape.
`zeros`, `ones`, `eye`, `diag`, `from_list`, `from_Matrix` and the
`from_dod`/`from_dok`/`from_flat_nz` constructors are available. Constructors
copy their inputs and convert entries to native elements of the chosen domain
(a Python `1` in a `QQ` matrix is stored as `QQ(1)`, where the reference keeps
the `int`). `to_Matrix`, `to_dod`, `to_dok`, `to_list`, `to_list_flat`,
`to_flat_nz`, `to_dense`, `to_sparse`, `copy`, `convert_to`, `to_field`,
`choose_domain` and `unify` expose conversion and views. As in SymPy, `eye`,
`zeros` and `from_Matrix` build sparse matrices, and a dense and a sparse
matrix compare unequal even with the same entries. `to_dfm()` raises
`NotImplementedError`, as in SymPy without python-flint.

Addition, subtraction, matrix multiplication, scalar multiplication, field
division and integer powers use native coefficient arithmetic. Operators
unify matrix domains; explicit `add`, `sub` and `matmul` require matching
domains. Negative powers require an invertible square matrix over a field.
Transpose, repeated-index extraction, slices, horizontal and vertical
concatenation preserve the domain. Scalar indexing returns a `DomainScalar`
with `.element`, `.domain` and `.to_sympy()` (a native symox expression).
`clear_denoms(convert=True)` returns a scalar denominator and numerator matrix
over the associated ring; without conversion the numerator retains its field.
`clear_denoms_rowwise`, `cancel_denom`, `content`, `primitive`, `scc`,
`eval_poly` and `eval_poly_mul` follow the reference.

## Elimination and solving

`rref()` returns a matrix over the associated field and a tuple of pivot
columns. `method='GJ'` uses Gauss-Jordan division, `method='FF'` exact
fraction-free elimination and `method='CD'` clears denominators first; the
`_dense` variants force dense storage. `rref_den()` keeps a common denominator
in the coefficient domain. The default `auto` method follows the reference
heuristics for choosing the algorithm and storage. Other algorithm names are
rejected explicitly. Sparse matrices use the sparse `sdm_irref` and
`sdm_rref_den` eliminations.

`rank`, `nullspace`, `nullspace_from_rref`, `rowspace` and `columnspace` are
available. Nullspace basis vectors are rows; the default need not normalize
or remove common factors. Over `ZZ`, nullspaces are fraction-free and exact
(the reference returns floats there). `divide_last=True` requires a field. Row
and column spaces require fields (`DMNotAField` otherwise); `rowspace` returns
the pivot rows of the transpose, which are independent, where the reference
can keep a dependent row. A different valid basis or common denominator may
be returned by SymPy.

`inv()` requires a field. `inv_den()` and `solve_den(b)` return a numerator
matrix and scalar denominator in the original domain, satisfying
`A * numerator == b * denominator`. Solving accepts square and consistent
overdetermined systems with a unique solution; singular, inconsistent and
underdetermined systems raise `DMNonInvertibleMatrixError`. The optional
`method='charpoly'` requires a square matrix and uses its adjugate.

`det()` uses Bareiss elimination; sparse matrices are converted to dense
storage first, as in the reference. `charpoly()` returns descending
coefficients using the division-free Samuelson-Berkowitz recurrence after
splitting block-triangular structure; `charpoly_factor_list`,
`charpoly_factor_blocks`, `charpoly_base` and `charpoly_berk` expose the
steps. `adj_det()`, `adjugate()` and `adj_poly_det()` return the actual adjugate
and determinant and support singular matrices. Empty square matrices have
determinant one and characteristic polynomial `[one]`. Negative powers return
`inv(M)**n`, where the reference raises `NotImplementedError`.

## Factorizations and normal forms

`lu()`, `lu_solve()` and `qr()` require a field and raise `DMNotAField` or
`DMDomainError` on rings, as in SymPy; `fflu()` returns the fraction-free `(P, L, D, U)` factorization over
rings. `symox.polys.matrices.normalforms` provides `smith_normal_form`,
`smith_normal_decomp` (with both unimodular transformations),
`invariant_factors` and `hermite_normal_form` (ordinary, modulo `D` and
`check_rank`). `lll()` and `lll_transform()` reduce integer row lattices; the
`Matrix` methods documented in [matrices](matrices.md) use the same kernel.

## Native kernels

For the stock domain objects and native elements, these kernels run in Rust
and reproduce the reference operation order, pivots and published objects:
DDM products and dense Gauss-Jordan over `ZZ`/`QQ`; sparse Gauss-Jordan over
`QQ`; dense and sparse Gauss-Jordan, LU and triangular solves over `GF(p)`
with `p < 2**32` and algebraic fields `QQ<a>` (Gauss-Jordan also over `QQ_I`);
LU and triangular solves over `QQ`; fraction-free `rref_den`, FFLU, Bareiss
determinants and dense/sparse Berkowitz over `ZZ`/`QQ`, `GF(p)`, `ZZ_I`/`QQ_I`,
`QQ<a>` and polynomial domains; Hermite normal forms over `ZZ` (ordinary and
modulo `D`); Smith decompositions over `ZZ`, `QQ`, `GF(p)`, `QQ<a>`, `QQ_I`,
`QQ[x]` and rational function fields; and exact integer LLL. Every cell must
be an exact native element of the matrix domain, otherwise the call declines
before writing and the ported Python algorithm runs. Composite `GF(n)`,
numerical and expression domains, and `QQ_I` LU, keep the Python algorithms.

Numerical fields retain their precision through the algorithms. Rank and
pivot selection use the domain's zero predicate; they do not apply a numerical
conditioning threshold. Numerical identities should be checked with the
domain's `almosteq`, allowing for finite precision. Over `EX` and `EXRAW`,
`rref` and `inv` use the native expression kernel, so zero pivots such as
`sin(x)**2 + cos(x)**2 - 1` are recognized; conditions on symbolic pivots are
not tracked. `sympy.polys.matrices.eigen` is not shipped (it needs
`agca.extensions.FiniteExtension`); see
[the parity gaps](algorithmic-parity.md).
