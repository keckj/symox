# Matrices over coefficient domains

[README](../README.md) · [Matrices](matrices.md) · [Domains](domains.md)

`DomainMatrix` stores native typed coefficients, with sparse internal storage
and independently selectable dense or sparse output views. Arithmetic does
not convert coefficients into scalar expression matrices or import SymPy.

```python
from symox import ZZ, QQ, GF
from symox.polys.matrices import DM, DomainMatrix

A = DM([[2, -1], [-1, 2]], ZZ)
numerator, denominator = A.inv_den()
assert A * numerator == DomainMatrix.eye(2, ZZ) * denominator
assert A.det() == ZZ(3)
assert A.charpoly() == [ZZ(1), ZZ(-4), ZZ(3)]

reduced, pivots = A.rref()
assert reduced.domain == QQ
assert pivots == (0, 1)

modular = A.convert_to(GF(3))
assert modular.rank() == 1
kernel = modular.nullspace()
assert modular * kernel.transpose() == DomainMatrix.zeros((2, 1), GF(3))
```

## Construction and arithmetic

Construction accepts rectangular rows or a dictionary of row dictionaries,
an explicit shape and a native domain. `DM(rows, domain)` infers the shape.
`zeros`, `ones`, `eye`, `diag`, `from_list` and `from_Matrix` are available.
Constructors copy their inputs and discard coefficients equal to zero in the
chosen domain. `to_Matrix`, `to_dod`, `to_dok`, `to_list`, `to_list_flat`,
`to_dense`, `to_sparse`, `copy`, `convert_to`, `to_field`, `choose_domain` and
`unify` expose conversion and views.

Addition, subtraction, matrix multiplication, scalar multiplication, field
division and integer powers use native coefficient arithmetic. Operators
unify matrix domains; explicit `add`, `sub` and `matmul` require matching
domains. Negative powers require an invertible square matrix over a field.
Transpose, repeated-index extraction, slices, horizontal and vertical
concatenation preserve the domain. Scalar indexing returns a `DomainScalar`
with `.element`, `.domain` and `.to_sympy()` (a native symox expression).
`clear_denoms(convert=True)` returns a scalar denominator and numerator matrix
over the associated ring; without conversion the numerator retains its field.

## Elimination and solving

`rref()` returns a matrix over the associated field and a tuple of pivot
columns. `method='GJ'` uses Gauss-Jordan division; `method='FF'` uses exact
fraction-free elimination before converting to the field. `rref_den()` keeps
a common denominator in the coefficient domain. Its default uses division
over fields and fraction-free elimination over rings. Other algorithm names
are rejected explicitly. Elimination currently materializes coefficient rows,
including for sparse input; there is no separate sparse elimination kernel.

`rank`, `nullspace`, `nullspace_from_rref`, `rowspace` and `columnspace` are
available. Nullspace basis vectors are rows; the default need not normalize
or remove common factors. `divide_last=True` requires a field. Row and column
spaces require fields and select independent original rows or columns. A
different valid basis or common denominator may be returned by SymPy.

`inv()` requires a field. `inv_den()` and `solve_den(b)` return a numerator
matrix and scalar denominator in the original domain, satisfying
`A * numerator == b * denominator`. Solving accepts square and consistent
overdetermined systems with a unique solution; singular, inconsistent and
underdetermined systems raise `DMNonInvertibleMatrixError`. The optional
`method='charpoly'` requires a square matrix and uses its adjugate.

`det()` uses Bareiss elimination. `charpoly()` returns descending coefficients
using the division-free Samuelson-Berkowitz recurrence, including over small
finite fields. `adj_det()` returns the actual adjugate and determinant and
supports singular matrices. Empty square matrices have determinant one and
characteristic polynomial `[one]`.

Numerical fields retain their precision through the algorithms. Rank and
pivot selection use the domain's zero predicate; they do not apply a numerical
conditioning threshold. Numerical identities should be checked with the
domain's `almosteq`, allowing for finite precision. Formal expression domains
use their existing simplification and zero conventions and do not track
conditions on symbolic pivots. Factorizations, normal forms, optimized sparse
elimination and the remaining full SymPy `DomainMatrix` API are not implemented.
