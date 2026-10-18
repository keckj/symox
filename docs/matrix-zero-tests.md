# Matrix zero tests and pivot simplification

[Matrices](matrices.md) · [README](../README.md)

`rref`, `echelon_form`, `rank`, `nullspace` and `rank_decomposition` accept an
`iszerofunc` callable. It receives native scalar entries and returns `True`
for zero, `False` for nonzero or `None` when undecided. Native boolean atoms
are accepted as well as Python booleans.
Boolean scalar containers exposing `item()`, including NumPy boolean scalars,
are accepted without converting numeric scalar items into boolean answers.

```python
from symox import Matrix

matrix = Matrix([[1e-12, 1], [0, 1]])
zero = lambda value: abs(complex(value)) < 1e-9
reduced, pivots = matrix.rref(iszerofunc=zero)
assert reduced == Matrix([[0, 1], [0, 0]])
assert pivots == (1,)
assert matrix.rank(iszerofunc=zero) == 1
```

The callable controls pivot selection and row elimination. Numerical columns
use partial pivoting, testing the largest absolute entry with the supplied
predicate. Arbitrary-precision entries are compared directly, including values
too small to be represented by `f64`. A tolerance therefore changes the rank
and subspaces according to that predicate.

These methods, plus `rowspace` and `columnspace`, also accept a `simplify`
callable. It is used when deciding an unknown pivot. Default boolean flags use
the native simplifier. If a pivot remains undecided, elimination assumes it is
nonzero and gives a generic result, with the corresponding denominators.

```python
from symox import Symbol

x = Symbol('x')
matrix = Matrix([[x, 1], [0, 1]])
simplify = lambda value: value.subs(x, 0)
assert matrix.rank(simplify=simplify) == 1
```

`rref` honors `pivots` and `normalize_last`; `echelon_form` honors
`with_pivots`. Dense and sparse adapters preserve the input and the matrix
mutability convention. Callback errors stop elimination and propagate the
original Python exception. All matrix arithmetic remains native; callbacks do
not cause a SymPy import.

## LU factorization and solves

`LUdecomposition` and `LUdecomposition_Simple` accept `iszerofunc`, `simpfunc`
and `rankcheck`. `LUsolve` accepts the same zero predicate for pivot selection
and consistency checks. A pivot proved zero by `simpfunc` is discarded before
division. An undecided pivot can still be used as a generic nonzero value.

```python
matrix = Matrix([[x, 1], [1, 0]])
zero = lambda value: value.is_zero
lower, upper, swaps = matrix.LUdecomposition(iszerofunc=zero)
assert lower*upper == matrix.permute_rows(swaps)
rhs = Matrix([1, 2])
assert matrix*matrix.LUsolve(rhs, iszerofunc=zero) == rhs
```

Default zero testing also respects symbols declared zero or nonzero, positive
powers of known-zero bases, and sums whose nonzero contributions can be counted
exactly. A lone symbol declared zero has matrix rank zero.

## Determinants and inverses

`det(method='bareiss', iszerofunc=...)` and `det(method='lu', iszerofunc=...)`
use the predicate in pivot selection, including for rational matrices that
would otherwise use an exact-domain fast path. Closed formulas for matrices
of order at most three, Berkowitz, Bird, Laplace and domain elimination compute
the algebraic determinant without consulting that predicate.

`inv`, `inverse_GE`, `inverse_LU`, `inverse_ADJ`, `inverse_CH`, `inverse_LDL`
and `inverse_QR` accept `iszerofunc`. GE and LU use it during elimination;
the other decomposition methods use it when their determinant verification is
undecided. Explicit `DM` and `DMNC` inversion uses the formal coefficient field
and does not consult a custom zero predicate. A supplied predicate makes the
default inversion use GE.

```python
matrix = Matrix([[2, 1], [1, 3]])
assert matrix*matrix.inv(iszerofunc=lambda value: value.is_zero) == Matrix.eye(2)
```

`try_block_diag=True` splits consecutive diagonal blocks before applying the
selected inversion method and predicate to each block. Every inverse preserves
the input and its mutability convention. Singular matrices raise
`NonInvertibleMatrixError`, and callback exceptions propagate unchanged.

LU inversion applies its predicate throughout factorization and solving, even
when the matrix has no symbols. This is stronger than SymPy 1.14's LU inverse,
which only uses that option in its preliminary symbolic singularity check.
As with row reduction, tolerance predicates define a different zero convention
from the exact algebraic determinant.

## Characteristic coefficients

`charpoly(x=None, simplify=None)` accepts a coefficient transformation callable.
It receives each coefficient in descending degree order, using native domain
arithmetic. General expression coefficients reach the callable before default
symbolic simplification. The result infers its coefficient domain again and
retains the generator's identity or the fresh name chosen to avoid capture.

```python
x = Symbol('x')
matrix = Matrix([[x, 1], [2, x]])
polynomial = matrix.charpoly(simplify=lambda value: value.subs(x, 0))
assert polynomial.all_coeffs() == [1, 0, -2]
```

Leading zero coefficients are removed; the zero polynomial has degree `-oo`.
Pure polynomial equality and hashing ignore generator names while retaining
the coefficient domain. A polynomial over `ZZ` and one over `QQ` are distinct
even when their coefficient lists coincide.
