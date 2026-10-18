# Matrix expression protocols

[README](../README.md) · [Matrices](matrices.md)

The `MatrixExpr` interface exposes all 105 method names measured against
SymPy 1.14. This is API presence, with the supported behaviors described below.
Matrix DAG nodes retain their types, dimensions and ordered matrix factors.
Calculations use native coefficients without importing SymPy.

## Structure and transformations

`atoms`, `find`, `count`, `has`, `has_free`, `has_xfree`, `is_same`, `compare`,
`dummy_eq`, `as_dummy`, `copy`, `replace`, `match` and `matches` traverse matrix
nodes and their scalar children. Structural matrix names are distinct from
scalar symbols. Bound variables in `Lambda`, Integral, Sum and Product are
excluded from free-argument queries. Matrix products match in factor order;
matrix sums allow term permutation. Matching has a 10,000-step search limit.

`fromiter` calls the concrete matrix constructor. `as_terms`,
`as_ordered_terms`, `as_coefficients_dict`, `args_cnc`, coefficient and content
views preserve noncommutative factors. Explicit monomial orders propagate
through scalar coefficients and nested presentation keys while retaining matrix
factor order. `count_ops(visual=True)` exposes operation
labels, including `MATMUL`, subtraction, negation and scalar operations.

Scalar simplifiers, rewriting, refinement and numerical evaluation transform
scalar children. Numerical evaluation preserves matrix indices, dimensions
and integer matrix exponents. `as_real_imag` uses matrix conjugation. Rational
content extraction handles matrix sums directly, including cases where SymPy
1.14 raises while rebuilding a scalar/matrix mixture.

`equals` checks canonical differences, then at most 256 explicit entries. It
returns `None` when it cannot prove a result. `as_poly` returns `None` because
the scalar `Poly` backend cannot represent noncommutative matrix coefficients.
Direct construction of a scalar `Poly` from a matrix raises `TypeError`.

## Scalar calculus

Scalar coefficients in matrix patterns accept `match(..., old=True)` and
`matches(..., old=True)`. The legacy scalar fallback can isolate one linear
wildcard while matrix factor order and structural shape matching are retained.

For concrete dimensions, differentiation with respect to an entire matrix or
array now returns a finite array Jacobian, with variable axes first. This
includes matrix products, transposes, inverses, traces and determinants; see
[arrays and tensors](arrays.md). Symbolic dimensions still require a future
lazy array-expression representation. Scalar calculus retains the matrix DAG.

Series, limits, Taylor terms, leading terms and asymptotic series use native
scalar coefficient kernels. Integration preserves constant ordered matrix
factors and handles FunctionMatrix bodies. More general dependent products
can integrate entrywise when dimensions are concrete, up to 256 entries.
Variable matrix dimensions do not define integration in a fixed vector space.

`collect` groups scalar powers with matrix coefficients. Scalar exponent views,
additive and multiplicative extraction, branch factors, meromorphic tests and
hypergeometric ratios operate on constant matrix monomials. A hypergeometric
matrix combination must have a common scalar consecutive-term ratio.
Unresolved local regularity returns `None`. Leading exponents must be orderable;
leading terms accept finite real and complex approach directions through the
scalar coefficient kernels, with `logx` substitution. `lseries` yields
nonzero terms lazily through order 256 and terminates for polynomials.

Modular inversion uses the explicit adjugate and the inverse of the determinant
modulo the supplied scalar polynomial. Rounding operates on explicit scalar
entries; symbolic matrix entries raise the corresponding scalar error.

## Index sums

```python
from symox import MatrixSymbol, MatrixExpr, Sum, symbols

i, j, k = symbols('i j k', integer=True)
A, B = MatrixSymbol('A', 2, 3), MatrixSymbol('B', 3, 4)
expression = Sum(A[i, j]*B[j, k], (j, 0, 2))
assert MatrixExpr.from_index_summation(expression) == A*B
```

Inference recognizes rectangular matrix chains, transpositions, closed traces,
linear combinations, disconnected scalar traces and parallel free edges as
Hadamard products. Every bound contraction index must occur twice and range
over its complete dimension. General higher-rank array contractions, partial
contractions and index-dependent scalar multipliers remain unsupported.

## Matrix coefficient series

`matrix_expr.fps(x)` returns a `FormalPowerSeries` subclass with exact matrix
coefficients. Products use ordered convolution. Inversion uses the recurrence
for a square matrix series with an invertible constant coefficient; a symbolic
constant matrix inverse remains formal. Finite-center inversion is supported.
Coefficients are available through index 255 and truncation through order 256.
Composition requires a scalar inner series. Matrix monomials must be independent
of the scalar expansion variable; Laurent and fractional series remain outside
the ordinary formal interface. The truncation remainder is represented by
`O((x-x0)**n)*OneMatrix(rows, cols)`, as an entrywise order bound.

`matrix_expr.fourier_series(...)` applies native scalar Fourier kernels to
constant matrix coefficients. Harmonics, exact mean, nonzero truncation,
Lanczos sigma approximation and scalar argument/value transformations retain
matrix shapes. A value shift requires a constant matrix of the same shape.
Known finite harmonic support terminates iteration; other searches retain the
scalar limit of harmonic 4096. Explicit SymPy conversion produces a matrix of
scalar series rather than a SymPy matrix-valued series class.
