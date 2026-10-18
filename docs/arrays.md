# Multidimensional arrays and tensor operations

[README](../README.md)

`Array` is an immutable dense array. Mutable and immutable dense and sparse
variants are exported at the package root and from `symox.tensor.array`.
Python containers own the shape and storage; entries use the native scalar
expression engine for arithmetic and calculus. Calculations do not import SymPy.

```python
from symox import Array, MutableSparseNDimArray, symbols, tensorproduct
from symox import tensorcontraction, tensordiagonal, permutedims, diff, shape

x, y = symbols('x y')
A = Array([[x, y], [1, 2]])
assert shape(A) == (2, 2)
assert A[:, 0].tolist() == [x, 1]
assert tensorcontraction(A, (0, 1)) == x + 2
assert tensordiagonal(A, (0, 1)).tolist() == [x, 2]
assert permutedims(A, (1, 0)) == A.T
assert tensorproduct(A, Array([1, 2])).shape == (2, 2, 2)
assert diff(x*x*y, Array([x, y])).tolist() == [2*x*y, x*x]

S = MutableSparseNDimArray({(1, 2): x}, (1000, 1000))
S[1, 2] = y
assert S[1, 2] == y
```

## Storage and indexing

Constructors accept rectangular nested iterables, a flat iterable with an
explicit shape, arrays, matrices and NumPy arrays. Sparse dictionaries require
a shape and accept flat or coordinate keys. Inputs are copied; mutable slice
assignment validates the complete right hand side before writing and handles
overlapping slices. Broadcasting is not implemented.

Negative indices, partial indices and Python slice steps are supported. Sparse
slices operate on stored entries. Symbolic indexing produces a finite native
`Piecewise`, including negative index normalization and integer/range guards.
Arrays also provide `reshape`, `tolist`, `applyfunc`, `subs`, `xreplace`,
`evalf`, `simplify`, `diff`, `conjugate`, `transpose`, `adjoint`, `tomatrix`,
mutable/immutable conversion, pickling and NumPy conversion. Transpose and
adjoint require rank two; use `permutedims` for other ranks.

Rank-zero and zero-length axes are supported. `len(A)` counts scalar positions;
iteration yields first-axis subarrays, following the SymPy interface. Equality
compares shape and entries across storage variants. Equal immutable arrays
have equal hashes; mutable arrays are unhashable. Explicit `to_sympy()` is an
optional conversion and is the only array operation that imports SymPy.

## Tensor operations and differentiation

`tensorproduct` concatenates input axes. Array multiplication and division
accept scalars; multiplying arrays requires an explicit tensor product.
`tensorcontraction` sums diagonals over disjoint groups of equal-length axes;
`tensordiagonal` retains one diagonal axis per group. `permutedims` accepts an
axis permutation or old/new index labels.

`derive_by_array(expression, variables)` places variable axes before expression
axes, matching SymPy's array convention. It accepts scalar, array and finite
matrix inputs. `diff` and scalar/matrix `.diff(...)` dispatch to this operation
when the variable is an array or matrix. Ordinary scalar derivatives remain
entrywise. Gradients and Hessians can be formed by repeated differentiation.

```python
from symox import MatrixSymbol, Trace, Determinant

B = MatrixSymbol('B', 2, 2)
assert B.diff(B).shape == (2, 2, 2, 2)
assert B.diff(B)[0, 1, 0, 1] == 1
assert Trace(B).diff(B).tolist() == [[1, 0], [0, 1]]
assert Determinant(B).diff(B)[0, 0] == B[1, 1]
```

Finite matrix Jacobians are computed from explicit native scalar entries.
They cover products, transposes, inverses, traces and determinants. Symbolic
dimensions and unevaluated array-variable differentiation are unsupported.
Scalar derivatives that cannot be evaluated can retain native `Derivative`
nodes, including derivatives of conjugates without suitable assumptions.

Public `simplify`, `expand`, `factor`, `trigsimp`, `powsimp`, `ratsimp`,
`radsimp`, `cancel`, `together`, `conjugate` and `N` dispatch entrywise.
`shape` supports registration for additional Python types through
`shape.register`; scalar inputs raise `symox.tensor.functions.NoShapeError`.

## Limits

Arrays require concrete nonnegative dimensions and rank at most 32. Dense
arrays and explicit materializations are limited to 100,000 positions or
nested containers; sparse arrays are limited to 100,000 stored entries.
Symbolic indexing is limited to 4,096 relevant cases. Tensor products and
Jacobians enforce their intermediate storage bounds before materializing.
Sparse tensors can have much larger shapes when operations preserve sparsity.
Nonfinite products respect native scalar rules, including `0*oo = nan`.

Lazy symbolic array-expression nodes, general indexed tensor algebra,
array-specific printers and code generation are not implemented by this layer.
