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
matrix inputs. `diff` and array/matrix `.diff(...)` dispatch to this operation
when the variable is an array or explicit matrix; a scalar differentiated by an
explicit matrix gives an immutable matrix of the same shape, as in SymPy.
Ordinary scalar derivatives remain entrywise. Gradients and Hessians can be
formed by repeated differentiation.

```python
from symox import MatrixSymbol, Trace, Determinant, Identity

B = MatrixSymbol('B', 2, 2)
assert B.diff(B).shape == (2, 2, 2, 2)
assert B.diff(B)[0, 1, 0, 1] == 1
assert Trace(B).diff(B) == Identity(2)
assert Determinant(B).diff(B)[0, 0] == B[1, 1]
```

Differentiation of a matrix expression by a `MatrixSymbol` uses SymPy's
matrix-line algorithm (see [matrix expressions](matrix-expressions.md#scalar-calculus)):
`Trace(B).diff(B)` is a matrix expression and `B.diff(B)` the rank-four array
expression `PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))`, for symbolic
dimensions too. Expressions outside that algorithm, such as determinants, use
a finite Jacobian of explicit native scalar entries, which requires concrete
dimensions. Symbolic derivative orders return a held `ArrayDerivative`.
Scalar derivatives that cannot be evaluated can retain native `Derivative`
nodes, including derivatives of conjugates without suitable assumptions.

Public `simplify`, `expand`, `factor`, `trigsimp`, `powsimp`, `ratsimp`,
`radsimp`, `cancel`, `together`, `conjugate` and `N` dispatch entrywise.
`shape` supports registration for additional Python types through
`shape.register`; scalar inputs raise `symox.tensor.functions.NoShapeError`.

## Deferred comprehensions and element kinds

`ArrayComprehension(expr, (index, lower, upper), ...)` describes a rectangular
array with inclusive integer limits. `shape`, `rank()`, `variables`,
`bound_symbols`, `function`, `limits` and `free_symbols` remain available with
symbolic bounds. `doit()` expands concrete bounds into an immutable dense
array; symbolic bounds retain the comprehension. `tolist()` and `tomatrix()`
require concrete bounds, with matrix conversion also requiring rank two.

```python
from symox import ArrayComprehension, ArrayComprehensionMap

i, j, n = symbols('i j n')
table = ArrayComprehension(10*i+j, (i, 1, n), (j, 1, 2))
assert table.subs(n, 2).tolist() == [[11, 12], [21, 22]]
assert ArrayComprehensionMap(lambda k: k*k, (i, 1, 3)).tolist() == [1, 4, 9]
```

Substitution preserves bound indices and renames them to avoid capturing
replacement symbols; literal `xreplace` can replace the indices themselves.
`ArrayComprehensionMap` accepts a Python lambda with zero or one positional
argument. A one-argument lambda receives the product of the index values,
following SymPy's multidimensional convention. Symbolic comprehensions defer
the callback until concrete expansion. The same rank and materialization
limits apply as for ordinary arrays.

`ArrayKind(element_kind)` is interned and records `NumberKind`, `BooleanKind`
or `UndefinedKind` element metadata. Dense and sparse arrays expose `.kind`;
mixed element kinds and empty arrays have `ArrayKind(UndefinedKind)`. Sparse
implicit zeros contribute `NumberKind` without materializing the full array.
This includes sparse metadata that SymPy 1.14 leaves as `UndefinedKind`.

## Array expressions

`symox.tensor.array.expressions` (also importable through SymPy's
`array_expressions`, `from_*` and `conv_*` module paths) provides lazy shaped
nodes: `ArraySymbol`, `ArrayElement`, `ZeroArray`, `OneArray`,
`ArrayTensorProduct`, `ArrayAdd`, `PermuteDims`, `ArrayContraction`,
`ArrayDiagonal`, `Reshape` and `ArrayElementwiseApplyFunc`. Shapes can be
symbolic; indexing evaluates one entry, and `as_explicit()` materializes
concrete shapes within the dense array limit.

```python
from symox.tensor.array.expressions import (
    convert_matrix_to_array, convert_array_to_matrix,
)

C, D = MatrixSymbol('C', n, n), MatrixSymbol('D', n, n)
array = convert_matrix_to_array(C*D)        # ArrayContraction(ArrayTensorProduct(C, D), (1, 2))
assert convert_array_to_matrix(array) == C*D
```

SymPy's canonical forms (`doit`, `nest_permutation` and the canonicalizing
constructors), `convert_matrix_to_array`, `convert_array_to_matrix` with its
passes and `_EditArrayContraction`, `array_derive`, `convert_indexed_to_array`
and `convert_array_to_indexed` run in Rust (`crates/symox-core/src/arrayexpr*.rs`).
Python encodes the expression tree and rebuilds the result. When the native
`convert_array_to_matrix` declines an input it does not model, the unmodified
SymPy 1.14 algorithm vendored in `symox._arrayexpr` runs on a structural
translation; the matrix derivative helper `_chain_array_to_matrix` is still
Python. Remaining differences from SymPy's own array-expression tests are
tracked in [the parity gaps](algorithmic-parity.md).

Dense and sparse arrays follow SymPy's subclass protocol: user subclasses keep
their type through arithmetic, slicing, `applyfunc`, `diff` and the tensor
operations. `latex` prints arrays as SymPy does, and `lambdify` evaluates
immutable arrays to NumPy arrays with the NumPy and SciPy modules; mutable
arrays and the `math`/`mpmath` modules are rejected, as in SymPy.

## Limits

Arrays require concrete nonnegative dimensions and rank at most 32. Dense
arrays and explicit materializations are limited to 100,000 positions or
nested containers; sparse arrays are limited to 100,000 stored entries.
Symbolic indexing is limited to 4,096 relevant cases. Tensor products and
Jacobians enforce their intermediate storage bounds before materializing.
Sparse tensors can have much larger shapes when operations preserve sparsity.
Nonfinite products respect native scalar rules, including `0*oo = nan`.

Typed abstract indices, Einstein contractions and slot symmetries are
provided separately by [indexed tensor algebra](indexed-tensors.md).
