# Permutation matrices

[README](../README.md) · [Matrices](matrices.md) · [Matrix expressions](matrix-expressions.md)

`PermutationMatrix` and `MatrixPermute` are native matrix expression nodes.
They store finite permutations without allocating a dense matrix. Arithmetic,
entry evaluation, scalar differentiation, substitutions and serialization run
without importing SymPy.

```python
from symox import Matrix, MatrixPermute, PermutationMatrix, Identity
from symox.combinatorics import Permutation

p = Permutation([1, 2, 0])
P = PermutationMatrix(p)
assert P * P.T == Identity(3)
assert P.det() == 1

A = Matrix([[1, 2], [3, 4], [5, 6]])
assert MatrixPermute(A, p, axis=0).as_explicit() == Matrix([[3, 4], [5, 6], [1, 2]])
```

The array gives the image of each index. Matrix entry `(i, p(i))` is one;
every other entry is zero. Left multiplication by `P` permutes rows;
right multiplication by `P.T` applies the same permutation to columns.
Permutation multiplication applies its left factor first, matching SymPy's
convention: `(p*q)(i) == q(p(i))`, so `P(p)*P(q) == P(p*q)`.

`PermutationMatrix` supports native composition, integer powers within the
native signed 64-bit range, inverse, transpose, adjoint, conjugation,
determinant and trace. Identity permutations remain held until `doit()`.
Its permutation is an atomic constructor argument; its images are not free
scalar symbols or independently substituted expression children.

`MatrixPermute(A, p, axis)` holds row (`axis=0`) or column (`axis=1`)
permutation. It accepts a permutation matrix too. The relevant axis must have
a finite integer size. Shorter permutations gain fixed points; longer ones
can shrink only by removing trailing fixed points. Scalar differentiation
preserves the permutation. Transpose and adjoint switch axes. `doit()` combines
nested permutations on the same axis and simplifies identity, zero and
all-ones matrices.

Concrete indices support negative indexing. Symbolic permutation application
uses finite native `Piecewise` expressions, so subsequent substitutions can
evaluate entries exactly. SymPy may represent these applications differently.
`rewrite(MatMul)` expresses a held permutation as matrix multiplication;
`PermutationMatrix.rewrite(BlockDiagMatrix)` splits contiguous invariant
intervals into diagonal permutation blocks.

The native matrix DAG codec preserves permutation payloads, held nodes and
shared parents. Decoding validates bijections, axis values and sizes. Native
and SymPy textual representations of permutation matrix expressions can also
be parsed. Text parsing bounds permutation indices below 100,000; the binary
codec avoids this text representation. Explicit matrix conversion retains the
usual native allocation limits.

The Python `Permutation` class accepts array images, cycle lists and the
positional cycle syntax, and carries the SymPy `Permutation` methods; see
[combinatorics](combinatorics.md). Native textual parsing additionally
understands the cycle syntax emitted by SymPy's `srepr` for these matrix
expressions.
