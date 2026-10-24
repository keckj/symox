"""Constructors and rewrites for native finite permutation matrix nodes."""
import operator

from . import _core
from ._matrixexpr import MatrixExpr, _register_classes


class PermutationMatrix(MatrixExpr):
    __slots__ = ()
    _kinds = ('PermutationMatrix',)

    def __new__(cls, perm):
        from .combinatorics import Permutation
        if not isinstance(perm, Permutation):
            raise ValueError('Expected a native Permutation')
        return _core._matrixexpr_node('PermutationMatrix', (perm,))


class MatrixPermute(MatrixExpr):
    __slots__ = ()
    _kinds = ('MatrixPermute',)

    def __new__(cls, mat, perm, axis=0):
        from .combinatorics import Permutation
        if not (isinstance(mat, (_core.Matrix, _core._MatrixExpr)) or
                getattr(mat, '_symox_sparse', False) or getattr(mat, 'is_Matrix', False)):
            raise ValueError('Expected a matrix')
        if isinstance(perm, PermutationMatrix):
            perm = perm.args[0]
        if not isinstance(perm, Permutation):
            raise ValueError('Expected a Permutation or PermutationMatrix')
        try:
            axis = operator.index(axis)
        except TypeError:
            raise ValueError('Permutation axis must be zero or one') from None
        if axis not in (0, 1):
            raise ValueError('Permutation axis must be zero or one')
        return _core._matrixexpr_node('MatrixPermute', (mat, perm, axis))


def rewrite_permutation(value, target, *, deep=True):
    """Return a matrix-specific rewrite, or None to use scalar rewrites."""
    from . import BlockDiagMatrix, MatMul
    from .combinatorics import Permutation
    name = target if isinstance(target, str) else getattr(target, '__name__', None)
    if isinstance(value, MatrixPermute) and name == 'MatMul':
        matrix, permutation, axis = value.args
        if deep:
            matrix = (matrix.rewrite(target, deep=True) if isinstance(matrix, MatrixExpr)
                      else matrix.applyfunc(lambda entry: entry.rewrite(target, deep=True)))
        if axis == 0:
            return MatMul(PermutationMatrix(permutation), matrix)
        return MatMul(matrix, PermutationMatrix(~permutation))
    if isinstance(value, PermutationMatrix) and name == 'BlockDiagMatrix':
        images = value.args[0].array_form
        blocks, start, end = [], 0, -1
        # A prefix closes exactly when all its images lie in that prefix.
        # Each closed interval gives an invariant diagonal block.
        for i, image in enumerate(images):
            end = max(end, image)
            if end == i:
                blocks.append(PermutationMatrix(Permutation(
                    [j-start for j in images[start:i+1]])))
                start = i+1
        return BlockDiagMatrix(*blocks)
    return None


if hasattr(_core, '_register_matrixexpr_type'):
    _register_classes(PermutationMatrix, MatrixPermute)
