"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'DeferredVector MatrixBase MatrixError NonSquareMatrixError ShapeError')


def a2idx(value,n=None):
    import operator
    index=operator.index(value)
    if n is not None:
        if index<0:index+=n
        if not 0<=index<n:raise IndexError('matrix index out of range')
    return index


def classof(A,B):
    if not hasattr(A,'shape') or not hasattr(B,'shape'):raise TypeError('matrix classes required')
    from symox import ImmutableMatrix,Matrix,ImmutableSparseMatrix,SparseMatrix
    immutable=isinstance(A,(ImmutableMatrix,ImmutableSparseMatrix)) or isinstance(B,(ImmutableMatrix,ImmutableSparseMatrix))
    sparse=isinstance(A,(SparseMatrix,ImmutableSparseMatrix)) and isinstance(B,(SparseMatrix,ImmutableSparseMatrix))
    return ImmutableSparseMatrix if immutable and sparse else ImmutableMatrix if immutable else SparseMatrix if sparse else Matrix

__all__ += ['a2idx','classof']
