"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'DeferredVector MatrixBase MatrixError NonSquareMatrixError ShapeError')


def a2idx(j,n=None):
    import operator
    index=operator.index(j)
    if n is not None:
        if index<0:index+=n
        if not 0<=index<n:raise IndexError('matrix index out of range')
    return index


def classof(A,B):
    if not hasattr(A,'shape') or not hasattr(B,'shape'):raise TypeError('matrix classes required')
    from symox import ImmutableMatrix,DenseMatrix as Matrix,ImmutableSparseMatrix,SparseMatrix
    if isinstance(A,Matrix) and isinstance(B,Matrix) and not A._immutable and not B._immutable:
        return type(A) if A._class_priority>B._class_priority else type(B)
    immutable=isinstance(A,(ImmutableMatrix,ImmutableSparseMatrix)) or isinstance(B,(ImmutableMatrix,ImmutableSparseMatrix))
    sparse=isinstance(A,(SparseMatrix,ImmutableSparseMatrix)) and isinstance(B,(SparseMatrix,ImmutableSparseMatrix))
    return ImmutableSparseMatrix if immutable and sparse else ImmutableMatrix if immutable else SparseMatrix if sparse else Matrix

__all__ += ['a2idx','classof']
