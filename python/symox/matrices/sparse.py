"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'MutableSparseMatrix SparseMatrix')
from symox import MutableSparseMatrix,SparseMatrix


class _SparseRepMeta(type):
    def __instancecheck__(cls,obj):
        from symox import ImmutableSparseMatrix
        return isinstance(obj,(MutableSparseMatrix,ImmutableSparseMatrix))
class SparseRepMatrix(metaclass=_SparseRepMeta):
    def __new__(cls,*args,**kwargs):return SparseMatrix(*args,**kwargs)
__all__ += ['SparseRepMatrix']
