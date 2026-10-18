"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'ImmutableDenseMatrix ImmutableMatrix ImmutableSparseMatrix')
from symox import ImmutableDenseMatrix,ImmutableMatrix,ImmutableSparseMatrix


class _ImmutableRepMeta(type):
    def __instancecheck__(cls,obj):return isinstance(obj,(ImmutableDenseMatrix,ImmutableSparseMatrix))
class ImmutableRepMatrix(metaclass=_ImmutableRepMeta):
    def __new__(cls,*args,**kwargs):return ImmutableMatrix(*args,**kwargs)


def sympify_matrix(arg):return ImmutableMatrix(arg)
def sympify_mpmath_matrix(arg):
    from symox import Float
    return ImmutableMatrix(arg.rows,arg.cols,[Float(str(arg[i,j])) for i in range(arg.rows) for j in range(arg.cols)])

__all__ += ['ImmutableRepMatrix','sympify_matrix','sympify_mpmath_matrix']
