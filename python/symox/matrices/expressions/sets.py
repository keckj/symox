"""Shape-checked matrices over a scalar set."""
from ..._set_extensions import _StructuralSet
from ... import Set, sympify, And, Eq, true, false

class MatrixSet(_StructuralSet):
    def __init__(self,n,m,set):
        self.n,self.m=sympify(n),sympify(m)
        if any(d.is_integer is False or d.is_nonnegative is False for d in (self.n,self.m)):
            raise ValueError('matrix dimensions must be nonnegative integers')
        if not isinstance(set,Set):raise TypeError('matrix entries require a scalar set')
        self.set=set
    args=property(lambda self:(self.n,self.m,self.set))
    shape=property(lambda self:(self.n,self.m))
    is_empty=False
    is_finite_set=False
    is_iterable=False
    @property
    def kind(self):
        from ..kind import MatrixKind
        from ...sets.sets import SetKind
        return SetKind(MatrixKind())
    def contains(self,matrix):
        if not hasattr(matrix,'shape'):raise TypeError('matrix membership requires a matrix')
        conditions=[Eq(a,b) for a,b in zip(self.shape,matrix.shape)]
        if false in conditions:return false
        if hasattr(matrix,'as_explicit'):
            try:matrix=matrix.as_explicit()
            except (TypeError,ValueError,NotImplementedError):return None
        try:return And(*conditions,*(self.set.contains(entry) for entry in matrix))
        except TypeError:return None
    def __contains__(self,matrix):
        answer=self.contains(matrix)
        if answer is None:raise TypeError('matrix membership is undecidable')
        return bool(answer)

__all__=['MatrixSet']
