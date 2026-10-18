"""Held statistical operators with executable exact counterparts."""
from functools import wraps
import symox as S
from ._base import scalar

HELD = {}


class StatisticalExpression:
    operation = None
    def __init__(self,*args,**kwargs):
        self.inputs = args;self.options = kwargs
        native = [S.MatrixExpr(arg) if isinstance(arg,S.Matrix) else scalar(arg) if arg is not None else S.true for arg in args]
        self._expr = S._core._apply(type(self).__name__,*native,evaluate=False)
        HELD[self._expr] = self
    @property
    def args(self):return self.inputs
    @property
    def func(self):return type(self)
    @property
    def free_symbols(self):return self._expr.free_symbols
    def _symox_(self):return self._expr
    def __repr__(self):return str(self._expr)
    def __hash__(self):return hash(self._expr)
    def __eq__(self,other):return isinstance(other,StatisticalExpression) and self._expr == other._expr
    def doit(self,**hints):
        from . import _queries
        return getattr(_queries,self.operation)(*self.inputs,**self.options)
    def evaluate_integral(self):return self.doit()
    def subs(self,*args,**kwargs):
        inputs = tuple((arg if isinstance(arg,S.Matrix) else scalar(arg)).subs(*args,**kwargs) if arg is not None else None for arg in self.inputs)
        return type(self)(*inputs,**self.options)
    def expand(self,**kwargs):return self
    def __add__(self,other):return self._expr+other
    def __radd__(self,other):return other+self._expr
    def __sub__(self,other):return self._expr-other
    def __rsub__(self,other):return other-self._expr
    def __mul__(self,other):return self._expr*other
    def __rmul__(self,other):return other*self._expr
    def __truediv__(self,other):return self._expr/other
    def __rtruediv__(self,other):return other/self._expr
    def __pow__(self,other):return self._expr**other
    def __neg__(self):return -self._expr


class Probability(StatisticalExpression):
    operation = 'P'
    def __init__(self,prob,condition=None,**kwargs):super().__init__(prob,condition,**kwargs)


class Expectation(StatisticalExpression):
    operation = 'E'
    def __init__(self,expr,condition=None,**kwargs):super().__init__(expr,condition,**kwargs)
    def expand(self,**kwargs):
        expr,condition = self.inputs;expr = scalar(expr)
        if expr.is_Add:return S.Add(*(Expectation(part,condition)._expr for part in expr.args))
        return self


class Variance(StatisticalExpression):
    operation = 'variance'
    def __init__(self,arg,condition=None,**kwargs):super().__init__(arg,condition,**kwargs)
    def rewrite(self,target,**kwargs):
        if target is Expectation:
            X,condition = self.inputs
            return Expectation(scalar(X)**2,condition)-Expectation(X,condition)**2
        return self


class Covariance(StatisticalExpression):
    operation = 'covariance'
    def __init__(self,arg1,arg2,condition=None,**kwargs):super().__init__(arg1,arg2,condition,**kwargs)
    def rewrite(self,target,**kwargs):
        if target is Expectation:
            X,Y,condition = self.inputs
            return Expectation(scalar(X)*scalar(Y),condition)-Expectation(X,condition)*Expectation(Y,condition)
        return self


class Moment(StatisticalExpression):
    operation = 'moment'
    def __init__(self,X,n,c=0,condition=None,**kwargs):super().__init__(X,n,c,condition,**kwargs)


class CentralMoment(StatisticalExpression):
    operation = 'cmoment'
    def __init__(self,X,n,condition=None,**kwargs):super().__init__(X,n,condition,**kwargs)


def _matrix_value(value):
    from ._matrix import MATRIX_REGISTRY
    try:return MATRIX_REGISTRY.get(value,value)
    except TypeError:return value


class ExpectationMatrix(StatisticalExpression):
    def __init__(self,expr,condition=None):super().__init__(expr,condition)
    def doit(self,**hints):
        from ._queries import E
        expr,condition = self.inputs
        return E(_matrix_value(expr),condition)


class VarianceMatrix(StatisticalExpression):
    def __init__(self,arg,condition=None):super().__init__(arg,condition)
    @property
    def shape(self):
        value = _matrix_value(self.inputs[0]);dimension = value.rows if value.cols == 1 else value.cols
        return dimension,dimension
    def doit(self,**hints):
        expr,condition = self.inputs;expr = _matrix_value(expr)
        if getattr(expr,'_stats_matrix',False):return expr.variance()
        from ._queries import covariance
        if not isinstance(expr,S.Matrix):raise TypeError('matrix variance requires a random matrix or an explicit matrix of random scalars')
        if expr.cols == 1:
            return S.Matrix([[covariance(expr[i,0],expr[j,0],condition) for j in range(expr.rows)] for i in range(expr.rows)])
        return S.Matrix([[S.Add(*(covariance(expr[i,j],expr[i,k],condition) for i in range(expr.rows))) for k in range(expr.cols)] for j in range(expr.cols)])


class CrossCovarianceMatrix(StatisticalExpression):
    def __init__(self,arg1,arg2,condition=None):super().__init__(arg1,arg2,condition)
    @property
    def shape(self):
        left,right = map(_matrix_value,self.inputs[:2])
        return (left.rows,right.rows) if left.cols == right.cols == 1 else (left.cols,right.cols)
    def doit(self,**hints):
        left,right,condition = self.inputs;left,right = _matrix_value(left),_matrix_value(right)
        column_vectors = left.cols == right.cols == 1
        if not column_vectors and left.rows != right.rows:raise ValueError('cross covariance requires equal row counts')
        if left is right:return VarianceMatrix(left,condition).doit()
        if getattr(left,'_stats_matrix',False) and getattr(right,'_stats_matrix',False):return S.zeros(*self.shape)
        from ._queries import covariance
        if column_vectors:return S.Matrix([[covariance(left[i,0],right[j,0],condition) for j in range(right.rows)] for i in range(left.rows)])
        return S.Matrix([[S.Add(*(covariance(left[i,j],right[i,k],condition) for i in range(left.rows))) for k in range(right.cols)] for j in range(left.cols)])


def install():
    if getattr(S.Expr,'_stats_doit_installed',False):return
    original = S.Expr.doit
    @wraps(original)
    def doit(self,*args,**kwargs):
        if self in HELD:return HELD[self].doit(**kwargs)
        if kwargs.get('deep',True):
            replacements,seen,stack = {},set(),[self]
            while stack:
                node = stack.pop()
                if not isinstance(node,S.Expr) or node in seen:continue
                seen.add(node)
                if node in HELD:replacements[node] = HELD[node].doit(**kwargs)
                else:stack.extend(node.args)
            if replacements:return original(self.xreplace(replacements),*args,**kwargs)
        return original(self,*args,**kwargs)
    S.Expr.doit = doit;S.Expr._stats_doit_installed = True


__all__ = ['Probability','Expectation','Variance','Covariance','Moment','CentralMoment',
 'ExpectationMatrix','VarianceMatrix','CrossCovarianceMatrix']
