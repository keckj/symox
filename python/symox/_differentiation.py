"""Evaluated and held scalar differentiation, also used entrywise by matrices."""
from . import _core

_native_diff=_core.diff
_native_scalar_diff=_core.Expr.diff
_native_matrix_diff=_core.Matrix.diff


def _options(evaluate, options):
    unknown=set(options)-{'simplify'}
    if unknown:raise TypeError('Unknown differentiation options: '+', '.join(sorted(unknown)))
    if evaluate and options.get('simplify',True) is not True:
        raise NotImplementedError('Evaluated derivatives use native canonical simplification')


def _matrix_variables(expression, symbols):
    if symbols:return symbols
    free=expression.free_symbols
    if len(free)>1:raise ValueError('Supply differentiation variables for a multivariate matrix')
    return tuple(free)


def scalar_diff(self, *symbols, evaluate=True, **options):
    _options(evaluate,options)
    if not evaluate:return _core.Derivative(self,*symbols,evaluate=False)
    return _native_scalar_diff(self,*symbols)


def matrix_diff(self, *symbols, evaluate=True, **options):
    from . import S
    _options(evaluate,options)
    if not evaluate:
        raise NotImplementedError('Held matrix derivatives require lazy array expressions')
    symbols=_matrix_variables(self,symbols)
    if not symbols:return self.applyfunc(lambda entry:S.Zero)
    if not isinstance(self,_core.Matrix):
        return self.applyfunc(lambda entry:_native_scalar_diff(entry,*symbols))
    return _native_matrix_diff(self,*symbols)


def diff(f, *symbols, evaluate=True, **options):
    """Differentiate scalar or matrix expressions, optionally holding evaluation."""
    _options(evaluate,options)
    if isinstance(f,_core.Matrix):
        return matrix_diff(f,*symbols,evaluate=evaluate,**options)
    if getattr(f,'_symox_sparse',False):
        return matrix_diff(f,*symbols,evaluate=evaluate,**options)
    if not evaluate:return _core.Derivative(f,*symbols,evaluate=False)
    return _native_diff(f,*symbols)


def install():
    from ._sparse import _SparseBase
    _core.Expr.diff=scalar_diff
    _core.Matrix.diff=matrix_diff
    _SparseBase.diff=matrix_diff
