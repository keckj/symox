"""Evaluated and held scalar differentiation, also used entrywise by matrices."""
from . import _core

_native_diff=_core.diff
_native_scalar_diff=_core.Expr.diff
_native_matrix_diff=_core.Matrix.diff


def _options(evaluate, options):
    unknown=set(options)-{'simplify'}
    if unknown:raise TypeError('Unknown differentiation options: '+', '.join(sorted(unknown)))
    return bool(evaluate) and not options.get('simplify',True)


def unsimplified_diff(expression, symbols):
    """SymPy's ``Derivative(..., evaluate=True, simplify=False)``.

    The reference differentiates one variable at a time with
    ``_eval_derivative`` and only skips the final
    ``factor_terms(signsimp(expr))`` of multiple derivatives; a native
    first derivative performs the same single step.
    """
    held=_core.Derivative(expression,*symbols,evaluate=False)
    if not getattr(held,'is_Derivative',False):return held
    pairs=held.variable_count
    if any(not (isinstance(count,int) or getattr(count,'is_Integer',False)) for variable,count in pairs):
        raise NotImplementedError('Unsimplified derivatives of symbolic order are not implemented')
    result=expression
    for variable,count in pairs:
        for _ in range(int(count)):
            result=_native_diff(result,variable)
            if result.is_zero is True:return result
    return result


def _derivative_count(expression, symbols):
    held=_core.Derivative(expression,*symbols,evaluate=False)
    if not getattr(held,'is_Derivative',False):return 0
    total=0
    for variable,count in held.variable_count:
        if not (isinstance(count,int) or getattr(count,'is_Integer',False)):return 0
        total+=int(count)
    return total


def _reference_simplify(result, expression, symbols):
    # Derivative.__new__ ends multiple derivatives with
    # factor_terms(signsimp(expr)) unless simplify=False.
    if len(symbols)==1 and not isinstance(symbols[0],(tuple,list)):return result
    if not isinstance(result,_core.Expr) or result.is_Atom:return result
    if _derivative_count(expression,symbols)<=1:return result
    from .core.exprtools import factor_terms
    from .simplify.simplify import signsimp
    return factor_terms(signsimp(result))


def _matrix_variables(expression, symbols):
    if symbols:return symbols
    free=expression.free_symbols
    if len(free)>1:raise ValueError('Supply differentiation variables for a multivariate matrix')
    return tuple(free)


def scalar_diff(self, *symbols, evaluate=True, **options):
    if _options(evaluate,options):return unsimplified_diff(self,symbols)
    if not evaluate:return _core.Derivative(self,*symbols,evaluate=False)
    return _reference_simplify(_native_scalar_diff(self,*symbols),self,symbols)


def matrix_diff(self, *symbols, evaluate=True, **options):
    from . import S
    if _options(evaluate,options):
        raise NotImplementedError('Unsimplified matrix derivatives are not implemented')
    if not evaluate:
        raise NotImplementedError('Held matrix derivatives require lazy array expressions')
    symbols=_matrix_variables(self,symbols)
    if not symbols:return self.applyfunc(lambda entry:S.Zero)
    if not isinstance(self,_core.Matrix):
        return self.applyfunc(lambda entry:_native_scalar_diff(entry,*symbols))
    return _native_matrix_diff(self,*symbols)


def diff(f, *symbols, evaluate=True, **options):
    """Differentiate scalar or matrix expressions, optionally holding evaluation."""
    unsimplified=_options(evaluate,options)
    if isinstance(f,_core.Matrix):
        return matrix_diff(f,*symbols,evaluate=evaluate,**options)
    if getattr(f,'_symox_sparse',False):
        return matrix_diff(f,*symbols,evaluate=evaluate,**options)
    if isinstance(f,_core.Expr) and type(f) is not _core.Expr:
        # The reference diff() calls f.diff: user host classes such as
        # Quaternion differentiate componentwise.
        from ._host import _defines
        if _defines(type(f),'diff'):
            return f.diff(*symbols,evaluate=evaluate,**options)
    if not evaluate:return _core.Derivative(f,*symbols,evaluate=False)
    if unsimplified:return unsimplified_diff(f,symbols)
    # Ordered operator polynomials (physics) differentiate themselves.
    hook=getattr(f,'_eval_symox_diff',None)
    if hook is not None:return hook(*symbols)
    return _reference_simplify(_native_diff(f,*symbols),f,symbols)


def install():
    from ._sparse import _SparseBase
    _core.Expr.diff=scalar_diff
    _core.Matrix.diff=matrix_diff
    _SparseBase.diff=matrix_diff
