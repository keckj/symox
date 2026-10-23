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


def _matrix_variables(expression, symbols):
    if symbols:return symbols
    free=expression.free_symbols
    if len(free)>1:raise ValueError('Supply differentiation variables for a multivariate matrix')
    return tuple(free)


def _ambiguous_variable(error,f):
    # SymPy's Derivative text names the expression.
    if not symbols_missing(error):raise error
    from ._public_tools import filldedent
    raise ValueError(filldedent('''
        Since there is more than one variable in the
        expression, the variable(s) of differentiation
        must be supplied to differentiate %s''' % f)) from None


def symbols_missing(error):
    return isinstance(error,ValueError) and str(error).startswith('Since there is more than one variable in the expression, the variable(s) of differentiation must be supplied')


def _derivative_order(expression, symbols):
    # Number of differentiations requested, None when it is not a known
    # integer (symbolic orders, or forms the held Derivative rejects).
    if not symbols:
        return 1
    try:
        held=_core.Derivative(expression,*symbols,evaluate=False)
    except Exception:
        return None
    if not getattr(held,'is_Derivative',False):
        return None
    total=0
    for variable,count in held.variable_count:
        if not (isinstance(count,int) or getattr(count,'is_Integer',False)):
            return None
        total+=int(count)
    return total


def _finish_multiple(result, expression, symbols):
    """The last step of the reference Derivative: after more than one
    differentiation, ``factor_terms(signsimp(expr))``."""
    if not isinstance(result,_core.Expr) or result.is_zero is True:
        return result
    order=_derivative_order(expression,symbols)
    if order is None or order <= 1:
        return result
    from .simplify import signsimp
    from . import factor_terms
    return factor_terms(signsimp(result))


def scalar_diff(self, *symbols, evaluate=True, **options):
    if _options(evaluate,options):return unsimplified_diff(self,symbols)
    if not evaluate:return _core.Derivative(self,*symbols,evaluate=False)
    try:result=_native_scalar_diff(self,*symbols)
    except ValueError as error:_ambiguous_variable(error,self)
    return _finish_multiple(result,self,symbols)


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
    try:result=_native_diff(f,*symbols)
    except ValueError as error:_ambiguous_variable(error,f)
    return _finish_multiple(result,f,symbols)


def _derivative_doit(native):
    # Derivative.doit re-creates Derivative(expr, *variable_count,
    # evaluate=True) in the reference, which ends with the multiple
    # derivative step.
    def doit(self, *args, **hints):
        deep=hints.get('deep',True) is not False
        if (deep and not getattr(self,'is_Derivative',False) and isinstance(self,_core.Expr)
                and self.has(_core.Derivative)):
            from . import Derivative
            held=[d for d in self.atoms(Derivative) if d.is_Derivative]
            replaced={d:doit(d,*args,**hints) for d in held}
            replaced={d:r for d,r in replaced.items() if r is not d and r != d}
            if replaced:
                self=self.xreplace(replaced)
        result=native(self,*args,**hints)
        if getattr(self,'is_Derivative',False) and deep:
            try:counts=self.variable_count
            except Exception:return result
            if result is not self and not getattr(result,'is_Derivative',False):
                return _finish_multiple(result,self.expr,tuple(counts))
        return result
    doit.__name__=getattr(native,'__name__','doit')
    doit.__doc__=getattr(native,'__doc__',None)
    doit._symox_native=getattr(native,'_symox_native',native)
    return doit


def install():
    from ._sparse import _SparseBase
    _core.Expr.doit=_derivative_doit(_core.Expr.doit)
    _core.Expr.diff=scalar_diff
    _core.Matrix.diff=matrix_diff
    _SparseBase.diff=matrix_diff
