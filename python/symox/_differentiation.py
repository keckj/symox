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
    closed=_core.Expr._integer_order_value(held)
    if closed is not None:return closed
    if any(not (isinstance(count,int) or getattr(count,'is_Integer',False)) for variable,count in pairs):
        raise NotImplementedError('Unsimplified derivatives of symbolic order are not implemented')
    result=expression
    for variable,count in pairs:
        for _ in range(int(count)):
            result=_native_diff(result,variable)
            if result.is_zero is True:return result
    return result


def _is_array_like(value):
    from .tensor.array import NDimArray
    return isinstance(value, (list, tuple, NDimArray, _core.Matrix)) or getattr(value, '_symox_sparse', False)


def _array_variable_spec(symbols):
    """SymPy's ``_derivative_dispatch``/``Derivative.__new__`` parsing of
    derive-by-array variables such as ``[[x, y]]``, ``([x, y], 2)`` or
    ``((x, y),)``. Returns None when no such variable is present."""
    if not any(isinstance(item, (tuple, list)) and item and _is_array_like(item[0]) for item in symbols):
        return None
    from .tensor.array import Array
    spec = []
    for item in symbols:
        if isinstance(item, (tuple, list)) and item and _is_array_like(item[0]):
            if len(item) == 1:
                variable, count = item[0], 1
            else:
                variable, count = item
            spec.append(('array', Array(variable), count))
        else:
            spec.append(('scalar', item, None))
    return spec


def _array_diff(expression, spec):
    """Evaluate ``ArrayDerivative`` one variable group at a time."""
    from .tensor.array import NDimArray, derive_by_array
    from . import sympify
    result = expression
    chunk = []

    def flush(value):
        if not chunk:
            return value
        symbols = tuple(chunk)
        chunk.clear()
        if isinstance(value, NDimArray):
            return value.applyfunc(lambda entry: diff(entry, *symbols)).as_immutable()
        return diff(value, *symbols)

    for kind, variable, count in spec:
        if kind == 'scalar':
            chunk.append(variable)
            continue
        result = flush(result)
        count = sympify(count)
        if not count.is_Integer or count < 0:
            raise ValueError('derivative order must be a nonnegative integer')
        for _ in range(int(count)):
            result = derive_by_array(result, variable)
    return flush(result)


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
    if not isinstance(result,_core.Expr) or result.is_zero is True or _core.Expr._has_deferred_numeric_power.__get__(result):
        return result
    order=_derivative_order(expression,symbols)
    if order is None or order <= 1:
        return result
    from .simplify import signsimp
    from . import factor_terms
    return factor_terms(signsimp(result))


def scalar_diff(self, *symbols, evaluate=True, **options):
    spec=_array_variable_spec(symbols)
    if spec is not None and evaluate and not options:
        return _array_diff(self,spec)
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
    spec=_array_variable_spec(symbols)
    if spec is not None:
        result=_array_diff(self,spec)
        # SymPy's MatrixCalculus.diff: a mutable matrix gives a mutable array.
        if not getattr(self,'_immutable',False) and hasattr(result,'as_mutable'):
            result=result.as_mutable()
        return result
    symbols=_matrix_variables(self,symbols)
    if not symbols:return self.applyfunc(lambda entry:S.Zero)
    if not isinstance(self,_core.Matrix):
        return self.applyfunc(lambda entry:_native_scalar_diff(entry,*symbols))
    result=_native_matrix_diff(self,*symbols)
    # SymPy's MatrixCalculus.diff: the derivative of a mutable matrix by a
    # matrix is a mutable array.
    if (not getattr(self,'_immutable',False) and not isinstance(result,_core.Matrix)
            and hasattr(result,'as_mutable')):
        result=result.as_mutable()
    return result


def diff(f, *symbols, evaluate=True, **options):
    """Differentiate scalar or matrix expressions, optionally holding evaluation."""
    unsimplified=_options(evaluate,options)
    if isinstance(f,_core.Matrix):
        return matrix_diff(f,*symbols,evaluate=evaluate,**options)
    if evaluate and not options and isinstance(f,_core.Expr):
        spec=_array_variable_spec(symbols)
        if spec is not None:
            return _array_diff(f,spec)
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


def _derivative_eval_derivative(self, v):
    """Reference ``Derivative._eval_derivative``."""
    from . import Derivative
    v=_core.sympify(v)
    if v not in self._wrt_variables:
        dedv=self.expr.diff(v)
        if getattr(dedv,'is_Derivative',False):
            return Derivative(dedv.expr,*(tuple(self.variable_count)+tuple(dedv.variable_count)))
        # self.variables raises TypeError for a symbolic count, as in SymPy.
        return Derivative(dedv,*self.variables,evaluate=True)
    return Derivative(self.expr,*(tuple(self.variable_count)+((v,1),)),evaluate=False)


class _DerivativeOnly:
    """A method that exists only on held Derivative nodes."""
    def __init__(self, function):
        self.function=function
        self.__name__=function.__name__

    def __get__(self, obj, owner=None):
        if obj is None or not getattr(obj,'is_Derivative',False):
            raise AttributeError(self.__name__)
        return self.function.__get__(obj,owner)


def install():
    from ._sparse import _SparseBase
    if '_eval_derivative' not in _core.Expr.__dict__:
        _derivative_eval_derivative.__name__='_eval_derivative'
        _core.Expr._eval_derivative=_DerivativeOnly(_derivative_eval_derivative)
    _core.Expr.doit=_derivative_doit(_core.Expr.doit)
    _core.Expr.diff=scalar_diff
    _core.Matrix.diff=matrix_diff
    _SparseBase.diff=matrix_diff
