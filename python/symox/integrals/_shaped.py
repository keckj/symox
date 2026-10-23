"""Held integrals of arrays, explicit matrices and matrix expressions.

SymPy's ``Integral`` stores any Basic integrand; a shaped integrand (an
N-dimensional array, an explicit matrix, a matrix expression) is not a
native scalar expression, so its held integral is this small node. It keeps
the integrand's kind and shape, as SymPy's ``Integral.kind`` and ``shape``
do, and integrates entrywise on ``doit``.
"""
from .. import _core


def is_shaped_integrand(value):
    from ..tensor.array import NDimArray
    from ..matrices import MatrixBase
    return isinstance(value, (NDimArray, MatrixBase, _core._MatrixExpr))


def _normalized_limits(limits):
    from .. import Dummy
    return _core._calculus("Integral", Dummy('f'), limits).limits


class ShapedIntegral:
    """``Integral(f, *limits)`` for a shaped integrand ``f``."""
    is_Integral = True
    is_commutative = False
    is_number = False

    def __init__(self, function, limits):
        from ..matrices import MatrixBase
        if isinstance(function, MatrixBase):
            from .. import ImmutableDenseMatrix
            function = ImmutableDenseMatrix(function)
        if not limits:
            free = function.free_symbols
            if len(free) != 1:
                raise ValueError("specify dummy variables for %s" % (function,))
            limits = (next(iter(free)),)
        self.function = function
        self.limits = _normalized_limits(limits)

    @property
    def args(self):
        return (self.function,) + self.limits

    @property
    def func(self):
        from .. import Integral
        return lambda function, *limits: Integral(function, *limits)

    @property
    def variables(self):
        return [limit[0] for limit in self.limits]

    @property
    def kind(self):
        return self.function.kind

    @property
    def shape(self):
        return tuple(self.function.shape)

    @property
    def free_symbols(self):
        free = set(self.function.free_symbols)
        for limit in reversed(self.limits):
            if len(limit) == 1:
                # an indefinite variable stays free, as in SymPy
                free.add(limit[0])
                continue
            free.discard(limit[0])
            for bound in limit[1:]:
                free |= bound.free_symbols
        return free

    def doit(self, **hints):
        from .. import integrate
        from ..tensor.array import NDimArray
        function = self.function
        if hints.get('deep', True) and hasattr(function, 'doit'):
            function = function.doit(**hints)
        if isinstance(function, _core._MatrixExpr):
            # a matrix expression constant in the integration variables
            # integrates to a scalar multiple of itself
            result = function
            for limit in self.limits:
                x = limit[0]
                if x in result.free_symbols:
                    return ShapedIntegral(function, self.limits)
                if len(limit) == 1:
                    result = x * result
                elif len(limit) == 3:
                    result = (limit[2] - limit[1]) * result
                else:
                    return self
            return result
        limits = self.limits

        def entry(value):
            return integrate(value, *limits)
        if isinstance(function, NDimArray):
            return function.applyfunc(entry)
        return function.applyfunc(entry)

    def subs(self, *args, **kwargs):
        from .. import Integral
        function = self.function.subs(*args, **kwargs)
        limits = [tuple(limit[:1]) + tuple(b.subs(*args, **kwargs) for b in limit[1:]) for limit in self.limits]
        return Integral(function, *limits)

    def _limit_str(self, limit):
        return str(limit[0]) if len(limit) == 1 else '(' + ', '.join(map(str, limit)) + ')'

    def __str__(self):
        return 'Integral(%s, %s)' % (self.function, ', '.join(self._limit_str(l) for l in self.limits))

    __repr__ = __str__

    def __eq__(self, other):
        return isinstance(other, ShapedIntegral) and self.args == other.args

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((ShapedIntegral, self.args))
