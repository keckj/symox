"""Lazy exact formal series over the native Laurent/Puiseux series engine.

Coefficients are computed on demand by the native recurrence engine. This
does not depend on SymPy or manufacture a closed coefficient formula when
the native engine cannot establish one.
"""
import operator

from .. import _core
from .._expr_analysis import _variable

_DEFAULT = object()
_MAX_ORDER = 256


def _order(value):
    value = operator.index(value)
    if not 0 <= value <= _MAX_ORDER:
        raise ValueError('Formal series order must be between 0 and 256')
    return value


class _Sequence:
    def __init__(self, series, powers=False):
        self._series, self._powers = series, powers

    def coeff(self, index):
        index = _order(index)
        if self._powers:
            return self._series._small**index
        return self._series.coefficient(index)

    __getitem__ = coeff


def _split_symbolic_power(f, x, x0, direction):
    """``(x**e, g)`` when ``f = x**e*g`` at ``x0 = 0`` from the right with a
    symbolic exponent ``e`` free of ``x`` (the reference ``_compute_fps``
    separates such a term as ``symb`` and multiplies it back), else
    ``(1, f)``."""
    from .. import S, Mul
    if x0 != 0 or direction != 1:
        return S.One, f
    symb, rest = [], []
    for factor in (f.args if f.is_Mul else (f,)):
        if (factor.is_Pow and factor.base == x and not factor.exp.is_number
                and not factor.exp.has(x)):
            symb.append(factor)
        else:
            rest.append(factor)
    if not symb:
        return S.One, f
    return Mul(*symb), Mul(*rest)


class FormalPowerSeries:
    """A lazily evaluated formal expansion, with exact finite coefficients."""
    def __init__(self, function, x, x0=0, dir=1):
        from .. import S
        self.function, self.x, self.x0 = map(_core.sympify, (function, x, x0))
        self.dir = _core.sympify(dir)
        self._symb, self._rest = _split_symbolic_power(self.function, self.x, self.x0, self.dir)
        self._default_order = 6
        self._cache = {}
        if self.x0 in (S.Infinity, S.NegativeInfinity):
            self._small = (1 if self.x0 == S.Infinity else -1)/self.x
        else:
            self._small = self.dir*(self.x-self.x0)
        self.ak, self.xk = _Sequence(self), _Sequence(self, powers=True)

    @property
    def free_symbols(self):
        return self.function.free_symbols | self.x0.free_symbols

    @property
    def args(self):
        return self.function, self.x, self.x0, self.dir

    def _series(self, n):
        from .. import O
        n = _order(n)
        if n not in self._cache:
            if self._symb != 1:
                from .. import Add, expand, powsimp
                value = expand(self._rest.series(self.x, self.x0, n, '+').removeO())
                terms = value.args if value.is_Add else (value,)
                self._cache[n] = (Add(*(powsimp(self._symb*t, combine='exp') for t in terms))
                                  + O(powsimp(self._symb*self._small**n, combine='exp'), (self.x, self.x0)))
                return self._cache[n]
            value = self.function.series(self.x, self.x0, n, '+' if self.dir == 1 else '-')
            self._cache[n] = value.removeO() + O(self._small**n, (self.x, self.x0))
        return self._cache[n]

    def _plain_series(self, n):
        """The expansion of the function without its symbolic power factor."""
        if self._symb == 1:
            return self._series(n)
        from .. import O
        return (self._rest.series(self.x, self.x0, n, '+').removeO()
                + O(self._small**n, (self.x, self.x0)))

    def truncate(self, n=_DEFAULT):
        if n is _DEFAULT:
            n = self._default_order
        if n is None:
            return iter(self)
        return self._series(n)

    def polynomial(self, n=6):
        return self._series(n).removeO()

    def coefficient(self, index):
        from .. import Dummy, S
        index = _order(index)
        slot = Dummy('t', positive=True)
        if self.x0 in (S.Infinity, S.NegativeInfinity):
            change = (1 if self.x0 == S.Infinity else -1)/slot
        else:
            change = self.x0 + self.dir*slot
        polynomial = self._plain_series(index+1).removeO().subs(self.x, change).expand()
        return polynomial.coeff(slot, index)

    def term(self, index):
        index = _order(index)
        if self._symb != 1:
            from .. import powsimp
            return powsimp(self.coefficient(index)*self._symb*self._small**index, combine='exp')
        return self.coefficient(index)*self._small**index

    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.stop is None:
                raise ValueError('A finite stop is required for formal-series slices')
            return [self.term(i) for i in range(index.start or 0, _order(index.stop), index.step or 1)]
        return self.term(index)

    def __iter__(self):
        index = 0
        while True:
            yield self.term(index)
            index += 1

    def _compatible(self, other, x):
        if not isinstance(other, FormalPowerSeries):
            raise ValueError('Expected a FormalPowerSeries')
        x = self.x if x is None else _core.sympify(x)
        if self.x != other.x or x != self.x or self.x0 != other.x0 or self.dir != other.dir:
            raise ValueError('Formal series must have the same variable, point and direction')

    def product(self, other, x=None, n=6):
        self._compatible(other, x)
        result = FormalPowerSeries(self.function*other.function, self.x, self.x0, self.dir)
        result._default_order = _order(n)
        return result

    def compose(self, other, x=None, n=6):
        self._compatible(other, x)
        if other.function.subs(self.x, self.x0) != self.x0:
            raise ValueError('The inner series must have zero constant term relative to the expansion point')
        result = FormalPowerSeries(self.function.subs(self.x, other.function), self.x, self.x0, self.dir)
        result._default_order = _order(n)
        return result

    def inverse(self, x=None, n=6):
        x = self.x if x is None else _core.sympify(x)
        if x != self.x or self.function.subs(self.x, self.x0).is_zero is not False:
            raise ValueError('Formal inverse requires a nonzero constant coefficient')
        result = FormalPowerSeries(1/self.function, self.x, self.x0, self.dir)
        result._default_order = _order(n)
        return result

    def diff(self, x=None):
        x = self.x if x is None else _core.sympify(x)
        return fps(self.function.diff(x), self.x, self.x0, self.dir)

    def integrate(self, x=None):
        from .. import integrate
        x = self.x if x is None else _core.sympify(x)
        function = integrate(self.function, x)
        return fps(function, self.x, self.x0, self.dir)

    def to_sympy(self):
        import sympy
        return sympy.fps(self.function.to_sympy(), self.x.to_sympy(),
                         self.x0.to_sympy(), dir=int(self.dir))

    def __repr__(self):
        return 'FormalPowerSeries(%s, %s, %s, %s)' % self.args


def fps(f, x=None, x0=0, dir=1, hyper=True, order=4, rational=True, full=False):
    """Construct an exact lazy formal expansion where native series succeed.

    Algorithm-selection hints are validated; either enabled strategy uses the
    native recurrence engine. With both strategies disabled the input is held.
    """
    from .. import S, Dummy
    if getattr(f, 'is_MatrixExpr', False):
        from .matrix import matrix_fps
        return matrix_fps(f, x, x0, dir, hyper, order, rational, full)
    f = _core.sympify(f)
    x = _variable(f, x)
    x0 = _core.sympify(x0)
    if dir in ('+', 1):
        direction = 1
    elif dir in ('-', -1):
        direction = -1
    else:
        raise ValueError("dir must be 1, -1, '+' or '-'")
    _order(order)
    if x not in f.free_symbols or not (hyper or rational):
        return f
    # A fractional power, logarithm of the local variable or essential
    # singularity does not have an ordinary integer-power formal expansion.
    probe = _split_symbolic_power(f, x, x0, direction)[1]
    if x not in probe.free_symbols:
        return f
    probe = probe.series(x, x0, 3, '+' if direction == 1 else '-')
    slot = Dummy('t', positive=True)
    change = (1/slot if x0 == S.Infinity else -1/slot) if x0 in (S.Infinity, S.NegativeInfinity) else x0+direction*slot
    polynomial = probe.removeO().subs(x, change).expand()
    for term in polynomial.args if polynomial.is_Add else (polynomial,):
        coefficient, exponent = term.as_coeff_exponent(slot)
        if slot in coefficient.free_symbols or not exponent.is_Integer or exponent.is_negative:
            return f
    return FormalPowerSeries(f, x, x0, direction)
