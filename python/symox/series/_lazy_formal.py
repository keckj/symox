"""Lazy exact expansions backing the matrix formal series (``matrix_fps``).

Scalar ``fps`` is the native port of the reference formal power series
(``formal.py``). This base keeps the coefficient-on-demand container used by
``MatrixFormalPowerSeries``, whose coefficients come from the native series
engine.
"""
import operator

from .. import _core

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


class LazyFormalSeries:
    """A lazily evaluated expansion with exact finite coefficients."""
    def __init__(self, function, x, x0=0, dir=1):
        from .. import S
        self.function, self.x, self.x0 = map(_core.sympify, (function, x, x0))
        self.dir = _core.sympify(dir)
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

    def truncate(self, n=_DEFAULT):
        if n is _DEFAULT:
            n = self._default_order
        if n is None:
            return iter(self)
        return self._series(n)

    def term(self, index):
        index = _order(index)
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
        from .formal import FormalPowerSeries
        if not isinstance(other, (LazyFormalSeries, FormalPowerSeries)):
            raise ValueError('Expected a FormalPowerSeries')
        x = self.x if x is None else _core.sympify(x)
        if self.x != other.x or x != self.x or self.x0 != other.x0 or self.dir != other.dir:
            raise ValueError('Formal series must have the same variable, point and direction')

    def __repr__(self):
        return 'FormalPowerSeries(%s, %s, %s, %s)' % self.args
