"""Formal and Fourier series with ordered, constant matrix coefficients."""
import copy
from itertools import islice

from .formal import FormalPowerSeries, _order
from .fourier import FourierSeries, _index, _MAX_HARMONIC
from .._matrixexpr_algebra import linear_terms, matrix_sum, coefficient
from .._matrixexpr_methods import _zero


def _variable(function, x):
    from .. import Symbol, MatrixSymbol, sympify
    if x is None:
        variables = {s for s in function.free_symbols if isinstance(s, Symbol)}
        dimensions = set()
        for matrix in function.atoms(MatrixSymbol):
            dimensions.update(matrix.rows.free_symbols | matrix.cols.free_symbols)
        coefficients = {s for scalar, _ in linear_terms(function) for s in scalar.free_symbols
                        if isinstance(s, Symbol)}
        variables = variables.difference(dimensions) | coefficients
        if not variables:
            return None
        if len(variables) != 1:
            raise ValueError('specify the scalar variable for a matrix series')
        x = next(iter(variables))
    x = sympify(x)
    if not isinstance(x, Symbol):
        raise ValueError('matrix series require a scalar Symbol')
    return x


class MatrixFormalPowerSeries(FormalPowerSeries):
    def __init__(self, function, x, x0=0, dir=1, _check=True):
        super().__init__(0, x, x0, dir)
        self.function = function
        if _check:
            list(linear_terms(function.expand(), self.x))

    def polynomial(self, n=6):
        n = _order(n)
        return self.function.series(self.x, self.x0, n,
                                    '+' if self.dir == 1 else '-').removeO().doit()

    def _series(self, n):
        from .. import O, OneMatrix
        n = _order(n)
        if n not in self._cache:
            remainder = O(self._small**n, (self.x, self.x0))*OneMatrix(*self.function.shape)
            self._cache[n] = self.polynomial(n)+remainder
        return self._cache[n]

    def coefficient(self, index):
        from .. import Dummy, S
        index = _order(index)
        if index == 256:
            raise ValueError('matrix coefficient extraction is limited to index 255')
        small = Dummy('t', positive=True)
        if self.x0 in (S.Infinity, S.NegativeInfinity):
            change = (1 if self.x0 == S.Infinity else -1)/small
        else:
            change = self.x0+self.dir*small
        polynomial = self.polynomial(index+1).subs(self.x, change).expand()
        return coefficient(polynomial, small, index)

    def product(self, other, x=None, n=6):
        self._compatible(other, x)
        result = _ProductMatrixFormalSeries(self, other)
        result._default_order = _order(n)
        return result

    def compose(self, other, x=None, n=6):
        from .. import MatrixExpr
        self._compatible(other, x)
        if isinstance(other.function, MatrixExpr):
            raise ValueError('matrix series composition requires a scalar inner series')
        if other.function.subs(self.x, self.x0) != self.x0:
            raise ValueError('the inner series must map the expansion point to itself')
        result = MatrixFormalPowerSeries(self.function.subs(self.x, other.function),
                                        self.x, self.x0, self.dir)
        result._default_order = _order(n)
        return result

    def inverse(self, x=None, n=6):
        if x is not None and x != self.x:
            raise ValueError('inverse variable must match the series variable')
        if self.function.rows != self.function.cols:
            raise ValueError('formal matrix inverse requires a square matrix')
        result = _InverseMatrixFormalSeries(self)
        result._default_order = _order(n)
        return result

    def diff(self, x=None):
        x = self.x if x is None else x
        return matrix_fps(self.function.diff(x), self.x, self.x0, self.dir)

    def integrate(self, x=None):
        x = self.x if x is None else x
        return matrix_fps(self.function.integrate(x), self.x, self.x0, self.dir)

    def to_sympy(self):
        import sympy
        entries = self.function.as_explicit().to_sympy()
        return entries.applyfunc(lambda e: sympy.fps(e, self.x.to_sympy(),
                                                     self.x0.to_sympy(), dir=int(self.dir)))


class _InverseMatrixFormalSeries(MatrixFormalPowerSeries):
    def __init__(self, original):
        from .. import S
        if original.x0 in (S.Infinity, S.NegativeInfinity):
            raise NotImplementedError('formal matrix inversion currently requires a finite center')
        super().__init__(original.function.inv(), original.x, original.x0, original.dir, _check=False)
        self._original = original
        constant = original.coefficient(0)
        if _zero(constant):
            raise ValueError('the constant matrix coefficient must be invertible')
        # For symbolic matrices the inverse remains formal. Native explicit
        # matrices reject a known singular constant coefficient.
        self._coefficients = [constant.inv()]

    def coefficient(self, index):
        index = _order(index)
        while len(self._coefficients) <= index:
            n = len(self._coefficients)
            convolution = matrix_sum([self._original.coefficient(k)*self._coefficients[n-k]
                                      for k in range(1, n+1)], self.function.shape)
            self._coefficients.append((-self._coefficients[0]*convolution).doit())
        return self._coefficients[index]

    def polynomial(self, n=6):
        n = _order(n)
        return matrix_sum([self.coefficient(k)*self._small**k for k in range(n)], self.function.shape)


class _ProductMatrixFormalSeries(MatrixFormalPowerSeries):
    def __init__(self, left, right):
        super().__init__(left.function*right.function, left.x, left.x0, left.dir, _check=False)
        self._left, self._right = left, right
        self._coefficients = {}

    def coefficient(self, index):
        index = _order(index)
        if index not in self._coefficients:
            self._coefficients[index] = matrix_sum([
                self._left.coefficient(k)*self._right.coefficient(index-k)
                for k in range(index+1)], self.function.shape)
        return self._coefficients[index]

    def polynomial(self, n=6):
        return matrix_sum([self.coefficient(k)*self._small**k for k in range(_order(n))], self.function.shape)


def matrix_fps(function, x=None, x0=0, dir=1, hyper=True, order=4, rational=True, full=False):
    from .. import S, Dummy
    x = _variable(function, x)
    if x is None:
        return function
    if dir not in (1, -1, '+', '-'):
        raise ValueError('formal series direction must be 1 or -1')
    direction = 1 if dir in (1, '+') else -1
    _order(order)
    if x not in function.free_symbols or not (hyper or rational):
        return function
    from .formal import fps
    for scalar, _ in linear_terms(function.expand(), x):
        # Keep Laurent, logarithmic, fractional and essential expansions
        # outside the ordinary formal-power-series interface.
        probe = fps(scalar, x, x0, direction, hyper=hyper, order=order, rational=rational, full=full)
        if x in scalar.free_symbols and not isinstance(probe, FormalPowerSeries):
            return function
    return MatrixFormalPowerSeries(function, x, x0, direction)


class MatrixFourierSeries(FourierSeries):
    def __init__(self, function, limits, finite=True):
        from .fourier import fourier_series
        super().__init__(0, limits)
        self.function = function
        self._components = []
        for scalar, matrix in linear_terms(function.expand(), self.x):
            series = fourier_series(scalar, limits, finite=finite)
            if not isinstance(series, FourierSeries):
                series = FourierSeries(series, limits, (series, {}, {}))
            self._components.append((series, matrix))
        self._harmonics = set()
        for series, _ in self._components:
            if series._finite_terms is None:
                self._harmonics = None
                break
            self._harmonics.update(series._finite_terms[1])
            self._harmonics.update(series._finite_terms[2])

    @property
    def a0(self):
        mean = matrix_sum([series.a0*matrix for series, matrix in self._components], self.function.shape)
        result = self._scale*mean
        return result if self._offset == 0 else self._offset+result

    @property
    def stop(self):
        from .. import S
        return S.Infinity if self._harmonics is None else S(max(self._harmonics, default=0)+2)

    def _coefficient(self, index, sine):
        return matrix_sum([series._coefficient(index, sine)*matrix
                           for series, matrix in self._components], self.function.shape)

    def _nonzero_terms(self):
        if not _zero(self.a0):
            yield self.a0
        if self._scale == 0:
            return
        harmonics = sorted(self._harmonics) if self._harmonics is not None else range(1, _MAX_HARMONIC+1)
        for harmonic in harmonics:
            term = self.term(harmonic)
            if not _zero(term):
                yield term
        if self._harmonics is None:
            raise NotImplementedError('matrix Fourier search exceeded harmonic 4096')

    def truncate(self, n=3):
        if n is None:
            return iter(self)
        return matrix_sum(list(islice(self._nonzero_terms(), _index(n))), self.function.shape)

    def sigma_approximation(self, n=3):
        from .. import pi, sinc
        n = _index(n)
        return matrix_sum([self.a0]+[self.term(k)*sinc(pi*k/n) for k in range(1, n)], self.function.shape)

    def shift(self, s):
        if s == 0:
            return copy.copy(self)
        if getattr(s, 'shape', None) != self.function.shape or self.x in s.free_symbols:
            raise ValueError('a matrix Fourier shift requires a constant matrix of the same shape')
        result = copy.copy(self)
        result.function += s
        result._offset = s if result._offset == 0 else result._offset+s
        return result

    def to_sympy(self):
        import sympy
        entries = self.function.as_explicit().to_sympy()
        limits = tuple(e.to_sympy() for e in (self.x, self._lower, self._upper))
        return entries.applyfunc(lambda e: sympy.fourier_series(e, limits))


def matrix_fourier_series(function, limits=None, finite=True):
    from .. import pi, sympify
    if limits is None:
        x = _variable(function, None)
        if x is None:
            return function
        limits = (x, -pi, pi)
    else:
        if len(limits) != 3:
            raise ValueError('Fourier limits must be (variable, start, stop)')
        limits = tuple(sympify(a) for a in limits)
        x = _variable(function, limits[0])
    from .fourier import fourier_series
    from .. import S
    fourier_series(S.Zero, limits, finite=finite)  # validate even a constant matrix's interval
    if x not in function.free_symbols:
        return function
    return MatrixFourierSeries(function, limits, finite=finite)
