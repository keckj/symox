"""Fourier coefficients from native exact definite integration.

Finite trigonometric polynomials are detected algebraically. Other series
compute harmonics lazily, retaining exact endpoints and symbolic parameters.
"""
import copy
import operator

from .. import _core
from .._expr_analysis import _variable

_MAX_HARMONIC = 4096


def _interval_pieces(expression, x, start, stop):
    """Split real affine Abs/sign/Piecewise boundaries before native integration."""
    from functools import cmp_to_key
    from .. import Add, S, integrate
    cuts, pending, seen = {start, stop}, [expression], set()
    while pending:
        node = pending.pop()
        if node in seen:
            continue
        seen.add(node)
        boundary = None
        if node.is_Function and not node.is_UndefinedFunction and node.name in ('Abs', 'sign', 'Heaviside'):
            boundary = node.args[0]
        elif node.is_Relational:
            boundary = node.args[0]-node.args[1]
        if boundary is not None and x in boundary.free_symbols:
            slope = boundary.diff(x)
            if x in slope.free_symbols or slope.is_zero is not False:
                raise NotImplementedError('Fourier branch splitting requires affine real boundaries')
            intercept = boundary.subs(x, 0)
            if slope.is_real is not True or intercept.is_real is not True:
                raise NotImplementedError('Fourier branch splitting requires real affine coefficients')
            root = -intercept/slope
            after_start, before_stop = (root-start).is_positive, (stop-root).is_positive
            if after_start is True and before_stop is True:
                cuts.add(root)
            elif root not in (start,stop) and after_start is not False and before_stop is not False:
                raise NotImplementedError('Cannot locate a symbolic Fourier branch boundary')
        pending.extend(node.args)
    if len(cuts) > 128:
        raise NotImplementedError('Fourier branch splitting limit exceeded')

    def compare(a, b):
        if a == b:
            return 0
        difference = a-b
        if difference.is_positive is True:
            return 1
        if difference.is_negative is True:
            return -1
        raise NotImplementedError('Cannot order symbolic Fourier interval boundaries')

    cuts = sorted(cuts, key=cmp_to_key(compare))
    results = []
    for left, right in zip(cuts, cuts[1:]):
        midpoint, memo = (left+right)/2, {}

        def branch(node):
            if node in memo:
                return memo[node]
            replacements = {arg: branch(arg) for arg in node.args if arg.args}
            replacements = {a: b for a, b in replacements.items() if a != b}
            value = node.xreplace(replacements) if replacements else node
            if value.is_Function and not value.is_UndefinedFunction and value.name in ('Abs', 'sign', 'Heaviside'):
                argument = value.args[0]
                sample = argument.subs(x, midpoint)
                if sample.is_positive is True:
                    value = argument if value.name == 'Abs' else S.One
                elif sample.is_negative is True:
                    value = -argument if value.name == 'Abs' else (-S.One if value.name == 'sign' else S.Zero)
            elif value.is_Piecewise:
                for term, condition in zip(value.args[::2], value.args[1::2]):
                    decision = condition.subs(x, midpoint)
                    if decision == S.true:
                        value = term
                        break
                    if decision != S.false:
                        raise NotImplementedError('Cannot decide a Fourier Piecewise branch')
            memo[node] = value
            return value

        results.append((branch(expression), left, right))
    return results


def _definite_integral(expression, x, start, stop):
    from .. import Add, integrate
    return Add(*(integrate(value, (x, left, right))
                 for value, left, right in _interval_pieces(expression, x, start, stop)))


def _index(value):
    value = operator.index(value)
    if not 0 <= value <= _MAX_HARMONIC:
        raise ValueError('Fourier harmonic must be between 0 and 4096')
    return value


def _finite_terms(function, x, omega):
    from .. import S, sin, cos, simplify
    mean, an, bn = S.Zero, {}, {}
    function = function.expand()
    for term in function.args if function.is_Add else (function,):
        if x not in term.free_symbols:
            mean += term
            continue
        coefficient, wave = term.as_independent(x, as_Add=False)
        if not wave.is_Function or wave.is_UndefinedFunction or wave.name not in ('sin', 'cos'):
            return None
        angle = wave.args[0]
        slope = angle.diff(x)
        if x in slope.free_symbols:
            return None
        harmonic = simplify(slope/omega)
        if not harmonic.is_Integer:
            return None
        harmonic = int(harmonic)
        if not harmonic:
            return None
        phase = simplify(angle-slope*x)
        sign, harmonic = (1 if harmonic > 0 else -1), abs(harmonic)
        if wave.name == 'sin':
            cosine, sine = coefficient*sin(phase), coefficient*sign*cos(phase)
        else:
            cosine, sine = coefficient*cos(phase), -coefficient*sign*sin(phase)
        an[harmonic] = an.get(harmonic, S.Zero)+cosine
        bn[harmonic] = bn.get(harmonic, S.Zero)+sine
    return mean, {k: v for k, v in an.items() if v != 0}, {k: v for k, v in bn.items() if v != 0}


class _FourierSequence:
    def __init__(self, owner, sine):
        self._owner, self._sine = owner, sine

    def coeff(self, index):
        from .. import cos, sin
        index = _index(index)
        coefficient = self._owner._coefficient(index, self._sine)
        wave = sin if self._sine else cos
        return self._owner._scale*coefficient*wave(index*self._owner._omega*self._owner._basis)

    __getitem__ = coeff


class FourierSeries:
    def __init__(self, function, limits, finite_terms=None):
        from .. import S, pi
        self.function = _core.sympify(function)
        self.x, self._lower, self._upper = limits
        self._integrand = self.function
        self._omega = 2*pi/(self._upper-self._lower)
        self._basis = self.x
        self._scale, self._offset = S.One, S.Zero
        self._finite_terms = finite_terms
        self._cache = {}

    @property
    def args(self):
        return self.function, (self.x, self._lower, self._upper)

    @property
    def interval(self):
        from .. import Interval
        return Interval(self.start, self.stop)

    @property
    def start(self):
        from .. import S
        return S.Zero

    @property
    def stop(self):
        from .. import S
        if not isinstance(self, FiniteFourierSeries):
            return S.Infinity
        harmonics = set(self._finite_terms[1]) | set(self._finite_terms[2])
        return S(max(harmonics, default=0)+1+(1 if self.a0 != 0 else 0))

    @property
    def period(self):
        return self._lower, self._upper

    @property
    def _length(self):
        return self._upper-self._lower

    @property
    def free_symbols(self):
        return self.function.free_symbols | self._lower.free_symbols | self._upper.free_symbols

    def _integrate(self, expression):
        from .. import integrate, simplify
        result = simplify(_definite_integral(expression, self.x, self._lower, self._upper))
        pending, seen = [result], set()
        while pending:
            node = pending.pop()
            if node in seen:
                continue
            seen.add(node)
            if node.is_Integral:
                raise NotImplementedError('The native integrator could not determine this Fourier coefficient')
            pending.extend(node.args)
        if result.is_finite is False:
            raise ValueError('The Fourier coefficient integral is not finite')
        return result

    @property
    def a0(self):
        if self._finite_terms is not None:
            mean = self._finite_terms[0]
        else:
            if 'mean' not in self._cache:
                self._cache['mean'] = self._integrate(self._integrand)/self._length
            mean = self._cache['mean']
        return self._offset + self._scale*mean

    def _coefficient(self, index, sine):
        from .. import sin, cos, S
        index = _index(index)
        if self._finite_terms is not None:
            return self._finite_terms[2 if sine else 1].get(index, S.Zero)
        key = index, sine
        if key not in self._cache:
            wave = sin if sine else cos
            self._cache[key] = 2*self._integrate(
                self._integrand*wave(index*self._omega*self.x))/self._length
        return self._cache[key]

    @property
    def an(self):
        return _FourierSequence(self, False)

    @property
    def bn(self):
        return _FourierSequence(self, True)

    def term(self, index):
        from .. import cos, sin
        index = _index(index)
        if isinstance(self, FiniteFourierSeries) and index > self.stop:
            raise IndexError('Fourier harmonic index out of bounds')
        if index == 0:
            return self.a0
        angle = index*self._omega*self._basis
        return self._scale*(self._coefficient(index, False)*cos(angle)
                            + self._coefficient(index, True)*sin(angle))

    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.stop is None:
                raise ValueError('A finite stop is required for Fourier slices')
            return [self.term(i) for i in range(index.start or 0, _index(index.stop), index.step or 1)]
        return self.term(index)

    def _nonzero_terms(self):
        from .. import S
        if self.a0 != 0:
            yield self.a0
        if self._scale == 0:
            return
        if self._finite_terms is not None:
            harmonics = sorted(set(self._finite_terms[1]) | set(self._finite_terms[2]))
            for harmonic in harmonics:
                term = self.term(harmonic)
                if term != 0:
                    yield term
            return
        for harmonic in range(1, _MAX_HARMONIC+1):
            term = self.term(harmonic)
            if term != S.Zero:
                yield term
        raise NotImplementedError('Fourier nonzero-term search exceeded harmonic 4096')

    def __iter__(self):
        return self._nonzero_terms()

    def truncate(self, n=3):
        from itertools import islice
        from .. import Add
        if n is None:
            return iter(self)
        return Add(*islice(self._nonzero_terms(), _index(n)))

    def sigma_approximation(self, n=3):
        from .. import Add, pi, sinc
        n = _index(n)
        if n == 0:
            return self.a0
        stop = min(n, int(self.stop)+1) if isinstance(self, FiniteFourierSeries) else n
        return self.a0+Add(*(self.term(k)*sinc(pi*k/n) for k in range(1, stop)))

    def shift(self, s):
        s = _core.sympify(s)
        if self.x in s.free_symbols:
            raise ValueError('Fourier shift must be independent of the series variable')
        result = copy.copy(self)
        result.function += s
        result._offset += s
        return result

    def shiftx(self, s):
        s = _core.sympify(s)
        if self.x in s.free_symbols:
            raise ValueError('Fourier argument shift must be independent of the series variable')
        result = copy.copy(self)
        result.function = self.function.subs(self.x, self.x+s)
        result._basis = self._basis.subs(self.x, self.x+s)
        return result

    def scale(self, s):
        s = _core.sympify(s)
        if self.x in s.free_symbols:
            raise ValueError('Fourier scale must be independent of the series variable')
        result = copy.copy(self)
        result.function *= s
        result._scale *= s
        result._offset *= s
        return result

    def scalex(self, s):
        s = _core.sympify(s)
        if self.x in s.free_symbols:
            raise ValueError('Fourier argument scale must be independent of the series variable')
        result = copy.copy(self)
        result.function = self.function.subs(self.x, s*self.x)
        result._basis = self._basis.subs(self.x, s*self.x)
        return result

    def to_sympy(self):
        import sympy
        from .. import simplify
        result = sympy.fourier_series(self._integrand.to_sympy(),
                                     (self.x.to_sympy(), self._lower.to_sympy(), self._upper.to_sympy()))
        slope = self._basis.diff(self.x)
        shift = simplify(self._basis-slope*self.x)
        if shift != 0:
            result = result.shiftx(shift.to_sympy())
        if slope != 1:
            result = result.scalex(slope.to_sympy())
        return result.scale(self._scale.to_sympy()).shift(self._offset.to_sympy())

    def __repr__(self):
        return 'FourierSeries(%s, (%s, %s, %s))' % (self.function, self.x, self._lower, self._upper)


class FiniteFourierSeries(FourierSeries):
    """A Fourier series whose harmonics are known to terminate."""


def fourier_series(f, limits=None, finite=True):
    from .. import pi, S
    if getattr(f, 'is_MatrixExpr', False):
        from .matrix import matrix_fourier_series
        return matrix_fourier_series(f, limits, finite)
    f = _core.sympify(f)
    if limits is None:
        x = _variable(f, None)
        start, stop = -pi, pi
    else:
        if len(limits) != 3:
            raise ValueError('Fourier limits must be (variable, start, stop)')
        x, start, stop = map(_core.sympify, limits)
        x = _variable(f, x)
    if x in start.free_symbols or x in stop.free_symbols:
        raise ValueError('Fourier endpoints must be independent of the series variable')
    length = stop-start
    if start.is_finite is False or stop.is_finite is False or length.is_positive is False or length == 0:
        raise ValueError('Fourier limits must describe a finite interval of positive length')
    if x not in f.free_symbols:
        return f
    terms = _finite_terms(f, x, 2*pi/length) if finite else None
    finite_expression = terms is not None
    if terms is None and finite:
        try:
            pieces = _interval_pieces(f, x, start, stop)
        except NotImplementedError:
            pieces = []
        if pieces and all(value == pieces[0][0] for value, _, _ in pieces):
            terms = _finite_terms(pieces[0][0], x, 2*pi/length)
    cls = FiniteFourierSeries if finite_expression else FourierSeries
    return cls(f, (x, start, stop), terms)
