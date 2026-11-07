"""Formal power series (reference: ``sympy/series/formal.py``).

The algorithms (``rational_algorithm``, ``rational_independent``,
``simpleDE``, ``exp_re``, ``hyper_re``, ``rsolve_hypergeometric``,
``solve_de``, ``hyper_algorithm``, ``compute_fps``) and the term arithmetic
of the series classes are native (``symox::series::formal``). This module
holds the API envelopes: argument admission, the class shells and the
dispatch to the native functions.

A series node stores its coefficient sequences as held ``Sum`` nodes,
``(function, x, x0, dir, (Sum(ak), Sum(xk), ind))``; ``ak`` and ``xk``
return the reference ``SeqFormula`` objects.
"""
import operator

from .. import _core
from .._expr_facade import Expr

__all__ = ['FormalPowerSeries', 'FiniteFormalPowerSeries', 'FormalPowerSeriesProduct',
           'FormalPowerSeriesCompose', 'FormalPowerSeriesInverse', 'Coeff', 'fps',
           'compute_fps', 'rational_algorithm', 'rational_independent', 'simpleDE',
           'exp_re', 'hyper_re', 'rsolve_hypergeometric', 'solve_de', 'hyper_algorithm']


def _sympify(value):
    from .. import sympify
    return sympify(value)


def _seq(node):
    """The ``SeqFormula`` encoded by a held ``Sum(formula, (k, start, stop))``."""
    from .sequences import SeqFormula
    formula, limits = node.args
    return SeqFormula(formula, tuple(limits.args))


def _held_sum(seq):
    from .. import Sum
    return Sum(seq.formula, (seq.variables[0], seq.start, seq.stop))


def _parts(value):
    """The fifth argument as the native ``Tuple(Sum(ak), Sum(xk), ind)``."""
    from .. import Tuple
    from .sequences import SeqBase
    if isinstance(value, (tuple, list)) or (isinstance(value, _core.Expr)
                                            and getattr(value.func, '__name__', '') == 'Tuple'):
        items = tuple(value.args if isinstance(value, _core.Expr) else value)
        if len(items) == 3 and any(isinstance(i, SeqBase) for i in items[:2]):
            return Tuple(*[_held_sum(i) if isinstance(i, SeqBase) else _sympify(i) for i in items])
    return _sympify(value)


def _rebuild(name, args):
    return globals()[name](*args)


class FormalPowerSeries(Expr):
    """A formal power series of a function (see the reference docstring)."""

    _symox_user_hooks = True
    __module__ = 'sympy.series.formal'

    def __new__(cls, *args):
        args = [_parts(a) if i == 4 else _sympify(a) for i, a in enumerate(args)]
        return Expr.__new__(cls, *args)

    def __reduce__(self):
        # The class reports the reference module; unpickle through symox's.
        return (_rebuild, (type(self).__name__, self.args))

    def _sympystr(self, printer):
        if len(self.args) != 5:
            return '%s(%s)' % (type(self).__name__, ', '.join(printer._print(a) for a in self.args))
        head = ', '.join(printer._print(a) for a in self.args[:4])
        seqs = []
        for node in self.args[4].args[:2]:
            formula, (k, start, stop) = node.args[0], node.args[1].args
            seqs.append('SeqFormula(%s, (%s, %s, %s))' % tuple(printer._print(e) for e in (formula, k, start, stop)))
        return 'FormalPowerSeries(%s, (%s, %s, %s))' % (head, seqs[0], seqs[1], printer._print(self.args[4].args[2]))

    def to_sympy(self):
        import sympy
        from sympy.series.formal import FormalPowerSeries as F
        conv = lambda e: e.to_sympy()
        return F(*map(conv, self.args[:4]), (self.ak.to_sympy(), self.xk.to_sympy(), conv(self.ind)))

    @property
    def function(self):
        return self.args[0]

    @property
    def x(self):
        return self.args[1]

    @property
    def x0(self):
        return self.args[2]

    @property
    def dir(self):
        return self.args[3]

    @property
    def ak(self):
        return _seq(self.args[4].args[0])

    @property
    def xk(self):
        return _seq(self.args[4].args[1])

    @property
    def ind(self):
        return self.args[4].args[2]

    @property
    def interval(self):
        from .. import Interval, oo
        return Interval(0, oo)

    @property
    def start(self):
        return self.interval.inf

    @property
    def stop(self):
        return self.interval.sup

    @property
    def length(self):
        from .. import oo
        return oo

    @property
    def variables(self):
        return ()

    @property
    def ak_seq(self):
        from .. import oo
        from .sequences import sequence
        ak = self.ak
        return sequence(ak.formula, (ak.variables[0], 1, oo))

    @property
    def fact_seq(self):
        from .. import factorial, oo
        from .sequences import sequence
        k = self.ak.variables[0]
        return sequence(factorial(k), (k, 1, oo))

    @property
    def bell_coeff_seq(self):
        return self.ak_seq*self.fact_seq

    @property
    def sign_seq(self):
        from .. import oo
        from .sequences import sequence
        return sequence((-1, 1), (self.ak.variables[0], 1, oo))

    def coeff_bell(self, n):
        from .. import Dummy, bell, oo
        from .sequences import sequence
        inner = [bell(n, j, tuple(self.bell_coeff_seq[:n-j+1])) for j in range(1, n+1)]
        return sequence(tuple(inner), (Dummy('k'), 1, oo))

    @property
    def infinite(self):
        return _core.fps_infinite(self)

    def _series(self, function, parts):
        return self.func(function, self.x, self.x0, self.dir, parts)

    def polynomial(self, n=6):
        return _core.fps_polynomial(self, _sympify(n))

    def truncate(self, n=6):
        if n is None:
            return iter(self)
        return _core.fps_truncate(self, _sympify(n))

    def zero_coeff(self):
        return self._eval_term(0)

    def _eval_term(self, pt):
        return _core.fps_term(self, _sympify(pt))

    def term(self, pt):
        if pt < self.start or pt > self.stop:
            raise IndexError("Index %s out of bounds %s" % (pt, self.interval))
        return self._eval_term(pt)

    def _ith_point(self, i):
        return self.start + i

    def __iter__(self):
        i = 0
        while i < self.length:
            yield self.term(self._ith_point(i))
            i += 1

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.term(self._ith_point(index))
        elif isinstance(index, slice):
            start, stop = index.start, index.stop
            if start is None:
                start = 0
            if stop is None:
                stop = self.length
            return [self.term(self._ith_point(i)) for i in range(start, stop, index.step or 1)]

    def _eval_subs(self, old, new):
        if old.has(self.x):
            return self

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        return _core.fps_leading_term(self)

    def _eval_derivative(self, x):
        return self._series(*_core.fps_derivative(self, x))

    def _eval_Integral(self, *limits, **kwargs):
        # The reference Integral.doit hands a formal power series to its
        # integrate method (one limit only).
        from .. import Tuple
        if not limits:
            free = self.free_symbols
            if len(free) != 1:
                raise ValueError('specify dummy variables for %s' % self)
            limits = tuple(free)
        if len(limits) > 1:
            raise NotImplementedError
        xab = limits[0]
        if isinstance(xab, (tuple, list, Tuple)):
            if len(xab) > 1:
                return self.integrate(tuple(xab), **kwargs)
            xab = xab[0]
        return self.integrate(xab, **kwargs)

    def integrate(self, x=None, **kwargs):
        from .. import integrate
        from ..utilities.iterables import iterable
        if x is None:
            x = self.x
        elif iterable(x):
            return integrate(self.function, x)
        return self._series(*_core.fps_integrate(self, _sympify(x)))

    def product(self, other, x=None, n=6):
        if n is None:
            return iter(self)
        other = _sympify(other)
        if not isinstance(other, FormalPowerSeries):
            raise ValueError('Both series should be an instance of FormalPowerSeries class.')
        _core.fps_check_pair(self, other, False)
        return FormalPowerSeriesProduct(self, other)

    def compose(self, other, x=None, n=6):
        if n is None:
            return iter(self)
        other = _sympify(other)
        if not isinstance(other, FormalPowerSeries):
            raise ValueError('Both series should be an instance of FormalPowerSeries class.')
        _core.fps_check_pair(self, other, True)
        return FormalPowerSeriesCompose(self, other)

    def inverse(self, x=None, n=6):
        if n is None:
            return iter(self)
        _core.fps_check_inverse(self)
        return FormalPowerSeriesInverse(self)

    def __add__(self, other):
        from .. import Add
        other = _sympify(other)
        if isinstance(other, FormalPowerSeries):
            r = _core.fps_add_fps(self, other)
            return r[0] if len(r) == 1 else self._series(*r)
        r = _core.fps_add_expr(self, other)
        return Add(self, other) if r is None else self._series(*r)

    def __radd__(self, other):
        # The reference reaches ``__radd__`` only from Python numbers; an
        # expression on the left builds ``Add`` (``Expr.__add__``).
        if isinstance(other, _core.Expr):
            from .. import Add
            return Add(other, self)
        return self.__add__(other)

    def __neg__(self):
        return self._series(*_core.fps_neg(self))

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        if isinstance(other, _core.Expr):
            from .. import Add
            return Add(other, -self)
        return (-self).__add__(other)

    def __mul__(self, other):
        from .. import Mul
        other = _sympify(other)
        r = _core.fps_mul_expr(self, other)
        return Mul(self, other) if r is None else self._series(*r)

    def __rmul__(self, other):
        if isinstance(other, _core.Expr):
            from .. import Mul
            return Mul(other, self)
        return self.__mul__(other)


class FiniteFormalPowerSeries(FormalPowerSeries):
    """Base of the series obtained by operations on formal power series."""

    _symox_user_hooks = True
    __module__ = 'sympy.series.formal'
    _kind = None

    def _sympystr(self, printer):
        return '%s(%s)' % (type(self).__name__, ', '.join(printer._print(a) for a in self.args))

    def to_sympy(self):
        import sympy.series.formal as F
        return getattr(F, type(self).__name__)(*(a.to_sympy() for a in self.args))

    @property
    def ffps(self):
        return self.args[0]

    @property
    def gfps(self):
        return self.args[1]

    @property
    def f(self):
        return self.ffps.function

    @property
    def g(self):
        return self.gfps.function

    @property
    def infinite(self):
        raise NotImplementedError('No infinite version for an object of FiniteFormalPowerSeries class.')

    def _eval_terms(self, n):
        if self._kind is None:
            raise NotImplementedError('(%s)._eval_terms()' % self)
        gfps = self.ffps if self._kind == 'inverse' else self.gfps
        return _core.fps_finite_terms(self._kind, self.ffps, gfps, operator.index(n))

    def _eval_term(self, pt):
        raise NotImplementedError('By the current logic, one can get termsupto a certain order, instead of getting term by term.')

    def polynomial(self, n):
        return self._eval_terms(n)

    def truncate(self, n=6):
        return _core.fps_finite_truncate(self.ffps, self.polynomial(n), _sympify(n))

    def _eval_derivative(self, x):
        raise NotImplementedError

    def integrate(self, x):
        raise NotImplementedError


class FormalPowerSeriesProduct(FiniteFormalPowerSeries):
    """The product of two formal power series, computed up to an order."""

    _symox_user_hooks = True
    __module__ = 'sympy.series.formal'
    _kind = 'product'

    @property
    def coeff1(self):
        from .. import oo
        from .sequences import sequence
        ak = self.ffps.ak
        return sequence(ak.formula, (ak.variables[0], 0, oo))

    @property
    def coeff2(self):
        from .. import oo
        from .sequences import sequence
        ak = self.gfps.ak
        return sequence(ak.formula, (ak.variables[0], 0, oo))

    @property
    def function(self):
        return self.f*self.g


class FormalPowerSeriesCompose(FiniteFormalPowerSeries):
    """The composition of two formal power series, computed up to an order."""

    _symox_user_hooks = True
    __module__ = 'sympy.series.formal'
    _kind = 'compose'

    @property
    def function(self):
        f, g, x = self.f, self.g, self.ffps.x
        return f.subs(x, g)


class FormalPowerSeriesInverse(FiniteFormalPowerSeries):
    """The multiplicative inverse of a formal power series, up to an order."""

    _symox_user_hooks = True
    __module__ = 'sympy.series.formal'
    _kind = 'inverse'

    @property
    def aux_seq(self):
        from .. import oo
        from .sequences import sequence
        ffps = self.ffps
        k = ffps.xk.variables[0]
        inv_seq = sequence(ffps.zero_coeff()**(-(k + 1)), (k, 1, oo))
        return ffps.sign_seq*ffps.fact_seq*inv_seq

    @property
    def function(self):
        return 1/self.f

    @property
    def g(self):
        raise ValueError('Only one function is considered while performing'
                         'inverse of a formal power series.')

    @property
    def gfps(self):
        raise ValueError('Only one function is considered while performing'
                         'inverse of a formal power series.')


def _coeff_base():
    from .. import Function
    return Function


class Coeff(_coeff_base()):
    """Coefficient of the ``n``-th power of ``x`` in the polynomial ``p``."""

    _symox_user_hooks = True
    __module__ = 'sympy.series.formal'

    @classmethod
    def eval(cls, p, x, n):
        return _core.fps_coeff_eval(p, x, n)

    def __reduce__(self):
        return (_rebuild, ('Coeff', self.args))


def rational_algorithm(f, x, k, order=4, full=False):
    """Rational algorithm for the formula of the coefficients (native)."""
    return _core.fps_rational_algorithm(_sympify(f), _sympify(x), _sympify(k),
                                        operator.index(order), bool(full))


def rational_independent(terms, x):
    """Group the terms that are rational multiples of each other (native)."""
    return _core.fps_rational_independent([_sympify(t) for t in terms], _sympify(x))


def simpleDE(f, x, g, order=4):
    """Generate simple differential equations satisfied by ``f`` (native, lazy)."""
    x = _sympify(x)
    return _core.FpsSimpleDe(_sympify(f), x, g(x), operator.index(order))


def exp_re(DE, r, k):
    """Recurrence for the coefficients of the exponential series (native)."""
    k = _sympify(k)
    return _core.fps_exp_re(_sympify(DE), r(k), k)


def hyper_re(DE, r, k):
    """Recurrence for the coefficients of the hypergeometric series (native)."""
    k = _sympify(k)
    return _core.fps_hyper_re(_sympify(DE), r(k), k)


def rsolve_hypergeometric(f, x, P, Q, k, m):
    """Solve the hypergeometric recurrence ``Q*r(k + m) = P*r(k)`` (native)."""
    return _core.fps_rsolve_hypergeometric(*map(_sympify, (f, x, P, Q, k, m)))


def solve_de(f, x, DE, order, g, k):
    """Solve the differential equation ``DE`` for the coefficients (native)."""
    x = _sympify(x)
    return _core.fps_solve_de(_sympify(f), x, _sympify(DE), operator.index(order), g(x), _sympify(k))


def hyper_algorithm(f, x, k, order=4):
    """Hypergeometric algorithm for the formula of the coefficients (native)."""
    return _core.fps_hyper_algorithm(_sympify(f), _sympify(x), _sympify(k), operator.index(order))


def _compute(f, x, x0, dir, hyper, order, rational, full):
    from .. import S
    f, x = _sympify(f), _sympify(x)
    if not f.has(x):
        return None
    x0 = _sympify(x0)
    if dir == '+':
        dir = S.One
    elif dir == '-':
        dir = -S.One
    elif dir not in [S.One, -S.One]:
        raise ValueError("Dir must be '+' or '-'")
    else:
        dir = _sympify(dir)
    return _core.fps_compute(f, x, x0, dir, bool(hyper), operator.index(order), bool(rational), bool(full))


def compute_fps(f, x, x0=0, dir=1, hyper=True, order=4, rational=True, full=False):
    """``(ak, xk, ind)`` of the formal power series of ``f``, or ``None`` (native)."""
    result = _compute(f, x, x0, dir, hyper, order, rational, full)
    if result is None:
        return None
    ak, xk, ind = result
    return _seq(ak), _seq(xk), ind


def fps(f, x=None, x0=0, dir=1, hyper=True, order=4, rational=True, full=False):
    """The formal power series of ``f`` about ``x0``, or ``f`` when none is found."""
    if getattr(f, 'is_MatrixExpr', False):
        from .matrix import matrix_fps
        return matrix_fps(f, x, x0, dir, hyper, order, rational, full)
    f = _sympify(f)
    if x is None:
        free = f.free_symbols
        if len(free) == 1:
            x = free.pop()
        elif not free:
            return f
        else:
            raise NotImplementedError('multivariate formal power series')
    result = _compute(f, x, x0, dir, hyper, order, rational, full)
    if result is None:
        return f
    return FormalPowerSeries(f, x, x0, dir, result)
