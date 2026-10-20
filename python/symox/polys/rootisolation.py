"""Isolation intervals of polynomial roots.

``RealInterval`` and ``ComplexInterval`` follow SymPy's classes: a real
interval stores a Mobius transform and the transformed polynomial over ZZ,
and refines by continued fractions; a complex interval stores an isolating
rectangle with its edge roots, quadrant sequences and edge polynomials, and
refines by bisection. Both states and their refinement are native and exact.
"""
# Module-level imports: the routing facade rebinds ``__spec__`` of aliased
# modules, so relative imports inside functions would warn.
from .. import _core
from .polyerrors import DomainError as _DomainError

_RealIsolation = _core.RealIsolation
_ComplexIsolation = _core.ComplexIsolation


def _pair(value):
    """``(numerator, denominator)`` of a rational domain element or number."""
    numerator = getattr(value, 'numerator', None)
    if numerator is None:
        numerator, denominator = value.p, value.q
    else:
        denominator = value.denominator
    return int(numerator), int(denominator)


def _ring(dom):
    return dom.get_ring() if dom.is_Field else dom


def _qq(dom, pair):
    return dom.get_field()(*pair)


class RealInterval:
    """A fully qualified representation of a real isolation interval. """

    def __init__(self, data, f, dom):
        """Initialize new real interval with complete information. """
        if len(data) == 2:
            s, t = data
            state = _RealIsolation.from_interval(
                [int(c) for c in f], _pair(s), _pair(t))
        else:
            state = _RealIsolation.from_mobius(
                [int(c) for c in f], tuple(int(c) for c in data[:-1]), bool(data[-1]))
        self._set(state, dom)

    @classmethod
    def _from_state(cls, state, dom):
        obj = cls.__new__(cls)
        obj._set(state, dom)
        return obj

    def _set(self, state, dom):
        self._state, self.dom = state, dom
        ring = _ring(dom)
        self.mobius = tuple(ring(c) for c in state.mobius)
        self.neg = state.neg
        self.f = [ring(c) for c in state.f]

    @property
    def func(self):
        return RealInterval

    @property
    def args(self):
        i = self
        return (i.mobius + (i.neg,), i.f, i.dom)

    def __eq__(self, other):
        if type(other) is not type(self):
            return False
        return self.args == other.args

    @property
    def a(self):
        """Return the position of the left end. """
        return _qq(self.dom, self._state.a())

    @property
    def b(self):
        """Return the position of the right end. """
        return _qq(self.dom, self._state.b())

    @property
    def dx(self):
        """Return width of the real isolating interval. """
        return self.b - self.a

    @property
    def center(self):
        """Return the center of the real isolating interval. """
        return (self.a + self.b)/2

    @property
    def max_denom(self):
        """Return the largest denominator occurring in either endpoint. """
        return max(self.a.denominator, self.b.denominator)

    def as_tuple(self):
        """Return tuple representation of real isolating interval. """
        return (self.a, self.b)

    def __repr__(self):
        return "(%s, %s)" % (self.a, self.b)

    def __contains__(self, item):
        """
        Say whether a complex number belongs to this real interval.

        Parameters
        ==========

        item : pair (re, im) or number re
            Either a pair giving the real and imaginary parts of the number,
            or else a real number.

        """
        if isinstance(item, tuple):
            re, im = item
        else:
            re, im = item, 0
        return im == 0 and self.a <= re <= self.b

    def is_disjoint(self, other):
        """Return ``True`` if two isolation intervals are disjoint. """
        if isinstance(other, RealInterval):
            return (self.b < other.a or other.b < self.a)
        assert isinstance(other, ComplexInterval)
        return (self.b < other.ax or other.bx < self.a
            or other.ay*other.by > 0)

    def _inner_refine(self):
        """Internal one step real root refinement procedure. """
        return RealInterval._from_state(self._state.refine(), self.dom)

    def refine_disjoint(self, other):
        """Refine an isolating interval until it is disjoint with another one. """
        expr = self
        while not expr.is_disjoint(other):
            expr, other = expr._inner_refine(), other._inner_refine()

        return expr, other

    def refine_size(self, dx):
        """Refine an isolating interval until it is of sufficiently small size. """
        expr = self
        while not (expr.dx < dx):
            expr = expr._inner_refine()

        return expr

    def refine_step(self, steps=1):
        """Perform several steps of real root refinement algorithm. """
        expr = self
        for _ in range(steps):
            expr = expr._inner_refine()

        return expr

    def refine(self):
        """Perform one step of real root refinement algorithm. """
        return self._inner_refine()


class ComplexInterval:
    """A fully qualified representation of a complex isolation interval.
    The printed form is shown as (ax, bx) x (ay, by) where (ax, ay)
    and (bx, by) are the coordinates of the southwest and northeast
    corners of the interval's rectangle, respectively.

    The rectangle is stored in the upper half plane; ``conj`` marks the
    conjugate root, whose rectangle is the mirror image.
    """

    def __init__(self, a, b, I, Q, F1, F2, f1, f2, dom, conj=False):
        """Initialize new complex interval with complete information. """
        if any(part is None for part in (I, Q, F1, F2, f1, f2)):
            # Incomplete data, as the reference stores it: the rectangle is
            # usable, refinement is not.
            self._state, self.dom, self.conj = None, dom, conj
            self.a, self.b = tuple(a), tuple(b)
            self.Q = Q
            self._parts = (I, F1, F2, f1, f2)
            return
        field = dom.get_field() if not dom.is_Field else dom
        edges = [[((_pair(s), _pair(t)), dict(indices), [int(c) for c in h])
                  for (s, t), indices, h in edge] for edge in I]
        quadrants = [[str(q) for q in qs] for qs in Q]
        dense = lambda polys: [[_pair(field.convert(c)) for c in g] for g in polys]
        f = [_pair(field.convert(c[-1])) if c else (0, 1) for c in f1]
        state = _ComplexIsolation.from_parts(
            (_pair(a[0]), _pair(a[1])), (_pair(b[0]), _pair(b[1])),
            edges, quadrants, dense(F1), dense(F2), f, bool(conj))
        self._set(state, dom)

    @classmethod
    def _from_state(cls, state, dom):
        obj = cls.__new__(cls)
        obj._set(state, dom)
        return obj

    def _set(self, state, dom):
        self._state, self.dom = state, dom
        self.conj = state.conj
        u, v, s, t = (_qq(dom, c) for c in state.corners())
        self.a, self.b = (u, v), (s, t)
        self.Q = [list(q) for q in state.quadrants()]
        self._parts = None

    def _rebuild(self):
        if self._parts is None:
            dom = self.dom
            field = dom.get_field() if not dom.is_Field else dom
            ring = field.get_ring()
            I = [[((_qq(dom, s), _qq(dom, t)), indices, [ring(c) for c in h])
                  for (s, t), indices, h in edge] for edge in self._state.edges()]
            F1, F2 = ([[_qq(dom, c) for c in g] for g in polys]
                      for polys in self._state.edge_polys())
            f1, f2 = ([[_qq(dom, c) for c in g] for g in polys]
                      for polys in self._state.real_imag())
            self._parts = (I, F1, F2, f1, f2)
        return self._parts

    @property
    def I(self):
        return self._rebuild()[0]

    @property
    def F1(self):
        return self._rebuild()[1]

    @property
    def F2(self):
        return self._rebuild()[2]

    @property
    def f1(self):
        return self._rebuild()[3]

    @property
    def f2(self):
        return self._rebuild()[4]

    @property
    def func(self):
        return ComplexInterval

    @property
    def args(self):
        i = self
        return (i.a, i.b, i.I, i.Q, i.F1, i.F2, i.f1, i.f2, i.dom, i.conj)

    def __eq__(self, other):
        if type(other) is not type(self):
            return False
        return self.args == other.args

    @property
    def ax(self):
        """Return ``x`` coordinate of south-western corner. """
        return self.a[0]

    @property
    def ay(self):
        """Return ``y`` coordinate of south-western corner. """
        if not self.conj:
            return self.a[1]
        else:
            return -self.b[1]

    @property
    def bx(self):
        """Return ``x`` coordinate of north-eastern corner. """
        return self.b[0]

    @property
    def by(self):
        """Return ``y`` coordinate of north-eastern corner. """
        if not self.conj:
            return self.b[1]
        else:
            return -self.a[1]

    @property
    def dx(self):
        """Return width of the complex isolating interval. """
        return self.b[0] - self.a[0]

    @property
    def dy(self):
        """Return height of the complex isolating interval. """
        return self.b[1] - self.a[1]

    @property
    def center(self):
        """Return the center of the complex isolating interval. """
        return ((self.ax + self.bx)/2, (self.ay + self.by)/2)

    @property
    def max_denom(self):
        """Return the largest denominator occurring in either endpoint. """
        return max(self.ax.denominator, self.bx.denominator,
                   self.ay.denominator, self.by.denominator)

    def as_tuple(self):
        """Return tuple representation of the complex isolating
        interval's SW and NE corners, respectively. """
        return ((self.ax, self.ay), (self.bx, self.by))

    def __repr__(self):
        return "(%s, %s) x (%s, %s)" % (self.ax, self.bx, self.ay, self.by)

    def conjugate(self):
        """This complex interval really is located in lower half-plane. """
        if self._state is None:
            return ComplexInterval(self.a, self.b, self.I, self.Q, self.F1,
                                   self.F2, self.f1, self.f2, self.dom, not self.conj)
        return ComplexInterval._from_state(self._state.conjugate(), self.dom)

    def __contains__(self, item):
        """
        Say whether a complex number belongs to this complex rectangular
        region.

        Parameters
        ==========

        item : pair (re, im) or number re
            Either a pair giving the real and imaginary parts of the number,
            or else a real number.

        """
        if isinstance(item, tuple):
            re, im = item
        else:
            re, im = item, 0
        return self.ax <= re <= self.bx and self.ay <= im <= self.by

    def is_disjoint(self, other):
        """Return ``True`` if two isolation intervals are disjoint. """
        if isinstance(other, RealInterval):
            return other.is_disjoint(self)
        if self.conj != other.conj:  # above and below real axis
            return True
        re_distinct = (self.bx < other.ax or other.bx < self.ax)
        if re_distinct:
            return True
        im_distinct = (self.by < other.ay or other.by < self.ay)
        return im_distinct

    def _inner_refine(self):
        """Internal one step complex root refinement procedure. """
        if self._state is None:
            raise TypeError("complex interval refinement needs its edge data")
        return ComplexInterval._from_state(self._state.refine(), self.dom)

    def refine_disjoint(self, other):
        """Refine an isolating interval until it is disjoint with another one. """
        expr = self
        while not expr.is_disjoint(other):
            expr, other = expr._inner_refine(), other._inner_refine()

        return expr, other

    def refine_size(self, dx, dy=None):
        """Refine an isolating interval until it is of sufficiently small size. """
        if dy is None:
            dy = dx
        expr = self
        while not (expr.dx < dx and expr.dy < dy):
            expr = expr._inner_refine()

        return expr

    def refine_step(self, steps=1):
        """Perform several steps of complex root refinement algorithm. """
        expr = self
        for _ in range(steps):
            expr = expr._inner_refine()

        return expr

    def refine(self):
        """Perform one step of complex root refinement algorithm. """
        return self._inner_refine()


def _bound(K, value):
    """``(numerator, denominator)`` of an optional bound, as the field of
    ``K`` converts it."""
    if value is None:
        return None
    from .. import sympify
    F = K.get_field() if not K.is_Field else K
    try:
        value = F.convert(value)
    except Exception:
        value = F.from_sympy(sympify(value))
    return _pair(value)


def _degree(f):
    k = 0
    for c in f:
        if c:
            break
        k += 1
    return len(f) - k - 1 if k < len(f) else -1


def _integer_coeffs(f, K):
    """The coefficients of ``f`` over ZZ, cleared of denominators over QQ."""
    if K.is_QQ:
        pairs = [_pair(c) for c in f]
        from math import lcm
        common = 1
        for _, q in pairs:
            common = lcm(common, q)
        return [p*(common//q) for p, q in pairs], K.get_ring()
    return [int(c) for c in f], K


def dup_isolate_real_roots_sqf(f, K, eps=None, inf=None, sup=None, fast=False, blackbox=False):
    """Isolate real roots of a square-free polynomial using the Vincent-Akritas-Strzebonski (VAS) CF approach. """
    if not (K.is_QQ or K.is_ZZ):
        raise _DomainError("isolation of real roots not supported over %s" % K)
    f, R = _integer_coeffs(f, K)
    if _degree(f) <= 0:
        return []
    F = R.get_field()
    states = _core.isolate_real_roots_sqf_state(
        f, _bound(R, eps), _bound(R, inf), _bound(R, sup), bool(fast))
    if not blackbox:
        return [(_qq(F, state.a()), _qq(F, state.b())) for state in states]
    return [RealInterval._from_state(state, R) for state in states]


def dup_isolate_complex_roots_sqf(f, K, eps=None, inf=None, sup=None, blackbox=False):
    """Isolate complex roots of a square-free polynomial using Collins-Krandick algorithm. """
    if not K.is_ZZ and not K.is_QQ:
        raise _DomainError("isolation of complex roots is not supported over %s" % K)
    F = K.get_field() if K.is_ZZ else K
    if _degree([F.convert(c, K) for c in f]) <= 0:
        return []
    coeffs = [_pair(F.convert(c, K)) for c in f]
    states = _core.isolate_complex_roots_sqf_state(
        coeffs, _bound(F, eps), _bound(F, inf), _bound(F, sup))
    roots = [ComplexInterval._from_state(state, F) for state in states]
    if blackbox:
        return roots
    return [r.as_tuple() for r in roots]


def dup_isolate_all_roots_sqf(f, K, eps=None, inf=None, sup=None, fast=False, blackbox=False):
    """Isolate real and complex roots of a square-free polynomial ``f``. """
    return (
        dup_isolate_real_roots_sqf( f, K, eps=eps, inf=inf, sup=sup, fast=fast, blackbox=blackbox),
        dup_isolate_complex_roots_sqf(f, K, eps=eps, inf=inf, sup=sup, blackbox=blackbox))


def dup_isolate_all_roots(f, K, eps=None, inf=None, sup=None, fast=False):
    """Isolate real and complex roots of a non-square-free polynomial ``f``. """
    if not K.is_ZZ and not K.is_QQ:
        raise _DomainError("isolation of real and complex roots is not supported over %s" % K)

    from .. import Dummy, Poly
    x = Dummy('x')
    _, factors = Poly(f, x, domain=K).sqf_list()

    if len(factors) == 1:
        (g, k), = factors

        real_part, complex_part = dup_isolate_all_roots_sqf(
            g.rep.to_list(), K, eps=eps, inf=inf, sup=sup, fast=fast)

        real_part = [ ((a, b), k) for (a, b) in real_part ]
        complex_part = [ ((a, b), k) for (a, b) in complex_part ]

        return real_part, complex_part
    else:
        raise NotImplementedError( "only trivial square-free polynomials are supported")


__all__ = ['RealInterval', 'ComplexInterval']
