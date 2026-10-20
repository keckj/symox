"""Isolation intervals of polynomial roots.

``RealInterval`` and ``ComplexInterval`` have SymPy's public interface:
endpoints in the fraction field of ``dom``, widths, centers, disjointness
tests and refinement. Endpoints come from native Sturm isolation, so they are
valid isolating intervals but are not always SymPy's continued-fraction
intervals, and ``args`` stores the endpoints instead of a Mobius transform.
Complex root isolation is not implemented natively; complex intervals can
be built and inspected but not refined.
"""
# Module-level imports: the routing facade rebinds ``__spec__`` of aliased
# modules, so relative imports inside functions would warn.
from .. import Dummy, Poly, Rational


def _field_element(field, value):
    value = Rational(value)
    return field(int(value.p), int(value.q))


def _rational(value):
    return Rational(int(value.numerator), int(value.denominator))


class RealInterval:
    """A fully qualified representation of a real isolation interval. """

    def __init__(self, data, f, dom):
        """Initialize new real interval with complete information. """
        if len(data) != 2:
            raise NotImplementedError('Mobius interval data is not supported; pass (s, t)')
        field = dom.get_field()
        s, t = (_field_element(field, value) for value in data)
        if t < s:
            raise ValueError('invalid isolation interval (%s, %s)' % (s, t))
        self.f, self.dom = list(f), dom
        self._bounds = (s, t)

    @property
    def func(self):
        return RealInterval

    @property
    def args(self):
        return (self._bounds, self.f, self.dom)

    def __eq__(self, other):
        if type(other) is not type(self):
            return False
        return self.args == other.args

    @property
    def a(self):
        """Return the position of the left end. """
        return self._bounds[0]

    @property
    def b(self):
        """Return the position of the right end. """
        return self._bounds[1]

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
        """Say whether a complex number, or a pair (re, im), belongs to this interval. """
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

    def _poly(self):
        x = Dummy('x')
        degree = len(self.f) - 1
        return Poly(sum((_rational(c)*x**(degree - i) for i, c in enumerate(self.f)), 0), x)

    def _inner_refine(self):
        """Internal one step real root refinement procedure. """
        if self.a == self.b:
            return self
        s, t = self._poly().refine_root(_rational(self.a), _rational(self.b), steps=1)
        return RealInterval((s, t), self.f, self.dom)

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

    The printed form is ``(ax, bx) x (ay, by)``, where ``(ax, ay)`` and
    ``(bx, by)`` are the southwest and northeast corners of the rectangle.
    ``conj`` marks the conjugate root, whose rectangle is mirrored.
    """

    def __init__(self, a, b, I, Q, F1, F2, f1, f2, dom, conj=False):
        """Initialize new complex interval with complete information. """
        self.a, self.b = a, b  # the southwest and northeast corner: (x1, y1), (x2, y2)
        self.I, self.Q = I, Q

        self.f1, self.F1 = f1, F1
        self.f2, self.F2 = f2, F2

        self.dom = dom
        self.conj = conj

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
        return ComplexInterval(self.a, self.b, self.I, self.Q,
            self.F1, self.F2, self.f1, self.f2, self.dom, conj=True)

    def __contains__(self, item):
        """Say whether a complex number, or a pair (re, im), belongs to this rectangle. """
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
        raise NotImplementedError('complex root isolation is not implemented natively')

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


__all__ = ['RealInterval', 'ComplexInterval']
