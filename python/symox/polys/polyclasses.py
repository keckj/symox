"""SymPy's dense polynomial representation ``DMP`` over native domains.

``Poly.rep`` is a ``DMP_Python`` holding the nonzero terms of a polynomial as
domain elements, as SymPy's ``Poly.rep`` does. Accessors, conversions and
ring arithmetic follow ``sympy.polys.polyclasses.DMP``; arithmetic runs on
native polynomials over placeholder generators, so the coefficients stay
exact and keep their domain. Algorithms specific to SymPy's dense
representation (root isolation, factorization, resultants...) are not part
of this view: use the corresponding ``Poly`` methods.
"""
from .orderings import monomial_key
from .polyerrors import PolynomialError, UnificationFailed

__all__ = ['DMP', 'DMP_Python']


def _zero_list(lev):
    rep = []
    for _ in range(lev):
        rep = [rep]
    return rep


def _depth(rep):
    lev = 0
    while isinstance(rep, list) and rep and isinstance(rep[0], list):
        rep, lev = rep[0], lev + 1
    return lev


def _dense_terms(rep, lev, convert):
    """Nonzero ``{monom: coeff}`` of a dense nested list."""
    terms = {}
    if lev == 0:
        degree = len(rep) - 1
        for i, coeff in enumerate(rep):
            coeff = convert(coeff)
            if coeff:
                terms[(degree - i,)] = coeff
        return terms
    degree = len(rep) - 1
    for i, sub in enumerate(rep):
        for monom, coeff in _dense_terms(sub, lev - 1, convert).items():
            terms[(degree - i,) + monom] = coeff
    return terms


def _dense_list(terms, lev, zero):
    if not terms:
        return _zero_list(lev)
    degree = max(monom[0] for monom in terms)
    if lev == 0:
        return [terms.get((n,), zero) for n in range(degree, -1, -1)]
    rows = [{} for _ in range(degree + 1)]
    for monom, coeff in terms.items():
        rows[degree - monom[0]][monom[1:]] = coeff
    return [_dense_list(row, lev - 1, zero) for row in rows]


def _tuple(rep):
    return tuple(_tuple(item) for item in rep) if isinstance(rep, list) else rep


class DMP:
    """Dense Multivariate Polynomials over `K`. """

    __slots__ = ('_terms', 'dom', 'lev')

    def __new__(cls, rep, dom, lev=None):
        if lev is None:
            lev = _depth(rep)
        return cls._new(_dense_terms(rep, lev, lambda c: c), dom, lev)

    @classmethod
    def _new(cls, terms, dom, lev):
        obj = object.__new__(DMP_Python if cls is DMP else cls)
        obj._terms, obj.dom, obj.lev = terms, dom, lev
        return obj

    @classmethod
    def from_list(cls, rep, lev, dom):
        """Create an instance of ``cls`` given a list of native coefficients. """
        return cls._new(_dense_terms(rep, lev, dom.convert), dom, lev)

    @classmethod
    def from_sympy_list(cls, rep, lev, dom):
        """Create an instance of ``cls`` given a list of SymPy coefficients. """
        return cls._new(_dense_terms(rep, lev, dom.from_sympy), dom, lev)

    @classmethod
    def from_dict(cls, rep, lev, dom):
        """Construct and instance of ``cls`` from a ``dict`` representation. """
        terms = {}
        for monom, coeff in dict(rep).items():
            coeff = dom.convert(coeff)
            if coeff:
                terms[tuple(monom)] = coeff
        return cls._new(terms, dom, lev)

    @classmethod
    def from_monoms_coeffs(cls, monoms, coeffs, lev, dom):
        return cls.from_dict(dict(zip(monoms, coeffs)), lev, dom)

    @classmethod
    def from_poly(cls, poly):
        """The representation of a native ``Poly`` (``Poly.rep``). """
        dom = poly.domain
        lev = len(poly.gens) - 1
        terms = {}
        for monom, coeff in poly.as_dict().items():
            coeff = dom.from_sympy(coeff)
            if coeff:
                terms[tuple(monom)] = coeff
        return cls._new(terms, dom, lev)

    def per(f, terms, dom=None, lev=None):
        return f._new(terms, f.dom if dom is None else dom, f.lev if lev is None else lev)

    @property
    def rep(f):
        return f.to_list()

    @property
    def zero(f):
        return f._new({}, f.dom, f.lev)

    @property
    def one(f):
        return f.ground_new(f.dom.one)

    def ground_new(f, coeff):
        """Construct a new ground instance of ``f``. """
        coeff = f.dom.convert(coeff)
        return f.per({(0,)*(f.lev + 1): coeff} if coeff else {})

    def __repr__(f):
        return "%s(%s, %s)" % (type(f).__name__, f.to_list(), f.dom)

    def __hash__(f):
        return hash((type(f).__name__, f.to_tuple(), f.lev, f.dom))

    def __reduce__(f):
        return DMP.from_list, (f.to_list(), f.lev, f.dom)

    def to_list(f):
        """Convert ``f`` to a list representation with native coefficients. """
        return _dense_list(f._terms, f.lev, f.dom.zero)

    def to_tuple(f):
        """Convert ``f`` to a tuple representation with native coefficients. """
        return _tuple(f.to_list())

    def to_sympy_list(f):
        """Convert ``f`` to a list representation with SymPy coefficients. """
        to_sympy = f.dom.to_sympy
        return _dense_list({m: to_sympy(c) for m, c in f._terms.items()}, f.lev, to_sympy(f.dom.zero))

    def to_dict(f, zero=False):
        """Convert ``f`` to a dict representation with native coefficients. """
        if not f._terms and zero:
            return {(0,)*(f.lev + 1): f.dom.zero}
        return dict(f._terms)

    def to_sympy_dict(f, zero=False):
        """Convert ``f`` to a dict representation with SymPy coefficients. """
        return {m: f.dom.to_sympy(c) for m, c in f.to_dict(zero=zero).items()}

    def convert(f, dom):
        """Convert ``f`` to a ``DMP`` over the new domain. """
        if f.dom == dom:
            return f
        return f.per({m: dom.convert_from(c, f.dom) for m, c in f._terms.items()}, dom=dom)

    def to_ring(f):
        """Make the ground domain a ring. """
        return f.convert(f.dom.get_ring())

    def to_field(f):
        """Make the ground domain a field. """
        return f.convert(f.dom.get_field())

    def to_exact(f):
        """Make the ground domain exact. """
        return f.convert(f.dom.get_exact())

    def unify_DMP(f, g):
        """Unify and return ``DMP`` instances of ``f`` and ``g``. """
        if not isinstance(g, DMP) or f.lev != g.lev:
            raise UnificationFailed("Cannot unify %s with %s" % (f, g))
        if f.dom != g.dom:
            dom = f.dom.unify(g.dom)
            f, g = f.convert(dom), g.convert(dom)
        return f, g

    def terms(f, order=None):
        """Returns all non-zero terms from ``f`` in lex order. """
        if not f._terms:
            return [((0,)*(f.lev + 1), f.dom.zero)]
        key = monomial_key(order)
        return sorted(f._terms.items(), key=lambda term: key(term[0]), reverse=True)

    def coeffs(f, order=None):
        """Returns all non-zero coefficients from ``f`` in lex order. """
        return [c for _, c in f.terms(order=order)]

    def monoms(f, order=None):
        """Returns all non-zero monomials from ``f`` in lex order. """
        return [m for m, _ in f.terms(order=order)]

    def _univariate(f):
        if f.lev:
            raise PolynomialError('multivariate polynomials not supported')

    def all_coeffs(f):
        """Returns all coefficients from ``f``. """
        f._univariate()
        return f.to_list() or [f.dom.zero]

    def all_monoms(f):
        """Returns all monomials from ``f``. """
        f._univariate()
        n = len(f.to_list()) - 1
        return [(n - i,) for i in range(n + 1)] if n >= 0 else [(0,)]

    def all_terms(f):
        """Returns all terms from a ``f``. """
        return list(zip(f.all_monoms(), f.all_coeffs()))

    def LC(f):
        """Returns the leading coefficient of ``f``. """
        return f.terms()[0][1]

    def TC(f):
        """Returns the trailing coefficient of ``f``. """
        return f._terms.get((0,)*(f.lev + 1), f.dom.zero)

    def nth(f, *N):
        """Returns the ``n``-th coefficient of ``f``. """
        if not all(isinstance(n, int) for n in N):
            raise TypeError("a sequence of integers expected")
        if len(N) != f.lev + 1:
            raise ValueError("expected %d exponents, got %d" % (f.lev + 1, len(N)))
        return f._terms.get(tuple(N), f.dom.zero)

    def degree(f, j=0):
        """Returns the leading degree of ``f`` in ``x_j``. """
        if not isinstance(j, int):
            raise TypeError("``int`` expected, got %s" % type(j))
        if not 0 <= j <= f.lev:
            raise IndexError("0 <= j <= %s expected, got %s" % (f.lev, j))
        return max((m[j] for m in f._terms), default=float('-inf'))

    def degree_list(f):
        """Returns a list of degrees of ``f``. """
        return tuple(f.degree(j) for j in range(f.lev + 1))

    def total_degree(f):
        """Returns the total degree of ``f``. """
        return max((sum(m) for m in f._terms), default=float('-inf'))

    @property
    def is_zero(f):
        """Returns ``True`` if ``f`` is a zero polynomial. """
        return not f._terms

    @property
    def is_one(f):
        """Returns ``True`` if ``f`` is a unit polynomial. """
        return f._terms == {(0,)*(f.lev + 1): f.dom.one}

    @property
    def is_ground(f):
        """Returns ``True`` if ``f`` is an element of the ground domain. """
        return all(not any(m) for m in f._terms)

    @property
    def is_monomial(f):
        """Returns ``True`` if ``f`` is zero or has only one term. """
        return len(f._terms) <= 1

    @property
    def is_monic(f):
        """Returns ``True`` if the leading coefficient of ``f`` is one. """
        return f.LC() == f.dom.one

    @property
    def is_linear(f):
        """Returns ``True`` if ``f`` is linear in all its variables. """
        return all(sum(m) <= 1 for m in f._terms)

    @property
    def is_quadratic(f):
        """Returns ``True`` if ``f`` is quadratic in all its variables. """
        return all(sum(m) <= 2 for m in f._terms)

    def _poly(f):
        from .. import Dummy, Poly
        gens = tuple(Dummy('x') for _ in range(f.lev + 1))
        return Poly(f.to_sympy_dict(), *gens, domain=f.dom)

    def _from_native(f, poly):
        return f.from_poly(poly)

    def _binary(f, g, operation):
        if not isinstance(g, DMP):
            g = f.ground_new(g)
        F, G = f.unify_DMP(g)
        from .. import Poly
        p = F._poly()
        q = Poly(G.to_sympy_dict(), *p.gens, domain=G.dom)
        return F, operation(p, q)

    def neg(f):
        """Negate all coefficients in ``f``. """
        return f.per({m: -c for m, c in f._terms.items()})

    def add(f, g):
        """Add two multivariate polynomials ``f`` and ``g``. """
        F, result = f._binary(g, lambda p, q: p.add(q))
        return F._from_native(result)

    def sub(f, g):
        """Subtract two multivariate polynomials ``f`` and ``g``. """
        F, result = f._binary(g, lambda p, q: p.sub(q))
        return F._from_native(result)

    def mul(f, g):
        """Multiply two multivariate polynomials ``f`` and ``g``. """
        F, result = f._binary(g, lambda p, q: p.mul(q))
        return F._from_native(result)

    def sqr(f):
        """Square a multivariate polynomial ``f``. """
        return f.mul(f)

    def pow(f, n):
        """Raise ``f`` to a non-negative power ``n``. """
        if not isinstance(n, int):
            raise TypeError("``int`` expected, got %s" % type(n))
        if n < 0:
            raise ValueError("Cannot raise polynomial to a negative power")
        return f._from_native(f._poly()**n)

    def div(f, g):
        """Polynomial division with remainder of ``f`` and ``g``. """
        F, (q, r) = f._binary(g, lambda p, q: p.div(q, auto=False))
        return F._from_native(q), F._from_native(r)

    def rem(f, g):
        """Computes polynomial remainder of ``f`` and ``g``. """
        return f.div(g)[1]

    def quo(f, g):
        """Computes polynomial quotient of ``f`` and ``g``. """
        return f.div(g)[0]

    def exquo(f, g):
        """Computes polynomial exact quotient of ``f`` and ``g``. """
        F, result = f._binary(g, lambda p, q: p.exquo(q, auto=False))
        return F._from_native(result)

    def gcd(f, g):
        """Returns polynomial GCD of ``f`` and ``g``. """
        F, result = f._binary(g, lambda p, q: p.gcd(q))
        return F._from_native(result)

    def lcm(f, g):
        """Returns polynomial LCM of ``f`` and ``g``. """
        F, result = f._binary(g, lambda p, q: p.lcm(q))
        return F._from_native(result)

    def eq(f, g, strict=False):
        if not strict:
            return f == g
        return isinstance(g, DMP) and f.lev == g.lev and f.dom == g.dom and f._terms == g._terms

    def ne(f, g, strict=False):
        return not f.eq(g, strict=strict)

    def __eq__(f, g):
        if f is g:
            return True
        if not isinstance(g, DMP):
            return NotImplemented
        try:
            F, G = f.unify_DMP(g)
        except UnificationFailed:
            return False
        return F._terms == G._terms

    def __ne__(f, g):
        result = f.__eq__(g)
        return result if result is NotImplemented else not result

    def __bool__(f):
        return not f.is_zero

    def __neg__(f):
        return f.neg()

    def __add__(f, g):
        return f.add(g)

    def __radd__(f, g):
        return f.add(g)

    def __sub__(f, g):
        return f.sub(g)

    def __rsub__(f, g):
        return f.neg().add(g)

    def __mul__(f, g):
        return f.mul(g)

    def __rmul__(f, g):
        return f.mul(g)

    def __pow__(f, n):
        return f.pow(n)

    def __divmod__(f, g):
        return f.div(g)

    def __mod__(f, g):
        return f.rem(g)

    def __floordiv__(f, g):
        return f.quo(g)


class DMP_Python(DMP):
    """Dense Multivariate Polynomials over `K`. """

    __slots__ = ()
