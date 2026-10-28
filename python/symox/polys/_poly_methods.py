"""SymPy ``Poly`` methods layered over the native polynomial type.

The native ``Poly`` stores exact sparse representations in Rust. The methods
installed here follow ``sympy.polys.polytools.Poly`` for the parts of its API
that are thin compositions of native operations: order-aware accessors,
ground-domain arithmetic, alternative constructors, generator injection, and
the ``Basic``-level protocol (``args``, ``atoms``, ``has``, ``compare``, ...).
Every result is rebuilt through the native constructor with an explicit
domain, so coefficients are converted exactly or the conversion fails with
the domain's exception.
"""
from math import ceil, log2

from .. import _core
from .polyerrors import (
    CoercionFailed, DomainError, ExactQuotientFailed, GeneratorsNeeded,
    MultivariatePolynomialError, NotReversible, PolynomialError,
    UnificationFailed,
)
from .monomials import Monomial
from .orderings import monomial_key

Poly = _core.Poly

_native = {name: getattr(Poly, name) for name in (
    "terms", "coeffs", "monoms", "LC", "LM", "as_list", "all_roots", "__eq__",
    "TC", "nth",
)}
_native_reduce=Poly.__reduce__
_native_as_dict=Poly.as_dict
_stock_poly_dict=Poly.__dict__
_native_ex_domain=_core.EX


def _stock_coefficient_view(f):
    if _stock_poly_dict.get('as_dict') is not _native_as_dict:return False
    if type(f) is Poly:return True
    from .. import _purepoly
    return type(f) is _purepoly.PurePoly and _purepoly.inspect.getattr_static(_purepoly.PurePoly,'as_dict') is _purepoly._stock_as_dict


_SAME = object()


def _restore_polynomial(cls,terms,gens,domain):return cls(terms,*gens,domain=domain)


def _reduce_polynomial(f):
    try:domain=f.domain
    except NotImplementedError:return _native_reduce(f)
    terms=f._coefficient_payload() if _stock_coefficient_view(f) and domain is _native_ex_domain else f.as_dict()
    return _restore_polynomial,(type(f),terms,f.gens,domain)


def _S():
    from .. import S
    return S


def _sympify(value):
    from .. import sympify
    return sympify(value)


def _is_poly(value):
    return isinstance(value, Poly)


def _domain(f):
    """The typed domain of ``f`` or ``None`` for legacy coefficient backends."""
    try:
        return f.domain
    except NotImplementedError:
        return None


def _domain_kwargs(f):
    if f.modulus is not None:
        return {"modulus": f.modulus}
    K = _domain(f)
    return {} if K is None else {"domain": K}


def _monomial_expr(gens, monom):
    from .. import Mul
    return Mul(*[g**e for g, e in zip(gens, monom)])


def _expr_from_dict(rep, gens):
    from .. import Add
    return Add(*[c*_monomial_expr(gens, m) for m, c in rep.items()])


def _from_dict(f, rep, gens=None, domain=_SAME):
    """A Poly with the terms ``rep`` over ``gens`` (default: those of ``f``)."""
    gens = tuple(f.gens if gens is None else gens)
    if domain is _SAME:
        kwargs = _domain_kwargs(f)
    elif domain is None:
        kwargs = {}
    else:
        kwargs = {"domain": domain}
    if all(g.is_Symbol for g in gens) and not any(
            getattr(c, 'has', None) is not None and not getattr(c, 'is_Rational', False)
            and c.has(*gens) for c in rep.values()):
        return Poly(_expr_from_dict(rep, gens), *gens, **kwargs)
    # Coefficients such as ``CRootOf(x**2 - 2, 0)`` hold a generator name
    # without depending on it: build over placeholders, as the reference
    # works on the representation.
    # Expression generators such as ``y + 1`` are not re-read from an
    # expression (``5*(y + 1)`` would distribute): build over placeholders.
    from .. import Dummy
    slots = tuple(Dummy('g') for _ in gens)
    return Poly(_expr_from_dict(rep, slots), *slots, **kwargs)._with_gens(*gens)


def _items(f):
    return [(m, c) for m, c in f.as_dict().items() if c != 0]


def _convert_ground(f, c):
    """Convert ``c`` into the ground domain of ``f`` (``CoercionFailed`` if impossible)."""
    c = _sympify(c)
    K = _domain(f)
    if K is None:
        if c.free_symbols & set(f.gens):
            raise CoercionFailed("cannot convert %s into the ground domain" % c)
        return c
    if K.is_RR:
        if not c.is_number or c.is_real is not True:
            raise CoercionFailed("expected a real number, got %s" % c)
        from .. import Float
        return Float(c)
    try:
        return K.to_sympy(K.convert(c))
    except CoercionFailed:
        if K.is_ZZ:
            raise CoercionFailed("expected an integer, got %s" % c) from None
        raise


# --- ordered accessors -------------------------------------------------------

def _rr_zero(f, c):
    """The reference keeps ``RR`` zeros as ``0.0`` (``Poly(0, x, domain=RR).LC()``)."""
    if c == 0 and not getattr(c, 'is_Float', False):
        K = _domain(f)
        if K is not None and K.is_RR:
            from .. import Float
            return Float(0, precision=K.precision)
    return c


def _terms(f, order=None):
    terms = _native["terms"](f)
    if len(terms) == 1 and terms[0][1] == 0:
        terms = [(terms[0][0], _rr_zero(f, terms[0][1]))]
    if order is None or (isinstance(order, str) and order == "lex"):
        return terms
    key = monomial_key(order)
    return sorted(terms, key=lambda t: key(t[0]), reverse=True)


def terms(f, order=None):
    """Nonzero terms of ``f`` as ``(monom, coeff)`` pairs sorted by ``order``."""
    return _terms(f, order)


def coeffs(f, order=None):
    """Nonzero coefficients of ``f`` sorted by ``order``."""
    return [c for _, c in _terms(f, order)]


def monoms(f, order=None):
    """Exponent tuples of the nonzero terms of ``f`` sorted by ``order``."""
    return [m for m, _ in _terms(f, order)]


def LC(f, order=None):
    """Leading coefficient with respect to ``order``."""
    if order is None:
        return _rr_zero(f, _native["LC"](f))
    return coeffs(f, order)[0]


def TC(f):
    """Trailing coefficient (``0.0`` for zero over ``RR``)."""
    return _rr_zero(f, _native["TC"](f))


def nth(f, *N):
    """Coefficient of the monomial with exponents ``N``."""
    return _rr_zero(f, _native["nth"](f, *N))


def EC(f, order=None):
    """Last nonzero coefficient with respect to ``order``."""
    return coeffs(f, order)[-1]


def LM(f, order=None):
    """Leading monomial as a ``Monomial`` with respect to ``order``."""
    return Monomial(monoms(f, order)[0], f.gens)


def EM(f, order=None):
    """Last monomial as a ``Monomial`` with respect to ``order``."""
    return Monomial(monoms(f, order)[-1], f.gens)


def LT(f, order=None):
    """Leading term ``(Monomial, coeff)`` with respect to ``order``."""
    monom, coeff = _terms(f, order)[0]
    return Monomial(monom, f.gens), coeff


def ET(f, order=None):
    """Last term ``(Monomial, coeff)`` with respect to ``order``."""
    monom, coeff = _terms(f, order)[-1]
    return Monomial(monom, f.gens), coeff


def as_list(f, native=False):
    """Nested dense coefficient lists, leading coefficient first."""
    rep = f._coefficient_payload() if _stock_coefficient_view(f) else dict(_items(f))
    n = len(f.gens)

    def build(level, prefix):
        if level == n:
            return rep.get(prefix, 0)
        tails = [m[level] for m in rep if m[:level] == prefix]
        if not tails:
            return [] if level == n - 1 else _zero_nest(n - level - 1)
        return [build(level + 1, prefix + (k,)) for k in range(max(tails), -1, -1)]

    def _zero_nest(depth):
        out = []
        for _ in range(depth):
            out = [out]
        return out

    if not rep:
        return _zero_nest(n - 1)
    return _strip(build(0, ()), n - 1)


def _strip(lst, depth):
    """Normalise inner zero polynomials as SymPy's dense representation does."""
    if depth == 0:
        return lst
    out = []
    for item in lst:
        out.append(_strip_zero(item, depth - 1))
    return out


def _strip_zero(item, depth):
    if depth == 0:
        while item and item[0] == 0:
            item = item[1:]
        return item
    stripped = [_strip_zero(i, depth - 1) for i in item]
    while stripped and _is_dense_zero(stripped[0]):
        stripped = stripped[1:]
    if not stripped:
        out = []
        for _ in range(depth):
            out = [out]
        return out
    return stripped


def _is_dense_zero(item):
    if isinstance(item, list):
        return all(_is_dense_zero(i) for i in item) if item else True
    return item == 0


# --- arithmetic --------------------------------------------------------------

def _ground_result(f, rep):
    return _from_dict(f, rep)


def add_ground(f, coeff):
    """Add an element of the ground domain to ``f``."""
    c = _convert_ground(f, coeff)
    return _ground_result(f, _add_constant(f, c))


def sub_ground(f, coeff):
    """Subtract an element of the ground domain from ``f``."""
    c = _convert_ground(f, coeff)
    return _ground_result(f, _add_constant(f, -c))


def _add_constant(f, c):
    rep = dict(_items(f))
    zero = (0,)*len(f.gens)
    rep[zero] = rep.get(zero, 0) + c
    return rep


def mul_ground(f, coeff):
    """Multiply ``f`` by an element of the ground domain."""
    c = _convert_ground(f, coeff)
    return _ground_result(f, {m: v*c for m, v in _items(f)})


def _quotient_str(K, a):
    """``a`` as the reference ``ExactQuotientFailed`` prints it: Gaussian
    elements (also inside ``ZZ_I[y]``) print as ``(x + y*I)``."""
    from .domains import ZZ_I, QQ_I
    if K == ZZ_I or K == QQ_I:
        a = K.convert(a)
        return '(%s + %s*I)' % (a.x, a.y)
    if K.is_Composite and K.dom in (ZZ_I, QQ_I):
        return str(K.convert(a))
    return str(a)


def _ground_div(f, coeff, exact):
    c = _convert_ground(f, coeff)
    if c == 0:
        raise ZeroDivisionError("polynomial division")
    K = _domain(f)
    out = {}
    for m, v in sorted(_items(f), reverse=True):
        if K is None:
            q = v/c
        elif K.is_RR:
            q = v/c
        elif exact:
            if K.is_Field and not K.is_Composite:
                q = K.to_sympy(K.quo(K.convert(v), K.convert(c)))
            else:
                try:
                    q = K.to_sympy(K.exquo(K.convert(v), K.convert(c)))
                except ExactQuotientFailed:
                    raise ExactQuotientFailed("%s does not divide %s in %s" % (
                        _quotient_str(K, c), _quotient_str(K, v), K)) from None
        else:
            q = K.to_sympy(K.quo(K.convert(v), K.convert(c)))
        out[m] = q
    return _ground_result(f, out)


def quo_ground(f, coeff):
    """Quotient of ``f`` by a ground element (floor division in a ring)."""
    return _ground_div(f, coeff, False)


def exquo_ground(f, coeff):
    """Exact quotient of ``f`` by a ground element."""
    return _ground_div(f, coeff, True)


def abs_(f):
    """Make all coefficients of ``f`` positive."""
    K = _domain(f)
    if K is not None and not K.is_Numerical and not K.is_FiniteField and not K.is_EX:
        raise TypeError("bad operand type for abs(): 'PolyElement'")
    return _ground_result(f, {m: abs(v) for m, v in _items(f)})


def neg(f):
    """Negate all coefficients of ``f``."""
    return -f


def add(f, g):
    """Sum of two polynomials (a non-Poly ``g`` is a ground element)."""
    g = _sympify(g) if not _is_poly(g) else g
    if not _is_poly(g):
        return f.add_ground(g)
    return f + g


def sub(f, g):
    """Difference of two polynomials (a non-Poly ``g`` is a ground element)."""
    g = _sympify(g) if not _is_poly(g) else g
    if not _is_poly(g):
        return f.sub_ground(g)
    return f - g


def mul(f, g):
    """Product of two polynomials (a non-Poly ``g`` is a ground element)."""
    g = _sympify(g) if not _is_poly(g) else g
    if not _is_poly(g):
        return f.mul_ground(g)
    return f*g


def sqr(f):
    """Square of ``f``."""
    return f*f


def pow_(f, n):
    """``f`` raised to a non-negative integer power ``n``."""
    n = int(n)
    if n < 0:
        raise ValueError("Cannot raise polynomial to a negative power")
    return f**n


def eq(f, g, strict=False):
    """``f == g``; with ``strict`` also compare the ground domains exactly."""
    if not strict:
        return f == g
    return _is_poly(g) and f.gens == g.gens and f == g


def ne(f, g, strict=False):
    return not f.eq(g, strict=strict)


def _abs(f):
    return f.abs()


# --- constructors ------------------------------------------------------------

def _options(gens, args):
    if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
        gens = tuple(gens[0])
    gens = tuple(_sympify(g) for g in gens)
    kwargs = {k: v for k, v in args.items() if k in ("domain", "modulus", "extension")}
    if "gens" in args:
        gens = tuple(args["gens"])
    return gens, kwargs


def from_dict(cls, rep, *gens, **args):
    """Construct a polynomial from a ``{monom: coeff}`` dictionary."""
    return cls(dict(rep),*gens,**args)


def from_list(cls, rep, *gens, **args):
    """Construct a univariate polynomial from dense coefficients."""
    return cls(list(rep),*gens,**args)


def from_poly(cls, rep, *gens, **args):
    """Construct a polynomial from another polynomial."""
    return cls(rep,*gens,**args)


def from_expr(cls, rep, *gens, **args):
    """Construct a polynomial from an expression."""
    return cls(_sympify(rep),*gens,**args)


def new(cls, rep, *gens):
    """Construct a polynomial from a Poly representation and generators."""
    if not _is_poly(rep):
        raise PolynomialError("invalid polynomial representation: %s" % (rep,))
    if len(rep.gens) != len(gens):
        raise PolynomialError("invalid number of generators for the representation")
    rename = dict(zip(rep.gens, gens))
    return _from_dict(rep, {m: c for m, c in _items(rep)}, gens=tuple(rename[g] for g in rep.gens))


def per(f, rep, gens=None, remove=None):
    """Create a polynomial over the generators of ``f`` from a representation.

    ``rep`` is a Poly with the same number of generators. With ``remove`` the
    generator at that index is dropped; a constant result is returned as a
    ground expression.
    """
    if gens is None:
        gens = f.gens
    gens = tuple(gens)
    from .polyclasses import DMP
    if isinstance(rep, DMP):
        # The reference ``f.__class__.new(rep, *gens)`` over ``rep.dom``.
        if remove is not None:
            gens = gens[:remove] + gens[remove + 1:]
        if len(gens) != rep.lev + 1:
            raise PolynomialError("invalid number of generators: expected %s, got %s"
                                  % (rep.lev + 1, len(gens)))
        return type(f).from_dict(rep.to_sympy_dict(), *gens, domain=rep.dom)
    if not _is_poly(rep):
        raise PolynomialError("invalid polynomial representation: %s" % (rep,))
    if len(rep.gens) != len(f.gens):
        raise PolynomialError("invalid number of generators for the representation")
    items = _items(rep)
    if remove is not None:
        gens = gens[:remove] + gens[remove + 1:]
        if any(m[remove] for m, _ in items):
            raise PolynomialError("cannot remove a generator that occurs in the representation")
        items = [(m[:remove] + m[remove + 1:], c) for m, c in items]
        if not gens:
            return sum((c for _, c in items), _S().Zero)
    return _from_dict(rep, dict(items), gens=gens)


def unify(f, g):
    """Make ``f`` and ``g`` belong to the same domain over the same generators."""
    if not _is_poly(g):
        g = _sympify(g)
        try:
            c = _convert_ground(f, g)
        except CoercionFailed:
            raise UnificationFailed("Cannot unify %s with %s" % (f, g)) from None
        return f, _from_dict(f, {(0,)*len(f.gens): c})
    s = f + g
    kwargs = _domain_kwargs(s)
    return (Poly(f.as_expr(), *s.gens, **kwargs), Poly(g.as_expr(), *s.gens, **kwargs))


def one(f):
    return _from_dict(f, {(0,)*len(f.gens): 1})


def zero(f):
    return _from_dict(f, {})


def free_symbols_in_domain(f):
    K = _domain(f)
    symbols = set()
    if K is None:
        for c in f.coeffs():
            symbols |= c.free_symbols
        return symbols
    if K.is_Composite:
        for gen in K.symbols:
            symbols |= gen.free_symbols
    elif K.is_FiniteField or K.is_Numerical or K.is_AlgebraicField:
        return set()
    return symbols


def coeff(f, x, n=1, right=False):
    raise NotImplementedError(
        "Either convert to Expr with `as_expr` method "
        "to use Expr's coeff method or else use the "
        "`coeff_monomial` method of Polys.")


def as_content_primitive(f, radical=False, clear=True):
    return _S().One, f


# --- domain changes ----------------------------------------------------------

def retract(f, field=None):
    """Recompute the ground domain of ``f`` from its coefficients."""
    if f.modulus is not None:
        return f
    g = Poly(f.as_expr(), *f.gens)
    if field:
        K = _domain(g)
        if K is not None and not K.is_Field:
            return g.to_field()
    return g


def to_exact(f):
    """Convert a polynomial over ``RR`` to one over ``QQ``."""
    K = _domain(f)
    if K is None or not K.is_RR:
        return f
    from fractions import Fraction
    from .. import Rational
    limit = int(1/(100*2**(1 - 53)))
    rep = {}
    for m, c in _items(f):
        q = Fraction(float(c)).limit_denominator(limit)
        rep[m] = Rational(q.numerator, q.denominator)
    return _from_dict(f, rep, domain="QQ")


def set_modulus(f, modulus):
    """Reduce the coefficients of ``f`` modulo a new modulus."""
    from .. import Integer
    modulus = int(modulus)
    old = f.modulus
    rep = {}
    for m, c in _items(f):
        if old is not None:
            c = Integer(int(c) % old)
        elif not getattr(c, "is_integer", False):
            K = _domain(f)
            raise CoercionFailed("Cannot convert %s from %s to GF(%s)" % (c, K, modulus))
        rep[m] = c
    return Poly(_expr_from_dict(rep, f.gens), *f.gens, modulus=modulus)


def inject(f, front=False):
    """Move the generators of a polynomial-ring domain into ``f.gens``."""
    K = _domain(f)
    if K is None:
        raise DomainError("Cannot inject generators over %s" % f.domain_descriptor)
    if K.is_Numerical or K.is_FiniteField:
        return f
    if not K.is_PolynomialRing:
        raise DomainError("Cannot inject generators over %s" % K)
    symbols = tuple(K.symbols)
    gens = symbols + tuple(f.gens) if front else tuple(f.gens) + symbols
    return Poly(f.as_expr(), *gens, domain=K.dom)


def eject(f, *gens):
    """Move the leading or trailing generators of ``f`` into the domain."""
    K = _domain(f)
    if K is None or not K.is_Numerical:
        raise DomainError("Cannot eject generators over %s" % (K if K is not None else f.domain_descriptor))
    k = len(gens)
    if f.gens[:k] == gens:
        rest = f.gens[k:]
    elif f.gens[-k:] == gens:
        rest = f.gens[:-k]
    else:
        raise NotImplementedError("can only eject front or back generators")
    return Poly(f.as_expr(), *rest, domain=K.poly_ring(*gens))


def rat_clear_denoms(f, g):
    """Clear denominators of the rational function ``f/g`` (keeping its value)."""
    if not _is_poly(g):
        raise UnificationFailed("Cannot unify %s with %s" % (f, g))
    f, g = f.unify(g)
    K = _domain(f)
    if K is None or not K.is_Field or not (K.is_QQ or K.is_FractionField):
        return f, g
    a, f = f.clear_denoms()
    b, g = g.clear_denoms()
    R = K.get_ring()
    f = Poly(f.as_expr(), *f.gens, domain=R)
    g = Poly(g.as_expr(), *g.gens, domain=R)
    return f.mul_ground(b), g.mul_ground(a)


# --- algebraic operations ----------------------------------------------------

def _gaussian(f):
    """True when ``f`` has Gaussian rational coefficients (legacy ``ZZ[I]``)."""
    from .. import I, im, re
    for c in f.coeffs():
        if not (re(c).is_rational and im(c).is_rational):
            return False
    return any(im(c) != 0 for c in f.coeffs())


def _norm(f):
    """Norm of ``f`` over an algebraic field (product of its conjugates)."""
    if len(f.gens) == 1:
        return f.norm()
    from .. import Dummy, minimal_polynomial
    K = _domain(f)
    n = Poly(minimal_polynomial(K.ext, Dummy("t")), domain="QQ").degree()
    degs = f.degree_list()
    base = [n*int(d) + 1 for d in degs]
    weights = []
    w = 1
    for b in base:
        weights.append(w)
        w *= b
    z = Dummy("z")
    rep = {}
    for m, c in _items(f):
        e = sum(k*wt for k, wt in zip(m, weights))
        rep[e] = rep.get(e, 0) + c
    from .. import Add
    uni = Poly(Add(*[c*z**e for e, c in rep.items()]), z, domain=K)
    r = uni.norm()
    out = {}
    for (e,), c in _items(r):
        monom = []
        for b in base:
            monom.append(e % b)
            e //= b
        out[tuple(monom)] = c
    return _from_dict(f, out, domain="QQ")


def lift(f):
    """Convert algebraic coefficients to rationals (the norm of ``f``)."""
    K = _domain(f)
    ring = False
    if K is not None and (getattr(K, "is_GaussianRing", False) or getattr(K, "is_GaussianField", False)):
        # The reference dmp_lift converts ZZ_I/QQ_I to QQ<I>; the norm stays
        # over ZZ for the ring.
        from .. import I, QQ
        ring = bool(getattr(K, "is_GaussianRing", False))
        f = Poly(f.as_expr(), *f.gens, domain=QQ.algebraic_field(I))
        K = f.domain
    elif K is None and _gaussian(f):
        from .. import I
        ring = all(c.is_integer for c in sum(((c.as_real_imag()) for c in f.coeffs()), ()))
        f = Poly(f.as_expr(), *f.gens, extension=I)
        K = f.domain
    if K is None or not K.is_AlgebraicField:
        raise DomainError("computation can be done only in an algebraic domain")
    r = _norm(f)
    return Poly(r.as_expr(), *r.gens, domain="ZZ" if ring else "QQ")


def sqf_norm(f):
    """``(s, g, r)`` with ``g = f(x - s*a)`` and ``r = norm(g)`` square-free."""
    K = _domain(f)
    if K is None or not K.is_AlgebraicField:
        raise DomainError("ground domain must be algebraic")
    a = K.to_sympy(K.unit)
    for s in _sqf_norm_shifts(f):
        rule = {x: x - k*a for x, k in zip(f.gens, s) if k}
        g = Poly(f.as_expr().xreplace(rule), *f.gens, domain=K) if rule else f
        r = _norm(g)
        if r.is_sqf:
            return list(s), g, r


def _sqf_norm_shifts(f):
    """Shift vectors in the order of SymPy's ``dup_sqf_norm``/``dmp_sqf_norm``."""
    n = len(f.gens)
    if n == 1:
        s = 0
        while True:
            yield (s,)
            s += 1
    yield (0,)*n
    degrees = f.degree_list()
    for _, i in sorted((int(d), i) for i, d in enumerate(degrees) if d > 0):
        s = [0]*n
        s[i] = 1
        yield tuple(s)
    j = 0
    while True:
        j += 1
        yield (j,)*n


def revert(f, n):
    """Series reversion: ``g`` with ``f*g = 1 mod x**n``."""
    if len(f.gens) != 1:
        raise ValueError("univariate polynomial expected")
    n = int(n)
    if n <= 0:
        raise ValueError("math domain error")
    K = _domain(f)
    tc = f.TC()
    field = f.modulus is not None or (K is not None and K.is_Field)
    if field:
        if tc == 0:
            raise NotReversible("zero is not reversible")
        inv = 1/tc
    else:
        if tc != 1 and tc != -1:
            raise NotReversible("only units are reversible in a ring")
        inv = tc
    x = f.gen
    kwargs = _domain_kwargs(f)
    g = Poly(inv, x, **kwargs)
    h = Poly(x**2, x, **kwargs)
    for _ in range(int(ceil(log2(n)))):
        g = (g*2 - f*g*g).rem(h)
        h = h*Poly(x**h.degree(), x, **kwargs)
    return g


def nth_power_roots_poly(f, n):
    """Polynomial whose roots are the ``n``-th powers of the roots of ``f``."""
    if f.is_multivariate:
        raise MultivariatePolynomialError("must be a univariate polynomial")
    N = _sympify(n)
    if N.is_Integer and N >= 1:
        n = int(N)
    else:
        raise ValueError("'n' must an integer and n >= 1, got %s" % n)
    from .. import Dummy
    x = f.gen
    t = Dummy("t")
    r = f.resultant(Poly(x**n - t, x, t))
    return Poly(r.as_expr().xreplace({t: x}), x, **_domain_kwargs(f))


def shift_list(f, a):
    """Shift each generator ``x_i`` of ``f`` to ``x_i + a[i]``."""
    rule = {g: g + _sympify(s) for g, s in zip(f.gens, a)}
    return Poly(f.as_expr().xreplace(rule), *f.gens, **_domain_kwargs(f))


def real_roots(f, multiple=True, radicals=True):
    """Real roots of ``f`` in increasing order (grouped when ``multiple`` is False)."""
    from .._public_roots import CRootOf
    reals = CRootOf.real_roots(f, radicals=radicals)
    return reals if multiple else _group(reals)


def all_roots(f, multiple=True, radicals=True):
    """All complex roots of ``f``: sorted real roots, then the non-real roots
    factor by factor (grouped when ``multiple`` is False)."""
    from .._public_roots import CRootOf
    roots = CRootOf.all_roots(f, radicals=radicals)
    return roots if multiple else _group(roots)


def _group(roots):
    out = []
    for r in roots:
        if out and out[-1][0] == r:
            out[-1] = (r, out[-1][1] + 1)
        else:
            out.append((r, 1))
    return out


def root(f, index, radicals=True):
    """The ``index``-th root of ``f`` in the order of ``all_roots``."""
    from .._public_roots import rootof
    return rootof(f, index, radicals=radicals)


def same_root(f, a, b):
    """Decide whether approximations ``a`` and ``b`` belong to the same root of ``f``."""
    from .. import Rational, N, re, im
    if f.is_multivariate:
        raise MultivariatePolynomialError("Must be a univariate polynomial")
    n = f.degree()
    if n < 2:
        raise PolynomialError("Polynomials of degree < 2 have no distinct roots.")
    D = abs(f.discriminant())
    l2sq = sum(c**2 for c in f.all_coeffs())
    delta_sq = Rational(3)*D/(Rational(n)**(n + 1)*l2sq**(n - 1))
    eps_sq = delta_sq/9
    inv = 1/eps_sq
    bits = (int(inv.p) // int(inv.q)).bit_length() + 1
    m = bits//2 + bits % 2
    digits = int(m*0.30103) + 15

    def ev(z):
        v = N(_sympify(z), digits)
        return Rational(re(v)), Rational(im(v))

    ar, ai = ev(a)
    br, bi = ev(b)
    return bool((ar - br)**2 + (ai - bi)**2 < eps_sq)


def dispersionset(f, g=None):
    from .dispersion import dispersionset as _ds
    return _ds(f, g)


def dispersion(f, g=None):
    from .dispersion import dispersion as _d
    return _d(f, g)


def trunc(f, p):
    """Reduce the integer coefficients of ``f`` modulo ``p`` (symmetric residues)."""
    K = _domain(f)
    p = _convert_ground(f, p)
    if f.modulus is None and K is not None and K.is_RR:
        # The reference reduces each element with the ``mpf`` floor
        # modulo (no symmetric residues) at the precision of the domain.
        import mpmath
        from .. import Float
        prec = K.precision
        with mpmath.workprec(prec):
            pm = mpmath.mpf(_sympify(p)._mpf_)
            rep = {}
            for m, c in _items(f):
                r = mpmath.mpf(_sympify(c)._mpf_) % pm
                if r:
                    rep[m] = Float(r, precision=prec)
        return _from_dict(f, rep)
    if f.modulus is None and (K is None or not K.is_ZZ):
        # The reference applies ``%`` to the domain elements; EX elements
        # are ``Expression`` objects.
        name = 'Expression' if K is not None and K.is_EX else K
        raise TypeError("unsupported operand type(s) for %%: '%s' and '%s'" % (name, name))
    from .. import Integer
    p = int(p)
    rep = {}
    for m, c in _items(f):
        c = int(c) % p
        if c > p // 2:
            c -= p
        rep[m] = Integer(c)
    return _from_dict(f, rep)


def is_primitive(f):
    """True when the content of ``f`` is one."""
    return f.content() == 1


def is_cyclotomic(f):
    """True when ``f`` is a cyclotomic polynomial (SymPy's ``dup_cyclotomic_p``,
    native in ``symox::poly::structure::cyclotomic_p``)."""
    from .._core import _is_cyclotomic
    return _is_cyclotomic(f)


def make_monic_over_integers_by_scaling_roots(f):
    """``(g, c)`` with ``g`` monic over ``ZZ`` whose roots are ``c`` times the
    roots of ``f`` (``f`` univariate over ``ZZ`` or ``QQ``)."""
    K = _domain(f)
    if not f.is_univariate or f.modulus is not None or K is None or not (K.is_ZZ or K.is_QQ):
        raise ValueError('Polynomial must be univariate over ZZ or QQ.')
    if K.is_ZZ and f.LC() == 1:
        return f, 1
    fm = f.monic()
    c, _ = fm.clear_denoms()
    return fm.transform(Poly(fm.gen), c).to_ring(), c


def _which_roots(f, candidates, num_roots):
    """SymPy's ``Poly._which_roots``: drop candidates ``r`` with ``f(r)``
    certifiably non-zero at doubling precisions until at most ``num_roots``
    remain (candidates are deduplicated in order)."""
    from .. import Abs, Integer, N
    counts = list(dict.fromkeys(_sympify(r) for r in candidates))
    expr, prec = f.as_expr(), 10
    while len(counts) > num_roots and prec <= 1280:
        keep = []
        for r in counts:
            value = expr.subs(f.gen, r)
            if value != 0:
                # A value whose magnitude survives a doubled working precision
                # above 10**-prec is not a rounding residue of a root.
                # Evaluate first: Abs of an unevaluated expression in
                # CRootOf nodes expands into conjugate products that N
                # leaves unevaluated.
                size = Abs(N(value, 2*prec))
                if size.is_zero is False and size > Integer(10)**(-prec):
                    continue
            keep.append(r)
        counts = keep
        prec *= 2
    return counts


def which_real_roots(f, candidates):
    """The candidates that are real roots of ``f`` (SymPy's heuristic)."""
    if f.is_multivariate:
        raise MultivariatePolynomialError("Must be a univariate polynomial")
    K = _domain(f)
    if f.modulus is not None or K is None or not (K.is_ZZ or K.is_QQ or K.is_AlgebraicField):
        raise NotImplementedError("root counting not supported over %s" % K)
    return _which_roots(f, candidates, f.count_roots())


def which_all_roots(f, candidates):
    """The candidates that are (complex) roots of ``f``."""
    if f.is_multivariate:
        raise MultivariatePolynomialError("Must be a univariate polynomial")
    return _which_roots(f, candidates, f.degree())


def transform(f, p, q):
    """``q**deg(f) * f(p/q)`` for univariate ``f``."""
    if f.is_multivariate:
        raise MultivariatePolynomialError("must be a univariate polynomial")
    p, q = f.unify(p)[1], f.unify(q)[1]
    n = f.degree()
    result = f.zero
    for (k,), c in _items(f):
        result = result + p**k*q**(n - k)*c
    return result


# --- Basic protocol ----------------------------------------------------------

def _args(f):
    return (f.as_expr(),) + tuple(f.gens)


def _func(f):
    return Poly


def class_key(cls):
    return 5, 0, cls.__name__


def sort_key(f, order=None):
    one = _S().One
    args = _args(f)
    inner = tuple(a.sort_key(order) if hasattr(a, "sort_key") else a for a in args)
    return Poly.class_key(), (len(args), inner), one.sort_key(), one


def _preorder(e):
    yield e
    if _is_poly(e):
        for a in _args(e):
            yield from _preorder(a)
        return
    for a in getattr(e, "args", ()):
        yield from _preorder(a)


def atoms(f, *types):
    out = set()
    if types:
        types = tuple(t if isinstance(t, type) else type(t) for t in types)
        if Poly in types and isinstance(f, types):
            out.add(f)
        for a in _args(f):
            out |= a.atoms(*types)
        return out
    for a in _args(f):
        out |= a.atoms()
    return out


def has(f, *patterns):
    for p in patterns:
        if isinstance(p, type):
            if p is Poly or any(_matches_type(e, p) for e in _preorder(f)):
                return True
            continue
        if _is_poly(p):
            if p == f and p.gens == f.gens:
                return True
            continue
        if any(a.has(p) for a in _args(f)):
            return True
    return False


def _matches_type(e, t):
    try:
        return isinstance(e, t)
    except TypeError:
        return False


def _is_wild(p):
    from .. import Symbol, srepr
    try:
        if srepr(p).startswith("Wild("):
            return True
        return bool(p.is_Symbol) and str(p) == p.name + "_" and p != Symbol(p.name)
    except Exception:
        return False


def match(f, pattern, old=False):
    pattern = _sympify(pattern)
    if _is_wild(pattern):
        return f.as_expr().match(pattern)
    return None


def find(f, query, group=False):
    if isinstance(query, type):
        results = [e for e in _preorder(f) if _matches_type(e, query)]
    elif _is_wild(query):
        results = [e for e in _preorder(f)
                   if (e.match(query) if not _is_poly(e) else match(e, query)) is not None]
    else:
        query = _sympify(query)
        results = [e for e in _preorder(f) if not _is_poly(e) and e == query]
    if not group:
        return set(results)
    counts = {}
    for r in results:
        counts[r] = counts.get(r, 0) + 1
    return counts


def count_ops(f, visual=False):
    from .. import Add, Symbol, Integer
    ops = {}
    stack = [f]
    while stack:
        a = stack.pop()
        args = _args(a) if _is_poly(a) else getattr(a, "args", ())
        if args:
            name = "POLY" if _is_poly(a) else a.func.__name__.upper()
            ops[name] = ops.get(name, 0) + 1
            stack.extend(args)
    if visual:
        return Add(*[Integer(k)*Symbol(name) for name, k in ops.items()])
    return sum(ops.values())


def _dense_key(f):
    return as_list(f)


def compare(f, other):
    if f is other:
        return 0
    if not _is_poly(other):
        return 1
    st, ot = len(f.gens) + 1, len(other.gens) + 1
    c = (st > ot) - (st < ot)
    if c:
        return c
    a, b = _dense_key(f), _dense_key(other)
    try:
        c = int(bool(a > b)) - int(bool(a < b))
    except TypeError:
        c = 0
    if c:
        return c
    for l, r in zip(f.gens, other.gens):
        c = l.compare(r) if hasattr(l, "compare") else _expr_compare(l, r)
        if c:
            return c
    return 0


def _expr_compare(l, r):
    kl, kr = l.sort_key(), r.sort_key()
    return (kl > kr) - (kl < kr)


def _query(query):
    """``Basic._make_find_query``: types, patterns (via ``match``) or callables."""
    if isinstance(query, type):
        return lambda e: _matches_type(e, query)
    if callable(query) and not hasattr(query, "free_symbols"):
        return query
    query = _sympify(query)
    if _is_wild(query):
        return lambda e: (e.match(query) if not _is_poly(e) else match(e, query)) is not None
    return lambda e: not _is_poly(e) and e == query


def count(f, query):
    """Number of subexpressions (the polynomial included) matching ``query``."""
    q = _query(query)
    return sum(bool(q(e)) for e in _preorder(f))


def has_xfree(f, s):
    """True when an element of the set ``s`` is a (free) subexpression."""
    if type(s) is not set:
        raise TypeError('expecting set argument')
    return any(e in s for e in _preorder(f))


def has_free(f, *patterns):
    """``has_xfree`` with types and patterns allowed, as in SymPy."""
    if not patterns:
        return False
    if len(patterns) == 1 and not hasattr(patterns[0], "free_symbols") \
            and not isinstance(patterns[0], type) and hasattr(patterns[0], "__iter__"):
        raise TypeError("\nExpecting 1 or more Basic args, not a single non-Basic "
                        "iterable. Don't\nforget to unpack iterables: "
                        "`eq.has_free(*patterns)`")
    plain = {p for p in patterns if not isinstance(p, type)}
    if has_xfree(f, plain):
        return True
    return any(any(_query(p)(e) for e in _preorder(f)) for p in patterns if isinstance(p, type))


def _postorder(e):
    for a in (_args(e) if _is_poly(e) else getattr(e, "args", ())):
        yield from _postorder(a)
    yield e


def is_same(a, b, approx=None):
    """Structural equality that also compares classes (``Basic.is_same``)."""
    from itertools import zip_longest
    for u, v in zip_longest(_postorder(a), _postorder(b)):
        if u is None or v is None:
            return False
        if getattr(u, "is_Number", False):
            if not getattr(v, "is_Number", False):
                return False
            if approx:
                return approx(u, v)
        if not (u == v and u.__class__ == v.__class__):
            return False
    return True


def dummy_eq(f, other, symbol=None):
    """Equality up to the renaming of a single dummy symbol."""
    from .. import Dummy
    s = f.as_dummy()
    o = other.as_dummy() if _is_poly(other) else _sympify(other)
    dummies = [i for i in s.free_symbols if getattr(i, "is_Dummy", False)]
    if len(dummies) != 1:
        return s == o
    dummy = dummies.pop()
    if symbol is None:
        symbols = o.free_symbols
        if len(symbols) != 1:
            return s == o
        symbol = symbols.pop()
    tmp = Dummy()
    return s.xreplace({dummy: tmp}) == o.xreplace({symbol: tmp})


def could_extract_minus_sign(f):
    return False


def fromiter(cls, args, **assumptions):
    return cls(*tuple(args), **assumptions)


def matches(f, expr, repl_dict=None, old=False):
    """``Basic.matches``: ``f`` is the pattern, matched argument-wise."""
    if not _is_poly(expr):
        return None
    d = {} if repl_dict is None else dict(repl_dict)
    if f == expr and f.gens == expr.gens:
        return d
    if len(_args(f)) != len(_args(expr)):
        return None
    for arg, other in zip(_args(f), _args(expr)):
        if arg == other:
            continue
        m = other.match(arg.xreplace(d))
        if m is None:
            return None
        d.update(m)
    return d


def rcall(f, *args):
    return f(*args)


def refine(f, assumption=True):
    """``refine`` of the expression, re-read over the same generators."""
    from .. import Add, Mul, refine as _refine
    # Refine the coefficients first: the expression-level ``refine`` does
    # not always descend into products.
    terms = [_refine(c, assumption)*Mul(*[g**e for g, e in zip(f.gens, m)])
             for m, c in f.terms()]
    return Poly(_refine(Add(*terms), assumption), *f.gens)


def simplify(f, **kwargs):
    """SymPy rebuilds a ``Poly`` from its arguments (the ground domain is
    re-inferred); coefficients are not simplified."""
    return Poly(f.as_expr(), *f.gens)


def is_hypergeometric(f, k):
    """True when ``f(k + 1)/f(k)`` is a rational function of ``k``."""
    from .. import simplify as _simplify
    k = _sympify(k)
    e = f.as_expr()
    if any(getattr(a, "is_Piecewise", False) for a in _preorder(e)):
        return None
    if not e.has(k) or e.is_polynomial(k):
        return True
    return bool(_simplify(e.subs(k, k + 1)/e).is_rational_function(k))


def xreplace(f, rule):
    if f in rule:
        return rule[f]
    if not rule:
        return f
    args = _args(f)
    new = tuple(a.xreplace(rule) for a in args)
    if all(n == a for n, a in zip(new, args)):
        return f
    return Poly(*new)


def doit(f, **hints):
    return Poly(*[a.doit(**hints) for a in _args(f)])


def rewrite(f, *args, **hints):
    if not args:
        return f
    new = tuple(a.rewrite(*args, **hints) if hasattr(a, "rewrite") else a for a in _args(f))
    return Poly(*new)


def copy(f):
    return Poly(*_args(f))


def as_dummy(f):
    return f


def replace(f, x, y=None, **_ignore):
    """Rename generator ``x`` of ``f`` to ``y``."""
    if y is None:
        if f.is_univariate:
            x, y = f.gen, x
        else:
            raise PolynomialError("syntax supported only in univariate case")
    x, y = _sympify(x), _sympify(y)
    if x == y or x not in f.gens:
        return f
    K = _domain(f)
    domain_symbols = ()
    if K is not None and K.is_Composite:
        domain_symbols = tuple(K.symbols)
    if y in f.gens or y in domain_symbols:
        raise PolynomialError("Cannot replace %s with %s in %s" % (x, y, f))
    gens = tuple(y if g == x else g for g in f.gens)
    return _from_dict(f, dict(_items(f)), gens=gens)


def _subs_one(f, old, new):
    if old in f.gens:
        if new.is_number:
            return f.eval(old, new)
        try:
            out = replace(f, old, new)
        except (PolynomialError, ValueError, TypeError):
            out = None
        if out is not None and set(out.gens) == {new if g == old else g for g in f.gens}:
            return out
    return f.as_expr().subs(old, new)


def subs(f, arg1, arg2=None, **kwargs):
    """Substitute like ``Basic.subs``: generators are evaluated or renamed."""
    if arg2 is None:
        if hasattr(arg1,'items'):
            from .. import default_sort_key
            from ..utilities.sorting import _nodes
            items=list(arg1.items())
            items.sort(key=lambda pair:(-_nodes(_sympify(pair[0])),default_sort_key(pair[0])))
        else:
            try:items=list(arg1)
            except TypeError:raise ValueError('subs requires a mapping or an iterable of pairs') from None
    else:
        items=[(arg1,arg2)]
    if kwargs.get('simultaneous',False):
        from .. import Dummy
        multiplier=Dummy('subs_m');replacements={multiplier:_S().One};rv=f
        for old,new in items:
            old,new=_sympify(old),_sympify(new)
            if old==new:continue
            commutative=getattr(new,'is_commutative',True)
            dummy=Dummy('subs_d',commutative=True if commutative is None else commutative)
            rv=_subs_one(rv,old,dummy*multiplier) if _is_poly(rv) else rv.subs(old,dummy*multiplier)
            replacements[dummy]=new
        return rv.xreplace(replacements)
    rv = f
    for old, new in items:
        old, new = _sympify(old), _sympify(new)
        rv = _subs_one(rv, old, new) if _is_poly(rv) else rv.subs(old, new, **kwargs)
    return rv


_FALSE_FLAGS = (
    "is_Add", "is_AlgebraicNumber", "is_Atom", "is_Boolean", "is_Derivative",
    "is_Dummy", "is_Equality", "is_Float", "is_Function", "is_Indexed",
    "is_Integer", "is_MatAdd", "is_MatMul", "is_Matrix", "is_Mul", "is_Not",
    "is_Number", "is_NumberSymbol", "is_Order", "is_Piecewise", "is_Point",
    "is_Pow", "is_Rational", "is_Relational", "is_Symbol", "is_Vector",
    "is_Wild", "is_comparable", "is_number", "is_symbol",
)

_NONE_FLAGS = (
    "is_algebraic", "is_antihermitian", "is_complex", "is_composite", "is_even",
    "is_extended_negative", "is_extended_nonnegative", "is_extended_nonpositive",
    "is_extended_nonzero", "is_extended_positive", "is_extended_real",
    "is_finite", "is_hermitian", "is_imaginary", "is_infinite", "is_integer",
    "is_irrational", "is_negative", "is_noninteger", "is_nonnegative",
    "is_nonpositive", "is_nonzero", "is_odd", "is_polar", "is_positive",
    "is_prime", "is_rational", "is_real", "is_transcendental",
)


def _rep(f):
    """SymPy's ``Poly.rep``: the dense representation over ``f.domain``."""
    from .polyclasses import DMP
    return DMP.from_poly(f)


def install():
    """Attach the SymPy-compatible methods to the native ``Poly`` class."""
    methods = {
        "terms": terms, "coeffs": coeffs, "monoms": monoms,
        "LC": LC, "EC": EC, "TC": TC, "nth": nth, "LM": LM, "EM": EM, "LT": LT, "ET": ET,
        "as_list": as_list,
        "__reduce__": _reduce_polynomial,
        "add_ground": add_ground, "sub_ground": sub_ground,
        "mul_ground": mul_ground, "quo_ground": quo_ground,
        "exquo_ground": exquo_ground, "abs": abs_, "__abs__": _abs, "neg": neg,
        "add": add, "sub": sub, "mul": mul, "sqr": sqr, "pow": pow_,
        "eq": eq, "ne": ne, "per": per, "unify": unify, "coeff": coeff,
        "as_content_primitive": as_content_primitive,
        "retract": retract, "to_exact": to_exact, "set_modulus": set_modulus,
        "inject": inject, "eject": eject, "rat_clear_denoms": rat_clear_denoms,
        "lift": lift, "sqf_norm": sqf_norm, "revert": revert,
        "nth_power_roots_poly": nth_power_roots_poly, "shift_list": shift_list,
        "real_roots": real_roots, "all_roots": all_roots, "root": root,
        "same_root": same_root, "dispersionset": dispersionset,
        "dispersion": dispersion, "trunc": trunc, "transform": transform,
        "sort_key": sort_key, "atoms": atoms, "has": has, "match": match,
        "find": find, "count_ops": count_ops, "compare": compare,
        "xreplace": xreplace, "doit": doit, "rewrite": rewrite, "copy": copy,
        "as_dummy": as_dummy, "replace": replace, "subs": subs,
        "make_monic_over_integers_by_scaling_roots": make_monic_over_integers_by_scaling_roots,
        "which_real_roots": which_real_roots, "which_all_roots": which_all_roots,
        "count": count, "has_free": has_free, "has_xfree": has_xfree,
        "is_same": is_same, "dummy_eq": dummy_eq,
        "could_extract_minus_sign": could_extract_minus_sign, "matches": matches,
        "rcall": rcall, "refine": refine, "simplify": simplify,
        "is_hypergeometric": is_hypergeometric,
    }
    for name, fn in methods.items():
        setattr(Poly, name, fn)
    native_as_poly = Poly.as_poly

    def as_poly(self, *gens, **args):
        """``Basic.as_poly``: the polynomial in new generators, or None when
        it is not a polynomial in them."""
        if not gens and not args:
            return native_as_poly(self)
        from .polyerrors import PolynomialError, GeneratorsNeeded
        try:
            poly = Poly(self, *gens, **args)
        except (PolynomialError, GeneratorsNeeded):
            return None
        return poly if poly.is_Poly else None
    Poly.as_poly = as_poly
    from .numberfields.galoisgroups import _poly_galois_group
    Poly.galois_group = _poly_galois_group
    for name, fn in {"from_dict": from_dict, "from_list": from_list,
                     "from_poly": from_poly, "from_expr": from_expr,
                     "new": new, "class_key": class_key, "fromiter": fromiter}.items():
        setattr(Poly, name, classmethod(fn))
    for name, fn in {"args": _args, "func": _func, "one": one, "zero": zero, "rep": _rep,
                     "free_symbols_in_domain": free_symbols_in_domain,
                     "is_primitive": is_primitive,
                     "is_cyclotomic": is_cyclotomic}.items():
        setattr(Poly, name, property(fn))
    for name in _FALSE_FLAGS:
        setattr(Poly, name, False)
    for name in _NONE_FLAGS:
        setattr(Poly, name, None)
    Poly.is_Poly = True
    Poly.is_commutative = True
