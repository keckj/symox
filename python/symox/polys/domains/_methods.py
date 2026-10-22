"""Domain protocols implemented with typed native coefficients.

Conversion hooks return None when an embedding does not exist; convert_from
raises CoercionFailed. No optional computer algebra system is used here.
"""
from math import isqrt

_native_convert = None
_native_invert = None
_native_from_sympy = None


def _gaussian(domain):
    return getattr(domain, 'is_GaussianRing', False) or getattr(domain, 'is_GaussianField', False)


def _expr(domain, value):
    return domain.to_sympy(domain.convert(value))


def _sign_coefficient(domain, value):
    value = domain.convert(value)
    if domain.is_FiniteField:
        return int(value)
    if domain.is_Composite:
        from ... import Poly
        numerator = domain.numer(value)
        ring = domain.get_ring()
        return Poly(ring.to_sympy(numerator), *domain.symbols,
                    domain=domain.dom).LC(order=domain.order)
    if domain.is_AlgebraicField:
        coefficients = value.to_list()
        return domain.dom.to_sympy(coefficients[0]) if coefficients else 0
    return domain.to_sympy(value)


def _gaussian_composite(domain):
    # PolynomialRing and FractionField delegate signs to the ground of the
    # leading coefficient, which has none over ZZ_I and QQ_I.
    return domain.is_Composite and _gaussian(domain.dom)


def _has_assoc(domain, method):
    if getattr(domain, 'is_EX', False) or getattr(domain, 'is_EXRAW', False):
        return method == 'get_field'
    if method == 'get_ring' and not domain.is_Exact:
        return False
    try:
        getattr(domain, method)()
    except Exception:
        return False
    return True


def is_positive(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self) or _gaussian_composite(self):
        self.convert(a)
        return False
    if getattr(self, 'is_EX', False):
        return _expr(self, a).as_coeff_mul()[0].is_positive
    return bool(_sign_coefficient(self, a) > 0)


def is_negative(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self) or _gaussian_composite(self):
        self.convert(a)
        return False
    if getattr(self, 'is_EX', False):
        return _expr(self, a).could_extract_minus_sign()
    return False if self.is_FiniteField else bool(_sign_coefficient(self, a) < 0)


def is_nonnegative(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self) or _gaussian_composite(self):
        self.convert(a)
        return False
    return not self.is_negative(a)


def is_nonpositive(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self) or _gaussian_composite(self):
        self.convert(a)
        return False
    return not self.is_positive(a)


def canonical_unit(self, a):
    if _gaussian(self):
        return self.units[-self.convert(a).quadrant()]
    if self.is_Composite and not self.is_Field and _gaussian(self.dom):
        # PolynomialRing.canonical_unit: the ground unit of the leading
        # coefficient, as a constant of the ring.
        from ... import Poly
        value = self.convert(a)
        lc = self.dom.zero
        if value:
            lc = self.dom.convert(Poly(self.to_sympy(value), *self.symbols,
                                       domain=self.dom).LC(order=self.order))
        return self.convert(self.dom.to_sympy(self.dom.canonical_unit(lc)))
    return -self.one if self.is_negative(a) else self.one


def abs_(self, a):
    if _gaussian_composite(self):
        self.convert(a)
        raise TypeError("bad operand type for abs(): '%s'"
                        % ('FracElement' if self.is_Field else 'PolyElement'))
    if getattr(self, 'is_CC', False):
        return abs(self.convert(a))
    if self.is_ZZ or self.is_QQ or self.is_Composite:
        a = self.convert(a)
        return -a if self.is_negative(a) else a
    from ... import Abs
    return self.convert(Abs(_expr(self, a)))


def pos(self, a):
    return self.convert(a)


def almosteq(self, a, b, tolerance=None):
    if not self.is_Exact:
        from ... import Abs, RealField, S, sympify
        a, b = _expr(self, a), _expr(self, b)
        tolerance = S(2)**(4-self.precision) if tolerance is None else sympify(tolerance)
        if tolerance.is_real is not True or tolerance.is_negative is True:
            raise ValueError('numerical tolerance must be nonnegative and real')
        real_field = RealField(self.precision)
        scale = max(real_field(Abs(a)), real_field(Abs(b)), real_field.one)
        return real_field(Abs(a-b)) <= real_field(tolerance)*scale
    return self.convert(a) == self.convert(b)


def cofactors(self, a, b):
    g = self.gcd(a, b)
    return g, self.quo(a, g), self.quo(b, g)


def half_gcdex(self, a, b):
    s, _, g = self.gcdex(a, b)
    return s, g


def revert(self, a):
    from ..polyerrors import NotReversible
    a = self.convert(a)
    if not self.is_unit(a):
        raise NotReversible('only units are reversible in this domain')
    return self.exquo(self.one, a)


def exsqrt(self, a):
    a = self.convert(a)
    if not self.is_Exact:
        from ... import sqrt, I
        expression = self.to_sympy(a)
        if self.is_RR and expression.is_negative is True:
            return None
        if self.is_CC and expression.is_real is True and expression.is_negative is True:
            return self.convert(I*sqrt(-expression))
        return self.convert(sqrt(expression))
    if _gaussian(self):
        from ... import QQ, I, CoercionFailed
        x, y = QQ.convert(a.x), QQ.convert(a.y)
        norm = QQ.exsqrt(x*x+y*y)
        if norm is None:
            return None
        u, v = QQ.exsqrt((norm+x)/QQ(2)), QQ.exsqrt((norm-x)/QQ(2))
        if u is None or v is None:
            return None
        if y < QQ.zero:
            v = -v
        try:
            return self.convert(QQ.to_sympy(u)+I*QQ.to_sympy(v))
        except CoercionFailed:
            return None
    if getattr(self, 'is_EX', False) or getattr(self, 'is_EXRAW', False):
        from ... import sqrt
        return self.convert(sqrt(self.to_sympy(a)))
    if self.is_ZZ or self.is_QQ:
        p, q = int(a.numerator), int(a.denominator)
        if p < 0:
            return None
        r, s = isqrt(p), isqrt(q)
        return self(r, s) if r*r == p and s*s == q else None
    if self.is_FiniteField:
        from ... import sqrt_mod
        r = sqrt_mod(int(a), self.characteristic())
        return None if r is None else self(r)
    if self.is_Composite:
        from ... import Poly
        ring = self.get_ring()
        result = self.one
        for part, inverse in ((self.numer(a), False), (self.denom(a), True)):
            coefficient, factors = Poly(ring.to_sympy(part), *self.symbols,
                                        domain=self.dom).factor_list()
            root = self.dom.exsqrt(self.dom.convert(coefficient))
            if root is None or any(power % 2 for _, power in factors):
                return None
            product = self.convert(root)
            for factor, power in factors:
                product *= self.convert(factor.as_expr()) ** (power // 2)
            result = result / product if inverse else result * product
        return result
    if self.is_AlgebraicField:
        from ... import sqrt, CoercionFailed
        try:
            r = self.convert(sqrt(self.to_sympy(a)))
        except CoercionFailed:
            return None
        return r if r*r == a else None
    raise NotImplementedError('exact square roots require an exact coefficient domain')


def is_square(self, a):
    return self.exsqrt(a) is not None


def sqrt(self, a):
    if self.is_ZZ:
        return self(isqrt(int(self.convert(a))))
    # QQ and general exact domains have no inexact square-root convention.
    # exsqrt is the operation that tests membership of a square root.
    raise NotImplementedError('use exsqrt for exact square roots outside ZZ')


def log(self, a, b):
    if not self.is_ZZ:
        raise NotImplementedError('integer logarithms require ZZ')
    from ... import integer_log
    a, b = int(self.convert(a)), int(self.convert(b))
    if a <= 0 or b <= 1:
        raise ValueError('log requires a > 0 and b > 1')
    return self(integer_log(a, b)[0])


def sum_(self, args):
    value = self.zero
    for a in args:
        value += self.convert(a)
    return value


def map_(self, seq):
    return [self.map(a) if isinstance(a, list) else self(a) for a in seq]


def new(self, *args):
    return self(*args)


def real(self, a):
    return self.convert(a)


def imag(self, a):
    self.convert(a)
    return self.zero


def evalf(self, a, prec=None, **options):
    expression = _expr(self, a)
    return expression.evalf(**options) if prec is None else expression.evalf(prec, **options)


def get_exact(self):
    if self.is_Composite and _inexact_ground(self):
        # The reference ``RR[x].get_exact()`` is ``QQ[x]``.
        ground = self.domain.get_exact()
        return ground.frac_field(*self.symbols) if self.is_FractionField else ground.poly_ring(*self.symbols)
    if self.is_Exact:
        return self
    from ... import QQ, QQ_I
    return QQ_I if getattr(self, 'is_CC', False) else QQ


_native_unify = None


def _inexact_ground(domain):
    ground = domain.domain if domain.is_Composite else domain
    return ground.is_RealField or ground.is_ComplexField


def unify(self, K1, symbols=None):
    # The reference composite rule when an inexact ground is involved:
    # ``ZZ[x].unify(RR)`` is ``RR[x]`` and ``RR[x].unify(CC)`` is ``CC[x]``.
    K0 = self
    if (symbols is None and (K0.is_Composite or K1.is_Composite)
            and (_inexact_ground(K0) or _inexact_ground(K1))
            and not (K0.is_EX or K1.is_EX)):
        g0 = K0.domain if K0.is_Composite else K0
        g1 = K1.domain if K1.is_Composite else K1
        s0 = tuple(K0.symbols) if K0.is_Composite else ()
        s1 = tuple(K1.symbols) if K1.is_Composite else ()
        if not s0 or not s1 or s0 == s1:
            gens = s0 or s1
            ground = g0.unify(g1)
            if K0.is_Composite and (not K1.is_Composite or K0.is_FractionField or K1.is_PolynomialRing):
                fraction = K0.is_FractionField
            else:
                fraction = K1.is_FractionField
            return ground.frac_field(*gens) if fraction else ground.poly_ring(*gens)
    if symbols is None:
        return _native_unify(self, K1)
    return _native_unify(self, K1, symbols)


def drop(self, *symbols):
    if not self.is_Composite:
        return self
    keep = tuple(s for s in self.symbols if s not in symbols)
    ground = self.dom.drop(*symbols)
    if not keep:
        return ground
    factory = ground.poly_ring if self.is_PolynomialRing else ground.frac_field
    return factory(*keep, order=self.order)


def inject(self, *gens):
    from ..polyerrors import GeneratorsError
    if not self.is_Composite:
        return self.poly_ring(*gens)
    if set(self.symbols) & set(gens):
        raise GeneratorsError('common generators in domain and injected symbols')
    factory = self.dom.poly_ring if self.is_PolynomialRing else self.dom.frac_field
    return factory(*(self.symbols + gens), order=self.order)


def unify_with_symbols(self, other, symbols):
    from ... import UnificationFailed
    symbols = set(symbols)
    if ((self.is_Composite and symbols.intersection(self.symbols)) or
            (other.is_Composite and symbols.intersection(other.symbols))):
        raise UnificationFailed('domain generators overlap polynomial generators')
    return self.unify(other)


def unify_composite(self, other):
    return self.unify(other)


def convert_from(self, element, base):
    if getattr(base, '_symox_monogenic_extension', False):
        from ... import CoercionFailed
        result = from_MonogenicFiniteExtension(self, element, base)
        if result is None:
            raise CoercionFailed('there is no canonical embedding of this extension element')
        return result
    return self.convert(element, base)


def _not_a_coeff(element):
    # The reference ``_not_a_coeff``: NaN and the infinities.
    import math
    if isinstance(element, float):
        return math.isnan(element) or math.isinf(element)
    from ... import Basic, S
    return isinstance(element, Basic) and any(
        element is value or element == value
        for value in (S.NaN, S.Infinity, S.NegativeInfinity, S.ComplexInfinity))


def _inexact_value(element):
    """``(complex, re, im)`` for a float, complex, mpf or mpc, as mpmath
    tuples rounded to the 53 bits of the reference ``RealField()`` and
    ``ComplexField()`` parents, or None for other objects."""
    from mpmath import libmp
    name = type(element).__name__
    if isinstance(element, float):
        return False, libmp.from_float(element), libmp.fzero
    if isinstance(element, complex):
        return True, libmp.from_float(element.real), libmp.from_float(element.imag)
    if name == 'mpf' and hasattr(element, '_mpf_'):
        return False, libmp.mpf_pos(element._mpf_, 53, 'n'), libmp.fzero
    if name == 'mpc' and hasattr(element, '_mpc_'):
        re, im = element._mpc_
        return True, libmp.mpf_pos(re, 53, 'n'), libmp.mpf_pos(im, 53, 'n')
    return None


def _inexact_failure(domain, value):
    import mpmath
    from mpmath import libmp
    from ... import CoercionFailed
    complex_, re, im = value
    if complex_:
        shown, kind, base = '(%s)' % libmp.mpc_to_str((re, im), 15), type(mpmath.mpc(0)), 'CC'
    else:
        shown, kind, base = libmp.to_str(re, 15), type(mpmath.mpf(0)), 'RR'
    return CoercionFailed('Cannot convert %s of type %s from %s to %s' % (shown, kind, base, domain))


def _exact(t):
    from mpmath import libmp
    from ... import Rational
    if t in (libmp.fnan, libmp.finf, libmp.fninf):
        return None
    return Rational(*libmp.to_rational(t))


def _limited(t):
    """The reference ``RealField.to_rational(a)`` (limited denominator)."""
    from ... import RR, Rational
    exact = _exact(t)
    if exact is None:
        return None
    return Rational(*RR.to_rational(RR.from_sympy(exact)))


def _as_float(t):
    from mpmath import libmp
    value = libmp.to_float(t)
    return value if libmp.from_float(value) == t else None


def _convert_inexact(domain, value):
    """The reference ``convert_from(parent(element), parent)`` for a
    Python or mpmath float or complex: the element is first rounded into
    ``RR`` or ``CC`` and each domain then applies its own ``from_RealField``
    or ``from_ComplexField`` rule (or has none, which raises)."""
    from mpmath import libmp
    from ... import I, Float
    complex_, re, im = value
    if getattr(domain, 'is_EX', False) or getattr(domain, 'is_EXRAW', False) or \
            getattr(domain, 'is_RealField', False) or getattr(domain, 'is_ComplexField', False):
        if getattr(domain, 'is_RealField', False) and complex_ and im != libmp.fzero:
            raise _inexact_failure(domain, value)
        fr, fi = _as_float(re), _as_float(im)
        if fr is not None and fi is not None:
            return _native_convert(domain, complex(fr, fi) if complex_ and im != libmp.fzero else fr, None)
        er, ei = _exact(re), _exact(im)
        if er is None or ei is None:
            raise _inexact_failure(domain, value)
        expression = Float(er, precision=53) + (I*Float(ei, precision=53) if complex_ and ei else 0)
        return _native_convert(domain, expression, None)
    if getattr(domain, 'is_PolynomialRing', False) or getattr(domain, 'is_FractionField', False):
        ground = domain.dom
        return _native_convert(domain, ground.to_sympy(_convert_inexact(ground, value)), None)
    if getattr(domain, 'is_ZZ_I', False):
        raise _inexact_failure(domain, value)
    if getattr(domain, 'is_QQ_I', False):
        if not complex_:
            raise _inexact_failure(domain, value)
        er, ei = _limited(re), _limited(im)
        if er is None or ei is None:
            raise _inexact_failure(domain, value)
        return _native_convert(domain, er + I*ei, None)
    if getattr(domain, 'is_ZZ', False) or getattr(domain, 'is_QQ', False) or \
            getattr(domain, 'is_FiniteField', False) or getattr(domain, 'is_AlgebraicField', False):
        if complex_:
            raise _inexact_failure(domain, value)
        r = _limited(re)
        if r is None or (not getattr(domain, 'is_Field', False) or domain.is_FiniteField) and not r.is_Integer:
            raise _inexact_failure(domain, value)
        return _native_convert(domain, r, None)
    return None


def convert(self, element, base=None):
    from ... import CoercionFailed, Expr, Basic, Set
    # EXRAW accepts every Expr before the coefficient check (``of_type``).
    raw = getattr(self, 'is_EXRAW', False) and base is None and isinstance(element, Expr)
    if not raw and _not_a_coeff(element):
        raise CoercionFailed('%s is not in any domain' % (element,))
    if base is None and isinstance(element, (set, frozenset, dict)) and (
            getattr(self, 'is_EX', False) or getattr(self, 'is_EXRAW', False)):
        # Not a sequence: strict sympify gives a FiniteSet or a Dict.
        from ... import sympify
        element = sympify(element, strict=True)
        if getattr(self, 'is_EXRAW', False):
            raise CoercionFailed('Expecting an Expr instance but found: %s' % type(element).__name__)
        return _native_convert(self, element, None)
    if base is None and isinstance(element, (set, frozenset, dict)) and _sympy_failure(self, None) is not None:
        from ... import sympify
        a = sympify(element, strict=True)
        raise _sympy_failure(self, a)
    if base is None and (element is None or isinstance(
            element, (str, bytes, bytearray, list, tuple, set, frozenset, dict))):
        # The reference reads only strict sympify input: no strings, and
        # no sequences.
        raise CoercionFailed("Cannot convert %s of type %s to %s" % (element, type(element), self))
    if getattr(base, '_symox_monogenic_extension', False):
        return convert_from(self, element, base)
    if base is None:
        if isinstance(element, bool):
            # ``isinstance(True, int)``: booleans convert as integers.
            element = int(element)
        else:
            value = _inexact_value(element)
            if value is not None:
                if value[1] in _INFINITIES:
                    raise CoercionFailed('%s is not in any domain' % (element,))
                result = _convert_inexact(self, value)
                if result is not None:
                    return result
        if getattr(self, 'is_EXRAW', False) and isinstance(element, (Basic, Set)) and not isinstance(element, Expr):
            raise CoercionFailed('Expecting an Expr instance but found: %s' % type(element).__name__)
        result = _float_from_sympy(self, element)
        if result is not None:
            return result
    try:
        return _native_convert(self, element, base)
    except CoercionFailed:
        from fractions import Fraction
        if base is None and isinstance(element, (Basic, Fraction)):
            if isinstance(element, Fraction):
                from ... import Rational
                element = Rational(element.numerator, element.denominator)
            error = _sympy_failure(self, element)
            if error is not None:
                raise error from None
        raise
    except TypeError:
        if base is not None or isinstance(element, Basic):
            raise
        raise CoercionFailed("Cannot convert %s of type %s to %s" % (element, type(element), self)) from None


# ``mpmath.libmp.finf`` and ``fninf``, spelled out so importing symox
# does not import mpmath.
_INFINITIES = ((0, 0, -456, -2), (1, 0, -789, -3))
_native_call = None


def _call(self, *args):
    """``K(a)`` builds ``K.dtype(a)``: ``ZZ(2.5)`` truncates like ``int``
    and ``QQ(2.5)`` is the exact binary fraction, as ``PythonMPQ``."""
    if len(args) == 1 and _inexact_value(args[0]) is not None:
        a = args[0]
        if getattr(self, 'is_ZZ', False):
            return _native_call(self, int(a))
        if getattr(self, 'is_QQ', False):
            if not isinstance(a, float):
                raise TypeError('PythonMPQ() requires numeric or string argument')
            if isinstance(a, float):
                from fractions import Fraction
                f = Fraction(a)
                return _native_call(self, f.numerator, f.denominator)
        if getattr(self, 'is_RealField', False) and (isinstance(a, complex) or type(a).__name__ == 'mpc'):
            raise TypeError('cannot create mpf from %r' % (a,))
        if getattr(self, 'is_FractionField', False) and not self.dom.is_Field:
            # FracField.ground_new retries through the field of the ground:
            # ZZ(y)(1.5) is 3/2 and ZZ_I(y)(1.5) reports QQ_I.
            from ... import CoercionFailed
            try:
                return self.convert(a)
            except CoercionFailed:
                field = self.dom.get_field()
                return self.convert(field.to_sympy(_convert_inexact(field, _inexact_value(a))))
        return self.convert(a)
    if len(args) == 1 and isinstance(args[0], bool):
        return _native_call(self, int(args[0]))
    return _native_call(self, *args)


def from_sympy(self, *args, **kwargs):
    if getattr(self, 'is_EXRAW', False) and len(args) + len(kwargs) == 1:
        # The reference EXRAW keeps any Expr unchanged and nothing else.
        from ... import CoercionFailed, Expr
        a = args[0] if args else next(iter(kwargs.values()))
        if not isinstance(a, Expr):
            raise CoercionFailed('Expecting an Expr instance but found: %s' % type(a).__name__)
    if len(args) + len(kwargs) == 1:
        result = _float_from_sympy(self, args[0] if args else next(iter(kwargs.values())))
        if result is not None:
            return result
    try:
        return _native_from_sympy(self, *args, **kwargs)
    except Exception as exc:
        from ... import CoercionFailed, Basic
        if not isinstance(exc, CoercionFailed) or len(args) + len(kwargs) != 1:
            raise
        a = args[0] if args else next(iter(kwargs.values()))
        error = _sympy_failure(self, a) if isinstance(a, Basic) else None
        if error is None:
            raise
        raise error from None


def _sympy_failure(domain, a):
    """The exception the reference ``from_sympy`` of ``ZZ``, ``QQ``,
    ``RR``, ``CC``, ``ZZ_I`` and ``QQ_I`` raises for a value it refuses
    (None for other domains; with ``a`` None, only tells whether the
    domain is one of them)."""
    from ... import CoercionFailed, I, Dict, FiniteSet
    if getattr(domain, 'is_ZZ_I', False) or getattr(domain, 'is_QQ_I', False):
        if a is None:
            return True
        if not hasattr(a, 'as_coeff_Add'):
            return AttributeError("'%s' object has no attribute 'as_coeff_Add'" % type(a).__name__)
        r, b = a.as_coeff_Add()
        ground = domain.dom
        try:
            ground.from_sympy(r)
        except CoercionFailed as exc:
            return exc
        if b:
            r, b = b.as_coeff_Mul()
            try:
                ground.from_sympy(r)
            except CoercionFailed as exc:
                return exc
        return CoercionFailed('%s is not Gaussian' % (a,))
    if getattr(domain, 'is_ZZ', False):
        message = 'expected an integer, got %s'
    elif getattr(domain, 'is_QQ', False):
        message = 'expected `Rational` object, got %s'
    elif getattr(domain, 'is_RealField', False):
        if a is not None and isinstance(a, Dict):
            return AttributeError("'Dict' object has no attribute 'evalf'")
        message = 'expected real number, got %s'
    elif getattr(domain, 'is_ComplexField', False):
        if a is not None and isinstance(a, Dict):
            return AttributeError("'Dict' object has no attribute 'evalf'")
        if a is not None and isinstance(a, FiniteSet):
            return AttributeError("'FiniteSet' object has no attribute 'as_real_imag'")
        message = 'expected complex number, got %s'
    else:
        return None
    return True if a is None else CoercionFailed(message % (a,))


def _float_ground(domain, r, a):
    """The reference ``ZZ.from_sympy`` / ``QQ.from_sympy`` of a number
    ``r``: an integer-valued Float is an integer, and ``QQ`` reads a Float
    through ``RR`` (limited denominator)."""
    from ... import CoercionFailed, Integer
    if not getattr(r, 'is_Float', False):
        return r if getattr(r, 'is_Rational', False) else None
    if domain.is_ZZ or domain.is_ZZ_I:
        from mpmath import libmp
        p, q = libmp.to_rational(r._mpf_)
        if q == 1:
            return Integer(p)
        raise CoercionFailed('expected an integer, got %s' % a if domain.is_ZZ else 'expected an integer, got %s' % r)
    from mpmath import libmp
    return _limited(libmp.mpf_pos(r._mpf_, 53, 'n'))


def _float_from_sympy(domain, a):
    """Float coefficients for ``ZZ``, ``QQ``, ``ZZ_I``, ``QQ_I`` and
    algebraic fields, as the reference ``from_sympy``; None elsewhere."""
    from ... import Basic, Float, I, CoercionFailed
    if not isinstance(a, Basic) or not getattr(a, 'is_number', False) or not a.has(Float):
        return None
    if domain.is_ZZ or domain.is_QQ or getattr(domain, 'is_AlgebraicField', False):
        if not a.is_Float:
            return None
        ground = domain if not getattr(domain, 'is_AlgebraicField', False) else domain.dom
        r = _float_ground(ground, a, a)
        return None if r is None else _native_convert(domain, r, None)
    if getattr(domain, 'is_ZZ_I', False) or getattr(domain, 'is_QQ_I', False):
        r, b = a.as_coeff_Add()
        x = _float_ground(domain, r, a)
        if x is None:
            return None
        if not b:
            return _native_convert(domain, x, None)
        r, b = b.as_coeff_Mul()
        y = _float_ground(domain, r, a)
        if y is None or b is not I:
            return None
        return _native_convert(domain, x + I*y, None)
    return None


def from_MonogenicFiniteExtension(self, a, K0):
    from ..agca.extensions import MonogenicFiniteExtension
    from ... import CoercionFailed
    if not isinstance(K0, MonogenicFiniteExtension):
        raise TypeError('the source must be a native monogenic finite extension')
    a = K0.convert(a)
    expression = K0.ring.to_sympy(a.data)
    if K0.symbol in expression.free_symbols:
        return None
    ground = K0.domain.convert(expression)
    try:
        return self.convert(ground, K0.domain)
    except CoercionFailed:
        return None


def invert(self, a, b):
    # The reference ``Ring.invert``: ``s % b`` when the gcd is exactly one.
    gaussian = getattr(self, 'is_ZZ_I', False) or getattr(self, 'is_QQ_I', False)
    complex_field = getattr(self, 'is_ComplexField', False)
    if not (self.is_ZZ or self.is_QQ or self.is_RR or self.is_EX or gaussian or complex_field):
        return _native_invert(self, a, b)
    from ... import NotInvertible
    from ... import _core
    # The reference ``s % b`` keeps ``b`` unconverted: mpc for CC elements.
    modulus = 'mpc' if isinstance(b, _core.DomainElement) else type(b).__name__
    a, b = self.convert(a), self.convert(b)
    s, _, h = self.gcdex(a, b)
    if not self.is_one(h):
        raise NotInvertible('zero divisor')
    if complex_field:
        # mpmath's mpc has no ``%``.
        raise TypeError("unsupported operand type(s) for %%: 'mpc' and '%s'" % modulus)
    if self.is_QQ or self.is_EX:
        # The reference ground types have no ``%`` (PythonMPQ, Expression).
        name = 'PythonMPQ' if self.is_QQ else 'Expression'
        raise TypeError("unsupported operand type(s) for %%: '%s' and '%s'" % (name, name))
    if self.is_RR:
        # The reference ``s % b`` is mpmath's ``mpf_mod``: it returns ``s``
        # unchanged when ``b`` (even zero) has the same sign and a larger
        # exponent, so ``RR.invert(4, 0)`` is 0.25 while ``RR.invert(-2, 0)``
        # divides by zero.
        import mpmath
        from mpmath.libmp import mpf_mod, round_nearest
        from ... import Float
        x, m = self.to_sympy(s)._mpf_, self.to_sympy(self.convert(b))._mpf_
        return self.convert(Float(mpmath.mpf(mpf_mod(x, m, self.precision, round_nearest)),
                                  precision=self.precision))
    if not b:
        raise ZeroDivisionError('integer modulo by zero' if self.is_ZZ else 'modulo by zero')
    if getattr(self, 'is_ZZ_I', False):
        return self.rem(s, b)
    if self.is_Field:
        # Gaussian rationals reduce to zero.
        return self.zero
    x, m = int(self.to_sympy(s)), int(self.to_sympy(b))
    return self.convert(x % m)


def alg_field_from_poly(self, poly, alias=None, root_index=-1):
    import operator
    from ... import Poly, QQ
    if not self.is_QQ:
        raise NotImplementedError('algebraic extensions currently require QQ')
    if not isinstance(poly, Poly):
        poly = Poly(poly)
    roots = poly.all_roots()
    root = roots[operator.index(root_index)]
    if root.is_Rational:
        return QQ
    return self.algebraic_field(root, alias=alias)


def cyclotomic_field(self, n, ss=False, alias='zeta', gen=None, root_index=-1):
    import operator
    from ... import Dummy, Poly, cyclotomic_poly
    n = operator.index(n)
    if n < 1:
        raise ValueError('cyclotomic order must be positive')
    if ss:
        alias += str(n)
    gen = Dummy('x') if gen is None else gen
    return self.alg_field_from_poly(Poly(cyclotomic_poly(n, gen), gen),
                                    alias=alias, root_index=root_index)


def free_module(self, rank):
    from ._modules import FreeModule
    return FreeModule(self, rank)


def ideal(self, *gens):
    from ._modules import Ideal
    return Ideal(self, gens)


def quotient_ring(self, e):
    from ._modules import Ideal, QuotientRing
    if not isinstance(e, Ideal):
        e = self.ideal(*e)
    return QuotientRing(self, e)


def old_poly_ring(self, *symbols, **kwargs):
    """Construct the native ring, including its module and ideal protocols."""
    return self.poly_ring(*symbols, **kwargs)


def old_frac_field(self, *symbols, **kwargs):
    return self.frac_field(*symbols, **kwargs)


def _conversion_hook(name):
    def from_domain(self, a, K0):
        from ... import CoercionFailed, Domain
        if not isinstance(K0, Domain):
            raise TypeError('the source must be a native Domain')
        # Fraction-field hooks do not specialize functions at a point.
        if K0.is_FractionField and (self.is_ZZ or self.is_QQ):
            return None
        try:
            return self.convert(a, K0)
        except CoercionFailed:
            return None
    from_domain.__name__ = name
    return from_domain


# SageMath's element conversion protocol. SymPy's ZZ and QQ elements are
# gmpy2 or Python integers and rationals, which SageMath's Integer and
# Rational constructors read directly; SageMath asks other objects for
# ``_integer_`` and ``_rational_``, and polynomial rings for ``_polynomial_``.
# SageMath is imported only when it makes these calls.

def _rational_parts(element):
    domain = element.parent()
    if not (domain.is_ZZ or domain.is_QQ):
        raise TypeError('unable to convert %r to a rational' % (element,))
    return int(element.numerator), int(element.denominator)


def _integer_(self, Z=None):
    p, q = _rational_parts(self)
    if q != 1:
        raise TypeError('unable to convert %r to an integer' % (self,))
    if Z is None:
        from sage.rings.integer_ring import ZZ as Z
    return Z(p)


def _rational_(self):
    p, q = _rational_parts(self)
    from sage.rings.rational_field import QQ
    return QQ((p, q))


def _polynomial_(self, R):
    domain = self.parent()
    if not domain.is_PolynomialRing:
        raise TypeError('unable to convert %r to a polynomial' % (self,))
    from ... import Poly
    poly = Poly(domain.to_sympy(self), *domain.symbols, domain=domain.domain)
    gens = [R(str(symbol)) for symbol in domain.symbols]
    base = R.base_ring()
    result = R.zero()
    for monom, coeff in poly.rep.terms():
        term = R(base(coeff))
        for gen, exponent in zip(gens, monom):
            term *= gen**exponent
        result += term
    return result


_ring_cache = {}


def _domain_ring(self):
    """`PolynomialRing.ring`: the `PolyRing` of the domain elements."""
    if not self.is_PolynomialRing:
        raise AttributeError("'%s' domain has no attribute 'ring'" % self)
    key = ('ring', self)
    if key not in _ring_cache:
        from ..rings import PolyRing
        _ring_cache[key] = PolyRing(self.symbols, self.domain, str(self.order))
    return _ring_cache[key]


def _domain_field(self):
    """`FractionField.field`: the `FracField` of the domain elements."""
    if not self.is_FractionField:
        raise AttributeError("'%s' domain has no attribute 'field'" % self)
    key = ('field', self)
    if key not in _ring_cache:
        from ..fields import FracField
        _ring_cache[key] = FracField(self.symbols, self.domain, str(self.order))
    return _ring_cache[key]


def install():
    from ... import Domain
    from ... import _core
    from . import _element_print
    _element_print.install()
    for function in (_integer_, _rational_, _polynomial_):
        setattr(_core.DomainElement, function.__name__, function)
    global _native_convert, _native_invert, _native_from_sympy, _native_unify
    if _native_unify is None:
        _native_unify = Domain.unify
        Domain.unify = unify
        native_exact = Domain.is_Exact
        Domain.ring = property(_domain_ring)
        Domain.field = property(_domain_field)
        Domain.has_assoc_Ring = property(lambda self: _has_assoc(self, 'get_ring'))
        Domain.has_assoc_Field = property(lambda self: _has_assoc(self, 'get_field'))
        Domain.is_Exact = property(
            lambda self: False if self.is_Composite and _inexact_ground(self)
            else native_exact.__get__(self))
    if _native_convert is None:
        _native_convert, _native_invert = Domain.convert, Domain.invert
        _native_from_sympy = Domain.from_sympy
    for name, function in {
        'is_positive': is_positive, 'is_negative': is_negative,
        'is_nonnegative': is_nonnegative, 'is_nonpositive': is_nonpositive,
        'canonical_unit': canonical_unit, 'abs': abs_, 'pos': pos,
        'almosteq': almosteq, 'cofactors': cofactors, 'half_gcdex': half_gcdex,
        'revert': revert, 'exsqrt': exsqrt, 'sqrt': sqrt, 'is_square': is_square,
        'log': log, 'sum': sum_, 'map': map_, 'new': new, 'normal': new, 'tp': new,
        'real': real, 'imag': imag, 'evalf': evalf, 'n': evalf,
        'get_exact': get_exact, 'drop': drop, 'inject': inject,
        'unify_with_symbols': unify_with_symbols, 'unify_composite': unify_composite,
        'convert_from': convert_from,
        'alg_field_from_poly': alg_field_from_poly, 'cyclotomic_field': cyclotomic_field,
        'free_module': free_module, 'ideal': ideal, 'quotient_ring': quotient_ring,
        'old_poly_ring': old_poly_ring, 'old_frac_field': old_frac_field,
        '__truediv__': quotient_ring,
    }.items():
        setattr(Domain, name, function)
    for source in ('AlgebraicField', 'FF', 'FF_gmpy', 'FF_python',
                   'FractionField', 'PolynomialRing',
                   'QQ', 'QQ_gmpy', 'QQ_python',
                   'ZZ', 'ZZ_gmpy', 'ZZ_python'):
        name = 'from_' + source
        setattr(Domain, name, _conversion_hook(name))
    global _native_call
    if _native_call is None:
        _native_call = Domain.__call__
        Domain.__call__ = _call
    if hasattr(_core, 'EX'):
        Domain.convert = convert
        Domain.invert = invert
        Domain.from_sympy = from_sympy
        Domain.from_MonogenicFiniteExtension = from_MonogenicFiniteExtension
        for source in ('ComplexField', 'RealField', 'ExpressionDomain', 'ExpressionRawDomain',
                       'GaussianRationalField', 'GeneralizedPolynomialRing', 'GlobalPolynomialRing'):
            name = 'from_'+source
            setattr(Domain, name, _conversion_hook(name))
