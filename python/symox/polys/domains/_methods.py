"""Domain protocols implemented with typed native coefficients.

Conversion hooks return None when an embedding does not exist; convert_from
raises CoercionFailed. No optional computer algebra system is used here.
"""
from math import isqrt

_native_convert = None
_native_invert = None


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


def is_positive(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self):
        self.convert(a)
        return False
    if getattr(self, 'is_EX', False):
        return _expr(self, a).as_coeff_mul()[0].is_positive
    return bool(_sign_coefficient(self, a) > 0)


def is_negative(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self):
        self.convert(a)
        return False
    if getattr(self, 'is_EX', False):
        return _expr(self, a).could_extract_minus_sign()
    return False if self.is_FiniteField else bool(_sign_coefficient(self, a) < 0)


def is_nonnegative(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self):
        self.convert(a)
        return False
    return not self.is_negative(a)


def is_nonpositive(self, a):
    if getattr(self, 'is_CC', False) or _gaussian(self):
        self.convert(a)
        return False
    return not self.is_positive(a)


def canonical_unit(self, a):
    if _gaussian(self):
        return self.units[-self.convert(a).quadrant()]
    return -self.one if self.is_negative(a) else self.one


def abs_(self, a):
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
    if self.is_Exact:
        return self
    from ... import QQ, QQ_I
    return QQ_I if getattr(self, 'is_CC', False) else QQ


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


def convert(self, element, base=None):
    if getattr(base, '_symox_monogenic_extension', False):
        return convert_from(self, element, base)
    return _native_convert(self, element, base)


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
    if not (self.is_ZZ or self.is_QQ or self.is_RR or self.is_EX or gaussian):
        return _native_invert(self, a, b)
    from ... import NotInvertible
    a, b = self.convert(a), self.convert(b)
    s, _, h = self.gcdex(a, b)
    if not self.is_one(h):
        raise NotInvertible('zero divisor')
    if self.is_QQ or self.is_EX:
        # The reference ground types have no ``%`` (PythonMPQ, Expression).
        name = 'PythonMPQ' if self.is_QQ else 'Expression'
        raise TypeError("unsupported operand type(s) for %%: '%s' and '%s'" % (name, name))
    if not b:
        raise ZeroDivisionError('integer modulo by zero' if self.is_ZZ else 'modulo by zero')
    if getattr(self, 'is_ZZ_I', False):
        return self.rem(s, b)
    if self.is_Field:
        # Gaussian rationals reduce to zero; real floats keep Python's
        # floored remainder.
        if not self.is_RR:
            return self.zero
        from ... import floor
        x, m = self.to_sympy(s), self.to_sympy(b)
        return self.convert(x - m*floor(x/m))
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


def install():
    from ... import Domain
    from ... import _core
    for function in (_integer_, _rational_, _polynomial_):
        setattr(_core.DomainElement, function.__name__, function)
    global _native_convert, _native_invert
    if _native_convert is None:
        _native_convert, _native_invert = Domain.convert, Domain.invert
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
    if hasattr(_core, 'EX'):
        Domain.convert = convert
        Domain.invert = invert
        Domain.from_MonogenicFiniteExtension = from_MonogenicFiniteExtension
        for source in ('ComplexField', 'RealField', 'ExpressionDomain', 'ExpressionRawDomain',
                       'GaussianRationalField', 'GeneralizedPolynomialRing', 'GlobalPolynomialRing'):
            name = 'from_'+source
            setattr(Domain, name, _conversion_hook(name))
