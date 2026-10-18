"""Exact Gaussian arithmetic, monic quotients and domain conversion hooks."""
import pickle
import random
import subprocess
import sys

import pytest


import sympy as sp
import symox as sx
from symox.polys.agca.extensions import MonogenicFiniteExtension, ExtensionElement
from symox.polys.domains import GeneralizedPolynomialRing, GlobalPolynomialRing


@pytest.mark.parametrize('integral',[False,True])
def test_gaussian_coefficients_coerce_into_scalar_expressions(integral):
    domain=sx.ZZ_I if integral else sx.QQ_I
    x=sx.Symbol('x')
    value=domain.convert(1+sx.I)
    scalar=value.as_expr()
    assert sx.sympify(value)==scalar
    for actual,expected in [(x+value,x+scalar),(value+x,scalar+x),
        (x-value,x-scalar),(value-x,scalar-x),(x*value,x*scalar),
        (value*x,scalar*x),(x/value,x/scalar),(value/x,scalar/x)]:
        assert sx.simplify(actual-expected)==0
    # Coercion into expressions does not erase typed structural identity.
    imaginary=domain.convert(sx.I)
    assert imaginary!=sx.I and sx.I!=imaginary


def compare(k, value, ref, expected):
    assert sp.simplify(k.to_sympy(value).to_sympy()-ref.to_sympy(expected)) == 0


@pytest.mark.parametrize('integral', [True, False])
@pytest.mark.parametrize('coordinates', [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),
                                        (2, 3), (-3, 4), (4, -5), (-6, -7)])
def test_gaussian_signs_coordinates_and_units(integral, coordinates):
    k, ref = (sx.ZZ_I, sp.ZZ_I) if integral else (sx.QQ_I, sp.QQ_I)
    x, y = coordinates
    a, b = k(x, y), ref(x, y)
    assert a.parent() == k and int(a.x) == x and int(a.y) == y
    assert a.quadrant() == b.quadrant()
    compare(k, k.canonical_unit(a), ref, ref.canonical_unit(b))
    for name in ['is_positive', 'is_negative', 'is_nonnegative', 'is_nonpositive']:
        assert getattr(k, name)(a) == getattr(ref, name)(b)
    assert k.real(a) == a and k.imag(a) == k.zero
    assert pickle.loads(pickle.dumps(a)) == a
    assert k.get_ring() == sx.ZZ_I and k.get_field() == sx.QQ_I
    assert k.is_PID


def test_seeded_gaussian_division_gcd_lcm_and_bezout():
    k, ref = sx.ZZ_I, sp.ZZ_I
    rng = random.Random(319)
    for _ in range(120):
        ar, ai, br, bi = [rng.randrange(-30, 31) for _ in range(4)]
        if br == bi == 0:
            br = 1
        a, b, A, B = k(ar, ai), k(br, bi), ref(ar, ai), ref(br, bi)
        q, r = k.div(a, b)
        Q, R = ref.div(A, B)
        compare(k, q, ref, Q)
        compare(k, r, ref, R)
        assert q*b+r == a
        assert r.x*r.x+r.y*r.y < b.x*b.x+b.y*b.y
        s, t, g = k.gcdex(a, b)
        assert s*a+t*b == g
        assert k.rem(a, g) == k.zero and k.rem(b, g) == k.zero
        compare(k, g, ref, ref.gcd(A, B))
        compare(k, k.lcm(a, b), ref, ref.lcm(A, B))


def test_rational_gaussian_content_and_fraction_parts():
    rng = random.Random(101)
    k, ref = sx.QQ_I, sp.QQ_I
    for _ in range(40):
        raw = [rng.randrange(-20, 21) for _ in range(4)]
        if not any(raw):
            raw[0] = 1
        a = k(sx.Rational(raw[0], 3), sx.Rational(raw[1], 5))
        b = k(sx.Rational(raw[2], 7), sx.Rational(raw[3], 11))
        A = ref(sp.QQ(raw[0], 3), sp.QQ(raw[1], 5))
        B = ref(sp.QQ(raw[2], 7), sp.QQ(raw[3], 11))
        compare(k, k.gcd(a, b), ref, ref.gcd(A, B))
        compare(k, k.lcm(a, b), ref, ref.lcm(A, B))
        n, d = k.numer(a), k.denom(a)
        assert n.parent() == sx.ZZ_I and d.parent() == sx.ZZ_I
        assert k.convert(n)/k.convert(d) == a
        if a:
            assert a*(a**-1) == k.one


@pytest.mark.parametrize('expression', ['0', '1', '-1', '3+4*I', '5+12*I',
                                       'I', '2', '(3+4*I)/4'])
def test_exact_gaussian_square_roots(expression):
    value = sx.sympify(expression)
    for k in [sx.ZZ_I, sx.QQ_I]:
        try:
            a = k.convert(value)
        except sx.CoercionFailed:
            continue
        root = k.exsqrt(a)
        if root is not None:
            assert root*root == a
            assert root.x >= 0
    assert sx.QQ_I.exsqrt(sx.QQ_I(3, 4)) == sx.QQ_I(2, 1)
    assert sx.QQ_I.exsqrt(sx.QQ_I(sx.Rational(3, 4), 1)) == sx.QQ_I(1, sx.Rational(1, 2))
    assert sx.ZZ_I.exsqrt(sx.ZZ_I(2)) is None


def test_gaussian_modular_inverse_and_ideals():
    ring = sx.ZZ_I
    modulus = ring(2, 1)
    for a in [ring(1), ring(2), ring(0, 1), ring(1, 1)]:
        inverse = ring.invert(a, modulus)
        assert ring.rem(a*inverse-ring.one, modulus) == ring.zero
    ideal = ring.ideal(modulus)
    assert ideal.contains(modulus*ring(2, -3))
    quotient = ring/ideal
    for a in [quotient(2), quotient(ring(0, 1))]:
        assert a**-1*a == quotient.one
    with pytest.raises(sx.NotInvertible):
        ring.invert(modulus, modulus)
    with pytest.raises(sx.CoercionFailed):
        ring.convert(sx.Rational(1, 2))


@pytest.mark.parametrize('domain', [sx.ZZ, sx.QQ])
@pytest.mark.parametrize('degree', [2, 3, 5, 8])
def test_monogenic_extension_relations_and_representation(domain, degree):
    x = sx.Symbol('x')
    modulus = sx.Poly(x**degree+1, x, domain=domain)
    ext = MonogenicFiniteExtension(modulus)
    t = ext.generator
    assert isinstance(t, ExtensionElement)
    assert t**degree == ext(-1)
    assert t**(2*degree) == ext.one
    assert t*(t**-1) == ext.one
    assert len(ext.basis) == degree
    for power, basis in enumerate(ext.basis):
        assert basis == t**power
    assert (t+1).rep == sx.Poly(x+1, x, domain=domain)
    assert pickle.loads(pickle.dumps(t)) == t
    assert sx.QQ.convert_from(ext(2), ext) == sx.QQ(2)
    assert sx.QQ.convert(ext(2), ext) == sx.QQ(2)
    assert sx.QQ.from_MonogenicFiniteExtension(t, ext) is None
    with pytest.raises(sx.CoercionFailed):
        sx.QQ.convert_from(t, ext)


def test_monogenic_nonfield_units_and_monic_normalization():
    x = sx.Symbol('x')
    ext = MonogenicFiniteExtension(sx.Poly(x*x, x, domain=sx.QQ))
    assert not ext.is_Field
    t = ext.generator
    assert t*t == ext.zero
    assert (ext.one+t)**-1 == ext.one-t
    with pytest.raises(sx.NotInvertible):
        t**-1
    normalized = MonogenicFiniteExtension(sx.Poly(2*x*x-2, x, domain=sx.ZZ))
    assert normalized.modulus == sx.Poly(x*x-1, x, domain=sx.ZZ)
    with pytest.raises((sx.ExactQuotientFailed, sx.CoercionFailed, ValueError)):
        MonogenicFiniteExtension(sx.Poly(2*x*x+1, x, domain=sx.ZZ))
    with pytest.raises(ValueError):
        MonogenicFiniteExtension(sx.Poly(0, x))


def test_all_rational_domain_source_hooks_have_working_inputs():
    x = sx.Symbol('x')
    glob, gen = GlobalPolynomialRing(sx.QQ, x), GeneralizedPolynomialRing(sx.QQ, x)
    ext = MonogenicFiniteExtension(sx.Poly(x*x-2, x, domain=sx.QQ))
    cases = [('ComplexField', sx.CC('0.5'), sx.CC, sx.QQ(1, 2)),
             ('RealField', sx.RR('0.5'), sx.RR, sx.QQ(1, 2)),
             ('ExpressionDomain', sx.EX(2), sx.EX, sx.QQ(2)),
             ('ExpressionRawDomain', sx.EXRAW(2), sx.EXRAW, sx.QQ(2)),
             ('GaussianRationalField', sx.QQ_I(sx.Rational(1, 2)), sx.QQ_I, sx.QQ(1, 2)),
             ('GlobalPolynomialRing', glob(2), glob, sx.QQ(2)),
             ('GeneralizedPolynomialRing', gen(2), gen, sx.QQ(2)),
             ('MonogenicFiniteExtension', ext(2), ext, sx.QQ(2))]
    for name, value, source, expected in cases:
        assert getattr(sx.QQ, 'from_'+name)(value, source) == expected
        assert sx.QQ.convert_from(value, source) == expected
    assert sx.QQ.from_ComplexField(sx.CC(1, 2), sx.CC) is None
    assert sx.QQ.from_GaussianRationalField(sx.QQ_I(1, 2), sx.QQ_I) is None
    assert sx.QQ.from_GlobalPolynomialRing(glob(x), glob) is None
    assert sx.QQ.from_ExpressionDomain(sx.EX(sx.sqrt(2)), sx.EX) is None


def test_gaussian_and_extension_domains_without_sympy():
    code = '''
import sys
sys.modules['sympy']=None
import symox as s
from symox.polys.agca.extensions import MonogenicFiniteExtension
a,b=s.ZZ_I(3,4),s.ZZ_I(1,2)
q,r=s.ZZ_I.div(a,b)
assert q*b+r==a
u,v,g=s.ZZ_I.gcdex(a,b)
assert u*a+v*b==g
x=s.Symbol('x'); ext=MonogenicFiniteExtension(s.Poly(x**5+1,x,domain=s.QQ))
assert ext.generator**5==ext(-1)
assert s.QQ.from_MonogenicFiniteExtension(ext(2),ext)==s.QQ(2)
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
