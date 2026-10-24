"""RR/CC and Gaussian composite domain parity (wave 16, rrcc).

Expected values are those of SymPy 1.14, except the documented
intentional correction of the ``ZZ_I`` square-free unit.
"""
from symox import (I, Float, Poly, RR, RealField, S, expand, gcdex,
                   half_gcdex, invert, sqf_list, factor_list, symbols)

x, y, z = symbols('x y z')
R100 = RealField(100)


def F30(s):
    return Float(s, 30)


# High precision gcdex / invert keep the precision of the floats.

def test_invert_gcdex_expression_prec30():
    assert str(invert(F30('1.5')*x, x**2 + 1)) == '-0.666666666666666666666666666667*x'
    assert str(gcdex(F30('1.5')*x + 1, x**2 + 1)) == (
        '(0.307692307692307692307692307692 - 0.461538461538461538461538461538*x, '
        '0.692307692307692307692307692308, 1.00000000000000000000000000000)')
    assert str(half_gcdex(F30('1.5')*x + 1, x**2 + 1)) == (
        '(0.307692307692307692307692307692 - 0.461538461538461538461538461538*x, '
        '1.00000000000000000000000000000)')


def test_invert_gcdex_poly_prec():
    p = Poly(F30('1.5')*x, x).invert(Poly(x**2 + 1, x))
    assert str(p) == "Poly(-0.666666666666666666666666666667*x, x, domain='RR')"
    s, t, h = Poly(F30('1.5')*x + 1, x).gcdex(Poly(x**2 + 1, x))
    assert str(s) == ("Poly(-0.461538461538461538461538461538*x + "
                      "0.307692307692307692307692307692, x, domain='RR')")
    assert str(t) == "Poly(0.692307692307692307692307692308, x, domain='RR')"
    assert str(h) == "Poly(1.0, x, domain='RR')"
    assert s.domain.precision == 103
    q = Poly(F30('1.5')*x, x, domain=R100).invert(Poly(x**2 + 1, x))
    assert str(q) == "Poly(-0.66666666666666666666666666667*x, x, domain='RR')"
    assert q.domain.precision == 100
    assert str(Poly(1.5*x, x).invert(Poly(x**2 + 1, x))) == \
        "Poly(-0.666666666666667*x, x, domain='RR')"


# trunc over RR: floor modulo of every coefficient.

def test_trunc_rr():
    p = Poly(F30('1.1')*x**2 + x, x).trunc(3)
    assert str(p) == "Poly(1.1*x**2 + 1.0*x, x, domain='RR')"
    assert [str(c) for c in p.coeffs()] == [
        '1.10000000000000000000000000000', '1.00000000000000000000000000000']
    assert str(Poly(1.1*x**2 + 5.0*x + 7, x).trunc(3)) == \
        "Poly(1.1*x**2 + 2.0*x + 1.0, x, domain='RR')"
    assert str(Poly(-1.5*x - 7.0, x).trunc(3)) == "Poly(1.5*x + 2.0, x, domain='RR')"
    assert str(Poly(1.5*x + 7.0, x).trunc(S(3)/2)) == "Poly(1.0, x, domain='RR')"


# Zero Float coefficients and the zero polynomial over RR.

def test_zero_float_domain():
    assert str(Poly(0.0, x)) == "Poly(0.0, x, domain='RR')"
    assert Poly(F30('0.0'), x, y).domain.precision == 103
    assert str(Poly([0.0, 1], x)) == "Poly(1.0, x, domain='RR')"
    assert str(Poly({(0,): 0.0, (1,): 1}, x)) == "Poly(1.0*x, x, domain='RR')"
    assert str(Poly(Float(0.0), x, domain='QQ')) == "Poly(0, x, domain='QQ')"


def test_zero_rr_accessors():
    p = Poly(0, x, domain=RR)
    assert str(p.LC()) == '0.0' and p.LC().is_Float
    assert [str(c) for c in p.coeffs()] == ['0.0']
    assert str(p.terms()) == '[((0,), 0.0)]'
    assert str(p.nth(0)) == '0.0'
    assert str(p.TC()) == '0.0'
    q = Poly(x**2 + 3*x, x)
    assert q.TC() == 0 and not q.TC().is_Float and q.nth(1) == 3


# numer / denom of RR(x) elements.

def test_rr_fraction_field_parts():
    F = RR.frac_field(x)
    g = F.field.gens[0]
    e = g/(g + 1)
    assert str(e) == 'x/(x + 1.0)'
    assert str(e.numer) == 'x'
    assert str(e.denom) == 'x + 1.0'
    v = F.convert(x)/F.convert(x + 1)
    assert str(F.numer(v)) == 'x'
    assert str(F.denom(v)) == 'x + 1.0'


# sqf_list of expressions with float factors: components on their own domain.

def test_sqf_list_float_components():
    assert str(sqf_list(2.0*(x - y)**2*(x + y))) == \
        '(2.00000000000000, [(x + y, 1), (x - y, 2)])'
    assert str(sqf_list(F30('1.1')*(x - 1)**2)) == \
        '(1.10000000000000000000000000000, [(x - 1, 2)])'
    assert str(sqf_list(1.0*x**2 - 2.0*x + 1.0)) == '(1.00000000000000, [(1.0*x - 1.0, 2)])'
    assert str(sqf_list(2.0*(x - 1)**2*(x + 1.5))) == \
        '(2.00000000000000, [(1.0*x + 1.5, 1), (x - 1, 2)])'
    assert str(sqf_list(1.5*(x - y)**3, x)) == '(1.50000000000000, [(x - y, 3)])'
    assert str(sqf_list(2.0*(x - y)**2*(x + y), x, y, polys=True)) == (
        "(2.00000000000000, [(Poly(x + y, x, y, domain='ZZ'), 1), "
        "(Poly(x - y, x, y, domain='ZZ'), 2)])")
    assert str(factor_list(2.0*(x - y)**2*(x + y))) == \
        '(2.00000000000000, [(x + y, 1), (x - y, 2)])'


# RR or RR[y] against ZZ[y]: the domains unify to RR[y].

def test_mixed_rr_composite_univariate():
    a, b = Poly(1.5*x + 1, x), Poly(x*y + 1, x)
    expected = "Poly(1.5*y*x**2 + (1.0*y + 1.5)*x + 1.0, x, domain='RR[y]')"
    assert str(a*b) == expected
    assert str(b*a) == expected
    assert str(a + b) == "Poly((1.0*y + 1.5)*x + 2.0, x, domain='RR[y]')"
    assert str(a.gcd(Poly(x*y + y, x))) == "Poly(1.0, x, domain='RR[y]')"
    assert str(Poly(x + 1.5, x).compose(Poly(x*y, x))) == "Poly(1.0*y*x + 1.5, x, domain='RR[y]')"
    c = Poly(F30('1.5')*x + 1, x)*b
    assert str(c) == expected


def test_mixed_rr_composite_multivariate():
    expected = ("Poly(1.5*y*x**2*z + 1.5*x*z**2 + 1.0*y**2*x + 1.0*y*z, x, z, "
                "domain='RR[y]')")
    assert str(Poly(1.5*x*z + y, x, z)*Poly(x*y + z, x, z)) == expected
    assert str(Poly(x*y + z, x, z)*Poly(1.5*x*z + y, x, z)) == expected
    assert str(Poly(1.5*x*z + y*x, x, z).rem(Poly(x + y, x, z))) == \
        "Poly(-1.5*y*z - 1.0*y**2, x, z, domain='RR[y]')"


# content / primitive over RR[y] in several generators: a monic gcd.

def test_rr_composite_multivariate_content():
    assert str(Poly(1.5*x*z + 1.5*y, x, z).content()) == '1.00000000000000'
    assert str(Poly(1.5*y*x*z + 3*y, x, z).content()) == '1.0*y'
    assert str(Poly(1.5*y*x*z + 3*y, x, z).primitive()) == \
        "(1.0*y, Poly(1.5*x*z + 3.0, x, z, domain='RR[y]'))"
    assert str(Poly(-1.5*y*x*z - 3*y**2, x, z).primitive()) == \
        "(1.0*y, Poly(-1.5*x*z - 3.0*y, x, z, domain='RR[y]'))"
    assert str(Poly(2.5*y**2*x*z + 5*y, x, z).primitive()) == \
        "(1.0*y, Poly(2.5*y*x*z + 5.0, x, z, domain='RR[y]'))"
    assert str(Poly(1.5*y*x*z + 3*y*z, x, z).primitive()) == \
        "(1.0*y, Poly(1.5*x*z + 3.0*z, x, z, domain='RR[y]'))"
    assert str(Poly(1.5*x*z + 1.5*y, x, z).primitive()) == \
        "(1.00000000000000, Poly(1.5*x*z + 1.5*y, x, z, domain='RR[y]'))"


# Intentional correction: the SymPy ZZ_I dmp_sqf_list loses the unit
# (it gives (1, [(-x + I, 2)]) for (I*x + 1)**2, whose product is
# -(I*x + 1)**2). symox keeps the SymPy factors and the exact unit.

def _product(unit, factors):
    r = unit
    for f, k in factors:
        r *= f.as_expr()**k
    return expand(r)


def test_gaussian_sqf_unit_correction():
    cases = [
        (Poly((I*x + 1)**2, x), "(-1, [(Poly(-x + I, x, domain='ZZ_I'), 2)])"),
        (Poly((I*x + y)**3, x, y), "(-1, [(Poly(-I*x - y, x, y, domain='ZZ_I'), 3)])"),
        (Poly((I*x + 2*y)**2, x, y), "(-1, [(Poly(-x + 2*I*y, x, y, domain='ZZ_I'), 2)])"),
        (Poly((2*I*x + 1)**2*(x + I), x),
         "(-1, [(Poly(x + I, x, domain='ZZ_I'), 1), (Poly(-2*x + I, x, domain='ZZ_I'), 2)])"),
    ]
    for p, expected in cases:
        unit, factors = p.sqf_list()
        assert str((unit, factors)) == expected
        assert _product(unit, factors) == expand(p.as_expr())
    # Unchanged where the reference is already exact.
    assert str(Poly((x + 1)**2*(x + I), x).sqf_list()) == (
        "(1, [(Poly(x + I, x, domain='ZZ_I'), 1), (Poly(x + 1, x, domain='ZZ_I'), 2)])")
    assert str(Poly((I*x + 1)**2, x).sqf_part()) == "Poly(-x + I, x, domain='ZZ_I')"
