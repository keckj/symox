"""Explicit and inferred CC, RR[y] gcd, Gaussian composites and the
half-integer power branch factor (wave 17, ccdom).

Expected values are those of SymPy 1.14; ``test_zz_i_sqf_list_unit`` checks
the documented intentional correction numerically.
"""
import pytest
from symox import *
del test

x, y, z, t = symbols('x y z t')
a, b = symbols('a b')
R100 = RealField(100)
C100 = ComplexField(100)


def F30(s):
    return Float(s, 30)


CASES = [
    ('Poly(1.5*x*z + y, x, z).gcd(Poly(x*y + z, x, z))',
     "Poly(0.5, x, z, domain='RR[y]')"),
    ('Poly(1.0*x + I, x).invert(Poly(x**2 + 2, x))',
     "Poly(-1.0*x + 1.0*I, x, domain='CC')"),
    ('Poly(1.0*x + I, x).gcdex(Poly(x**2 + 2, x))',
     "(Poly(-1.0*x + 1.0*I, x, domain='CC'), Poly(1.0, x, domain='CC'), Poly(1.0, x, domain='CC'))"),
    ('Poly(x, x, domain=CC)',
     "Poly(1.0*x, x, domain='CC')"),
    ('Poly(x + 1, x, domain=CC)',
     "Poly(1.0*x + 1.0, x, domain='CC')"),
    ('Poly(x + 2*I, x, domain=CC)',
     "Poly(1.0*x + 2.0*I, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC).factor_list()',
     "(1, [(Poly(1.0*x**2 + 1.0, x, domain='CC'), 1)])"),
    ('Poly(x**2 + 1, x, domain=CC).nroots()',
     '[-1.0*I, 1.0*I]'),
    ('Poly(x**2 + 1, x, domain=CC).degree()',
     '2'),
    ('Poly(x + 1, x, domain=CC) * Poly(x + I, x)',
     "Poly(1.0*x**2 + (1.0 + 1.0*I)*x + 1.0*I, x, domain='CC')"),
    ('Poly(1.5*x + I, x, domain=CC).gcd(Poly(x + 2, x, domain=CC))',
     "Poly(1.0, x, domain='CC')"),
    ('Poly(x**2 + 2*x + 1, x, domain=CC).sqf_list()',
     "(1.00000000000000, [(Poly(1.0*x + 1.0, x, domain='CC'), 2)])"),
    ('Poly(x**2 - 1, x, domain=CC).div(Poly(x - 1, x, domain=CC))',
     "(Poly(1.0*x + 1.0, x, domain='CC'), Poly(0, x, domain='CC'))"),
    ('Poly(x**2 - 1, x, domain=CC).monic()',
     "Poly(1.0*x**2 - 1.0, x, domain='CC')"),
    ('Poly(2*x**2 - 2, x, domain=CC).LC()',
     '2.00000000000000'),
    ('Poly(x, x, domain=CC).domain',
     'CC'),
    ('Poly(x + 1, x, domain=C100)',
     "Poly(1.0*x + 1.0, x, domain='CC')"),
    ('Poly(x + I, x, domain=C100)',
     "Poly(1.0*x + 1.0*I, x, domain='CC')"),
    ('factor(x**2 + 1, domain=CC)',
     '1.0*x**2 + 1.0'),
    ('factor(x**2 - 2, domain=RR)',
     '2.0*(0.5*x**2 - 1.0)'),
    ('gcd(x**2 - 1, x - 1, domain=CC)',
     '1.0*x - 1.0'),
    ('div(x**2 + 1, x + I, domain=CC)',
     '(1.0*x - 1.0*I, 0)'),
    ('Poly(x + I/2, x).gcd(Poly(x**2 + 1/4, x))',
     "Poly(1.0*x + 0.5*I, x, domain='CC')"),
    ('Poly(x**2/3 - 1, x, domain=C100).coeffs()',
     '[0.33333333333333333333333333333, -1.0000000000000000000000000000]'),
    ("Poly(F30('1.1')*x**2 + 1, x).set_domain(CC)",
     "Poly(1.1*x**2 + 1.0, x, domain='CC')"),
    ("Poly(F30('1.1')*x**2 + 1, x).unify(Poly(x, x, domain=CC))",
     "(Poly(1.1*x**2 + 1.0, x, domain='CC'), Poly(1.0*x, x, domain='CC'))"),
    ('Poly(x**2 + 1, x, domain=CC)',
     "Poly(1.0*x**2 + 1.0, x, domain='CC')"),
    ('Poly(x**2 + I, x, domain=CC).LC()',
     '1.00000000000000'),
    ('Poly(x + 1, x, domain=CC).domain',
     'CC'),
    ('Poly(x**2 + 1, x, domain=C100)',
     "Poly(1.0*x**2 + 1.0, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC) * Poly(x + I, x, domain=CC)',
     "Poly(1.0*x**3 + 1.0*I*x**2 + 1.0*x + 1.0*I, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC) + Poly(x, x)',
     "Poly(1.0*x**2 + 1.0*x + 1.0, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC) - Poly(x + 1.5, x)',
     "Poly(1.0*x**2 - 1.0*x - 0.5, x, domain='CC')"),
    ('Poly(x**2 - 1, x, domain=CC).gcd(Poly(x - 1, x, domain=CC))',
     "Poly(1.0*x - 1.0, x, domain='CC')"),
    ('Poly(x**2 - 1, x, domain=CC).lcm(Poly(x + 1, x, domain=CC))',
     "Poly(1.0*x**2 - 1.0, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC).div(Poly(x + I, x, domain=CC))',
     "(Poly(1.0*x - 1.0*I, x, domain='CC'), Poly(0, x, domain='CC'))"),
    ('Poly(x**2 + 1, x, domain=CC).rem(Poly(x + 2, x, domain=CC))',
     "Poly(5.0, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC).quo(Poly(2*x + 2, x, domain=CC))',
     "Poly(0.5*x - 0.5, x, domain='CC')"),
    ('Poly(x + I, x, domain=CC).invert(Poly(x**2 + 2, x, domain=CC))',
     "Poly(-1.0*x + 1.0*I, x, domain='CC')"),
    ('Poly(x - 1, x, domain=CC).invert(Poly(x**2 + 1, x, domain=CC))',
     "Poly(-0.5*x - 0.5, x, domain='CC')"),
    ('Poly(x + I, x, domain=CC).gcdex(Poly(x**2 + 2, x, domain=CC))',
     "(Poly(-1.0*x + 1.0*I, x, domain='CC'), Poly(1.0, x, domain='CC'), Poly(1.0, x, domain='CC'))"),
    ('Poly(x + I, x, domain=CC).half_gcdex(Poly(x**2 + 2, x, domain=CC))',
     "(Poly(-1.0*x + 1.0*I, x, domain='CC'), Poly(1.0, x, domain='CC'))"),
    ('Poly(x**2 + 2*x + 1, x, domain=CC).sqf_list()',
     "(1.00000000000000, [(Poly(1.0*x + 1.0, x, domain='CC'), 2)])"),
    ('Poly(x**2 + 1, x, domain=CC).factor_list()',
     "(1, [(Poly(1.0*x**2 + 1.0, x, domain='CC'), 1)])"),
    ('Poly(3, x, domain=CC).factor_list()',
     '(3.00000000000000, [])'),
    ('Poly(x**2 + 2*I*x - 1, x, domain=CC).sqf_part()',
     "Poly(1.0*x + 1.0*I, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC).monic()',
     "Poly(1.0*x**2 + 1.0, x, domain='CC')"),
    ('Poly(2*x**2 + 1, x, domain=CC).monic()',
     "Poly(1.0*x**2 + 0.5, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC).diff(x)',
     "Poly(2.0*x, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=CC).eval(I)',
     '0'),
    ('Poly(x**2 + 1, x, domain=CC).nroots()',
     '[-1.0*I, 1.0*I]'),
    ('Poly(x**2 + 1, x, domain=CC).as_expr()',
     '1.0*x**2 + 1.0'),
    ('Poly(x**2 + 1, x, domain=CC)**2',
     "Poly(1.0*x**4 + 2.0*x**2 + 1.0, x, domain='CC')"),
    ('Poly(x**2 + 1, x, domain=C100).gcd(Poly(x + I, x, domain=C100))',
     "Poly(1.0*x + 1.0*I, x, domain='CC')"),
    ('Poly(x/3 + I, x, domain=C100).LC()',
     '0.33333333333333333333333333333'),
    ('Poly(x**2 + I, x, domain=QQ_I).set_domain(CC)',
     "Poly(1.0*x**2 + 1.0*I, x, domain='CC')"),
    ("Poly(x**2 + 1, x, domain='CC')",
     "Poly(1.0*x**2 + 1.0, x, domain='CC')"),
    ('Poly(1.0*x + I, x).invert(Poly(x**2 + 2, x))',
     "Poly(-1.0*x + 1.0*I, x, domain='CC')"),
    ('Poly(x + 0.5*I, x).gcdex(Poly(x**2 + 1, x))',
     "(Poly(-1.33333333333333*x + 0.666666666666667*I, x, domain='CC'), Poly(1.33333333333333, x, domain='CC'), Poly(1.0, x, domain='CC'))"),
    ('Poly(x + 0.5*I, x).invert(Poly(x**2 + 2, x))',
     "Poly(-0.571428571428571*x + 0.285714285714286*I, x, domain='CC')"),
    ('Poly(1.0*x + I, x).gcd(Poly(x**2 + 1, x))',
     "Poly(1.0*x + 1.0*I, x, domain='CC')"),
    ('Poly(1.0*x**2 + 2.0*I*x - 1.0, x).sqf_list()',
     "(1.00000000000000, [(Poly(1.0*x + 1.0*I, x, domain='CC'), 2)])"),
    ('Poly(x + 1, x, domain=CC) * Poly(x + I, x)',
     "Poly(1.0*x**2 + (1.0 + 1.0*I)*x + 1.0*I, x, domain='CC')"),
    ('Poly(1.0*x + I, x).lcm(Poly(x + I, x))',
     "Poly(1.0*x + 1.0*I, x, domain='CC')"),
    ('gcd(x**2 - 1, x - 1, domain=CC)',
     '1.0*x - 1.0'),
    ('gcd(x**2 - 1, x - 1, domain=RR)',
     '1.0*x - 1.0'),
    ('gcd(x**2 - 1, 2*x - 2, domain=RR)',
     '1.0*x - 1.0'),
    ('lcm(x**2 - 1, x - 1, domain=RR)',
     '1.0*x**2 - 1.0'),
    ('lcm(x**2 - 1, x + 1, domain=CC)',
     '1.0*x**2 - 1.0'),
    ('factor(x**2 + 1, domain=CC)',
     '1.0*x**2 + 1.0'),
    ('factor(x**2 - 2, domain=RR)',
     '2.0*(0.5*x**2 - 1.0)'),
    ('factor(2*x**2 - 4, domain=RR)',
     '4.0*(0.5*x**2 - 1.0)'),
    ('factor((x - 1)**2*(x + 2), domain=RR)',
     '2.0*(0.5*x + 1.0)*(1.0*x - 1.0)**2'),
    ('factor_list(x**2 + 1, domain=CC)',
     '(1, [(1.0*x**2 + 1.0, 1)])'),
    ('sqf(x**2 + 2*x + 1, domain=RR)',
     '1.0*(1.0*x + 1.0)**2'),
    ('sqf_list(x**2 + 2*x + 1, domain=CC)',
     '(1.00000000000000, [(1.0*x + 1.0, 2)])'),
    ('div(x**2 + 1, x + I, domain=CC)',
     '(1.0*x - 1.0*I, 0)'),
    ('gcdex(x**2 + 1, x - 1, domain=RR)',
     '(0.500000000000000, -0.5*x - 0.5, 1.00000000000000)'),
    ('invert(x - 1, x**2 + 1, domain=RR)',
     '-0.5*x - 0.5'),
    ('invert(x + I, x**2 + 2, domain=CC)',
     '-1.0*x + 1.0*I'),
    ('gcdex(x + I, x**2 + 2, domain=CC)',
     '(-1.0*x + 1.0*I, 1.00000000000000, 1.00000000000000)'),
    ('gcd(x**2 - 1, x - 1, domain=C100)',
     '1.0*x - 1.0'),
    ('factor(x**2 - 2, domain=R100)',
     '2.0*(0.5*x**2 - 1.0)'),
    ('gcd(x*y - y, x - 1, domain=RR)',
     '1.0*x - 1.0'),
    ('factor(x*y + y, domain=RR)',
     '1.0*y*(1.0*x + 1.0)'),
    ('factor(x**2 + 1, x, domain=CC)',
     '1.0*x**2 + 1.0'),
    ('Poly(1.5*x*y + 1.5*y, x, z).gcd(Poly(x*y + y, x, z))',
     "Poly(0.5*y*x + 0.5*y, x, z, domain='RR[y]')"),
    ('sqrt((-1 - 2*I)**2)',
     '1 + 2*I'),
    ('sqrt((1 + 2*I)**2)',
     '1 + 2*I'),
    ('sqrt((-1 + 2*I)**2)',
     '1 - 2*I'),
    ('((1 + I)**2)**Rational(3, 2)',
     '(1 + I)**3'),
    ('((-1 - I)**2)**Rational(3, 2)',
     '-(-1 - I)**3'),
    ('sqrt(CRootOf(x**3 - 2, 1)**2)',
     '-CRootOf(x**3 - 2, 1)'),
    ('sqrt(CRootOf(x**3 - 2, 2)**2)',
     '-CRootOf(x**3 - 2, 2)'),
    ('sqrt(CRootOf(x**5 - x + 1, 1)**2)',
     '-CRootOf(x**5 - x + 1, 1)'),
    ('sqrt(CRootOf(x**5 - x + 1, 3)**2)',
     'CRootOf(x**5 - x + 1, 3)'),
    ('sqrt(CRootOf(x**4 + 1, 0)**2)',
     '-CRootOf(x**4 + 1, 0)'),
    ('sqrt(CRootOf(x**4 + 1, 3)**2)',
     'CRootOf(x**4 + 1, 3)'),
    ('sqrt((1 - sqrt(2))**3)',
     '-(1 - sqrt(2))**(3/2)'),
    ('((1 - sqrt(2))**3)**Rational(3, 2)',
     '-(1 - sqrt(2))**(9/2)'),
    ('sqrt((-1 + sqrt(2)*I)**3)',
     '-(-1 + sqrt(2)*I)**(3/2)'),
    ('((-2 + I)**3)**Rational(-1, 2)',
     '-1/(-2 + I)**(3/2)'),
    ('((-2 - I)**4)**Rational(1, 2)',
     '(-2 - I)**2'),
    ('sqrt((pi*I - 1)**2)',
     '1 - I*pi'),
    ('sqrt((-1 - 2*I)**3)',
     '-(-1 - 2*I)**(3/2)'),
]


@pytest.mark.parametrize('case, expected', CASES)
def test_ccdom_case_matches_reference(case, expected):
    assert str(eval(case)) == expected


def test_zz_i_sqf_list_unit():
    # Intentional correction: SymPy 1.14 returns the unit 1 here, but
    # (I*x + y)**2 == -(-x + I*y)**2, so the unit is -1.
    f = Poly((I*x + y)**2, x, y)
    c, factors = f.sqf_list()
    assert str((c, factors)) == "(-1, [(Poly(-x + I*y, x, y, domain='ZZ_I'), 2)])"
    product = c*Mul(*[g.as_expr()**k for g, k in factors])
    assert expand(product - f.as_expr()) == 0
    point = {x: Rational(3, 7), y: Rational(-2, 5)}
    assert expand(product.subs(point) - f.as_expr().subs(point)) == 0
    assert expand(Mul(*[g.as_expr()**k for g, k in factors]).subs(point)
                  - f.as_expr().subs(point)) != 0


def test_multivariate_composite_sqf_list():
    f = Poly((x + I*y)**2*(x - z)*z, x, z)
    c, factors = f.sqf_list()
    assert expand(c*Mul(*[g.as_expr()**k for g, k in factors]) - f.as_expr()) == 0
    assert [k for _, k in factors] == [1, 2]


def test_half_branch_numeric():
    # (b**a)**(p/2) == s*b**(a*p/2) with the reference branch factor s.
    for base, a, e in [(-1 - 2*I, 2, S.Half), (-1 + sqrt(2)*I, 3, S.Half),
                       (-2 + I, 3, Rational(-1, 2)), (1 - sqrt(2), 3, Rational(3, 2)),
                       (CRootOf(x**5 - x + 1, 1), 2, S.Half)]:
        lhs = Pow(Pow(base, a, evaluate=False), e, evaluate=False)
        rhs = Pow(base, a)**e
        assert abs(N(lhs - rhs, 30)) < 1e-25


def test_inexact_options_unknown_keys_use_native_path():
    assert str(gcd(x**2 - 1, x - 1, domain=QQ_I)) == 'x - 1'
    assert str(gcd(2*x + 4*I, x**2 + 4, domain=QQ_I)) == 'x + 2*I'
