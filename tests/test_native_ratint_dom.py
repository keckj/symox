"""ratint over Gaussian grounds and EX coefficient growth.

The expected strings are SymPy 1.14 outputs (hard-coded); each value is also
checked by differentiating it back to the integrand at exact rational points.
"""
import mpmath
import pytest

import symox as sx
from symox.integrals.rationaltools import ratint, ratint_logpart, ratint_ratpart

x = sx.Symbol('x')
a, b = sx.symbols('a b')


def assert_primitive(f, r, points=('-7/3', '1/5', '3/2', '11/4')):
    d = r.diff(x)
    with mpmath.workdps(40):
        for p in points:
            v = sx.Rational(p)
            err = abs(sx.N(d.subs(x, v) - f.subs(x, v), 40))
            assert err < 1e-30, (p, err)


# SymPy 1.14 `str(ratint(f, x))`; the reference crashed in symox before
# (`incompatible coefficient fields` on `ZZ[a...]` against `QQ_I`).
GAUSSIAN = [
    ('(x**2 + I)**(-2)',
     '-I*x/(2*x**2 + 2*I) + sqrt(-I)*log(x - sqrt(-I))/4 - sqrt(-I)*log(x + sqrt(-I))/4'),
    ('x/(x**2 + 2*I*x - 3)',
     '(1/2 - sqrt(2)*I/4)*log(x - 4*I*(1/2 - sqrt(2)*I/4) + 3*I) + (1/2 + sqrt(2)*I/4)*log(x - 4*I*(1/2 + sqrt(2)*I/4) + 3*I)'),
    ('(x + I)/(x**3 - I)',
     'sqrt(3)*log(2*x - sqrt(3) - I)/3 - sqrt(3)*log(2*x + sqrt(3) - I)/3'),
    ('(x - I)**(-3)', '-1/(2*x**2 - 4*I*x - 2)'),
    ('I/(x**2 + 1)**2', 'I*x/(2*x**2 + 2) + log(x - I)/4 - log(x + I)/4'),
    ('1/((x + 2)*(x - I))',
     'I*(1 + 2*I)*log(2*x + 2 - I + I*(-3 - 4*I)*(1 + 2*I)/5)/5 - I*(1 + 2*I)*log(2*x + 2 - I*(-3 - 4*I)*(1 + 2*I)/5 - I)/5'),
]


@pytest.mark.parametrize('source,expected', GAUSSIAN)
def test_gaussian_ratint_matches_reference(source, expected):
    f = sx.sympify(source, locals={'x': x})
    r = ratint(f, x)
    assert str(r) == expected
    assert_primitive(f, r)


def test_gaussian_quadratic_value_and_logpart():
    # SymPy: -sqrt(I)*log(x - sqrt(I)*I)/2 + sqrt(I)*log(x + sqrt(I)*I)/2;
    # symox's core Mul folds `sqrt(I)*I` into `I**(3/2)` (documented gap).
    f = 1/(x**2 + sx.I)
    r = ratint(f, x)
    s = sx.sqrt(sx.I)
    assert sx.expand(r - (-s*sx.log(x - s*sx.I)/2 + s*sx.log(x + s*sx.I)/2)) == 0
    assert_primitive(f, r)
    assert str(ratint_logpart(sx.Poly(1, x, domain='QQ'), sx.Poly(x**2 + sx.I, x, domain='QQ_I'), x)) == \
        "[(Poly(x + 2*_t*I, x, domain='ZZ_I[_t]'), Poly(4*I*_t**2 + 1, _t, domain='ZZ_I'))]"


def test_gaussian_ratpart():
    g, h = ratint_ratpart(sx.Poly(1, x, domain='QQ'), sx.Poly((x**2 + sx.I)**2, x, domain='QQ_I'), x)
    assert str(g) == '-I*x/(2*x**2 + 2*I)'
    assert str(h) == '-I/(2*x**2 + 2*I)'


@pytest.mark.parametrize('expr,expected', [
    ('A*v', "Poly(a*x**3 + b*x**2 + I*a*x + I*b, x, domain='QQ_I[a,b]')"),
    ('v*A', "Poly(a*x**3 + b*x**2 + I*a*x + I*b, x, domain='QQ_I[a,b]')"),
    ('A+v', "Poly(x**2 + a*x + b + I, x, domain='QQ_I[a,b]')"),
    ('v-A', "Poly(x**2 - a*x - b + I, x, domain='QQ_I[a,b]')"),
    ('A*w', "Poly(2*a*x**2 + (I*a + 2*b)*x + I*b, x, domain='ZZ_I[a,b]')"),
    ('w-A', "Poly((2 - a)*x - b + I, x, domain='ZZ_I[a,b]')"),
    ('R*v', "Poly(x**3 + a*x**2 + I*x + I*a, x, domain='QQ_I(a)')"),
    ('v+R', "Poly(x**2 + x + a + I, x, domain='QQ_I(a)')"),
])
def test_gaussian_against_composite_ground_unifies(expr, expected):
    A = sx.Poly(a*x + b, x, domain=sx.ZZ[a, b])
    v = sx.Poly(x**2 + sx.I, x, domain='QQ_I')
    w = sx.Poly(2*x + sx.I, x, domain='ZZ_I')
    R = sx.Poly(x + a, x, domain='QQ(a)')
    assert repr(eval(expr, {'A': A, 'v': v, 'w': w, 'R': R})) == expected


# SymPy 1.14 `len(str(ratint(f, x)))`; symox was 82704 for the first case
# before the EX `simplify` applied the reference `signsimp`/`factor_terms`.
# The third case is 9261 since the EX quadratic roots of `log_to_real` are
# simplified as in SymPy: its r_v now has SymPy's form, but the native EX
# cancel of A**2 + B**2 does not reach SymPy's normal form of the 2**(1/3)
# denominators (log term 1549 characters, SymPy 1424). Tracked in
# docs/algorithmic-parity.md (exroots).
GROWTH = [('1/(x**3+x+1)', 16653), ('1/(x**3-x+1)', 16611),
          ('1/((x**2+1)*(x**3+2))', 9261)]


@pytest.mark.parametrize('source,reference_length', GROWTH)
def test_cubic_log_to_real_does_not_grow(source, reference_length):
    f = sx.sympify(source, locals={'x': x})
    r = ratint(f, x)
    assert len(str(r)) <= reference_length
    assert r.has(sx.atan) and r.has(sx.log)
    assert_primitive(f, r)



# SymPy 1.14 returns 0 for 1/(x**8 + 1) and drops the complex-pair terms of
# 1/(x**6 - 2) (its derivative misses 0.3175... at x = 1/3): `_get_real_roots`
# trusts an incomplete `roots(..., filter='R')` when `count_roots` raises
# DomainError over EX. symox declines that step, so log_to_real gives up and
# ratint keeps the exact RootSum (the SymPy real=False value).
@pytest.mark.parametrize('source,expected', [
    ('1/(x**8+1)', 'RootSum(16777216*_t**8 + 1, Lambda(_t, _t*log(8*_t + x)))'),
    ('1/(x**6-2)', 'RootSum(1492992*_t**6 - 1, Lambda(_t, _t*log(-12*_t + x)))'),
])
def test_incomplete_real_roots_keep_exact_rootsum(source, expected):
    f = sx.sympify(source, locals={'x': x})
    r = ratint(f, x)
    assert str(r) == expected
    assert str(ratint(f, x, real=False)) == expected
    g = sx.lambdify(x, f, 'mpmath')
    with mpmath.workdps(30):
        v = sx.N(r.subs(x, sx.Rational(1, 2)), 30) - sx.N(r.subs(x, sx.Rational(1, 10)), 30)
        assert abs(v - mpmath.quad(g, [mpmath.mpf(1)/10, mpmath.mpf(1)/2])) < 1e-25


def test_parametric_real_roots_keep_reference_behaviour():
    ar = sx.Symbol('a', real=True)
    assert sx.integrate(1/(x**2 + ar**2), (x, 0, 1)) == sx.atan(1/ar)/ar
