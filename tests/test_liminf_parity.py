"""Limits of special functions at infinite arguments (tractable rewrites
and asymptotic series), their series at finite points, and the integrals
that depend on them. Expected values are those of SymPy 1.14, except the
intentional corrections checked numerically below."""

import pytest

from symox import *
del test

x = Symbol('x')
a = Symbol('a', positive=True)

SAME_AS_SYMPY = [
    ('limit(Si(x), x, oo)', 'pi/2'),
    ('limit(Si(x), x, -oo)', '-pi/2'),
    ('limit(Ci(x), x, oo)', '0'),
    ('limit(Shi(x), x, oo)', 'oo'),
    ('limit(Shi(x), x, -oo)', '-oo'),
    ('limit(Chi(x), x, oo)', 'oo'),
    ('limit(Ei(x), x, oo)', 'oo'),
    ('limit(Ei(x), x, -oo)', '0'),
    ('limit(Ei(-x), x, oo)', '0'),
    ('limit(erf(x), x, oo)', '1'),
    ('limit(erf(x), x, -oo)', '-1'),
    ('limit(erfc(x), x, oo)', '0'),
    ('limit(erfc(x), x, -oo)', '2'),
    ('limit(erfi(x), x, oo)', 'oo'),
    ('limit(li(x), x, oo)', 'oo'),
    ('limit(fresnels(x), x, oo)', '1/2'),
    ('limit(fresnelc(x), x, oo)', '1/2'),
    ('limit(fresnels(x), x, -oo)', '-1/2'),
    ('limit(fresnelc(x), x, -oo)', '-1/2'),
    ('limit(uppergamma(2, x), x, oo)', '0'),
    ('limit(uppergamma(a, x), x, oo)', '0'),
    ('limit(lowergamma(3, x), x, oo)', '2'),
    ('limit(Li(x), x, oo)', 'oo'),
    ('limit(x*exp(-x)*Ei(x), x, oo)', '1'),
    ('limit(x*exp(x)*Ei(-x), x, oo)', '-1'),
    ('limit(x*exp(x**2)*erfc(x), x, oo)', '1/sqrt(pi)'),
    ('limit(x*exp(x**2)*(erf(x) - 1), x, oo)', '-1/sqrt(pi)'),
    ('limit(li(x)*log(x)/x, x, oo)', '1'),
    ('limit(uppergamma(3, x)*exp(x)/x**2, x, oo)', '1'),
    ('limit(exp(x**2)*erfc(x), x, oo)', '0'),
    ('limit(Ei(x)/exp(x), x, oo)', '0'),
    ('limit(erf(x)**x, x, oo)', '1'),
    ('limit(x*(1 - erf(x)), x, oo)', '0'),
    ('limit(erfi(x)*exp(-x**2)*x, x, oo)', '1/sqrt(pi)'),
    ('limit(Ei(-log(x)), x, oo)', '0'),
    ('limit(Ei(-log(x)) - log(log(x))/x, x, oo)', '0'),
    ("limit(erfc(sqrt(log(1/x))), x, 0, '+')", '0'),
    ("limit(li(x) - li(2*x), x, 0, '+')", '0'),
    ("limit(li(x), x, 0, '+')", '0'),
    ("limit(Si(1/x), x, 0, '+')", 'pi/2'),
    ("limit(Si(1/x), x, 0, '-')", '-pi/2'),
    ("limit(erf(1/x), x, 0, '+')", '1'),
    ("limit(erf(1/x), x, 0, '-')", '-1'),
    ("limit(Ei(1/x), x, 0, '-')", '0'),
    ('limit(Si(x)/x, x, 0)', '1'),
    ("limit(Ci(x) - log(x), x, 0, '+')", 'EulerGamma'),
    ('limit(erf(x)/x, x, 0)', '2/sqrt(pi)'),
    ("limit(Ei(x) - log(x), x, 0, '+')", 'EulerGamma'),
    ('limit(Si(x) - pi/2, x, oo)', '0'),
    ('limit(Si(x**2), x, oo)', 'pi/2'),
    ('limit(Si(exp(x)), x, oo)', 'pi/2'),
    ('limit(erf(log(x)), x, oo)', '1'),
    ('limit(Ei(2*x)/Ei(x), x, oo)', 'oo'),
    ('limit(Si(x) + erf(x), x, oo)', '1 + pi/2'),
    ('limit(Ei(x) - Ei(x - 1), x, oo)', 'oo'),
    ('limit(li(x) - li(x/2), x, oo)', 'oo'),
    ('limit(1/Ei(x), x, oo)', '0'),
    ('limit(erfc(x)**(1/x**2), x, oo)', 'exp(-1)'),
    ('limit(fresnels(x)*x, x, oo)', 'oo'),
    ('limit(erf(x)*x - x, x, oo)', '0'),
    ('limit(Si(x)*Ci(x), x, oo)', '0'),
]


@pytest.mark.parametrize("case,expected", SAME_AS_SYMPY)
def test_limit_matches_sympy(case, expected):
    assert str(eval(case)) == expected


# SymPy 1.14 is wrong on these (or leaves them unevaluated); the values
# below are the true limits, checked with mpmath.
CORRECTIONS = [
    ("limit(expint(1, x), x, oo)", "0"),          # SymPy: expint(1, oo)
    ("limit(expint(2, x), x, oo)", "0"),          # SymPy: -oo*I
    ("limit(expint(3, x), x, oo)", "0"),          # SymPy: oo*I
    ("limit(Chi(x)*x*exp(-x), x, oo)", "1/2"),    # SymPy: 0
    ("limit(Shi(x)*x*exp(-x), x, oo)", "1/2"),    # SymPy: 0
    ("limit(expint(a, x)*x*exp(x), x, oo)", "1"),  # SymPy: oo*sign(...)
    ("limit(x*exp(x)*expint(1, x), x, oo)", "1"),  # SymPy: unevaluated
    ("limit(erfc(x)/erfc(2*x), x, oo)", "oo"),    # SymPy: unevaluated
    ("limit(Shi(x) - Chi(x), x, oo)", "0"),       # SymPy: unevaluated
]


@pytest.mark.parametrize("case,expected", CORRECTIONS)
def test_limit_corrections(case, expected):
    assert str(eval(case)) == expected


def test_corrections_numerically():
    import mpmath
    with mpmath.workdps(30):
        X = mpmath.mpf(200)
        assert abs(mpmath.chi(X)*X*mpmath.exp(-X) - mpmath.mpf(1)/2) < 0.01
        assert abs(mpmath.shi(X)*X*mpmath.exp(-X) - mpmath.mpf(1)/2) < 0.01
        assert abs(X*mpmath.exp(X)*mpmath.expint(1, X) - 1) < 0.01
        assert abs(X*mpmath.exp(X)*mpmath.expint(mpmath.mpf(5)/2, X) - 1) < 0.02
        assert mpmath.expint(2, X) < mpmath.mpf(10)**-80


SERIES = [
    ("series(Si(x),x,0,8)", "x - x**3/18 + x**5/600 - x**7/35280 + O(x**8)"),
    ("series(Shi(x),x,0,8)", "x + x**3/18 + x**5/600 + x**7/35280 + O(x**8)"),
    ("series(erfi(x),x,0,8)", "2*x/sqrt(pi) + 2*x**3/(3*sqrt(pi)) + x**5/(5*sqrt(pi)) + x**7/(21*sqrt(pi)) + O(x**8)"),
    ("series(fresnels(x),x,0,12)", "pi*x**3/6 - pi**3*x**7/336 + pi**5*x**11/42240 + O(x**12)"),
    ("series(fresnelc(x),x,0,10)", "x - pi**2*x**5/40 + pi**4*x**9/3456 + O(x**10)"),
    ("series(Si(x+1),x,0,3)", "Si(1) + x*sin(1) + x**2*(-sin(1)/2 + cos(1)/2) + O(x**3)"),
    ("series(Si(2*x**2),x,0,8)", "2*x**2 - 4*x**6/9 + O(x**8)"),
    ("series(erfi(x),x,1,3)", "erfi(1) + 2*E*(x - 1)/sqrt(pi) + 2*E*(x - 1)**2/sqrt(pi) + O((x - 1)**3, (x, 1))"),
    ("limit(Si(x)/x,x,0)", "1"),
    ("limit(Shi(x)/x,x,0)", "1"),
    ("limit(fresnels(x)/x**3,x,0)", "pi/6"),
    ("limit(fresnelc(x)/x,x,0)", "1"),
    ("limit(erfi(x)/x,x,0)", "2/sqrt(pi)"),
    ("limit((Si(x)-x)/x**3,x,0)", "-1/18"),
]


@pytest.mark.parametrize("case,expected", SERIES)
def test_series_matches_sympy(case, expected):
    assert str(eval(case)) == expected


def test_integrals_through_special_limits():
    assert integrate(log(log(x))/x**2, (x, 1, oo)) == -EulerGamma
    assert integrate(sqrt(log(1/x)), (x, 0, 1)) == sqrt(pi)/2
    assert limit(li(x) - li(2*x), x, 0) == 0


def test_fraction_field_powers_of_exp_generator():
    t = Symbol('t')
    p = Poly(exp(S(1)/2)/(1 + E + 2*exp(S(1)/2)), t, domain='ZZ(exp(1/2))')
    assert str(p) == "Poly(exp(1/2)/(1 + E + 2*exp(1/2)), t, domain='ZZ(exp(1/2))')"
    p = Poly(E*t + exp(-S(3)/2), t, domain='QQ(exp(1/2))')
    assert str(p) == "Poly(E*t + exp(-3/2), t, domain='QQ(exp(1/2))')"
    p = Poly(E*t + exp(S(1)/3), t, domain='QQ(exp(1/6))')
    assert str(p) == "Poly(E*t + exp(1/3), t, domain='QQ(exp(1/6))')"


def test_risch_exp_half_and_E_one_transcendental():
    f = 1/(1 + exp(x - S(1)/2)/(1 + exp(x)))
    r = integrate(f, x)
    assert str(r) == "x - exp(1/2)*log(exp(x) + exp(1/2)/(1 + exp(1/2)))/(exp(1/2) + E)"
    assert simplify(diff(r, x) - f) == 0
