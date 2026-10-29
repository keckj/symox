"""ratint growth family: SymPy 1.14 forms, cubic size, structural declines.

Expected strings are hard-coded from SymPy 1.14 `ratint(f, x)`.
"""
import importlib

import mpmath
import pytest
import symox as sx

rt = importlib.import_module('symox.integrals.rationaltools')

SYMPY = {
    '1/(x**6+8)': '-sqrt(6)*log(x**2 - sqrt(6)*x + 2)/96 + sqrt(6)*log(x**2 + sqrt(6)*x + 2)/96 + sqrt(2)*atan(sqrt(2)*x/2)/24 + sqrt(2)*atan(sqrt(2)*x - sqrt(3))/48 + sqrt(2)*atan(sqrt(2)*x + sqrt(3))/48',
    '1/(x**3-3*x+1)': 'RootSum(81*_t**3 - 9*_t + 1, Lambda(_t, _t*log(-54*_t**2 - 9*_t + x + 4)))',
    '1/(x**4+x+1)': 'RootSum(229*_t**4 + 18*_t**2 + 8*_t + 1, Lambda(_t, _t*log(2061*_t**3/64 - 687*_t**2/64 + 391*_t/64 + x + 27/64)))',
    '1/(x**4-4*x**2+x+1)': 'RootSum(1957*_t**4 - 174*_t**2 + 8*_t + 1, Lambda(_t, _t*log(22557*_t**3/64 + 721*_t**2/64 - 32235*_t/1216 + x + 705/1216)))',
    'x/(x**3-3*x+1)': 'RootSum(81*_t**3 - 9*_t - 1, Lambda(_t, _t*log(27*_t**2 - 9*_t + x - 2)))',
    '1/(x**3-4*x+1)': 'RootSum(229*_t**3 - 12*_t + 1, Lambda(_t, _t*log(-1832*_t**2/9 - 229*_t/9 + x + 64/9)))',
    '1/((x-1)**4+1)': '-sqrt(2)*log(x**2 + x*(-2 - sqrt(2)) + sqrt(2) + 2)/8 + sqrt(2)*log(x**2 + x*(-2 + sqrt(2)) - sqrt(2) + 2)/8 - sqrt(2)*atan(-sqrt(2)*x + 1 + sqrt(2))/4 + sqrt(2)*atan(sqrt(2)*x - sqrt(2) + 1)/4',
    '1/(x**4-2*x**2+4)': '-sqrt(6)*log(x**2 - sqrt(6)*x + 2)/48 + sqrt(6)*log(x**2 + sqrt(6)*x + 2)/48 + sqrt(2)*atan(sqrt(2)*x - sqrt(3))/8 + sqrt(2)*atan(sqrt(2)*x + sqrt(3))/8',
    '1/(x**6+27)': '-log(x**2 - 3*x + 3)/108 + log(x**2 + 3*x + 3)/108 + sqrt(3)*atan(sqrt(3)*x/3)/81 + sqrt(3)*atan(2*sqrt(3)*x/3 - sqrt(3))/162 + sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3))/162',
    '1/(x**4+2)': '-2**(3/4)*log(x**2 - 2**(3/4)*x + sqrt(2))/16 + 2**(3/4)*log(x**2 + 2**(3/4)*x + sqrt(2))/16 + 2**(3/4)*atan(2**(1/4)*x - 1)/8 + 2**(3/4)*atan(2**(1/4)*x + 1)/8',
    'x**2/(x**6+8)': 'sqrt(2)*atan(sqrt(2)*x**3/4)/12',
    '1/(x**4+x**2+1)': '-log(x**2 - x + 1)/4 + log(x**2 + x + 1)/4 + sqrt(3)*atan(2*sqrt(3)*x/3 - sqrt(3)/3)/6 + sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3)/3)/6',
}


def groups(source):
    x, t = sx.symbols('x t')
    f = sx.sympify(source, locals={'x': x})
    n, d = sx.fraction(sx.cancel(f))
    _, r = sx.Poly(n, x).div(sx.Poly(d, x))
    pairs = rt.ratint_logpart(r, sx.Poly(d, x), x, t)
    return x, t, [(h, q) for h, q in pairs if q.degree() >= 3]


@pytest.mark.parametrize('source', sorted(SYMPY))
def test_public_ratint_prints_the_sympy_form(source):
    x = sx.Symbol('x')
    f = sx.sympify(source, locals={'x': x})
    assert str(rt.ratint(f, x)) == SYMPY[source]


# Non-real roots outside every QQ(i, sqrt(k)): a quartic factor whose root
# field does not contain i, a factor of degree >= 5, or a quartic with real
# roots. These decline before the resultant (same result, Python path).
@pytest.mark.parametrize('source', ['1/(x**6+8)', '1/(x**4-2*x**2+4)', '1/(x**4+2)',
                                    '1/(x**8+1)', '1/(x**4+x+1)'])
def test_nonreal_field_groups_decline(source):
    x, t, pairs = groups(source)
    assert pairs
    for h, q in pairs:
        assert sx._core._log_to_real(h, q, x, t) is None


# Fields QQ(i, sqrt(k)) and quadratic factors stay native.
@pytest.mark.parametrize('source', ['1/((x-1)**4+1)', '1/(x**4-x**2+1)', '1/(x**6+27)',
                                    '1/(x**4+1)', '1/(x**6+1)'])
def test_gaussian_biquadratic_groups_stay_native(source):
    x, t, pairs = groups(source)
    assert pairs
    for h, q in pairs:
        assert sx._core._log_to_real(h, q, x, t) is not None


def test_x8_plus_16_keeps_the_exact_rootsum():
    # Intentional correction: SymPy 1.14 returns 0 (incomplete real roots).
    x = sx.Symbol('x')
    f = 1/(x**8 + 16)
    value = rt.ratint(f, x)
    assert value.has(sx.RootSum)
    derivative = sx.diff(value, x) - f
    for point in (sx.Rational(1, 3), sx.Rational(-7, 5)):
        assert abs(complex(derivative.subs(x, point).evalf(30))) < 1e-20


def test_cubic_casus_form_is_smaller_and_equal():
    x = sx.Symbol('x')
    f = 1/(x**3 + x + 1)
    value = rt.ratint(f, x)
    assert sx.integrate(f, x) == value
    # SymPy 1.14 prints 16653 characters; the native EX cancel gives 7892.
    assert len(str(value)) < 16653
    assert not value.has(sx.RootSum, sx.Integral)
    logs = sorted(value.atoms(sx.log), key=lambda e: sx.Poly(e.args[0], x).degree())
    assert [sx.Poly(e.args[0], x).degree() for e in logs] == [1, 2]
    assert len(value.atoms(sx.atan)) == 1
    with mpmath.workdps(60):
        # log(x**2 + b*x + c) is |x - r|**2 for the complex roots r: b is the real root.
        quadratic = sx.Poly(logs[1].args[0], x)
        b = quadratic.nth(1)/quadratic.LC()
        assert abs(sx.N(b, 50) - sx.Float('-0.6823278038280193273694837397110482568912', 50)) < 1e-38
        definite = sx.N(value.subs(x, 2) - value.subs(x, 1), 50)
        assert abs(definite - sx.Float('0.1841354443888388763172177529151723566604', 50)) < 1e-38
        derivative = sx.diff(value, x).subs(x, sx.Rational(3, 7)) - f.subs(x, sx.Rational(3, 7))
        assert abs(complex(sx.N(derivative, 50))) < 1e-40


def test_atan_parametric_definite_matches_sympy():
    x, a = sx.symbols('x a')
    assert str(sx.integrate(sx.atan(a*x)/x**2, (x, 1, sx.oo))) == \
        '-a*log(a**2)/2 + a*log(a**2 + 1)/2 + atan(a)'
    b = sx.Symbol('b', positive=True)
    value = sx.integrate(sx.atan(b*x)/x**2, (x, 1, sx.oo))
    # SymPy 1.14: -b*log(b) + b*log(b**2 + 1)/2 + atan(b) (equal for b > 0).
    expected = -b*sx.log(b) + b*sx.log(b**2 + 1)/2 + sx.atan(b)
    for point in (sx.Rational(1, 3), 2, 7):
        assert abs(sx.N((value - expected).subs(b, point), 40)) < 1e-35
