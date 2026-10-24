"""Limits with AccumBounds, symbolic oscillation, complex infinities and the
Piecewise/Heaviside leftovers (wave 11, limacc).

Expected values are SymPy 1.14's, unless a test says otherwise.
"""
import pytest

from symox import *  # noqa: F401,F403
from symox import AccumBounds

del test  # noqa: F821

x, y, z = symbols('x y z')
p = Symbol('p', positive=True)
s = Symbol('s', nonzero=True, real=True)
r = Symbol('r', real=True)


@pytest.mark.parametrize('expr, point, expected', [
    (x*Ci(x), oo, 'AccumBounds(-1, 1)'),
    (x*(Si(x) - pi/2), oo, 'AccumBounds(-1, 1)'),
    (x*(fresnels(x) - S(1)/2), oo, 'AccumBounds(-1, 1)/pi'),
    (x*(fresnelc(x) - S(1)/2), oo, 'AccumBounds(-1, 1)/pi'),
    (x*(Si(x) + pi/2), -oo, 'AccumBounds(-1, 1)'),
    (x*(fresnels(x) + S(1)/2), -oo, 'AccumBounds(-1, 1)/pi'),
    (x*Ci(x) + sin(x), oo, 'AccumBounds(-2, 2)'),
    (x*Ci(x)*2 + 1, oo, 'AccumBounds(-1, 3)'),
    (x*Ci(2*x), oo, 'AccumBounds(-1/2, 1/2)'),
    (exp(x)*Ci(exp(x)), oo, 'AccumBounds(-1, 1)'),
    (x*Ci(x) + 1/x, oo, 'AccumBounds(-1, 1)'),
    (p*x*Ci(x), oo, 'p*AccumBounds(-1, 1)'),
    (x*Ci(x)*(1 + 1/x), oo, 'AccumBounds(-1, 1)'),
    (sqrt(x)*(fresnelc(x) - S(1)/2), oo, '0'),
    (sin(x) + cos(x), oo, 'AccumBounds(-2, 2)'),
])
def test_accumbounds_affine_limits(expr, point, expected):
    assert str(limit(expr, x, point)) == expected


def test_accumbounds_at_zero():
    assert str(limit(Ci(1/x)/x, x, 0, '+')) == 'AccumBounds(-1, 1)'


@pytest.mark.parametrize('expr, point, expected', [
    (sin(s*x), oo, 'sin(oo*sign(s))'),
    (cos(s*x), oo, 'AccumBounds(-1, 1)'),
    (x*sin(s*x), oo, 'oo*sin(oo*sign(s))'),
    (sin(s*x) + 1, oo, 'sin(oo*sign(s)) + 1'),
    (sin(s*x)**2, oo, 'sin(oo*sign(s))**2'),
    (sin(s*x)/x, oo, '0'),
    (tan(s*x), oo, 'tan(oo*sign(s))'),
    (sin(s*x), -oo, '-sin(oo*sign(s))'),
    (sin(s*x)*sin(x), oo, 'sin(oo*sign(s))*AccumBounds(-1, 1)'),
])
def test_symbolic_frequency_limits(expr, point, expected):
    assert str(limit(expr, x, point)) == expected


def test_symbolic_frequency_unevaluated():
    # s may be zero for r, and exp(sin(oo*sign(s))) is not a value
    assert isinstance(limit(sin(r*x), x, oo), Limit)
    assert isinstance(limit(exp(sin(s*x)), x, oo), Limit)
    assert str(limit(sin(s/x), x, 0)) == 'sin(oo*sign(s))'


def test_symbolic_frequency_integrals():
    assert str(integrate(cos(s*x), (x, 0, oo))) == 'sin(oo*sign(s))/s'
    assert str(integrate(cos(s*x), (x, -oo, 0))) == 'sin(oo*sign(s))/s'
    assert str(integrate(sin(s*x), (x, 0, oo))) == 'AccumBounds(-1, 1)/s + 1/s'


@pytest.mark.parametrize('expr, args, expected', [
    (I*x, (oo,), 'oo*I'),
    ((1 + I)*x, (oo,), 'oo*sign(1 + I)'),
    (I/x, (0, '-'), '-oo*I'),
    (I/x, (0, '+'), 'oo*I'),
    (I*x**2 + x, (oo,), 'oo*I'),
    (exp(x)*I, (oo,), 'oo*I'),
    (I*x, (-oo,), '-oo*I'),
    (-I*pi*x, (oo,), '-oo*I'),
    (x*Ci(x), (-oo,), '-oo*I'),
    (I*log(x), (oo,), 'oo*I'),
    (sqrt(-x), (oo,), 'oo*I'),
    (I*x/(x + 1)*x, (oo,), 'oo*I'),
    (x*(1 + 2*I) + 1, (oo,), 'oo*sign(1 + 2*I)'),
    (I/x**2, (0,), 'oo*I'),
    (x*Shi(x)*I, (oo,), 'oo*I'),
])
def test_complex_infinite_limits(expr, args, expected):
    assert str(limit(expr, x, *args)) == expected


@pytest.mark.parametrize('expr, n, expected', [
    (Heaviside(x)/x, 6, '1/x'),
    (Heaviside(x)/x**2, 6, 'x**(-2)'),
    (Heaviside(x)*x, 6, 'x'),
    (Heaviside(x)/x + 1, 6, '1 + 1/x'),
    (Heaviside(x)*sin(x)/x**2, 6, '1/x - x/6 + x**3/120 - x**5/5040 + O(x**6)'),
    (Heaviside(x)/x, 3, '1/x'),
    (Heaviside(-x)/x, 6, '0'),
    (Heaviside(x)*exp(x)/x, 3, '1/x + 1 + x/2 + x**2/6 + O(x**3)'),
    (Heaviside(x), 6, '1'),
    (1/x + Heaviside(x), 6, '1 + 1/x'),
])
def test_heaviside_series_exact(expr, n, expected):
    assert str(series(expr, x, 0, n)) == expected


@pytest.mark.parametrize('expr, subs, expected', [
    (DiracDelta(x), {x: 0}, 'DiracDelta(x)'),
    (DiracDelta(x + 1), {x: -1}, 'DiracDelta(x + 1)'),
    # SymPy floats the numbers around DiracDelta(0); symox's evalf keeps an
    # expression holding DiracDelta(0) exact (equivalent form).
    pytest.param(2*DiracDelta(x), {x: 0}, '2.0*DiracDelta(0)', marks=pytest.mark.xfail(strict=True)),
    pytest.param(DiracDelta(x) + 1, {x: 0}, 'DiracDelta(0) + 1.0', marks=pytest.mark.xfail(strict=True)),
    (DiracDelta(x, 1), {x: 0}, 'DiracDelta(x, 1)'),
    (DiracDelta(x), {x: 0.5}, '0'),
    (DiracDelta(x*y), {x: 0}, 'DiracDelta(x*y)'),
    (DiracDelta(x), {y: 0}, 'DiracDelta(x)'),
    (Heaviside(x), {x: 0}, '0.500000000000000'),
    (sin(x), {x: 0}, '0'),
])
def test_dirac_delta_evalf_subs(expr, subs, expected):
    assert str(expr.evalf(subs=subs)) == expected
    assert str(DiracDelta(x).n(subs={x: 0})) == 'DiracDelta(x)'


@pytest.mark.parametrize('expr, old, new, names', [
    (O(x*y), y, x, '[x, x]'),
    (O(x*y), x, y, '[y, y]'),
    (O(x + y), y, x, '[x, x]'),
    (O(x*y) + x, y, x, '[x, x]'),
    (O(x*y, (x, 1), (y, 1)), y, x, '[x, x]'),
])
def test_order_subs_repeated_variables(expr, old, new, names):
    with pytest.raises(ValueError, match=r'Variables are supposed to be unique symbols, got \[%s\]' % names[1:-1]):
        expr.subs(old, new)


def test_order_subs_distinct_variables():
    assert str(O(x*y).subs(y, z)) == 'O(x*z, x, z)'
    assert str(O(x).subs(x, y)) == 'O(y)'
    assert str(O(x).subs(y, x)) == 'O(x)'


def test_periodic_inequality_one_period():
    # SymPy restricts periodic inequalities to one period [0, period)
    assert solveset(sin(x) + cos(x) > 0, x, S.Reals) == Union(
        Interval.Ropen(0, 3*pi/4), Interval.open(7*pi/4, 2*pi))


def test_accumbounds_intentional():
    # SymPy gives the loose enclosure AccumBounds(-2, 2); the limit is 0
    # (x*(Si(x) - pi/2) ~ -cos(x) at oo).
    assert limit(x*(Si(x) - pi/2) + cos(x), x, oo) == 0
    # SymPy raises NotImplementedError; x*Ci(x)**2 ~ sin(x)**2/x -> 0
    assert limit(x*Ci(x)**2, x, oo) == 0
    assert limit(x*Ci(x)/log(x), x, oo) == 0


@pytest.mark.parametrize('expr, args, expected', [
    (log(x - I) - log(-x - I), (oo,), 'I*pi'),
    (log(x - I) - log(-x - I), (-oo,), '-I*pi'),
    (log(-x - I) - log(x), (oo,), '-I*pi'),
    (log(-x + I) - log(x), (oo,), 'I*pi'),
    (sqrt(-x - I)/sqrt(x), (oo,), '-I'),
    (log(-1 - I/x), (oo,), '-I*pi'),
    (log(I*x - 1), (0, '-'), '-I*pi'),
    (log(I*x - 1), (0, '+'), 'I*pi'),
    (asin(2 + I*x), (0, '-'), 'asin(2)'),
    (acos(2 + I*x), (0, '+'), '-acos(2)'),
    (atan(2*I + x), (0, '-'), '-pi + I*atanh(2)'),
])
def test_branch_cut_limits(expr, args, expected):
    assert str(limit(expr, x, *args)) == expected
