"""Limits, series and evaluation of discontinuous functions at their jumps.

The expected strings are the outputs of SymPy 1.14 on the same inputs,
except in the ``CORRECTIONS`` table, where SymPy 1.14 returns a wrong value
and the mathematically correct one is kept (see docs/algorithmic-parity.md).
"""
import pytest

import symox as sx
from symox import (Abs, DiracDelta, Heaviside, Limit, Piecewise, S, Symbol, ceiling, cos, exp,
                   floor, frac, integrate, limit, log, oo, series, sign, sin, symbols)

x, y, a = symbols('x y a')
p = Symbol('p', positive=True)
xr = Symbol('xr', real=True)
H, D = Heaviside, DiracDelta

DIFFER = 'ValueError'

# (expression, point, {dir: expected}); DIFFER means SymPy raises ValueError
LIMITS = [
    (lambda: H(x), 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x, 0), 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x, 1), 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x, S.Half), 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x, y), 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x - 1), 1, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(-x), 0, {'+': '0', '-': '1', '+-': DIFFER}),
    (lambda: H(x**2), 0, {'+': '1', '-': '1', '+-': '1'}),
    (lambda: H(x)*x, 0, {'+': '0', '-': '0', '+-': '0'}),
    (lambda: H(x) + H(-x), 0, {'+': '1', '-': '1', '+-': '1'}),
    (lambda: H(x)/x, 0, {'+': 'oo', '-': '0', '+-': DIFFER}),
    (lambda: H(x)*sin(x)/x, 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x)*exp(x), 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x)**2, 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: H(x)*H(-x), 0, {'+-': '0'}),
    (lambda: H(x) - H(x, 0), 0, {'+-': '0'}),
    (lambda: H(p*x), 0, {'+': '1', '-': '0'}),
    (lambda: H(x), 1, {'+': '1'}),
    (lambda: H(x), oo, {'-': '1'}),
    (lambda: H(x), -oo, {'+': '0'}),
    (lambda: H(1/x), oo, {'-': '1'}),
    (lambda: H(-1/x), oo, {'-': '0'}),
    (lambda: H(x)/x, oo, {'-': '0'}),
    (lambda: H(x, 0)*cos(x), 0, {'-': '0'}),
    (lambda: x*H(x)/sin(x), 0, {'+': '1'}),
    (lambda: H(x)*log(x), 0, {'+': '-oo'}),
    (lambda: H(x)*log(-x), 0, {'-': '0'}),
    (lambda: H(-x)*exp(1/x), 0, {'-': '0'}),
    (lambda: H(x)*exp(-1/x), 0, {'+': '0'}),
    (lambda: H(x)*x**x, 0, {'+': '1'}),
    (lambda: H(x - y), 0, {'+': 'Heaviside(-y)'}),
    (lambda: sign(x), 0, {'+': '1', '-': '-1', '+-': DIFFER}),
    (lambda: sign(xr), 0, {'+': '1', '-': '-1', '+-': DIFFER}),
    (lambda: sign(x)*x, 0, {'+': '0', '-': '0', '+-': '0'}),
    (lambda: sign(x - 1), 1, {'-': '-1'}),
    (lambda: sign(sin(x)), 0, {'+': '1'}),
    (lambda: Abs(x)/x, 0, {'+': '1', '-': '-1', '+-': DIFFER}),
    (lambda: Abs(x - 1)/(x - 1), 1, {'-': '-1'}),
    (lambda: floor(x), 0, {'+': '0', '-': '-1', '+-': DIFFER}),
    (lambda: ceiling(x), 0, {'+': '1', '-': '0', '+-': DIFFER}),
    (lambda: frac(x), 0, {'+': '0', '-': '1', '+-': DIFFER}),
    (lambda: floor(x - 1), 1, {'+': '0', '-': '-1', '+-': DIFFER}),
    (lambda: floor(x), S.Half, {'+': '0'}),
    (lambda: floor(x), 1, {'-': '0'}),
    (lambda: ceiling(x), 1, {'-': '1'}),
    (lambda: frac(x), 1, {'-': '1'}),
    (lambda: D(x), 0, {'+': '0', '-': '0', '+-': '0'}),
    (lambda: D(x)*x, 0, {'+': '0', '-': '0', '+-': '0'}),
    (lambda: D(x - 1), 1, {'+': '0'}),
    (lambda: D(x, 1), 0, {'+': '0'}),
    (lambda: D(x), 1, {'+': '0'}),
]


def _limit_str(f, point, d):
    try:
        return str(limit(f(), x if 'xr' not in str(f()) else xr, point, d))
    except ValueError:
        return DIFFER


@pytest.mark.parametrize('f,point,d,expected', [
    (f, point, d, expected) for f, point, table in LIMITS for d, expected in table.items()])
def test_one_sided_limits_at_jumps(f, point, d, expected):
    assert _limit_str(f, point, d) == expected


def test_heaviside_limit_messages_and_held_forms():
    with pytest.raises(ValueError, match='left hand limit = 0 and right hand limit = 1'):
        limit(H(x), x, 0, '+-')
    assert str(Limit(H(x), x, 0).doit()) == '1'
    assert str(H(x).limit(x, 0, '-')) == '0'


# SymPy 1.14 is wrong on these (or never finishes); the correct value is kept.
CORRECTIONS = [
    # SymPy: 1/2 on every side (gruntz evaluates Heaviside(sin(0)))
    (lambda: limit(H(sin(x)), x, 0, '+'), '1'),
    (lambda: limit(H(sin(x)), x, 0, '-'), '0'),
    (lambda: limit(H(x**3 - x), x, 0, '+'), '0'),
    (lambda: limit(H(x**3 - x), x, 0, '-'), '1'),
    (lambda: limit(H(exp(x) - 1), x, 0, '+'), '1'),
    (lambda: limit(H(exp(x) - 1), x, 0, '-'), '0'),
    # SymPy: 1/2, but the sign of a is unknown
    (lambda: limit(H(a*x), x, 0), "Limit(Heaviside(a*x), x, 0, dir='+')"),
    # SymPy never finishes
    (lambda: limit(H(x), x, -1), '0'),
    (lambda: limit(H(x - 1), x, 0), '0'),
    # SymPy: 1 (the right branch) for the left limit of a Piecewise
    (lambda: limit(Piecewise((0, x < 0), (1, True)), x, 0, '-'), '0'),
    (lambda: limit(Piecewise((x - 1, x < 0), (x + 1, x >= 0)), x, 0, '-'), '-1'),
    # SymPy: 1/2 + x*DiracDelta(0) + ... (a Taylor expansion through the jump)
    (lambda: series(H(sin(x)), x, 0, 3), '1'),
    (lambda: H(x).nseries(x, 0, 3), '1'),
    (lambda: H(x, y).nseries(x, 0, 3), '1'),
    (lambda: DiracDelta(x).nseries(x, 0, 3), '0'),
    (lambda: (H(x)*exp(x)).nseries(x, 0, 3), '1 + x + x**2/2 + O(x**3)'),
    # SymPy raises NotImplementedError (no _eval_as_leading_term routine)
    (lambda: H(x).as_leading_term(x), '1'),
    (lambda: (H(x)*x).as_leading_term(x), 'x'),
    (lambda: H(-x).as_leading_term(x), '0'),
    (lambda: sign(x).as_leading_term(x), '1'),
]


@pytest.mark.parametrize('f,expected', CORRECTIONS)
def test_documented_corrections(f, expected):
    assert str(f()) == expected


@pytest.mark.parametrize('f,expected', [
    (lambda: limit(Piecewise((0, x < 0), (1, True)), x, 0, '+-'), DIFFER),
    (lambda: limit(H(sin(x)), x, 0, '+-'), DIFFER),
])
def test_two_sided_corrections_raise(f, expected):
    with pytest.raises(ValueError):
        f()


SERIES = [
    (lambda: series(H(x), x, 0, 3), '1'),
    (lambda: series(H(x, 0), x, 0, 3), '1'),
    (lambda: series(H(x, y), x, 0, 3), '1'),
    (lambda: series(H(x - 1), x, 1, 3), '1'),
    (lambda: series(H(-x), x, 0, 3), '0'),
    (lambda: series(H(x**2), x, 0, 3), '1'),
    (lambda: series(H(x)*x, x, 0, 3), 'x'),
    (lambda: series(H(x) + H(-x), x, 0, 3), '1'),
    (lambda: series(H(x)*sin(x)/x, x, 0, 3), '1 - x**2/6 + O(x**3)'),
    (lambda: series(H(x)*exp(x), x, 0, 3), '1 + x + x**2/2 + O(x**3)'),
    (lambda: series(H(x)**2, x, 0, 3), '1'),
    (lambda: series(H(x), x, 1, 3), '1'),
    (lambda: series(H(x), x, 0, 3, dir='-'), '0'),
    (lambda: series(H(x, 0), x, 0, 3, dir='-'), '0'),
    (lambda: series(H(x)*cos(x), x, 0, 4), '1 - x**2/2 + O(x**4)'),
    (lambda: series(H(x - 1)*x, x, 1, 3), 'x'),
    (lambda: series(H(x), x, oo, 3), '1'),
    (lambda: series(H(x + 1), x, 0, 3), '1'),
    (lambda: series(D(x), x, 0, 3), '0'),
    (lambda: series(D(x + 1), x, 0, 3), '0'),
    (lambda: series(D(x)*x, x, 0, 3), '0'),
    (lambda: series(sign(x), x, 0, 3, dir='-'), '-1'),
    (lambda: series(floor(x), x, 0, 3), '0'),
    (lambda: series(ceiling(x), x, 0, 3), '1'),
    (lambda: series(frac(x), x, 0, 3), 'x'),
    (lambda: series(floor(x - 1), x, 1, 3), '0'),
    (lambda: floor(x).as_leading_term(x), '0'),
    (lambda: ceiling(x).as_leading_term(x), '1'),
    (lambda: frac(x).as_leading_term(x), 'x'),
    (lambda: H(p).as_leading_term(p), '1'),
    (lambda: H(p).series(p, 0, 3), '1'),
]


@pytest.mark.parametrize('f,expected', SERIES)
def test_series_at_jumps(f, expected):
    assert str(f()) == expected


EVALF = [
    (lambda: H(x).evalf(subs={x: -1}), '0'),
    (lambda: H(x).evalf(subs={x: 0}), '0.500000000000000'),
    (lambda: H(x).evalf(subs={x: 2}), '1.00000000000000'),
    (lambda: H(x).evalf(subs={x: S(1)/3}), '1.00000000000000'),
    (lambda: H(x, 0).evalf(subs={x: -1}), '0'),
    (lambda: H(x, 0).evalf(subs={x: 0}), '0'),
    (lambda: H(x, 0).evalf(subs={x: 2}), '1.00000000000000'),
    (lambda: H(x, y).evalf(subs={x: 0, y: 0.25}), '0.250000000000000'),
    (lambda: H(x, y).evalf(subs={x: -1, y: 0.25}), '0'),
    (lambda: H(x, y).evalf(subs={x: 2, y: 0.25}), '1.00000000000000'),
    (lambda: H(x).subs(x, 0).evalf(), '0.500000000000000'),
    (lambda: H(x, 0).subs(x, 0).evalf(), '0'),
    (lambda: (H(x) + x).evalf(subs={x: 0}), '0.500000000000000'),
    (lambda: H(x - S.Half).evalf(subs={x: S.Half}), '0.500000000000000'),
    (lambda: (2*H(x)).evalf(subs={x: 0}), '1.00000000000000'),
    (lambda: H(x, 1).evalf(5, subs={x: 0}), '1.0000'),
    (lambda: D(x).evalf(subs={x: 1}), '0'),
    (lambda: H(x).evalf(subs={x: 0.1}), '1.00000000000000'),
    (lambda: H(x).evalf(subs={x: -0.1}), '0'),
    (lambda: sign(x).evalf(subs={x: 0}), '0'),
    (lambda: sx.N(H(0)), '0.500000000000000'),
    (lambda: sx.N(H(0, 0)), '0'),
    (lambda: H(sx.sqrt(2) - x).evalf(subs={x: 1}), '1.00000000000000'),
]


@pytest.mark.parametrize('f,expected', EVALF)
def test_evalf_subs(f, expected):
    assert str(f()) == expected


@pytest.mark.parametrize('f,expected', [
    (lambda: sx.diff(H(x, y), y), 'Derivative(Heaviside(x, y), y)'),
    (lambda: sx.diff(H(x, y), x), 'DiracDelta(x)'),
    (lambda: sx.diff(H(x, y**2), y), '2*y*Subs(Derivative(Heaviside(x, _xi_2), _xi_2), _xi_2, y**2)'),
    (lambda: sx.diff(H(x*y, y), y),
     'x*DiracDelta(x*y) + Subs(Derivative(Heaviside(x*y, _xi_2), _xi_2), _xi_2, y)'),
])
def test_derivative_through_h0(f, expected):
    assert str(f()) == expected


@pytest.mark.parametrize('f,expected', [
    (lambda: integrate(H(x)*exp(-x), (x, -oo, oo)), '1'),
    (lambda: integrate(H(x, 0), (x, -1, 1)), '1'),
    (lambda: integrate(H(x)*x, (x, -1, 1)), '1/2'),
    (lambda: integrate(D(x)*cos(x), (x, -1, 1)), '1'),
    (lambda: integrate(H(x)*exp(-p*x), (x, 0, oo)), '1/p'),
])
def test_definite_integrals_with_jumps(f, expected):
    assert str(f()) == expected


@pytest.mark.parametrize('f,expected', [
    # the delta derivative lies inside the interval or at infinity, where
    # the limit of the antiderivative is now 0 (the reference returns
    # DiracDelta(oo) or leaves the support undecided)
    (lambda: integrate(D(x - p, 1), (x, 0, oo)), '0'),
    (lambda: integrate(D(x - 1, 1), (x, 0, oo)), '0'),
    (lambda: integrate(D(x, 1), (x, -oo, oo)), '0'),
    (lambda: integrate(D(x - a), (x, 0, oo)), '1 - Heaviside(-a)'),
])
def test_delta_derivatives_away_from_the_bounds(f, expected):
    assert str(f()) == expected


@pytest.mark.parametrize('f,lo,hi', [
    (D(x, 1), 0, oo), (D(x, 1), -oo, 0), (D(x, 1), 0, 1),
    (x*D(x, 2), 0, oo), (D(x, 1)*H(x), 0, oo), (D(x - a, 1), 0, oo),
])
def test_delta_derivative_on_a_bound_stays_unevaluated(f, lo, hi):
    # the reference returns DiracDelta(0), which is not a value
    assert isinstance(integrate(f, (x, lo, hi)), sx.Integral)
