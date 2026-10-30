"""The public checksol (native port of SymPy 1.14's checksol): every value
below was produced by SymPy 1.14."""
import warnings

import pytest

from symox import (Abs, And, Derivative, Dict, Eq, Function, Integral, Interval, LambertW,
                   Matrix, Ne, Piecewise, Poly, Rational, S, Symbol, Tuple, cos, exp, log,
                   nan, oo, pi, sin, sqrt, symbols, zoo)
from symox.core.sympify import SympifyError
from symox.solvers.solvers import checksol, unrad

x, y, z = symbols('x y z')
f = Function('f')


def _run(*args, **flags):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result = checksol(*args, **flags)
    return result, [str(w.message) for w in caught if 'could not verify' in str(w.message)]


@pytest.mark.parametrize('args, flags, expected', [
    ((Eq(x, 1), x, 1), {}, True),
    ((Eq(x**2, 4), x, -2), {}, True),
    ((Eq(x, True), x, True), {}, True),
    ((Eq(True, x), x, True), {}, True),
    ((Eq(x, True), x, 1), {}, False),
    ((Eq(x, False), x, False), {}, True),
    ((Eq(x, True), x, y), {}, None),
    ((Eq(x, True), x, y), {'minimal': True}, None),
    ((True, x, 1), {}, True),
    ((False, x, 1), {}, False),
    ((S.true, {}), {}, True),
    ((0, {}), {}, True),
    ((x - x, x, 1), {}, True),
    (([x - 1, y - 2], {x: 1, y: 2}), {}, True),
    (([x - 1, y - 2], {x: 1, y: 3}), {}, False),
    (([x - 2, Ne(x, 2)], {x: 1}), {}, False),
    (((x - 1,), x, 1), {}, True),
    ((Poly(x**2 - 1, x), x, 1), {}, True),
    ((x - y, x, y), {}, True),
    ((x + 1, x, x), {}, False),
    ((x - x**2, x, x), {}, False),
    ((sin(x)**2 + cos(x)**2 - 1, {}), {}, True),
    ((sin(x)**2 + cos(x)**2 - 1, {}), {'minimal': True}, None),
    ((x - 1, x, oo), {}, False),
    ((x - 1, x, zoo), {}, False),
    ((1/x, x, 0), {}, False),
    ((nan, x, 1), {}, None),
    ((x - 1, x, 1), {'numerical': False}, True),
    ((x - sqrt(2), x, sqrt(2)), {'simplify': False}, True),
    ((x - 1, x, 1), {'rational': True}, True),
    ((x - 1, x, 1), {'foo': 3}, True),
    ((exp(x) - y, x, log(y)), {}, True),
    ((exp(x) - y, x, log(y)), {'numerical': False}, True),
    ((sqrt(x) + 1, x, 1), {}, False),
    ((sqrt(x) + 1, x, 1), {'force': False}, False),
    ((f(x) - 1, f(x), 1), {}, True),
    ((f(x) - x, {f(x): x}), {}, True),
    ((x - pi, x, 3.14159265358979), {}, True),
    ((x - 1, {x: 1.0}), {}, True),
    ((Matrix([x - 1, y - 2]), {x: 1, y: 2}), {}, True),
    ((Tuple(x - 1), x, 1), {}, True),
    ((Dict({x: 1}), y, 1), {}, False),
    ((x - 1, 'x', 1), {}, True),
    ((x - 1, {'x': 1}), {}, True),
    ((x - y, {x: 1}), {}, False),
    ((x**2 - y, {x: sqrt(y)}), {}, True),
    ((Abs(x) - x, x, y), {}, True),
    ((x*exp(x) - 1, x, LambertW(1)), {}, True),
    ((x*exp(x) - 1, x, LambertW(1)), {'numerical': False}, None),
    ((Eq(sin(x), 0), x, pi), {}, True),
    ((x**2 - 2, x, sqrt(2)), {'minimal': True}, True),
    ((x**y - 1, x, 1), {}, True),
    ((x - y, y, x + 1), {}, False),
    ((x + y - 1, {x: 1 - y}), {}, True),
    ((Piecewise((x, x > 0), (0, True)) - 1, x, 1), {}, True),
    ((z**x - y, {x: log(y)/log(z)}), {'minimal': True}, True),
    ((x - 1, {x: Integral(1, (y, 0, 1))}), {'numerical': False}, True),
    ((Derivative(f(x), x) - 1, {f(x): x}), {}, True),
    ((x - 2, {x: 2*sin(y)**2 + 2*cos(y)**2}), {'simplify': False}, True),
    ((x - 2, {x: 2*sin(y)**2 + 2*cos(y)**2}), {'minimal': True}, None),
])
def test_checksol_decisions(args, flags, expected):
    assert _run(*args, **flags)[0] is expected


@pytest.mark.parametrize('args, error, message', [
    ((Ne(x, 1), x, 2), AttributeError, "'BooleanTrue' object has no attribute 'is_constant'"),
    ((Ne(x, 1), x, 1), AttributeError, "'BooleanFalse' object has no attribute 'is_constant'"),
    ((x < y, x, 2), AttributeError, "'StrictLessThan' object has no attribute 'is_constant'"),
    ((And(x > 0, x < 2), x, 1), AttributeError, "'BooleanTrue' object has no attribute 'is_constant'"),
    ((Ne(x, y), {x: 1}), AttributeError, "'Unequality' object has no attribute 'is_constant'"),
    (([x - 1, Ne(x, 2)], {x: 1}), AttributeError, "'BooleanTrue' object has no attribute 'is_constant'"),
    (({x - 1: 1}, {x: 1}), AttributeError, "'Dict' object has no attribute 'is_constant'"),
    ((Interval(x, 2), x, 0), AttributeError, "'Interval' object has no attribute 'is_constant'"),
    (([], x, 1), ValueError, 'no functions to check'),
    ((x - 1, 1), ValueError, 'Expecting (sym, val) or ({sym: val}, None) but got (1, None)'),
    ((x - 1, x), ValueError, 'Expecting (sym, val) or ({sym: val}, None) but got (x, None)'),
])
def test_checksol_errors(args, error, message):
    with pytest.raises(error) as info:
        checksol(*args)
    assert str(info.value) == message


def test_checksol_string_is_not_sympified():
    with pytest.raises(SympifyError):
        checksol('x-1', x, 1)


def test_checksol_simplifies_the_mapping_in_place():
    sol = {x: 2*sin(y)**2 + 2*cos(y)**2}
    assert checksol(x - 2, sol) is True
    assert sol == {x: 2}
    sol = {x: 2*sin(y)**2 + 2*cos(y)**2}
    assert checksol(x - 2, sol, simplify=False) is True
    assert sol == {x: 2*sin(y)**2 + 2*cos(y)**2}


def test_checksol_warnings():
    w = LambertW(y)
    assert _run(x*exp(x) - y, {x: w*(sin(y)**2 + cos(y)**2)}, warn=True) == \
        (None, ['\n\tWarning: could not verify solution {x: LambertW(y)}.'])
    assert _run(x*exp(x) - 1, x, LambertW(1), numerical=False, warn=True) == \
        (None, ['\n\tWarning: could not verify solution {x: LambertW(1)}.'])
    # the mapping is shared by the equations of a list
    sol = {x: 2*sin(y)**2 + 2*cos(y)**2}
    assert _run([x - 2, x*exp(x) - y], sol, warn=True) == \
        (None, ['\n\tWarning: could not verify solution {x: 2}.'])
    assert sol == {x: 2}
    assert _run([x*exp(x) - y, x - 2], {x: w*(sin(y)**2 + cos(y)**2)}, warn=True) == \
        (None, ['\n\tWarning: could not verify solution {x: LambertW(y)}.'] * 2)
    # no warning without the flag
    assert _run(x*exp(x) - y, x, w) == (None, [])


def test_pow_recognizes_base_as_e():
    # Pow.__new__: b**(c*num/log(b)) is E**(c*num)
    assert z**(log(y)/log(z)) == y
    assert 2**(x/log(2)) == exp(x)
    assert z**((log(y) + 1)/log(z)) == S.Exp1*y
    assert z**(log(y)/log(z) + 1) == z**(log(y)/log(z) + 1)


def test_unrad_cov_updated_before_raise():
    for eq in (x**Rational(2, 3) + (x + 1)**Rational(1, 3) + 1,
               x**Rational(2, 3) + y + (x + 1)**Rational(1, 3)):
        cov = []
        with pytest.raises(NotImplementedError, match='no successful change of variable found'):
            unrad(eq, x, cov=cov)
        assert len(cov) == 2 and cov[1] == cov[0]**3 - x
    cov = []
    with pytest.raises(NotImplementedError, match='Cannot remove all radicals'):
        unrad(sqrt(x) + (x + 1)**Rational(1, 3) + (x + 2)**Rational(1, 5) + 1, x, cov=cov)
    assert cov == []
    with pytest.raises(TypeError, match="'tuple' object does not support item assignment"):
        unrad(x**Rational(2, 3) + (x + 1)**Rational(1, 3) + 1, x, cov=())
