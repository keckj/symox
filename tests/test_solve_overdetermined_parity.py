"""Overdetermined solve systems: SymPy 1.14 outputs hardcoded.

SymPy solves a system polynomial in the unknowns over the coefficient domain
of its parameters, so a residual equation that is a nonzero parametric
expression makes the system generically inconsistent. Equations with free
symbols but no unknown are dropped before solving.
"""
import pytest

import symox as sx

x, y, z, a, b = sx.symbols('x y z a b')
NS = dict(x=x, y=y, z=z, a=a, b=b, sin=sx.sin, cos=sx.cos, exp=sx.exp,
          sqrt=sx.sqrt, pi=sx.pi)


def _run(source, flags):
    args = eval('(' + source + ')', dict(NS))
    return sx.solve(*args, **flags)


def _expected(source):
    return eval(source, dict(NS))


CASES = [
    # floats, the reported case
    ("[1.5*x*y - 2, x + y - 1], [x]", {}, "[]"),
    ("[1.5*x*y - 2, x + y - 1], [x]", {'dict': True}, "[]"),
    ("[1.5*x*y - 2, x + y - 1], [x]", {'rational': False}, "[]"),
    ("[0.5*x - 1, x - 3], [x]", {}, "[]"),
    ("[1.5*x - 3, 2.0*x - 5], [x]", {'dict': True}, "[]"),
    # exact, generically inconsistent in the parameters
    ("[x*y - 2, x + y - 1], [x]", {}, "[]"),
    ("[x*y - 2, x + y - 1], [x]", {'dict': True}, "[]"),
    ("[x*y - 2, x + y - 1], [x]", {'set': True}, "([x], set())"),
    ("[x - y, x - 1], [x]", {}, "[]"),
    ("[x - y, x - 1], [x]", {'set': True}, "([x], set())"),
    ("[x + y, x - y], [x]", {'dict': True}, "[]"),
    ("[x - a, x - b], [x]", {}, "[]"),
    ("[a*x - b, x - 1], [x]", {'dict': True}, "[]"),
    ("[x**2 - 4, x - a], [x]", {}, "[]"),
    ("[x**2 - y, x - 1], [x]", {'dict': True}, "[]"),
    ("[x**2 + y**2 - 1, x - y], [x]", {}, "[]"),
    ("[x - y, x**2 - 1], [x]", {}, "[]"),
    ("[x**2 - 1, x - sin(y)], [x]", {'dict': True}, "[]"),
    ("[x - 1, x - sin(y)], [x]", {}, "[]"),
    ("[x - 1, x - sqrt(y)], [x]", {'dict': True}, "[]"),
    ("[x - 1, x*(sin(y)**2 + cos(y)**2) - 1], [x]", {}, "[]"),
    ("[x + y - 1, x - y - 1, x - a], [x, y]", {}, "[]"),
    ("[x*y - 2, x + y - 1, y - 1], [x]", {}, "[]"),
    ("[x/y - 1, x - 2], [x]", {}, "[]"),
    # non-polynomial: a residual without function or radical is rejected
    ("[exp(x) + x - y, x - 1], [x]", {}, "[]"),
    ("[x*exp(x) - y, x - 1], [x]", {'dict': True}, "[]"),
    # inconsistent constants
    ("[x - 1, x - 2], x", {}, "[]"),
    ("[x - 1, x - 2], [x]", {'set': True}, "([x], set())"),
    ("[x + y - 1, x - y - 1, x - 3], [x, y]", {}, "[]"),
    ("[x - 1, x - 2, y], [x, y]", {}, "[]"),
    # consistent
    ("[x - 1, 2*x - 2], [x]", {}, "{x: 1}"),
    ("[x**2 - 1, x - 1], [x]", {}, "[(1,)]"),
    ("[x**2 - 1, x + 1], [x]", {'dict': True}, "[{x: -1}]"),
    ("[x**2 - 1, x**3 - 1], [x]", {}, "[(1,)]"),
    ("[sin(x), x - pi], [x]", {}, "[(pi,)]"),
    ("[x + y - 1, x - y - 1, 2*x - 2], [x, y]", {}, "{x: 1, y: 0}"),
    ("[x + y - 2, x - y, x*y - 1], [x, y]", {}, "[(1, 1)]"),
    ("[x - 1, y*x - y], [x]", {}, "{x: 1}"),
    ("[x*y, x], [x]", {'dict': True}, "[{x: 0}]"),
    ("[x**2 - a**2, x - a], [x]", {}, "[(a,)]"),
    ("[2*x - 2*a, x - a], [x]", {'set': True}, "([x], {(a,)})"),
    ("[x - 1, a*(x - 1)], [x]", {}, "{x: 1}"),
    ("[x/y - 1, x - y], [x]", {}, "{x: y}"),
    # equations without the unknowns are dropped
    ("[x*y - 2, y - 1], [x]", {}, "{x: 2/y}"),
    ("[x - 1, y - 1], [x]", {'dict': True}, "[{x: 1}]"),
    ("[x + z - 1, x - y], [x]", {}, "[]"),
    # floats are restored only when the kept equations have floats
    ("[x - 1, 1.5*y - 2], [x]", {'dict': True}, "[{x: 1}]"),
    ("[0.5*x - 1, x - 2], [x]", {}, "{x: 2.0}"),
    ("[0.5*x - 1, x - 2], [x]", {'rational': False}, "{x: 2.0}"),
]


@pytest.mark.parametrize('source,flags,expected', CASES)
def test_overdetermined_matches_reference(source, flags, expected):
    assert _run(source, flags) == _expected(expected)


def test_float_restoration_ignores_dropped_equations():
    result = sx.solve([x - 1, 1.5*y - 2], [x], dict=True)
    assert result == [{x: 1}]
    assert result[0][x].is_Integer


def test_proved_nonzero_parameter_still_rules_out_the_system():
    # Kept from test_solve_parameter_constraints: SymPy 1.14 drops p and
    # returns {x: 1}, although no x satisfies p = 0 when p is nonzero.
    p = sx.Symbol('p', nonzero=True)
    assert sx.solve([x - 1, p], x) == []
