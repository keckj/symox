"""The radical branches of SymPy 1.14's ``_solve`` are native.

``symox._core.solve`` is the Rust solver itself: these equations used to
raise NotImplementedError there, after which the Python wrapper reached the
``_tsolve`` step (which returned ``[]`` for the ODE reproducer) or the
Python unrad fallback. The results below are SymPy 1.14's.
"""
import pytest

from symox import (Eq, I, Rational, Symbol, cbrt, root, solve, sqrt, symbols)
from symox import _core
import symox._solve_relations as relations

x, y, C1 = symbols('x y C1')
p = Symbol('p', positive=True)


@pytest.fixture
def no_python_fallback(monkeypatch):
    """Fail if the Python unrad/inversion fallbacks are reached."""
    def forbidden(*args, **kwargs):
        raise AssertionError('Python fallback reached')
    monkeypatch.setattr(relations, '_unrad_solution', forbidden)
    monkeypatch.setattr(relations, '_inverted_solution', forbidden)
    monkeypatch.setattr(relations, '_masked_solution', forbidden)
    import symox.solvers._tsolve_steps as steps
    monkeypatch.setattr(steps, 'tsolve_step', forbidden)


ODE = -2*y/sqrt(y**3) - C1 - x
ODE_ROOT = 4/(C1**2 + 2*C1*x + x**2)


def test_ode_reproducer_is_native(no_python_fallback):
    # several generators (y, sqrt(y**3)): the unrad fallback of `_solve`
    assert _core.solve(ODE, y) == [ODE_ROOT]
    # the reference leaves the unrad solution unsimplified
    assert solve(ODE, y) == [ODE_ROOT]
    assert solve(ODE, y, dict=True) == [{y: ODE_ROOT}]
    assert solve(Eq(-2*y/sqrt(y**3), C1 + x), y) == [ODE_ROOT]


def test_ode_reproducer_check_false(no_python_fallback):
    # check=False keeps the root y = 0 of the denominator
    assert solve(ODE, y, check=False) == [0, ODE_ROOT]


def test_ode_reproducer_simplify_flag():
    # SymPy 1.14's simplify keeps the expanded denominator
    assert solve(ODE, y, simplify=True) == [ODE_ROOT]


def test_single_generator_roots_are_native(no_python_fallback):
    # roots of t**3 - x for t = sqrt(y), then y = t**2 for every root
    w = -Rational(1, 2) + sqrt(3)*I/2
    found = _core.solve(y**Rational(3, 2) - x, y)
    assert len(found) == 3
    assert x**Rational(2, 3) in found
    assert solve(y**Rational(3, 2) - x, y) == [
        x**Rational(2, 3), (-cbrt(x)/2 - sqrt(3)*I*cbrt(x)/2)**2,
        (-cbrt(x)/2 + sqrt(3)*I*cbrt(x)/2)**2]
    assert w != 0
    assert solve(y**Rational(2, 3) - x, y) == [-x**Rational(3, 2), x**Rational(3, 2)]
    assert solve(sqrt(y**3) - y, y) == [0, 1]


def test_radicals_with_one_generator(no_python_fallback):
    assert solve(sqrt(y) + y - x, y) == [
        x - sqrt(4*x + 1)/2 + Rational(1, 2), x + sqrt(4*x + 1)/2 + Rational(1, 2)]
    r = (sqrt(69)/18 + Rational(1, 2))**Rational(1, 3)
    assert solve(x**Rational(1, 3) - x + 1, x) == [1/(3*r) + r + 1]


def test_radicals_with_several_generators(no_python_fallback):
    assert _core.solve(sqrt(x) + sqrt(x + 1) - 2, x) == [Rational(9, 16)]
    assert solve(sqrt(x) + sqrt(x + 1) - 2, x) == [Rational(9, 16)]
    assert solve(sqrt(2*x + 1)/sqrt(x) - sqrt(2 + 1/x) + 6*I, x) == [-Rational(1, 11)]
    assert solve(Eq(sqrt(5*x + 6) - 2, x), x) == [-1, 2]
    assert solve(root(x, 3) + sqrt(x) - 2, x) == [1]


def test_unrad_flag_is_accepted():
    eq = (-512*y**3 + 1344*(x + 2)**Rational(1, 3)*y**2
          - 1176*(x + 2)**Rational(2, 3)*y - 169*x + 686)
    assert len(solve(eq, y, _unrad=False)) == 3


def test_undecided_radical_is_not_empty():
    # SymPy raises: three radicals of the unknown and no closed form;
    # the native branch never turns an undecided case into []
    found = None
    try:
        found = solve(sqrt(x) + root(x, 3) + root(x + 1, 5) - 2, x)
    except NotImplementedError:
        return
    assert found != []


def test_positive_unknown_check():
    # the candidate 4/(C1 + x)**2 does not satisfy the equation for every
    # positive value: the reference check removes it
    assert solve(-2*p/sqrt(p**3) - C1 - x, p) == []


def test_nested_solve_simplifies_like_reference(no_python_fallback):
    # the nested `_vsolve` of the radical-free equation simplifies its
    # linear solution (the reference `_solve` default), the unrad call
    # itself does not
    a = Symbol('a')
    assert solve(sqrt(y + 1)/sqrt(y - 1) - a, y) == [(a**2 + 1)/(a**2 - 1)]


def test_solveset_keeps_radical_step():
    # `_solveset` reaches `_solve_radical`, not `_solve`'s unrad branch:
    # the undecided candidate stays in a ConditionSet
    from symox import Complement, ConditionSet, FiniteSet, S, Union, solveset
    xr, yr = symbols('xr yr', real=True)
    r = solveset(-yr + xr/sqrt(-xr**2 + 1), xr, S.Reals)
    assert isinstance(r, Union)
    assert any(isinstance(a, Complement) and isinstance(a.args[0], ConditionSet) for a in r.args)
    assert solve(-yr + xr/sqrt(-xr**2 + 1), xr) == [-yr/sqrt(yr**2 + 1), yr/sqrt(yr**2 + 1)]
    assert FiniteSet(0) != r
