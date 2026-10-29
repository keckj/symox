"""Regressions from SymPy 1.14's own solver tests (wave 25, solvesuite).

Expected values are hard-coded from SymPy 1.14 (sympy/solvers/tests/
test_solvers.py and direct reference runs).
"""
from symox import Eq, I, log, oo, pi, solve, sqrt, symbols
from symox.abc import x


def test_eq_denominators_cancelling_in_difference():
    # test_issue_8666: the reference keeps the denominators of both sides
    assert solve(Eq(x**2 - 1/(x**2 - 4), 4 - 1/(x**2 - 4)), x) == []
    assert solve(Eq(x + 1/x, 1/x), x) == []
    # without Eq the difference cancels at construction, as in SymPy
    assert solve(x**2 - 1/(x**2 - 4) - 4 + 1/(x**2 - 4), x) == [-2, 2]
    assert solve(x + 1/x - 1/x, x) == [0]
    assert solve(Eq(x**2, 4), x) == [-2, 2]


def test_infinite_coefficients():
    # test_inf
    assert solve(1 - oo*x) == []
    assert solve(1 - oo*x, x) == []
    assert solve(oo*x, x) == []
    assert solve(oo*x - oo, x) == []


def test_imaginary_symbol_powers():
    # test_real_imag_splitting
    A, B = symbols('a b', imaginary=True)
    assert solve(sqrt(A**2 + B**2) - 3, A) == []
    assert solve(A**2 + B**2 - 9, A) == []
    assert (B**2).is_negative is True
    assert (B**2).is_real is True
    assert (B**4).is_positive is True
    assert (9 - B**2).is_positive is True
    assert sqrt(9 - B**2).is_imaginary is False
    assert (B**3).is_imaginary is True
    a2, b2 = symbols('a b', real=True)
    assert solve(sqrt(a2**2 + b2**2) - 3, a2) == [-sqrt(9 - b2**2),
                                                  sqrt(9 - b2**2)]


def test_power_of_power_tsolve():
    # test_issue_17452
    assert solve((7**x)**x + pi, x) == [
        -sqrt(log(pi) + I*pi)/sqrt(log(7)), sqrt(log(pi) + I*pi)/sqrt(log(7))]
    assert solve((7**x)**x - 2, x) == [-sqrt(log(2))/sqrt(log(7)),
                                       sqrt(log(2))/sqrt(log(7))]
