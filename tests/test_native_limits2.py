"""limits2: limits and series against SymPy 1.14 -- the reference
`heuristics`, limits through a positive dummy at infinity, Order
construction, containment and substitution (orders inside the expression,
several variables, noncommutative factors, conjugates), the `Expr.series`
loop of held Subs/Derivative/Integral, the generic Taylor expansion of
undefined functions, acoth expansions and leading terms with real symbolic
exponents."""
import pytest

from symox import *  # noqa: F401,F403
from symox import _core
from symox.calculus.accumulationbounds import AccumBounds
from symox.series.limits import heuristics

del test  # noqa: F821

x, y, z, a, b, w, n = symbols('x y z a b w n')
A, B = symbols('A B', commutative=False)
f = Function('f')


def test_heuristics():
    xr = Symbol('x', real=True)
    assert heuristics(sin(1/xr) + atan(xr), xr, 0, '+') == AccumBounds(-1, 1)
    assert limit(log(2 + sqrt(atan(x))*sqrt(sin(1/x))), x, 0) == log(2)


def test_limit_of_zero_symbol_and_real_part():
    p = Symbol('p', zero=True)
    assert limit(n + p, n, 0) == 0
    assert limit(re(z), z, 0) == 0
    assert limit(2*pi*re(z)/5, z, 0) == 0
    assert re((1 - exp(-2*I*pi*z))/(1 - exp(-2*I*pi*z/5))).limit(z, 0) == 5


def test_acoth_limit_at_infinite_argument():
    assert limit(acoth(exp(1/x)), x, 0, '+') == 0


@pytest.mark.parametrize('expr, expected', [
    (acoth(x).nseries(x, n=6), -I*pi/2 + x + x**3/3 + x**5/5 + O(x**6)),
    (acoth(x + 1).nseries(x, n=3), log(2)/2 - log(x)/2 + x/4 - x**2/16 + O(x**3)),
    (acoth(x + S(1)/2).series(x, n=3), -I*pi/2 + log(3)/2 + 4*x/3 + 8*x**2/9 + O(x**3)),
])
def test_acoth_series(expr, expected):
    assert expr == expected


def test_order_api():
    assert Order(x).as_expr_variables(((x, 0),)) == (x, ((x, 0),))
    assert Order(y).as_expr_variables(((x, 0),)) == (y, ((x, 0), (y, 0)))
    assert O(x**2)._eval_as_leading_term(x, None, 1) == O(x**2)


def test_order_sums_and_powers():
    assert Order(a + b, a, b) + 1 + O(1, a, b) == O(1, a, b)
    assert O(1)**O(1) == O(1)
    assert (Order(x + y) + x).expr == x + y
    assert Order(x + y) + x**2 == Order(x + y)
    assert Order(x*(x + y), x, y).expr == x**2 + x*y


def test_order_symbols():
    e = x*y*sin(x)*Integral(x, (x, 1, 2))
    assert O(e) == O(x**2*y, x, y)
    assert O(e, x) == O(x**2)
    assert O(x, x, y).variables == (x, y)
    assert Order(exp(x)*exp(1/y)).expr == exp(x)*exp(1/y)


def test_order_noncommutative_and_conjugate():
    assert O(x)*(x*A + A) == O(x)
    assert x*A**2 + O(x**2)*A**2 == x*A**2 + O(x**2)
    xr = Symbol('x', real=True)
    xi = Symbol('x', imaginary=True)
    assert conjugate(O(xr)) == O(xr)
    assert conjugate(O(xi)) == O(xi)


def test_order_subs():
    assert (x + O(x**2)).subs(x, x + O(x**2)) == x + O(x**2)
    assert O(x**2).subs(x, x + O(x**2)) == O(x**2)
    assert (2 + O(w)).subs(w, 0) == 2
    assert sin(x).series(n=8).subs(x, sin(x).series(n=8)).expand() == \
        x - x**3/3 + x**5/10 - 8*x**7/315 + O(x**8)
    assert cos(x).series(n=8).subs(x, sin(x).series(n=8)).expand() == \
        1 - x**2/2 + 5*x**4/24 - 37*x**6/720 + O(x**8)


def test_series_of_undefined_function():
    assert f(x).series(x, 0, 3) == \
        f(0) + x*Subs(Derivative(f(x), x), x, 0) + x**2*Subs(Derivative(f(x), x, x), x, 0)/2 + O(x**3)
    assert f(x**2).series(x, 0, 3) == f(0) + x**2*Subs(Derivative(f(x), x), x, 0) + O(x**3)


def test_series_of_held_operations():
    assert Subs(sin(x), x, y).series(y) == \
        Subs(x**5/120, x, y) + Subs(-x**3/6, x, y) + Subs(x, x, y) + O(y**6)
    assert Derivative(sin(x), x).series(x, n=3).doit() == cos(x).series(x, n=3)
    assert Integral(cos(x), x).series(x) == sin(x).series(x)
    assert Integral(x, (x, 1, 3), (y, 1, x)).series(x) == 4*x - 4
    assert next(S.Zero.lseries(x)) == 0
    assert next(Derivative(cos(x), x).lseries()) == Derivative(1, x)


@pytest.mark.parametrize('expr, expected', [
    (Derivative(x**3, x), 3*x**2),
    (Derivative(x**3, y), 0),
    (Derivative(sin(x), x), 1),
    (Derivative(cos(x), x), -x),
    (Derivative(1, x), Derivative(1, x)),
])
def test_derivative_leading_term(expr, expected):
    assert expr.as_leading_term(x) == expected


def test_leadterm_real_symbolic_exponent():
    assert (3 + 2*x**(log(3)/log(2) - 1)).leadterm(x) == (3, 0)


def test_limit_free_symbols_non_symbol_variable():
    assert Limit(x + y, x + y, a).free_symbols == {a}
    assert Limit(-x**2 + y, x**2, a).free_symbols == {y, a}


def test_native_entry_points():
    xr = Symbol('x', real=True)
    assert _core._limit_heuristics(sin(1/xr) + atan(xr), xr, S(0), '+') == AccumBounds(-1, 1)
    assert _core._limit_heuristics(x**2 + 1, x, S(0), '+') == 1
