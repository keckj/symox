"""limheld (wave 24): limits left unevaluated after limrest3 and seriesrest.

Expected values are those of SymPy 1.14 (sympy/series/tests/test_limits.py
and test_order.py), except the one-sided Piecewise cases marked as
intentional corrections, where SymPy 1.14 evaluates the conditions at the
limit point itself (or raises TypeError at infinity).
"""

from symox import *

del test

x, y, z = symbols("x y z")
a = Symbol("a", real=True)
n = Symbol("n", integer=True, positive=True)
k = Symbol("k", positive=True)


def test_parametric_sign_of_leading_coefficient():
    # sign(c0, x) of a product: one symbolic sign per x-free factor
    assert limit(y/x/log(x), x, 0) == -oo*sign(y)
    assert limit(y*x/log(x), x, oo) == oo*sign(y)
    assert limit(-y*z*exp(x)/x, x, oo) == -oo*sign(y)*sign(z)
    assert limit(a*log(x)/x**2, x, 0) == -oo*sign(a)
    assert limit(y*x**2*log(x)/(1 + x), x, oo) == oo*sign(y)


def test_symbolic_power_differences():
    # the exponent difference -1 + (n - 5)/n has the sign of -5/n
    xr = Symbol("x", real=True)
    assert limit(xr**n - xr**(n - 0), xr, oo) == 0
    assert limit(xr**n - xr**(n - 5), xr, oo) == oo
    assert limit(xr**n - xr**(n - 2.5), xr, oo) == oo
    assert limit(xr**n - xr**(n - k - 1), xr, oo) == oo
    xp = Symbol("x", positive=True)
    assert limit(xp**n - xp**(n - 1), xp, oo) == oo
    assert limit(xp**n - xp**(n + 2), xp, oo) == -oo
    assert limit(xp**(n + 1) - xp**(n + 3)/2, xp, oo) == -oo


def test_piecewise_branch_near_the_point():
    func1 = 2*sqrt(x)*Piecewise(((4*x - 2)/Abs(sqrt(4 - 4*(2*x - 1)**2)), 4*x - 2 >= 0),
                                ((2 - 4*x)/Abs(sqrt(4 - 4*(2*x - 1)**2)), True))
    func2 = Piecewise((x**2/2, x <= 0.5), (x/2 - 0.125, True))
    func3 = Piecewise(((x - 9)/5, x < -1), ((x - 9)/5, x > 4), (sqrt(Abs(x - 3)), True))
    assert limit(func1, x, 0) == 1
    assert limit(func2, x, 0) == 0
    assert limit(func3, x, -1) == 2
    assert limit(Piecewise((sin(x)/x, x < 2), (x, True)), x, 0) == 1
    assert limit(Piecewise((cos(x), x < 1), (x, True)), x, 3) == 3
    assert limit(Piecewise((cos(x), x < 1), (x, True)), x, 3, "-") == 3


def test_piecewise_one_sided_intentional_corrections():
    # SymPy 1.14 gives 2 and 1 (conditions evaluated at x = 0 itself); the
    # one-sided neighbourhoods select the other branch.
    assert limit(Piecewise((1, x < 0), (2, True)), x, 0, "-") == 1
    assert limit(Piecewise((1, x <= 0), (2, True)), x, 0, "+") == 2
    # left limit of func3 at -1 uses the branch x < -1 (SymPy: 2)
    func3 = Piecewise(((x - 9)/5, x < -1), ((x - 9)/5, x > 4), (sqrt(Abs(x - 3)), True))
    assert limit(func3, x, -1, "-") == -2
    # SymPy 1.14 raises TypeError (comparison with zoo) at infinity
    assert limit(Piecewise((1/x, x < 1), (x, True)), x, oo) is oo
    # a condition on another symbol stays undecided (SymPy 1.14 gives 2)
    assert isinstance(limit(Piecewise((x, y < 1), (2, True)), x, 0), Limit)


def test_order_terms_in_limits():
    assert (1 + Order(x)).limit(x, 0) == 1
    assert ((x + Order(x**2))/x).limit(x, 0) == 1
    assert limit(cos(x) + O(x**3), x, 0) == 1
    # SymPy 1.14: Order(0, (x, 0)), which is not O(0)
    assert limit(O(x**2), x, 0) == Order(0, (x, 0))


def test_nan_limit_point():
    assert limit(Order(2)*x, x, S.NaN) is S.NaN
    assert limit(x, x, S.NaN) is S.NaN
