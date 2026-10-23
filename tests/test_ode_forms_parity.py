"""ODE corpus cases where Symox raised while SymPy 1.14 answers, and the
polynomial helpers they go through. Expected values are SymPy 1.14's."""
from symox import (CRootOf, Eq, Function, I, Symbol, cancel, collect,
                   cos, cse, dsolve, exp, im, re, simplify, sin, sqrt,
                   symbols,
                   terms_gcd)
from symox.solvers.ode import checkodesol, classify_ode

x = Symbol('x')
_x = Symbol('_x')
f = Function('f')
C1, C2, C3, C4, C5, C6 = symbols('C1:7')


def test_collect_zero_keeps_rootof():
    # nonhomogeneous._get_simplified_sol calls collect(sol, 0) for a real
    # root; ComplexRootOf._eval_subs returns self.
    r = CRootOf(x**5 + 11*x - 2, 0)
    assert collect(x + r*x, 0) == x*r + x
    assert r.subs(x, 2) == r


LCHOM_33 = (
    "C5*exp(x*CRootOf(_x**5 + 11*_x - 2, 0)) + (C1*sin(x*im(CRootOf(_x**5 "
    "+ 11*_x - 2, 1))) + C2*cos(x*im(CRootOf(_x**5 + 11*_x - 2, 1))))"
    "*exp(x*re(CRootOf(_x**5 + 11*_x - 2, 1))) + (C3*sin(x*im(CRootOf(_x**5 "
    "+ 11*_x - 2, 3))) + C4*cos(x*im(CRootOf(_x**5 + 11*_x - 2, 3))))"
    "*exp(x*re(CRootOf(_x**5 + 11*_x - 2, 3)))")


def test_lin_const_coeff_hom_rootof():
    # SymPy 1.14 proves im(r1) < 0 numerically, so sin(x*Abs(im(r1)))
    # becomes -sin(x*im(r1)) and the sign goes into C1; the CRootOf sign
    # is decided natively too (it raised in collect(sol, 0) before).
    eq = f(x).diff(x, 5) + 11*f(x).diff(x) - 2*f(x)
    sol = dsolve(eq, f(x), hint='nth_linear_constant_coeff_homogeneous')
    assert str(sol.rhs) == LCHOM_33


def test_euler_hom_rootof():
    eq = x**6*f(x).diff(x, 6) - x*f(x).diff(x) + f(x)
    sol = dsolve(eq, f(x), hint='nth_linear_euler_eq_homogeneous')
    p = "_x**5 - 14*_x**4 + 71*_x**3 - 154*_x**2 + 120*_x - 1"
    assert str(sol) == "Eq(f(x), C1*x + " + " + ".join(
        "C%d*x**CRootOf(%s, %d)" % (i + 2, p, i) for i in range(5)) + ")"


def test_terms_gcd_gaussian_composite_generator():
    e = C3 + x*(C5 + I*x/8) + 5*I*exp(2*I*x)/16
    assert str(terms_gcd(e, clear=False, expand=False)) == \
        "C3 + x*(C5 + I*x/8) + 5*(I*exp(2*I*x))/16"


def test_variation_of_parameters_gaussian():
    eq = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - exp(I*x)
    sol = dsolve(eq, f(x),
                 hint='nth_linear_constant_coeff_variation_of_parameters')
    assert sol == Eq(f(x), C1
                     + (C2 + x*(C3 - x/8) + 5*exp(2*I*x)/16)*sin(x)
                     + (C4 + x*(C5 + I*x/8) + 5*I*exp(2*I*x)/16)*cos(x)
                     - I*exp(I*x))


def test_cancel_extension_exponential_generator():
    e = (exp(2*x) - 2)/(exp(x) - sqrt(2))
    assert cancel(e, extension=True) == exp(x) + sqrt(2)


def test_rational_riccati_algebraic_roots():
    eq = Eq(f(x)**2 + f(x).diff(x) - 2, 0)
    sol = dsolve(eq, f(x), hint='1st_rational_riccati')
    assert sol == Eq(f(x), sqrt(2)*(-C1 - exp(2*sqrt(2)*x))
                     / (C1 - exp(2*sqrt(2)*x)))
    assert checkodesol(eq, sol) == (True, 0)


def test_rational_riccati_mixed_fields():
    eq = (f(x).diff(x) + (3*x**2 + 1)*f(x)**2/x
          + (6*x**2 - x + 3)*f(x)/(x*(x - 1))
          + (3*x**2 - 2*x + 2)/(x*(x - 1)**2))
    sol = dsolve(eq, f(x), hint='1st_rational_riccati')
    assert sol == Eq(f(x), (-C1 - x**3 + x**2 - 2*x + 1)
                     / (C1*x - C1 + x**4 - x**3 + x**2 - 2*x + 1))
    assert checkodesol(eq, sol) == (True, 0)


def test_cse_skips_rootof():
    # tree_cse's _find_repeated neither counts nor enters a RootOf.
    p = _x**5 + 11*_x - 2
    e = (C1*exp(x*CRootOf(p, 0))
         + C2*exp(x*re(CRootOf(p, 1)))*sin(x*im(CRootOf(p, 1))))
    assert cse(e) == ([], [e])


def test_simplify_count_ops_tie_keeps_sum():
    # SymPy's shorter steps undo factor_terms here (count_ops tie), and the
    # factorable classification of dsolve depends on that form.
    e = -2*x*f(x) + x*exp(-x**2)
    assert simplify(e) == e
    y = Symbol('y')
    assert simplify(x - x*y) == x*(1 - y)
    assert simplify(2*x*f(x) + 2*x) == 2*x*(f(x) + 1)


def test_factorable_classification_follows_solve_form():
    eq = 2*x*f(x) - x*exp(-x**2) + f(x).diff(x)
    assert classify_ode(eq, f(x)) == (
        '1st_exact', '1st_linear', 'Bernoulli', 'almost_linear',
        '1st_power_series', 'lie_group', '1st_exact_Integral',
        '1st_linear_Integral', 'Bernoulli_Integral',
        'almost_linear_Integral')
