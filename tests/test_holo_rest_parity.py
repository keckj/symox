"""Regression tests for the holonomic, multiple-derivative,
solve_triangulated(extension=True) and parametric nonlinsolve items
(holo-rest). Expected values are those of SymPy 1.14."""
import pytest

import symox as sp
from symox import (CRootOf, Derivative, Function, I, S, QQ, besselj, cos, csc,
                   exp, hyperexpand, log, nonlinsolve, sec, sin, solve, sqrt,
                   symbols)
from symox.holonomic import (DifferentialOperators, HolonomicFunction,
                             expr_to_holonomic)
from symox.solvers.polysys import solve_triangulated

x, y, z, a, b = symbols('x y z a b')
R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')


def _sqrt_x():
    return HolonomicFunction(-S(1)/2 + x*Dx, x, 0, {S(1)/2: [1]})


def test_frobenius_initial_conditions():
    h = _sqrt_x()
    assert str(h) == 'HolonomicFunction((-1/2) + (x)*Dx, x, 0, {1/2: [1]})'
    assert h.is_singularics() is True
    assert HolonomicFunction(Dx - 1, x, 0, [1]).is_singularics() is False
    assert str(h.to_sequence()) == '[(HolonomicSequence((n), n), u(0) = 1, 1/2, 1)]'
    assert h.to_expr() == sqrt(x)
    assert str(h.series()) == 'sqrt(x) + O(x**6)'
    assert h.to_hyper() == sqrt(x)
    assert str(HolonomicFunction(x*Dx - 2, x, 0, {2: [1]}).series()) == 'x**2 + O(x**6)'
    assert str(HolonomicFunction(4*x**2*Dx**2 + 1, x, 0, {S(1)/2: [1]}).series()) == 'sqrt(x) + O(x**6)'
    assert str(HolonomicFunction(x**2*Dx**2 + x*Dx - 1, x, 0, {1: [1], -1: [1]}).series()) == \
        '[1/x + C_2*x + C_1 + O(x**6), x + O(x**6)]'
    assert str(HolonomicFunction(2*x*Dx**2 + Dx + x, x, 0).series()) == \
        '[C_0 - C_0*x**2/6 + C_0*x**4/168 + O(x**6), C_0*sqrt(x) - C_0*x**(5/2)/10 + C_0*x**(9/2)/360 + O(x**6)]'
    assert str(HolonomicFunction(x*Dx**2 + Dx - x, x, 0, {0: [1]}).to_sequence()) == \
        '[(HolonomicSequence((-1) + (n**2 + 4*n + 4)Sn**2, n), u(0) = 1, u(1) = 0, 0, 0)]'
    # the initial conditions add; the annihilator is a scalar multiple of
    # the reference one
    assert str((h + h).y0) == '{1/2: [2]}'


def test_to_expr_besselj():
    bessel1 = HolonomicFunction(x**2*Dx**2 + x*Dx + (x**2 - 1), x, 0, [0, S(1)/2])
    assert bessel1.to_expr() == besselj(1, x)
    assert str(bessel1.to_hyper()) == 'x*hyper((), (2,), -x**2/4)/2'
    assert str(bessel1.series()) == 'x/2 - x**3/16 + x**5/384 + O(x**6)'
    assert str(bessel1.to_sequence()) == \
        '[(HolonomicSequence((1) + (n**2 + 6*n + 8)Sn**2, n), u(0) = 1/2, u(1) = 0, 1, 0)]'
    assert HolonomicFunction(x**2*Dx**2 + x*Dx + (x**2 - 4), x, 0, [0, 0, S(1)/8]).to_expr() == besselj(2, x)/2
    assert HolonomicFunction(x**2*Dx**2 + x*Dx + x**2, x, 0, [1, 0]).to_expr() == besselj(0, x)
    assert str(HolonomicFunction((x + 1)*Dx**3 + Dx**2, x, 0, [1, 2, 3]).to_expr()) == \
        '3*x*log(x + 1) - x + 3*log(x + 1) + 1'
    assert HolonomicFunction(x*Dx**2 + Dx, x, 1, [0, 1]).to_expr() == log(x)


def test_to_meijerg():
    h = expr_to_holonomic(cos(x) + sin(x)).to_meijerg()
    assert str(h) == ('sqrt(pi)*x*meijerg(((), ()), ((0,), (-1/2,)), x**2/4)/2'
                      ' + sqrt(pi)*meijerg(((), ()), ((0,), (1/2,)), x**2/4)')
    assert hyperexpand(h) == sin(x) + cos(x)
    assert str(expr_to_holonomic(exp(x)).to_meijerg()) == 'meijerg(((), ()), ((0,), ()), -x)'
    assert str(HolonomicFunction(Dx**2 + 1, x, 0, [0, 1]).to_meijerg()) == \
        'sqrt(pi)*x*meijerg(((), ()), ((0,), (-1/2,)), x**2/4)/2'


def test_to_hyper_reference_forms():
    assert str(HolonomicFunction((x + 1)*Dx**3 + Dx**2, x, 0, [1, 2, 3]).to_hyper()) == \
        '-x**3*hyper((1, 2), (4,), -x)/2 + 3*x**2/2 + 2*x + 1'
    assert str(HolonomicFunction(Dx**2 - 1, x, 0, [1, 1]).to_hyper(as_list=True)) == \
        '[(1, hyper((), (1/2,), x**2/4)), (x, hyper((), (3/2,), x**2/4))]'
    assert str(expr_to_holonomic(log(1 + x)).to_hyper()) == '-x**2*hyper((1, 2), (3,), -x)/2 + x'


def test_multiple_derivative_final_step():
    f = Function('f')
    d5 = sp.diff(sqrt((x + 1)**2 + x), (x, 5))
    assert d5.count_ops() == 30
    assert str(d5) == ('15*(2*x + 3)*(48 - 40*(2*x + 3)**2/(x + (x + 1)**2) + 7*(2*x + 3)**4/'
                       '(x + (x + 1)**2)**2)/(32*(x + (x + 1)**2)**(5/2))')
    assert str(sp.diff(exp(x)*sin(x), x, 3)) == '2*(-sin(x) + cos(x))*exp(x)'
    assert str(sp.diff(x**3*exp(2*x), x, 2)) == '2*x*(2*x**2 + 6*x + 3)*exp(2*x)'
    assert str((x**4*exp(-x)).diff(x, 2)) == 'x**2*(x**2 - 8*x + 12)*exp(-x)'
    assert str(sp.diff(exp(-x**2), x, 3)) == '4*x*(3 - 2*x**2)*exp(-x**2)'
    assert str(sp.diff(log(x)/x, x, 2)) == '(2*log(x) - 3)/x**3'
    assert str(sp.diff(f(x)*x**2, x, 2)) == \
        'x**2*Derivative(f(x), (x, 2)) + 4*x*Derivative(f(x), x) + 2*f(x)'
    assert str(Derivative(x**2*exp(-x), (x, 2)).doit()) == '(x**2 - 4*x + 2)*exp(-x)'
    assert str(Derivative(x**2*exp(-x), x, 2, evaluate=True)) == '(x**2 - 4*x + 2)*exp(-x)'
    assert str((Derivative(x**2*exp(-x), (x, 2)) + 1).doit()) == '(x**2 - 4*x + 2)*exp(-x) + 1'
    # a single derivative is not refactored
    assert str(sp.diff(x**2*exp(-x), x)) == '-x**2*exp(-x) + 2*x*exp(-x)'


def test_expand_mul_distributes_held_terms():
    held = sp.Mul(-1, a*sin(x) + b*cos(x), evaluate=False)
    assert sp.expand_mul(held + x) == x - a*sin(x) - b*cos(x)


def test_solve_triangulated_extension():
    r = solve_triangulated([x**2 + y + z - 1, x + y**2 + z - 1, x + y + z**2 - 1], x, y, z,
                           extension=True)
    c0, c1 = CRootOf(x**2 + 2*x - 1, 0), CRootOf(x**2 + 2*x - 1, 1)
    assert r == [(0, 0, 1), (0, 1, 0), (1, 0, 0), (c0, c0, c0), (c1, c1, c1)]
    r = solve_triangulated([x**2 - 2, y**2 - x], x, y, extension=True)
    assert r == [(CRootOf(x**2 - 2, 0), CRootOf(x**4 - 2, 2)), (CRootOf(x**2 - 2, 0), CRootOf(x**4 - 2, 3)),
                 (CRootOf(x**2 - 2, 1), CRootOf(x**4 - 2, 0)), (CRootOf(x**2 - 2, 1), CRootOf(x**4 - 2, 1))]
    r = solve_triangulated([x*y - 1, x**2 - 3], x, y, extension=True)
    assert r == [(CRootOf(x**2 - 3, 0), CRootOf(3*x**2 - 1, 0)), (CRootOf(x**2 - 3, 1), CRootOf(3*x**2 - 1, 1))]
    q = x**6 - 4*x**4 - x**3 + 4*x**2 - 1
    r = solve_triangulated([x**3 - x - 1, y**2 - x*y - 1], x, y, extension=True)
    assert r == [(CRootOf(x**3 - x - 1, i), CRootOf(q, j))
                 for i, j in [(0, 0), (0, 1), (1, 2), (1, 4), (2, 3), (2, 5)]]
    # the ground-domain solutions without extension
    assert solve_triangulated([x**2 - 2, y - x], x, y) == []


def test_crootof_coefficient_named_like_a_generator():
    root = CRootOf(x**2 - 2, 0)
    p = sp.Poly(y**2 - root, y, extension=True)
    assert p.all_roots(multiple=False, radicals=False) == \
        [(CRootOf(x**4 - 2, 2), 1), (CRootOf(x**4 - 2, 3), 1)]


def test_sec_csc_multiple_argument():
    assert str(solve(sec(2*x) - 3, x)) == (
        '[-I*log(-sqrt(1/3 - 2*sqrt(2)*I/3)), -I*log(1/3 - 2*sqrt(2)*I/3)/2, '
        '-I*log(-sqrt(1/3 + 2*sqrt(2)*I/3)), -I*log(1/3 + 2*sqrt(2)*I/3)/2]')
    assert str(solve(2*sec(2*x) - 5, x)) == (
        '[-I*log(-sqrt(2/5 - sqrt(21)*I/5)), -I*log(2/5 - sqrt(21)*I/5)/2, '
        '-I*log(-sqrt(2/5 + sqrt(21)*I/5)), -I*log(2/5 + sqrt(21)*I/5)/2]')
    assert str(solve(csc(3*x) + 4, x)) == (
        '[-I*log(-(1 - sqrt(3)*I)*(-sqrt(15)/4 - I/4)**(1/3)/2), '
        '-I*log(-(1 - sqrt(3)*I)*(sqrt(15)/4 - I/4)**(1/3)/2), '
        '-I*log(-(1 + sqrt(3)*I)*(-sqrt(15)/4 - I/4)**(1/3)/2), '
        '-I*log(-(1 + sqrt(3)*I)*(sqrt(15)/4 - I/4)**(1/3)/2), '
        '-I*log(-sqrt(15)/4 - I/4)/3, -I*log(sqrt(15)/4 - I/4)/3]')
    for equation in (csc(3*x) + 4, sec(2*x) - 3, sec(2*x) - a):
        for root in solve(equation, x):
            residual = equation.subs(a, S(7)/3).subs(x, root.subs(a, S(7)/3))
            assert abs(complex(residual.evalf())) < 1e-12
    with pytest.raises(NotImplementedError):
        solve(sec(x/2) - 3, x)


def test_nonlinsolve_generic_parameters():
    assert nonlinsolve([a*x**2 - 1, y - x], [x, y]) == \
        sp.FiniteSet((-sqrt(1/a), -sqrt(1/a)), (sqrt(1/a), sqrt(1/a)))
    s = nonlinsolve([a*x + y - 1, x*y - b], [x, y])
    expected = [(-(-sqrt(-4*a*b + 1)/2 - S(1)/2)/a, S(1)/2 - sqrt(-4*a*b + 1)/2),
                (-(sqrt(-4*a*b + 1)/2 - S(1)/2)/a, sqrt(-4*a*b + 1)/2 + S(1)/2)]
    assert len(s) == 2
    for got in s:
        assert any(all(sp.simplify(g - e) == 0 for g, e in zip(got, want)) for want in expected)


def test_nonlinsolve_drops_denominator_zeros():
    # SymPy also returns (0, 0) and (0, 1), where a denominator vanishes.
    assert nonlinsolve([1/x + 1/y - 1, x - 2*y], [x, y]) == sp.FiniteSet((3, S(3)/2))
    assert nonlinsolve([x/(y - 1) - 1, y**2 - 1], [x, y]) == sp.FiniteSet((-2, -1))
