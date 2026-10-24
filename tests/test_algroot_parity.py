"""Root objects and algebraic fields: values hard-coded from SymPy 1.14."""
import mpmath

from symox import (AlgebraicNumber, CRootOf, Lambda, Poly, Rational, RootSum, N, factor,
                   integrate, latex, log, minimal_polynomial, sqrt, symbols, exp, I)

x, y, a = symbols('x y a')


def test_factor_over_sixth_root_extensions():
    # The l1-norm test of dup_zz_zassenhaus before every trial division,
    # and the product cost of each recombination candidate.
    assert str(factor(x**6 + 2, extension=(-2)**Rational(1, 6))) == (
        '(x - (-2)**(1/6))*(x + (-2)**(1/6))*(x**2 - (-2)**(1/6)*x + (-2)**(1/3))'
        '*(x**2 + (-2)**(1/6)*x + (-2)**(1/3))')
    assert str(factor(x**6 - 2, extension=2**Rational(1, 6))) == (
        '(x - 2**(1/6))*(x + 2**(1/6))*(x**2 - 2**(1/6)*x + 2**(1/3))'
        '*(x**2 + 2**(1/6)*x + 2**(1/3))')
    assert str(factor(x**4 + 1, extension=sqrt(2))) == '(x**2 - sqrt(2)*x + 1)*(x**2 + sqrt(2)*x + 1)'
    assert str(factor(minimal_polynomial(sqrt(2) + sqrt(3) + sqrt(5), x))) == (
        'x**8 - 40*x**6 + 352*x**4 - 960*x**2 + 576')


def test_definite_rational_integral_without_radical_poles():
    # No real root of x**5 - x - 1 on [2, 3]: the RootSum antiderivative
    # is evaluated at the bounds although roots() finds no radicals.
    r = integrate(1/(x**5 - x - 1), (x, 2, 3))
    assert r.has(RootSum)
    assert str(r.evalf()) == '0.0132656737262677'
    assert str(r.evalf(30)) == '0.0132656737262677493098747616090'
    assert str(N(r)) == '0.0132656737262677'
    assert str(N(r, 15, chop=True)) == '0.0132656737262677'
    assert str(integrate(1/(x**5 - x - 1), (x, 2, 3)).evalf()) == '0.0132656737262677'


def test_root_sum_evalf_guard_bits_and_real_part():
    r = integrate(1/(x**4 + x + 1), (x, 0, 1))
    assert str(r.evalf()) == '0.644752916573651'
    assert str(N(r, 15, chop=True)) == '0.644752916573651'
    e = RootSum(x**5 - x - 1, Lambda(y, exp(y)))
    assert str(e.evalf()) == '5.20845872973884'
    assert str(N(e)) == '5.20845872973884'
    assert str((2*e + 1).evalf()) == '11.4169174594777'
    assert str(RootSum(x**3 - 2, Lambda(y, log(y))).evalf(30)) == '0.693147180559945309417232121458'
    assert str(RootSum(x**3 - 2, Lambda(y, I*exp(y))).evalf()) == '4.01673290451468*I'


def test_root_sum_subs_keeps_generator_bound():
    rs = RootSum(x**5 - x - 1, Lambda(y, log(x - y)))
    value = rs.subs(x, 2)
    assert str(value) == 'RootSum(x**5 - x - 1, Lambda(y, log(2 - y)))'
    # Intentional correction: SymPy 1.14 substitutes into the generator
    # too and gives 0. The sum of log(2 - r) over the roots is log p(2).
    with mpmath.workdps(30):
        assert abs(float(value.evalf(20)) - float(log(29).evalf(20))) < 1e-14
    assert str(value.evalf()) == '3.36729582998647'
    assert rs.subs(y, 2) == rs
    parametric = RootSum(x**3 - a*x - 1, Lambda(y, log(x - y)), x)
    assert str(parametric.subs(a, 3)) == 'RootSum(x**3 - 3*x - 1, Lambda(y, log(x - y)))'
    assert parametric.subs(a, 2) == log(x + 1) + RootSum(x**2 - x - 1, Lambda(y, log(x - y)))


def test_root_sum_latex():
    assert latex(RootSum(x**3 - x - 1, Lambda(y, log(x - y)))) == (
        r'\operatorname{RootSum} {\left(x^{3} - x - 1, \left( y \mapsto \log{\left(x - y \right)} \right)\right)}')


def test_algebraic_number_to_root_and_is_aliased():
    r0, r1 = CRootOf(x**3 - 2, 0), CRootOf(x**3 - 2, 1)
    assert str(AlgebraicNumber(r0, [1, 2, 3]).to_root()) == 'CRootOf(x**3 - 9*x**2 + 15*x - 11, 0)'
    assert str(AlgebraicNumber(r0, [1, 2, 3]).to_root(radicals=False)) == 'CRootOf(x**3 - 9*x**2 + 15*x - 11, 0)'
    assert str(AlgebraicNumber(sqrt(2), [1, 1]).to_root()) == '1 + sqrt(2)'
    assert str(AlgebraicNumber(sqrt(2), [-1, 1]).to_root()) == '1 - sqrt(2)'
    assert str(AlgebraicNumber(sqrt(2)).to_root()) == 'sqrt(2)'
    assert str(AlgebraicNumber(sqrt(2), [2]).to_root()) == '2'
    assert str(AlgebraicNumber(sqrt(2) + sqrt(3), [1, 0, 0]).to_root()) == '2*sqrt(6) + 5'
    assert str(AlgebraicNumber(r1, [1, 0]).to_root()) == 'CRootOf(x**3 - 2, 1)'
    assert str(AlgebraicNumber(r1, [1, 0]).to_root(minpoly=Poly(x**3 - 4))) == 'CRootOf(x**3 - 2, 1)'
    assert str(AlgebraicNumber(sqrt(-2), [1, 1]).to_root()) == '1 + sqrt(2)*I'
    assert AlgebraicNumber(sqrt(2), alias='t').is_aliased is True
    assert AlgebraicNumber(sqrt(2)).is_aliased is False
