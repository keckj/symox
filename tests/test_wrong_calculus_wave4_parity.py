"""Reference values (SymPy 1.14) for wave-4 wrong results in special
numbers, rising/falling factorials of polynomials, Ei at polar arguments,
sums and products of Indexed terms, ConditionSet substitution and limit_seq
of nested sums."""
from symox import (ConditionSet, Ei, FiniteSet, Function, IndexedBase,
                   Interval, Poly, Product, QQ, Rational, S, Sum, Symbol,
                   Union, euler, exp_polar, ff, genocchi, I, limit_seq, oo,
                   pi, polar_lift, product, rf, summation, symbols, digamma)
from symox.abc import x, y, z, n, m, k
from symox.testing.pytest import raises

f = IndexedBase('f')
g = Function('g')


def test_genocchi_polynomials_and_values():
    assert genocchi(2, x) == 1 - 2*x
    assert genocchi(4, x) == -4*x**3 + 6*x**2 - 1
    assert genocchi(6, x) == -6*x**5 + 15*x**4 - 15*x**2 + 3
    assert genocchi(8, x) == -8*x**7 + 28*x**6 - 70*x**4 + 84*x**2 - 17
    assert genocchi(12, x).subs(x, 0) == -2073
    assert [genocchi(nn) for nn in range(2, 13)] == [
        -1, 0, 1, 0, -3, 0, 17, 0, -155, 0, 2073]
    assert [genocchi(nn, 2) for nn in range(2, 13)] == [
        -3, -6, -9, -10, -9, -14, -33, -18, 135, -22, -2097]
    assert [genocchi(nn, -3) for nn in range(2, 13)] == [
        7, -36, 161, -660, 2541, -9324, 32977, -113508, 383285, -1276572,
        4204473]
    assert genocchi(4, Rational(1, 3)) == Rational(-13, 27)
    assert genocchi(10, Rational(1, 3)) == Rational(1525355, 19683)
    assert genocchi(7, Rational(-5, 2)) == Rational(-208985, 64)
    assert genocchi(12, Rational(-5, 2)) == Rational(145952937, 256)
    assert genocchi(2, y + 1) == -2*y - 1
    assert genocchi(4, x).diff(x) == -12*x**2 + 12*x
    assert genocchi(29, 2) == -58


def test_euler_negative_one():
    assert euler(-1) == pi/2
    assert euler(0) == 1
    assert euler(0, x) == 1
    assert euler(-1, x) == digamma((x + 1)/2) - digamma(x/2)


def test_rf_ff_infinity():
    assert ff(-oo, 2) is oo
    assert rf(-oo, 2) is oo
    assert ff(-oo, 3) is -oo
    assert rf(-oo, 3) is -oo


def test_rf_ff_poly():
    assert rf(Poly(x**3, x), 2) == Poly(x**6 + 3*x**5 + 3*x**4 + x**3, x)
    assert ff(Poly(x**3, x), 2) == Poly(x**6 - 3*x**5 + 3*x**4 - x**3, x)
    assert rf(Poly(x**3, x), 1) == Poly(x**3, x)
    assert rf(Poly(2*x + 1, x), 3) == Poly(8*x**3 + 36*x**2 + 46*x + 15, x)
    assert ff(Poly(2*x + 1, x), 3) == Poly(8*x**3 - 12*x**2 - 2*x + 3, x)
    assert ff(Poly(2*x + 1, x), 2) == Poly(4*x**2 - 1, x)
    assert rf(Poly(x**2 + y, x), 2) == Poly(
        x**4 + 2*x**3 + (2*y + 1)*x**2 + 2*y*x + y**2 + y, x)
    p = rf(Poly(x, x, domain=QQ), 3)
    assert p == Poly(x**3 + 3*x**2 + 2*x, x, domain=QQ)
    assert p.domain == QQ
    assert rf(Poly(2*x + 1, x), -1) == 1/(2*x - 1)
    assert ff(Poly(2*x + 1, x), -1) == 1/(2*x + 3)
    assert rf(Poly(x**2, x), -2) == 1/(x**4 - 6*x**3 + 13*x**2 - 12*x + 4)
    assert ff(Poly(x**2, x), -2) == 1/(x**4 + 6*x**3 + 13*x**2 + 12*x + 4)
    assert ff(Poly(x**3, x), -1) == 1/(x**3 + 3*x**2 + 3*x + 1)
    raises(ValueError, lambda: rf(Poly(x*y, x, y), 2))
    raises(ValueError, lambda: ff(Poly(x*y, x, y), 2))


def _close(a, b, tol=1e-12):
    return abs(complex(a) - complex(b)) < tol


def test_Ei_polar_evalf():
    v = Ei(polar_lift(-1)).n()
    assert str(v) == '-0.21938393439552 + 3.14159265358979*I'
    assert str(Ei(exp_polar(I*pi)).n()) == \
        '-0.21938393439552 + 3.14159265358979*I'
    assert str(Ei(exp_polar(-I*pi)).n()) == \
        '-0.21938393439552 - 3.14159265358979*I'
    assert str(Ei(polar_lift(-2)).n(5)) == '-0.0489 + 3.1416*I'
    assert str(Ei(polar_lift(-1)).evalf(30)) == (
        '-0.21938393439552027367716377546 + '
        '3.14159265358979323846264338328*I')
    assert _close(Ei(polar_lift(-3)).n(), -0.013048381094197 + 3.14159265358979j)


def test_sum_product_indexed():
    assert Sum(f[n], (n, 0, 3)).doit() == f[0] + f[1] + f[2] + f[3]
    assert Sum(f[n]**2, (n, 0, 3)).doit() == \
        f[0]**2 + f[1]**2 + f[2]**2 + f[3]**2
    assert Sum(f[n]*n, (n, 1, 3)).doit() == f[1] + 2*f[2] + 3*f[3]
    assert Sum(f[n + 1] - f[n], (n, 0, 3)).doit() == f[4] - f[0]
    assert Sum(f[2*n], (n, 0, 2)).doit() == f[0] + f[2] + f[4]
    assert Sum(f[n]**2, (n, 3, 0)).doit() == -f[1]**2 - f[2]**2
    assert Sum(f[n, m], (n, 0, 1), (m, 0, 1)).doit() == \
        f[0, 0] + f[0, 1] + f[1, 0] + f[1, 1]
    assert Sum(f[n]*g(n), (n, 0, 1)).doit() == g(0)*f[0] + g(1)*f[1]
    assert summation(f[n], (n, 0, 2)) == f[0] + f[1] + f[2]
    assert summation(f[n]**2, (n, 1, 2)) == f[1]**2 + f[2]**2
    assert Product(f[n], (n, 1, 3)).doit() == f[1]*f[2]*f[3]
    assert Product(f[n] + 1, (n, 0, 1)).doit() == (f[0] + 1)*(f[1] + 1)
    assert product(f[n]**2, (n, 1, 2)) == f[1]**2*f[2]**2
    s = Sum(f[n], (n, 0, k))
    assert s.doit() == s


def test_conditionset_subs_bound_symbol():
    assert ConditionSet(x, x < 1, FiniteSet(x, z)).subs(x, 0) == \
        Union(ConditionSet(x, x < 1, FiniteSet(z)), FiniteSet(0))
    assert ConditionSet(x, x < 1, FiniteSet(x, z)).subs({x: 0, z: 5}) == \
        FiniteSet(0)
    assert ConditionSet(x, x < 1, FiniteSet(x, z)).subs(x, 2) == \
        ConditionSet(x, x < 1, FiniteSet(z))
    assert ConditionSet(x, x < 1, S.Reals).subs(x, y) == \
        ConditionSet(x, x < 1, S.Reals)
    assert ConditionSet(x, x < y, S.Reals).subs(x, y) == \
        ConditionSet(x, x < y, S.Reals)
    assert ConditionSet(x, x < y, Interval(0, z)).subs(x, z) == \
        ConditionSet(x, x < y, Interval(0, z))


def test_conditionset_subs_free_symbol_intentional_divergence():
    # SymPy 1.14 ignores a substitution of a free (non-bound) symbol in the
    # condition; Symox applies it, which is the mathematically correct set.
    assert ConditionSet(x, x < y, S.Reals).subs(y, 1) == \
        ConditionSet(x, x < 1, S.Reals)


def test_limit_seq_sum_over_exponential():
    # The unevaluated Sum has no known mrv set, so gruntz refuses it (the
    # reference raises) and limit_seq takes the Stolz-Cesaro path.
    assert limit_seq(Sum(2**m/m, (m, 1, n))/(2**n/n), n) == 2
    nested = limit_seq(
        Sum(k**2*Sum(2**m/m, (m, 1, k)), (k, 1, n))/(2**n*n), n)
    assert nested == 4
