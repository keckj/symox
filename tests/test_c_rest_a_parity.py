"""Parity regressions for doc batches c1 and c2 (values checked against SymPy 1.14)."""

import pytest

from symox import (
    Derivative, Eq, Function, Lambda, O, Piecewise, Product, Rational, S, Symbol,
    cos, cosh, cot, exp, fps, gamma, limit, log, oo, simplify, sin, sinh, symbols,
    tan, trigsimp, pi, sqrt, I,
)
from symox.simplify.hyperexpand import (
    G_Function, MeijerShiftA, MeijerShiftB, MeijerUnShiftB, devise_plan_meijer,
)
from symox.stats import BernoulliProcess, Binomial, Die, E, P, density

x, y, z = symbols('x y z')
n = Symbol('n', integer=True)


# --- products of numeric-base powers -------------------------------------------

def test_numeric_bases_left_alone_when_unsafe():
    # (-2)**x*3**x != (-6)**x for non-integer x.
    assert str((-2)**x*3**x) == '(-2)**x*3**x'
    assert str(2**x*3**y) == '2**x*3**y'
    assert str(I**x*2**x) == '2**x*I**x'


def test_wallis_product():
    k = Symbol('k', integer=True)
    m = Symbol('m', integer=True)
    W = Product(2*k/(2*k - 1)*2*k/(2*k + 1), (k, 1, m))
    assert limit(simplify(W.rewrite(gamma)), m, oo) == pi/2


# --- devise_plan_meijer -----------------------------------------------------

def test_devise_plan_meijer_docstring():
    plan = devise_plan_meijer(G_Function([0], [], [0], []), G_Function([1], [], [1], []), z)
    assert [type(o) for o in plan] == [MeijerUnShiftB, MeijerShiftA]
    plan = devise_plan_meijer(G_Function([0], [], [0], []), G_Function([-1], [], [1], []), z)
    assert [type(o) for o in plan] == [MeijerShiftA, MeijerShiftB]
    assert devise_plan_meijer(G_Function([1], [2], [3], [4]), G_Function([1], [2], [3], [4]), z) == []


def test_devise_plan_meijer_unreachable():
    with pytest.raises(NotImplementedError):
        devise_plan_meijer(G_Function([0], [], [0], []), G_Function([S.Half], [], [1], []), z)


# --- series of x**y ----------------------------------------------------------

def test_symbolic_power_series():
    logx = Symbol('logx')
    assert (x**y).nseries(x, 0, 2) == x**y
    assert (x**y).nseries(x, 0, 2, logx=logx) == exp(logx*y)
    assert (x**y).series(x, 0, 2, logx=logx) == exp(logx*y)
    assert (x**y).series(x, 0, 2) == x**y
    assert (x**y).series(x, 0) == x**y
    assert (x**sin(y)).series(x, 0, 2) == x**sin(y)
    assert (x**y).series(x, 1, 2) == 1 + y*(x - 1) + O((x - 1)**2, (x, 1))
    assert str((x**y).series(x, oo, 2)) == '(1/x)**(-y)'


def test_limit_of_split_symbolic_powers():
    pp = Symbol('pp', positive=True)
    q = Symbol('q')
    assert limit(x**pp*x**(-pp - 1), x, oo) == 0
    assert limit(x**q*x**(-q - 1), x, oo) == 0
    assert limit(x**pp*x**(-pp - 1), x, 0, '+') == oo
    assert limit(x**q*x**(-q)*sin(x)/x, x, 0) == 1


# --- trigsimp(old=True) ---------------------------------------------------------

def test_trigsimp_old():
    assert trigsimp(2*sin(x)**2 + 2*cos(x)**2, old=True) == 2
    assert trigsimp(log(2*sin(x)**2 + 2*cos(x)**2), deep=True, old=True) == log(2)
    e = (-sin(x) + 1)/cos(x) + cos(x)/(-sin(x) + 1)
    assert trigsimp(e, old=True) == e
    assert trigsimp(e, method="groebner", old=True) == 2/cos(x)
    assert trigsimp(1/cot(x)**2, old=True) == cot(x)**(-2)
    assert trigsimp(sin(x)/cos(x), old=True) == tan(x)
    assert trigsimp(sin(x)*cos(y) + cos(x)*sin(y), old=True) == sin(x + y)
    assert trigsimp(sinh(x)/cosh(x), old=True).func.__name__ == 'tanh'
    assert trigsimp(tan(x)**2 + 1, old=True) == cos(x)**(-2)
    assert trigsimp(cosh(x)**2 - sinh(x)**2, old=True) == 1
    assert trigsimp((1 + cos(x))*(cos(x) - 1), old=True) == -sin(x)**2
    assert trigsimp(cos(x)**2 + sin(x)**2 + y*sin(z)**2 + y*cos(z)**2, old=True) == y + 1


# --- formal power series with a symbolic power prefactor -----------------------

def test_fps_symbolic_prefactor():
    m = Symbol('n')
    assert fps(x**m*sin(x**2), x).truncate(8) == \
        -x**(m + 6)/6 + x**(m + 2) + O(x**(m + 8), (x, 0))
    assert fps(x**m*exp(x), x).truncate(4) == \
        x**(m + 3)/6 + x**(m + 2)/2 + x**(m + 1) + x**m + O(x**(m + 4), (x, 0))
    assert fps(x**m*cos(x), x).truncate(5) == \
        x**(m + 4)/24 - x**(m + 2)/2 + x**m + O(x**(m + 5), (x, 0))
    assert fps(x**m, x) == x**m
    s = fps(x**m*exp(x), x)
    assert (s[0], s[1], s[2]) == (x**m, x**(m + 1), x**(m + 2)/2)


# --- stats: symbolic Die/Binomial and BernoulliProcess --------------------------

def test_symbolic_finite_density_dict():
    k = Symbol('k', positive=True, integer=True)
    X = Binomial('X', k, S.Half)
    assert density(X).dict.subs(k, 4).doit() == {
        0: Rational(1, 16), 1: Rational(1, 4), 2: Rational(3, 8),
        3: Rational(1, 4), 4: Rational(1, 16)}
    D = Die('D', k)
    assert density(D).dict.subs(k, 4).doit() == {i: Rational(1, 4) for i in range(1, 5)}
    assert str(density(D).dict) == 'Density(DieDistribution(k))'


def test_bernoulli_process_joint_distribution():
    B = BernoulliProcess("B", p=0.7, success=1, failure=0)
    assert str(B[1]) == 'B[1]'
    jd = B.joint_distribution(B[1], B[2])
    assert str(jd) == ('JointDistributionHandmade(Lambda((B[1], B[2]), '
                       'Piecewise((0.7, Eq(B[1], 1)), (0.3, Eq(B[1], 0)), (0, True))*'
                       'Piecewise((0.7, Eq(B[2], 1)), (0.3, Eq(B[2], 0)), (0, True))))')
    assert B.joint_distribution(1, 2) == jd
    assert abs(jd(1, 0) - 0.21) < 1e-12
    assert abs(jd.pdf(1, 1) - 0.49) < 1e-12
    assert abs(P(Eq(B[1] + B[2] + B[3], 2)) - 0.441) < 1e-12
    B2 = BernoulliProcess('B2', p=S(1)/3, success=2, failure=5)
    assert B2.joint_distribution(B2[1], B2[3])(2, 5) == Rational(2, 9)


def test_symbolic_finite_cdf():
    from symox.stats import cdf
    k = Symbol('k', positive=True, integer=True)
    c = cdf(Die('D', k))
    assert str(c) == 'Lambda(_k, Sum(Piecewise((1/k, (_ki >= 1) & (_ki <= k)), (0, True)), (_ki, 1, _k)))'
    assert c(3).subs(k, 6).doit() == Rational(1, 2)
    assert cdf(Binomial('X', k, S.Half))(2).subs(k, 4).doit() == Rational(11, 16)
    assert cdf(Binomial('Y', k, Rational(1, 3)))(1).subs(k, 3).doit() == Rational(20, 27)
    # The reference cdf of a shifted symbolic Binomial ignores success and
    # failure (wrong support); it stays refused.
    with pytest.raises(NotImplementedError):
        cdf(Binomial('Z', k, S.Half, 2, 1))
