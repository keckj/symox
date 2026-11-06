"""Reference regressions for the category-c (missing domain or feature)
callables of batch c1: integrals, Sum and Product, series, calculus,
holonomic functions and stats. Expected values were computed with SymPy
1.14."""
import pytest

from symox import (
    Abs, Integral, Matrix, Ne, O, Piecewise, Rational, S, Sum, binomial, cos,
    exp, factorial, limit_seq, log, oo, periodicity, pi, simplify, sin, sqrt,
    summation, symbols, sympify, QQ, ZZ,
)
from symox.calculus.finite_diff import apply_finite_diff
from symox.holonomic import DifferentialOperators, HolonomicFunction
from symox.stats import (
    DiscreteMarkovChain, DiscreteRV, Die, E, MultivariateNormal, P, density,
)

x, y, t, u, a, b, n, k = symbols('x y t u a b n k')


# Integral.as_sum, transform and principal_value

def test_as_sum_methods():
    assert Integral(sin(x), (x, 3, 7)).as_sum(2, 'midpoint') == 2*sin(4) + 2*sin(6)
    assert Integral(x**2, (x, 0, 2)).as_sum(4, 'trapezoid') == Rational(11, 4)
    assert Integral(x**3, (x, 0, 2)).as_sum(2, 'midpoint') == Rational(7, 2)
    # Reference form: sqrt(3)*(1 + sqrt(2) + sqrt(3))/9.
    right = Integral(sqrt(x), (x, 0, 1)).as_sum(3, 'right')
    assert simplify(right - sqrt(3)*(1 + sqrt(2) + sqrt(3))/9) == 0


def test_integral_transform():
    assert Integral(x, (x, 2, 5)).transform(x, 2*u + 1) == \
        Integral(4*u + 2, (u, Rational(1, 2), 2))
    assert Integral(x*exp(x**2), (x, 0, 1)).transform(x**2, u) == \
        Integral(exp(u)/2, (u, 0, 1))


def test_principal_value():
    assert Integral(x + 1, (x, -oo, oo)).principal_value() == oo
    assert Integral(1/x, (x, -1, 1)).principal_value() == 0
    assert Integral(1/(x**2 - 1), (x, -oo, oo)).principal_value() == 0


# Order.contains and multivariate normalization

def test_order_contains():
    assert (O(1) in O(1, x)) is True
    assert O(x).contains(x**2) is True
    assert O(x**2).contains(x) is False
    assert O(x, (x, oo)).contains(1) is True
    assert O(x**2, x).contains(O(x**3, x)) is True


def test_order_multivariate_leading_terms():
    assert str(O(1 + x*y)) == 'O(1, x, y)'
    assert str(O(x + y + x*y)) == 'O(x + y, x, y)'
    assert str(O(1 + x*y, (x, oo), (y, oo))) == 'O(x*y, (x, oo), (y, oo))'
    assert str(O(x**2*y + x*y)) == 'O(x*y, x, y)'


# Hypergeometric summation stage

def test_sum_hyper_closed_forms():
    assert Sum(x**k/factorial(k), (k, 0, oo)).doit() == exp(x)
    assert summation((-1)**k*x**(2*k)/factorial(2*k), (k, 0, oo)) == cos(x)
    assert summation(x**k/k, (k, 1, oo)) == Piecewise(
        (-log(1 - x), Ne(x, 1) & (Abs(x) <= 1)),
        (Sum(x**k/k, (k, 1, oo)), True))
    assert summation(binomial(2*k, k)*x**k, (k, 0, oo)) == Piecewise(
        (1/sqrt(1 - 4*x), Ne(x, Rational(1, 4)) & (Abs(x) <= Rational(1, 4))),
        (Sum(x**k*binomial(2*k, k), (k, 0, oo)), True))


# limit_seq with a symbolic exponent

def test_limit_seq_symbolic_exponent_undecided():
    # The limit depends on the sign of b: the reference returns None.
    assert limit_seq((n**2 + a)**(n + b)/(n + a)**(2*n), n) is None
    assert limit_seq((2*n**2 + 1)/(n**2 + 3), n) == 2


def test_limit_symbolic_power_needs_sign_of_exponent():
    from symox import limit
    # These used to be 0, which is wrong for b > 1. The reference gives
    # oo**(b - 1) for x**b/x and raises for exp(b*x)*exp(-x) (powsimp merges
    # it into exp(x*(b - 1)), whose rewrite depends on the sign of b - 1).
    assert limit(x**b/x, x, oo) == oo**(b - 1)
    with pytest.raises(NotImplementedError, match="Result depends on the sign"):
        limit(exp(b*x)*exp(-x), x, oo)
    assert limit(x**b*exp(-x), x, oo) == 0
    assert limit(x**b*exp(x), x, oo) == oo
    assert limit(log(x)**b/x, x, oo) == 0
    assert limit(exp(x)/x**a, x, oo) == oo
    assert limit(exp(-x)*exp(b*sqrt(x)), x, oo) == 0


# Calculus

def test_periodicity_power_with_varying_exponent():
    assert periodicity(sin(4*x)**cos(2*x), x) == pi
    assert periodicity(cos(x)**sin(2*x), x) == 2*pi


def test_apply_finite_diff_indexed_points():
    from symox import IndexedBase, Idx
    assert apply_finite_diff(1, [-1, 0, 1], [y, 2*y, 3*y], 0) == y
    xb, yb = IndexedBase('x'), IndexedBase('y')
    i = Idx('i')
    result = apply_finite_diff(1, [xb[i - 1], xb[i], xb[i + 1]],
                               [yb[i - 1], yb[i], yb[i + 1]], xb[i])
    expected = (((xb[i + 1] - xb[i])/(-xb[i - 1] + xb[i]) - 1)*yb[i]/(xb[i + 1] - xb[i])
                - (xb[i + 1] - xb[i])*yb[i - 1]/((xb[i + 1] - xb[i - 1])*(-xb[i - 1] + xb[i]))
                + (-xb[i - 1] + xb[i])*yb[i + 1]/((xb[i + 1] - xb[i - 1])*(xb[i + 1] - xb[i])))
    assert simplify(result - expected) == 0


# Holonomic functions

def test_holonomic_printing_and_sequences():
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    assert str(HolonomicFunction(Dx**3 - Dx**2 + Dx - 1, x, 0, [1, 2, 1])) == \
        'HolonomicFunction((-1) + (1)*Dx + (-1)*Dx**2 + (1)*Dx**3, x, 0, [1, 2, 1])'
    assert str(HolonomicFunction(Dx - 1, x, 0, [1]).to_sequence()) == \
        '[(HolonomicSequence((-1) + (n + 1)Sn, n), u(0) = 1, 0)]'
    assert str(HolonomicFunction((1 + x)*Dx**2 + Dx, x, 0, [0, 1]).to_sequence()) == \
        '[(HolonomicSequence((n**2) + (n**2 + n)Sn, n), u(0) = 0, u(1) = 1, u(2) = -1/2, 2)]'
    assert str(DifferentialOperators(ZZ.old_poly_ring(x), 'Dx')[0]) == \
        'Univariate Differential Operator Algebra in intermediate Dx over the base ring ZZ[x]'


# Stats

def test_discrete_rv_with_python_set_support():
    X = DiscreteRV(x, Rational(1, 10), set=set(range(10)))
    assert E(X) == Rational(9, 2)
    assert P(X > 3) == Rational(3, 5)
    with pytest.raises(NotImplementedError):
        DiscreteRV(x, Rational(1, 3), set={1, 2, 5})


def test_multivariate_normal_density_form():
    X = MultivariateNormal('X', [3, 4], [[2, 1], [1, 2]])
    assert density(X)(y, t) == sympify(
        'sqrt(3)*exp(-t**2/3 + t*y/3 + 5*t/3 - y**2/3 + 2*y/3 - 13/3)/(6*pi)',
        locals={'t': t, 'y': y})


def test_finite_density_dict():
    assert density(2*Die('D', 4)).dict == {2: Rational(1, 4), 4: Rational(1, 4),
                                          6: Rational(1, 4), 8: Rational(1, 4)}


def test_discrete_markov_chain_structure():
    Y = DiscreteMarkovChain('Y', [0, 1, 2], Matrix([[0, 1, 0], [1, 0, 0], [1, 0, 0]]))
    assert Y.communication_classes() == [([0, 1], True, 2), ([2], False, 1)]
    T = Matrix([[S(1)/2, S(1)/2, 0], [0, 1, 0], [S(1)/3, S(1)/3, S(1)/3]])
    X = DiscreteMarkovChain('X', [0, 1, 2], T)
    states, C = X.canonical_form()
    assert states == [1, 0, 2]
    assert C == Matrix([[1, 0, 0], [S(1)/2, S(1)/2, 0], [S(1)/3, S(1)/3, S(1)/3]])
    states, A, B, Q = X.decompose()
    assert states == [1, 0, 2]
    assert A == Matrix([[1]])
    assert B == Matrix([[S(1)/2], [S(1)/3]])
    assert Q == Matrix([[S(1)/2, 0], [S(1)/3, S(1)/3]])
