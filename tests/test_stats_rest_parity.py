"""Regression tests for the remaining sympy.stats gaps (batch stats-rest).

Expected values are those of SymPy 1.14 (simpler equal forms where SymPy
leaves erfc/Integral residue, see docs/algorithmic-parity.md)."""
import pytest

from symox import (Contains, Eq, Gt, Interval, Le, Matrix, MatrixSymbol, Not,
                   Rational, S, Tuple, exp, latex, pi, simplify, sqrt, symbols,
                   besselk, erf)
from symox.stats import (E, P, ContinuousMarkovChain, DiscreteMarkovChain, Die,
                         GammaProcess, GaussianSymplecticEnsemble, Laplace,
                         MultivariateNormal, Normal, WienerProcess, density,
                         dependent, given)


def test_random_symbols_print_as_their_names():
    N = Normal('N', 0, 1)
    assert str(N + N**2) == 'N**2 + N'
    assert latex(N + N**2) == 'N^{2} + N'
    assert E(N + N**2) == 1
    assert (N > 1).is_Relational


def test_given_tuple_dependent_finite():
    D1, D2 = Die('D1', 6), Die('D2', 6)
    s, d = given(Tuple(D1 + D2, D1), D1 > 4)
    assert E(s) == 9 and E(d) == Rational(11, 2)
    assert P(s > 9) == Rational(5, 12)
    a, b = given(Tuple(D1, D2), Eq(D1 + D2, 4))
    assert E(a) == 2 and E(a*b) == Rational(10, 3)
    assert dependent(a, b)


def test_given_tuple_dependent_continuous():
    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    A, B = given(Tuple(X, Y), Eq(X + Y, 3))
    assert E(A) == Rational(3, 2)
    assert E(A*B) == Rational(7, 4)
    assert simplify(density(A)(1)) == exp(Rational(-1, 4))/sqrt(pi)
    A, B = given(Tuple(X, Y), Eq(2*X + Y, 1))
    # SymPy: 2/5 behind erfc residue and an unevaluated Integral for E(B).
    assert E(A) == Rational(2, 5) and E(B) == Rational(1, 5)
    with pytest.raises(NotImplementedError, match='Multivariate Inequalities'):
        given(Tuple(X, Y), X + Y > 1)


def test_stationary_distribution_reference_forms():
    T = Matrix([[S(1)/2, S(1)/2, 0], [S(1)/4, S(1)/2, S(1)/4], [0, S(1)/2, S(1)/2]])
    assert str(DiscreteMarkovChain('Z', [0, 1, 2], T).stationary_distribution()) == 'Matrix([[1/4, 1/2, 1/4]])'
    assert str(DiscreteMarkovChain('Z', [0, 1], Matrix([[1, 0], [0, 1]])).stationary_distribution()) == 'Matrix([[1 - tau0, tau0]])'
    W = DiscreteMarkovChain('W', [0, 1, 2], MatrixSymbol('T', 3, 3))
    assert str(W.stationary_distribution(condition_set=True)) == 'ConditionSet(wm, Eq(wm*T, wm))'
    assert str(W.stationary_distribution()) == 'Lambda((wm, T), Eq(wm*T, wm))'
    U = DiscreteMarkovChain('U')
    assert str(U.stationary_distribution()) == 'Lambda((wm, _T), Eq(wm*_T, wm))'


def test_multivariate_laplace_dispatch():
    L = Laplace('L', [2, 4], [[3, 1], [1, 3]])
    assert simplify(density(L)(1, 2) - sqrt(2)*exp(Rational(11, 4))*besselk(0, sqrt(165)/4)/(4*pi)) == 0


def test_multivariate_normal_symbolic_dimension_is_unsupported():
    n = symbols('n', natural=True)
    with pytest.raises(NotImplementedError, match='symbolic dimension'):
        MultivariateNormal('X', MatrixSymbol('mu', n, 1), MatrixSymbol('Sg', n, n))


def test_gse_density_at_matrix_symbol():
    G = GaussianSymplecticEnsemble('G', 2)
    X = MatrixSymbol('X', 2, 2)
    assert str(density(G)(X)) == '32*exp(-2*Trace(X**2))/pi**3'


def _ctmc():
    G = Matrix([[-S(1), Rational(1, 10), Rational(9, 10)],
                [Rational(2, 5), -S(1), Rational(3, 5)],
                [Rational(1, 2), Rational(1, 2), -S(1)]])
    return ContinuousMarkovChain('C', state_space=[0, 1, 2], gen_mat=G)


def test_ctmc_float_times():
    C = _ctmc()
    assert str(P(Eq(C(3.65), 2), Eq(C(1.78), 1)).round(10)) == "0.4002723175"
    assert str(P(Gt(C(43.2), 0), Eq(C(3.29), 2)).round(10)) == "0.6832579186"


def test_markov_symbolic_times():
    a, b, c, d = symbols('a b c d')
    C = _ctmc()
    query = P(Eq(C(a), b), Eq(C(c), d))
    assert str(query.subs({a: 3.65, b: 2, c: 1.78, d: 1}).evalf().round(10)) == "0.4002723175"
    T = Matrix([[Rational(1, 10), Rational(4, 10), Rational(5, 10)],
                [Rational(3, 10), Rational(4, 10), Rational(3, 10)],
                [Rational(7, 10), Rational(2, 10), Rational(1, 10)]])
    Y = DiscreteMarkovChain('Y', [0, 1, 2], T)
    assert str(P(Eq(Y[a], b), Eq(Y[c], d)).subs({a: 10, b: 2, c: 5, d: 1}).round(4)) == "0.3096"
    assert str(P(Gt(Y[a], b), Eq(Y[c], d)).subs({a: 21, b: 0, c: 5, d: 0}).evalf().round(5)) == "0.64705"
    assert str(P(Le(Y[a], b), Eq(Y[c], d)).subs({a: 7, b: 1, c: 5, d: 0}).evalf().round(5)) == "0.78000"


def test_wiener_and_gamma_contains_conditions():
    t1, t2, t, d = symbols('t1 t2 t d', positive=True)
    W = WienerProcess('W')
    assert E(W(t1) + 2*W(t2), Contains(t1, Interval.Lopen(0, 1)) & Contains(t2, Interval.Lopen(1, 2))) == 0
    assert simplify(P(W(t1) > 1, Contains(t1, Interval(1, 3)))) == S(1)/2 - erf(S(1)/2)/2
    # SymPy adds the probabilities of the two (overlapping) events of an Or
    # and returns -erf(1)/2 + erf(2)/2 + 1; the event is certain.
    assert simplify(P((W(t1) > 2) | (W(t1) < 4), Contains(t1, Interval.Ropen(2, 4)))) == 1
    X = GammaProcess('X', 1, 2)
    value = P(Not((X(t) < 5) & (X(d) > 3)), Contains(t, Interval.Ropen(2, 4)) & Contains(d, Interval.Lopen(7, 8)))
    assert simplify(value) == -4*exp(-3) + Rational(472, 3)*exp(-8) + 1
