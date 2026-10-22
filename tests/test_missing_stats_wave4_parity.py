"""Wave-4 statistics parity: joint events over independent spaces and the
missing symbolic-probability features. Expected values come from SymPy 1.14,
except where SymPy is wrong (noted inline and in docs/algorithmic-parity.md)."""
import pytest
import symox as S
from symox.stats import (P, E, Poisson, Normal, Uniform, Die, Geometric, Exponential, variance)

Q = S.Rational
e = S.E


def same(a, b):
    d = S.sympify(a) - S.sympify(b)
    assert S.simplify(d) == 0, (a, b)


@pytest.fixture(scope='module')
def rv():
    return dict(A=Poisson('A', 1), B=Poisson('B', 2), N=Normal('N', 0, 1), U=Uniform('U', 0, 1),
                D=Die('D', 6), G=Geometric('G', Q(1, 3)), X=Exponential('X', 1))


def test_independent_discrete_and(rv):
    A, B = rv['A'], rv['B']
    same(P(S.Eq(A, 2) & S.Eq(B, 3)), 2*S.exp(-3)/3)
    same(P(S.Eq(rv['G'], 2) & S.Eq(A, 1)), 2*S.exp(-1)/9)
    same(P(S.Eq(A, 2) & S.Eq(B, 3) & S.Eq(rv['G'], 1)), 2*S.exp(-3)/9)
    same(P(S.Eq(rv['D'], 3) & S.Eq(A, 1)), S.exp(-1)/6)
    same(P((rv['G'] > 2) & (A < 2)), 8*S.exp(-1)/9)
    same(P((A > 1) & (B < 3)), 5*(e - 2)*S.exp(-3))


def test_independent_or_not(rv):
    A, B = rv['A'], rv['B']
    # SymPy adds the probabilities of non-disjoint events (4*exp(-2)/3 + exp(-1)/2);
    # inclusion-exclusion subtracts P(A = 2) P(B = 3).
    same(P(S.Eq(A, 2) | S.Eq(B, 3)), S.exp(-1)/2 + 4*S.exp(-2)/3 - 2*S.exp(-3)/3)
    # SymPy: 1/6 + exp(-1), which double counts the joint outcome
    same(P(S.Eq(rv['D'], 2) | S.Eq(A, 0)), Q(1, 6) + 5*S.exp(-1)/6)
    same(P((A < 2) | (B < 1)), (2*S.exp(2) + e - 2)*S.exp(-3))
    same(P(S.Not((A > 1) & (B > 1))), (-6 + 3*e + 2*S.exp(2))*S.exp(-3))


def test_discrete_continuous_mix(rv):
    A, N = rv['A'], rv['N']
    same(P((A > 1) & (N > 0)), Q(1, 2) - S.exp(-1))
    # SymPy: (2 - erfc(sqrt(2)/2))*exp(-1)/2
    same(P(S.Eq(A, 1) & (N < 1)), (1 + S.erf(S.sqrt(2)/2))*S.exp(-1)/2)
    same(P((A < 2) | (N > 0)), Q(1, 2) + S.exp(-1))
    same(P((rv['D'] > 3) & (N > 0)), Q(1, 4))
    same(P((N > 0) & (rv['U'] < Q(1, 2))), Q(1, 4))
    same(P((N > 0) | (rv['X'] < 1)), 1 - S.exp(-1)/2)
    same(P((A > 1) & (rv['B'] > 1) & (N > 0)), (S.exp(2) - 3)*(e - 2)*S.exp(-3)/2)


def test_conditional_independent(rv):
    A, B, N = rv['A'], rv['B'], rv['N']
    same(P(S.Eq(A, 2), S.Eq(B, 3)), S.exp(-1)/2)
    same(P(S.Eq(A, 2), B > 1), S.exp(-1)/2)
    same(P(S.Eq(A, 2) & S.Eq(B, 3), A > 1), 2*S.exp(-2)/(3*(e - 2)))
    same(P((A > 1) & (B > 1), S.Eq(A, 2)), 1 - 3*S.exp(-2))
    same(P(S.Eq(A, 2) & (B > 1), B < 4), 5*S.exp(-1)/19)
    same(P((rv['X'] > 1) & (N < 0), A > 0), S.exp(-1)/2)
    same(P(N > 0, (A > 1) & (B > 1)), Q(1, 2))
    assert E(A, B > 1) == 1
    assert E(N**2, A < 2) == 1
    same(E(rv['D'], A > 1), Q(7, 2))
    assert variance(A, B > 1) == 1


def test_sum_piecewise_binder_not_folded():
    x = S.Symbol('x', integer=True)
    s = S.Sum(S.Piecewise((1/S.factorial(x), x > 1), (0, True)), (x, 0, S.oo))
    assert S.simplify(s) == s
    assert not S.simplify(2*s).is_Piecewise
    y = S.Symbol('y')
    i = S.Integral(S.Piecewise((y, y > 1), (0, True)), (y, 0, 3))
    assert not S.piecewise_fold(i).is_Piecewise
    # a condition free of the bound variable may still leave the binder
    z = S.Symbol('z')
    assert S.piecewise_fold(S.Sum(S.Piecewise((x, z > 1), (0, True)), (x, 0, 3))).is_Piecewise


def test_rewrite_and_expand_symbolic():
    from symox.stats import Expectation, Variance, Covariance, Probability
    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    a, b = S.symbols('a b')
    assert Expectation(X + Y).expand() == Expectation(X) + Expectation(Y)
    assert Expectation(a*X + b).expand() == a*Expectation(X) + b
    assert Variance(a*X + b*Y).expand() == a**2*Variance(X) + b**2*Variance(Y) + 2*a*b*Covariance(X, Y)
    assert Covariance(X + Y, X).expand() == Covariance(X, X) + Covariance(X, Y)
    assert Covariance(X, X).expand() == Variance(X)
    assert Expectation(X**2).rewrite(S.Integral).doit() == 1
    assert Variance(X).rewrite(S.Integral).doit() == 1
    same(Probability(X > 1).rewrite(S.Integral).doit(), Q(1, 2) - S.erf(S.sqrt(2)/2)/2)
    A = Poisson('A', 3)
    # the native summation of k*3**k/k! is still a gap (SymPy gives 3)
    assert getattr(Expectation(A).rewrite(S.Sum), "is_Sum", False)
    assert Expectation(X).doit() == 0


def test_observed_random_symbol_condition():
    from symox.stats import given, density
    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    z = S.Symbol('z')
    assert E(X + Y, Y) == Y.symbol
    assert E(Y**2, Y) == Y.symbol**2
    assert variance(X + Y, Y) == 1
    same(P(X > 1, Y), Q(1, 2) - S.erf(S.sqrt(2)/2)/2)
    same(P(Y > 1, X), Q(1, 2) - S.erf(S.sqrt(2)/2)/2)
    # SymPy returns Probability(X > Y, Y) unevaluated
    same(P(X > Y, Y), Q(1, 2) - S.erf(S.sqrt(2)*Y.symbol/2)/2)
    assert given(X + Y, Y) == S.sympify(X) + Y.symbol
    assert given(X, Y) == X
    same(density(X + Y, Y)(z), S.sqrt(2)*S.exp(-(z - Y.symbol)**2/2)/(2*S.sqrt(S.pi)))
    same(density(X, Y)(z), S.sqrt(2)*S.exp(-z**2/2)/(2*S.sqrt(S.pi)))


def test_given_tuple_and_dependence():
    from symox.stats import given, dependent, independent
    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    A, B = Poisson('A', 1), Poisson('B', 2)
    T = given(S.Tuple(X, Y), X > 0)
    assert len(T) == 2
    same(E(T[0]), S.sqrt(2)/S.sqrt(S.pi))
    assert E(T[1]) == 0
    # SymPy reports the entries as dependent through the shared conditional
    # space; conditioning on X > 0 alone leaves them independent.
    assert independent(T[0], T[1])
    assert given(S.Tuple(X, Y), A > 1) == S.Tuple(X, Y)
    with pytest.raises(NotImplementedError):
        given(S.Tuple(X, Y), S.Eq(X + Y, 3))
    assert dependent(X, X + Y) and independent(X, Y**2) and independent(A, X)
    assert dependent(A, A + B)


def test_symbolic_soliton():
    from symox.stats import IdealSoliton, RobustSoliton, density
    k = S.Symbol('k', positive=True, integer=True)
    c = S.Symbol('c', positive=True)
    sol = IdealSoliton('sol', k)
    assert str(density(sol).dict) == 'Density(IdealSolitonDistribution(k))'
    assert density(sol).dict.subs(k, 10).doit() == {1: Q(1, 10), 2: Q(1, 2), 3: Q(1, 6), 4: Q(1, 12), 5: Q(1, 20),
                                                    6: Q(1, 30), 7: Q(1, 42), 8: Q(1, 56), 9: Q(1, 72), 10: Q(1, 90)}
    assert E(sol.subs(k, 10)) == Q(7381, 2520)
    assert P(sol.subs(k, 4) > 2) == Q(1, 4)
    assert density(sol)(3) == S.Piecewise((Q(1, 6), k >= 3), (0, True))
    with pytest.raises(NotImplementedError):
        E(sol)
    rob = RobustSoliton('rob', k, Q(1, 2), c)
    assert str(density(rob).dict).startswith('Density(RobustSolitonDistribution(k, 1/2, c')
    table = density(rob).dict.subs(k, 4).subs(c, Q(1, 10)).doit()
    assert sorted(table) == [1, 2, 3, 4]
    assert abs(float(sum(table.values())) - 1) < 1e-12
    m = S.Symbol('m')
    assert E(Normal('W', m, 1).subs(m, 2)) == 2


def test_coskewness_conditional_fast():
    import time
    from symox.stats import coskewness, Exponential, skewness
    p = S.symbols('p', positive=True)
    X, Y = Exponential('X', p), Exponential('Y', 2*p)
    start = time.time()
    assert coskewness(X + 2*Y, Y + X, Y + 2*X, X > 3) == 9*S.sqrt(170)/85
    assert time.time() - start < 20
    assert coskewness(X, Y + X, Y + 2*X) == 16*S.sqrt(85)/85
    assert coskewness(Y, Y, Y) == skewness(Y)
    same(coskewness(X, Y + p*X, Y + 2*p*X), 4/(S.sqrt(1 + 1/(4*p**2))*S.sqrt(4 + 1/(4*p**2))))


def test_markov_chain_float_symbolic():
    from symox.stats import DiscreteMarkovChain, ContinuousMarkovChain, TransitionMatrixOf
    T = S.Matrix([[0.5, 0.2, 0.3], [0.2, 0.5, 0.3], [0.2, 0.3, 0.5]])
    Y = DiscreteMarkovChain('Y', [0, 1, 2], T)
    assert abs(float(P(S.Eq(Y[3], 2), S.Eq(Y[1], 1))) - 0.36) < 1e-12
    st = Y.stationary_distribution()
    assert all(abs(float(u) - v) < 1e-12 for u, v in zip(st, [0.285714285714286, 0.339285714285714, 0.375]))
    assert Y.limiting_distribution == st
    sunny, cloudy, rainy = S.symbols('Sunny Cloudy Rainy')
    W = DiscreteMarkovChain('W', [sunny, cloudy, rainy], T)
    assert abs(float(P(S.Eq(W[3], rainy), S.Eq(W[1], cloudy))) - 0.36) < 1e-12
    e = E(W[3], S.Eq(W[1], cloudy))
    assert all(abs(float(e.coeff(s)) - v) < 1e-12 for s, v in [(cloudy, 0.38), (rainy, 0.36), (sunny, 0.26)])
    YS = DiscreteMarkovChain('YS')
    TS = S.MatrixSymbol('T', 3, 3)
    res = P(S.Eq(YS[3], 2), S.Eq(YS[1], 1) & TransitionMatrixOf(YS, TS))
    assert S.expand(res - (TS[0, 2]*TS[1, 0] + TS[1, 1]*TS[1, 2] + TS[1, 2]*TS[2, 2])) == 0
    wm = S.MatrixSymbol('wm', 1, 3)
    assert DiscreteMarkovChain('Z', [0, 1, 2], TS).stationary_distribution(condition_set=True) == S.ConditionSet(wm, S.Eq(wm*TS, wm))
    n = S.Symbol('n', nonnegative=True, integer=True)
    TN = S.MatrixSymbol('T', n, n)
    YN = DiscreteMarkovChain('YN', trans_probs=TN)
    assert YN.state_space == S.Range(0, n)
    assert P(S.Eq(YN[10], 2), S.Eq(YN[5], 1)) == (TN**5)[1, 2]
    a, c = S.symbols('a c')
    with pytest.raises(NotImplementedError):
        P(S.Eq(Y[a], 2), S.Eq(Y[c], 1))
    C = ContinuousMarkovChain('C', [0, 1], S.Matrix([[-1, 1], [1, -1]]))
    assert C.limiting_distribution() == S.ImmutableMatrix([[Q(1, 2), Q(1, 2)]])
