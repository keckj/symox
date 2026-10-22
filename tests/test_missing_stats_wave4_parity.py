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
