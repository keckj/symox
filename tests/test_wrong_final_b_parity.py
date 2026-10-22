"""Regressions for wrong results fixed against SymPy 1.14.

Expected values are hard-coded from the reference SymPy output.
"""
import pytest

import symox as S
from symox.physics.quantum import Dagger, TensorProduct
from symox.stats import (CrossCovarianceMatrix, Exponential, Expectation,
                         ExpectationMatrix, Normal, VarianceMatrix)
from symox.stats.rv import RandomMatrixSymbol


A, B, C, D = S.symbols('A B C D', commutative=False)
x, y = S.symbols('x y')


def _s(value):
    return str(value).replace('_', '')


# --- Dagger / adjoint of noncommutative symbols ---------------------------

def test_dagger_symbol_is_held():
    assert str(Dagger(A)) == 'Dagger(A)'
    assert str(S.adjoint(A)) == 'Dagger(A)'
    assert Dagger(Dagger(A)) == A


def test_dagger_reverses_products():
    assert Dagger(A*B) == Dagger(B)*Dagger(A)
    assert str(Dagger(A*B)) == 'Dagger(B)*Dagger(A)'
    assert str(Dagger(A*B*C)) == 'Dagger(C)*Dagger(B)*Dagger(A)'
    assert str(Dagger(x*A*B)) == 'conjugate(x)*Dagger(B)*Dagger(A)'
    assert str(S.adjoint(A*B)) == 'Dagger(B)*Dagger(A)'
    assert str(Dagger(A*B).doit()) == 'Dagger(B)*Dagger(A)'


def test_dagger_add_pow_numbers():
    assert str(Dagger(A + B)) == 'Dagger(A) + Dagger(B)'
    assert str(Dagger(A**2)) == 'Dagger(A)**2'
    assert str(Dagger(2*S.I*A)) == '-2*I*Dagger(A)'
    assert Dagger(S.I) == -S.I
    assert Dagger(x) == S.conjugate(x)


def test_dagger_matrix():
    assert Dagger(S.Matrix([[1, S.I], [2, 3]])) == S.Matrix([[1, 2], [-S.I, 3]])


def test_conjugate_transpose_nc_products():
    assert str(S.conjugate(A*B)) == 'conjugate(A)*conjugate(B)'
    assert str(S.transpose(A*B)) == 'transpose(B)*transpose(A)'


def test_tensor_product_dagger_and_expand():
    assert str(Dagger(TensorProduct(A, B))) == 'Dagger(A)xDagger(B)'
    assert str(TensorProduct(A + B, C).expand(tensorproduct=True)) == 'AxC + BxC'
    e = TensorProduct(A, B)*TensorProduct(C, D)
    assert str(e) == '(A*C)x(B*D)'
    assert str(Dagger(e)) == '(Dagger(C)*Dagger(A))x(Dagger(D)*Dagger(B))'
    m1 = S.Matrix([[1, 2], [3, 4]])
    m2 = S.Matrix([[1, 0], [0, 1]])
    assert TensorProduct(m1, m2) == S.Matrix(
        [[1, 0, 2, 0], [0, 1, 0, 2], [3, 0, 4, 0], [0, 3, 0, 4]])


def test_dagger_latex():
    assert S.latex(Dagger(A)) == 'A^{\\dagger}'


# --- Matrix expressions with random matrix symbols ------------------------

k = S.Symbol('k')
MA = S.MatrixSymbol('A', k, k)
MB = S.MatrixSymbol('B', k, k)
MC = S.MatrixSymbol('C', k, k)
MD = S.MatrixSymbol('D', k, k)
X = RandomMatrixSymbol('X', k, 1)
Y = RandomMatrixSymbol('Y', k, 1)
Z = RandomMatrixSymbol('Z', k, 1)
W = RandomMatrixSymbol('W', k, 1)
a, b = S.symbols('a b')


def test_matrix_times_random_matrix_keeps_order():
    assert str(MA*X) == 'A*X'
    assert (MA*X).shape == (k, 1)
    assert str(X.T*MA) == 'X.T*A'
    assert (X.T*MA).shape == (1, k)
    assert str(MA*X + MB*Y) == 'A*X + B*Y'
    assert (MA*X + MB*Y).shape == (k, 1)


def test_expectation_matrix():
    assert str(ExpectationMatrix(X)) == 'ExpectationMatrix(X)'
    assert ExpectationMatrix(X).shape == (k, 1)
    assert ExpectationMatrix(MA*X).shape == (k, 1)
    assert str(ExpectationMatrix(MA*X).expand()) == 'A*ExpectationMatrix(X)'
    assert str(ExpectationMatrix(MA*X + MB*Y).expand()) == \
        'A*ExpectationMatrix(X) + B*ExpectationMatrix(Y)'
    assert str(ExpectationMatrix(X + Y).expand()) == \
        'ExpectationMatrix(X) + ExpectationMatrix(Y)'
    assert str(ExpectationMatrix(MA*X + MB).expand()) == 'A*ExpectationMatrix(X) + B'


def test_variance_matrix():
    assert VarianceMatrix(X).shape == (k, k)
    assert str(VarianceMatrix(MA*X)) == 'VarianceMatrix(A*X)'
    assert str(VarianceMatrix(MA*X).expand()) == 'A*VarianceMatrix(X)*A.T'
    assert str(VarianceMatrix(MA*X + MB*Y).expand()) == (
        '2*A*CrossCovarianceMatrix(X, Y)*B.T + A*VarianceMatrix(X)*A.T'
        ' + B*VarianceMatrix(Y)*B.T')
    assert str(VarianceMatrix(X + Y).expand()) == \
        '2*CrossCovarianceMatrix(X, Y) + VarianceMatrix(X) + VarianceMatrix(Y)'


def test_cross_covariance_matrix():
    assert CrossCovarianceMatrix(X, Y).shape == (k, k)
    assert CrossCovarianceMatrix(MA*X, Y).shape == (k, k)
    assert str(CrossCovarianceMatrix(X + Y, Z).expand()) == \
        'CrossCovarianceMatrix(X, Z) + CrossCovarianceMatrix(Y, Z)'
    assert str(CrossCovarianceMatrix(MA*X, Y).expand()) == 'A*CrossCovarianceMatrix(X, Y)'
    assert str(CrossCovarianceMatrix(MA*X, MB*Y).expand()) == \
        'A*CrossCovarianceMatrix(X, Y)*B.T'
    assert str(CrossCovarianceMatrix(MA*X + MB*Y, MC*Z + MD*W).expand()) == (
        'A*CrossCovarianceMatrix(X, W)*D.T + A*CrossCovarianceMatrix(X, Z)*C.T'
        ' + B*CrossCovarianceMatrix(Y, W)*D.T + B*CrossCovarianceMatrix(Y, Z)*C.T')


def test_scalar_times_random_matrix():
    assert (a*X).shape == (k, 1)
    assert str(VarianceMatrix(a*X).expand()) == 'a**2*VarianceMatrix(X)'
    # SymPy prints the unflattened scalar as '(a*b)*CrossCovarianceMatrix(X, Y)'.
    r = CrossCovarianceMatrix(a*X, b*Y).expand()
    assert r == a*b*CrossCovarianceMatrix(X, Y)
    assert r.shape == (k, k)


# --- Nested Expectation doit ----------------------------------------------

mu = S.Symbol('mu')
NX = Normal('X', mu, 1)
EY = Exponential('Y', 1)


def test_expectation_nested_doit():
    assert Expectation(NX + EY).doit() == mu + 1
    assert Expectation(NX + Expectation(EY + Expectation(2*NX))).doit() == 3*mu + 1
    assert Expectation(Expectation(EY)).doit() == 1
    assert Expectation(Expectation(NX) + EY).doit() == mu + 1
    assert Expectation(a*NX + EY).doit() == a*mu + 1
    assert Expectation(NX*Expectation(EY)).doit() == mu
    assert Expectation(Expectation(NX)**2).doit() == mu**2


def test_expectation_doit_not_deep():
    r = Expectation(NX + Expectation(EY)).doit(deep=False)
    assert _s(r) == 'mu + Expectation(Expectation(Y))'
    assert r.doit() == mu + 1
    assert _s(Expectation(Expectation(EY)).doit(deep=False)) == 'Expectation(Y)'


# --- nonlinsolve with rational equations ----------------------------------

def test_nonlinsolve_rational_doc_example():
    r = S.nonlinsolve([x**2 + 2/y - 2, x + y - 3], [x, y])
    assert r == S.FiniteSet(
        S.Tuple(1, 2),
        S.Tuple(1 - S.sqrt(5), 2 + S.sqrt(5)),
        S.Tuple(1 + S.sqrt(5), 2 - S.sqrt(5)))


@pytest.mark.parametrize('eqs, expected', [
    ([x + 1/y - 2, x - y], {(1, 1)}),
    ([x*y - 1, x + 1/y - 2], {(1, 1)}),
    ([x**2/y - 1, x + y - 2], {(-2, 4), (1, 1)}),
    ([x/(y - 1) - 1, x + y - 3], {(1, 2)}),
    # SymPy also returns the spurious (0, 0), where 1/x is undefined;
    # symox drops roots that zero a denominator (intentional divergence).
    ([1/x + 1/y - 1, x - 2*y], {(3, S.Rational(3, 2))}),
])
def test_nonlinsolve_rational_variants(eqs, expected):
    r = S.nonlinsolve(eqs, [x, y])
    assert r == S.FiniteSet(*[S.Tuple(*t) for t in expected])


def test_matrix_statistics_doc_examples():
    r = ExpectationMatrix((X + Y)*(X - Y).T).expand()
    assert r == (ExpectationMatrix(X*X.T) - ExpectationMatrix(X*Y.T)
                 + ExpectationMatrix(Y*X.T) - ExpectationMatrix(Y*Y.T))
    assert r.shape == (k, k)
    assert str(CrossCovarianceMatrix(MA*X, MB.T*Y).expand()) == \
        'A*CrossCovarianceMatrix(X, Y)*B'
    assert str(CrossCovarianceMatrix(MA*X + MB*Y, MC.T*Z + MD.T*W).expand()) == (
        'A*CrossCovarianceMatrix(X, W)*D + A*CrossCovarianceMatrix(X, Z)*C'
        ' + B*CrossCovarianceMatrix(Y, W)*D + B*CrossCovarianceMatrix(Y, Z)*C')
