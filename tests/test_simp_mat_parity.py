"""Simplification and matrix parity regressions (wave 6, simp-mat).

Expected strings are SymPy 1.14 outputs.
"""
import pytest

import symox as sp

x, y = sp.symbols('x y')
j = sp.Symbol('j', integer=True)
q = (x**2 - 1)/(x - 1)
K = sp.KroneckerDelta

CASES = [
    ('sp.nsimplify(1 - sp.GoldenRatio)',
     '1/2 - sqrt(5)/2'),
    ('sp.nsimplify((1 + sp.sqrt(5))/4, [sp.GoldenRatio])',
     'GoldenRatio/2'),
    ('sp.nsimplify(.2, tolerance=0)',
     '1/5'),
    ('sp.nsimplify(-.2222, tolerance=0)',
     '-1111/5000'),
    ('sp.nsimplify(sp.S(2e-8))',
     '1/50000000'),
    ('sp.nsimplify(1/.3 + x, rational=True)',
     'x + 10/3'),
    ('sp.simplify(sp.sign(x)*sp.Abs(x))',
     'x'),
    ('sp.simplify(sp.Abs(x)/sp.Abs(x**2))',
     '1/Abs(x)'),
    ('sp.simplify(5*sp.Abs(q))',
     '5*Abs(x + 1)'),
    ('sp.simplify(sp.sin(q))',
     'sin(x + 1)'),
    ('sp.simplify(sp.sign(q))',
     'sign((x**2 - 1)/(x - 1))'),
    ('sp.simplify(sp.gamma(q))',
     'gamma(x + 1)'),
    ('sp.simplify(sp.floor(q))',
     'floor(x) + 1'),
    ('sp.simplify(sp.hyper([], [], x))',
     'exp(x)'),
    ('sp.simplify(K(0, j)*K(1, j))',
     '0'),
    ('(3 - (-1)**sp.Rational(2, 3) + (-1)**sp.Rational(1, 3)).equals(4)',
     'True'),
    ('sp.radsimp(1/(3 - (-1)**sp.Rational(2, 3) + (-1)**sp.Rational(1, 3)))',
     '1/4'),
    ('sp.Matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]]).is_positive_definite',
     'True'),
    ('sp.Matrix([[1, 2], [2, 1]]).is_indefinite',
     'True'),
    ('sp.Matrix([[1, 1], [1, 1]]).is_positive_semidefinite',
     'True'),
    ('sp.Matrix([[1, 1], [1, 1]]).is_positive_definite',
     'False'),
    ('sp.Matrix([[-2, 1], [1, -2]]).is_negative_definite',
     'True'),
    ('sp.Matrix([[x, 0], [0, 1]]).is_positive_definite',
     'None'),
    ('sp.Matrix([[1, 2], [3, 4]]).is_positive_definite',
     'False'),
    ('sp.Matrix([[1, 2, 3], [0, 0, 4], [0, 0, 0]]).is_echelon',
     'True'),
    ('sp.Matrix([[0, 1], [1, 0]]).is_echelon',
     'False'),
    ('sp.Matrix([[3, -1], [1, 4]]).is_strongly_diagonally_dominant',
     'True'),
    ('sp.Matrix([[1.5, 2], [3, 4]]).eigenvals()',
     '{-2.12902801124989e-32: 1, 5.50000000000000: 1}'),
    ('sp.Matrix([[1.5, 2], [3, 4]]).eigenvects()',
     '[(-2.12902801124989e-32, 1, [Matrix([ [ 0.8], [-0.6]])]), (5.50000000000000, 1, [Matrix([ [0.454545454545455], [0.909090909090909]])])]'),
    ('sp.Matrix([[0, -1.0], [1, 0]]).eigenvals(multiple=True)',
     '[-1.0*I, 1.0*I]'),
    ('sp.Matrix([[1, 2], [3, 4]]).eigenvals(simplify=True)',
     '{5/2 - sqrt(33)/2: 1, 5/2 + sqrt(33)/2: 1}'),
    ('sp.Matrix([[-0.5, 0], [0, 0.5]]).diagonalize()',
     '(Matrix([ [1.0,   0], [  0, 1.0]]), Matrix([ [-0.5,   0], [   0, 0.5]]))'),
]


@pytest.mark.parametrize("code, expected", CASES)
def test_simp_mat_matches_sympy(code, expected):
    assert str(eval(code)).replace('\n', ' ') == expected


def test_float_eigenvals_precision_exhausted_like_sympy():
    M = sp.Matrix([[sp.Float('0.1', 30), 1], [1, 2]])
    with pytest.raises(sp.PrecisionExhausted):
        M.eigenvals()
