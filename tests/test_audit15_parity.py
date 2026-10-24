"""Audit 15 parity regressions against SymPy 1.14 (expected values hard-coded
from the reference)."""
import pytest

from symox import (
    Add, FiniteSet, ImageSet, Intersection, Interval, Lambda, MatAdd, MatMul,
    MatPow, Matrix, MatrixSymbol, PermutationMatrix, MatrixPermute, S, Union,
    beta, latex, pretty, symbols,
)
from symox.combinatorics import Permutation
from symox.matrices.exceptions import ShapeError
from symox.tensor.array.expressions import ArraySymbol, Reshape

x, y, z, a, b, n = symbols('x y z a b n')
A, B = MatrixSymbol('A', 2, 2), MatrixSymbol('B', 2, 2)
C = MatrixSymbol('C', 2, 3)


def test_subs_promotes_scalar_operations_to_matrix_expressions():
    e = Add(x, y).subs({x: A, y: B})
    assert isinstance(e, MatAdd) and str(e) == 'A + B'
    e = (x*y).subs({x: A, y: B})
    assert isinstance(e, MatMul) and str(e) == 'B*A'
    e = (2*x).subs(x, A)
    assert isinstance(e, MatMul) and str(e) == '2*A'
    e = (x*y).subs(x, A)
    assert isinstance(e, MatMul) and str(e) == 'y*A'
    e = (x**2).subs(x, A)
    assert isinstance(e, MatPow) and str(e) == 'A**2'
    e = (x*y + z).subs({x: A, y: B, z: A})
    assert isinstance(e, MatAdd) and str(e) == 'A + B*A'
    assert str((x + y).subs(x, A)) == 'y + A'
    assert x.subs(x, A) == A


def test_xreplace_promotes_scalar_operations_to_matrix_expressions():
    e = (x + y).xreplace({x: A, y: B})
    assert isinstance(e, MatAdd) and str(e) == 'A + B'
    e = (x*y).xreplace({x: A, y: B})
    assert isinstance(e, MatMul) and str(e) == 'A*B'


def test_subs_matrix_shape_mismatch_raises():
    with pytest.raises(ShapeError):
        (x*y).subs({x: A, y: C})
    with pytest.raises(ShapeError):
        (x + y).subs({x: A, y: C})


def test_array_symbol_elements_print_like_reference():
    r = Reshape(A, (4,))
    assert str(r.as_explicit()) == '[A[0, 0], A[0, 1], A[1, 0], A[1, 1]]'
    Bs = ArraySymbol('B', (2, 3))
    assert str(Bs.as_explicit()) == (
        '[[B[0, 0], B[0, 1], B[0, 2]], [B[1, 0], B[1, 1], B[1, 2]]]')
    assert latex(Bs.as_explicit()[1, 2]) == '{{B}_{1, 2}}'
    As = ArraySymbol('A', (6,))
    assert str(As.as_explicit()[0]) == 'A[0]'


def test_permutation_matrix_pretty():
    P = PermutationMatrix(Permutation([0, 3, 8, 1, 6, 7, 2, 4, 5]))
    assert pretty(P) == 'PermutationMatrix((1 3)(2 8 5 7 4 6))'
    M = MatrixSymbol('M', 3, 3)
    assert pretty(MatrixPermute(M, Permutation([1, 2, 0]), 0)) == (
        'MatrixPermute(M, (0 1 2), 0)')


def test_beta_pretty_and_latex():
    assert pretty(beta(a, b), use_unicode=False) == 'B(a, b)'
    assert pretty(beta(a, b), use_unicode=True) == 'Β(a, b)'
    assert latex(beta(a, b)) == r'\operatorname{B}\left(a, b\right)'
    assert latex(beta(a, b)**2) == r'\operatorname{B}^{2}\left(a, b\right)'
    assert latex(beta(a)) == r'\operatorname{B}\left(a, a\right)'


def test_pinv_solve_arbitrary_dummy_names():
    M = Matrix([[1, 2], [2, 4]])
    r = M.pinv_solve(Matrix([3, 6]))
    assert str(r) == ('Matrix([[4*_w0_0/5 - 2*_w1_0/5 + 3/5], '
                      '[-2*_w0_0/5 + _w1_0/5 + 6/5]])')
    assert sorted(str(s) for s in r.free_symbols) == ['_w0_0', '_w1_0']


def test_finite_set_meets_image_set():
    sq = ImageSet(Lambda(x, x**2), S.Naturals)
    assert FiniteSet(0, 1, 2, 3, 4, 5, 6, 7, 9, 10).intersect(sq) == \
        FiniteSet(1, 4, 9)
    assert str(Intersection(FiniteSet(1, 2, x), sq)) == \
        'Intersection({1, x}, ImageSet(Lambda(x, x**2), Naturals))'
    assert str(Intersection(FiniteSet(1, 2, x, y), sq)) == \
        'Intersection({1, x, y}, ImageSet(Lambda(x, x**2), Naturals))'
    assert Intersection(FiniteSet(2, 3, 5), sq) == S.EmptySet
    assert Intersection(FiniteSet(1, 4), sq, Interval(0, 2)) == FiniteSet(1)
    assert Intersection(FiniteSet(1, 4, 9), sq, FiniteSet(4, 9, 16)) == \
        FiniteSet(4, 9)
    sqz = ImageSet(Lambda(n, n**2), S.Integers)
    assert FiniteSet(-4, 4, S.ImaginaryUnit).intersect(sqz) == FiniteSet(4)
    assert str(Union(FiniteSet(1, 2), Interval(5, 6)).intersect(sq)) == (
        'Union({1}, Intersection(ImageSet(Lambda(x, x**2), Naturals), '
        'Interval(5, 6)))')
