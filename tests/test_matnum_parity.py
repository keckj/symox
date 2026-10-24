"""Float matrix functions, inexact inverses and pivot zero tests (SymPy 1.14).

Expected values are hard-coded from SymPy 1.14.
"""
import pytest

from symox import *
del test

x = Symbol('x')


def _F():
    return Matrix([[Float('0.1', 30), Float('0.2', 30)],
                   [Float('0.3', 30), Float('0.4', 30)]])


def test_float_exp_via_jordan_form():
    assert str(Matrix([[1.0, 2.0], [3.0, 4.0]]).exp()) == (
        'Matrix([[51.9689561987050, 74.7365645670032], '
        '[112.104846850505, 164.073803049210]])')
    assert str(Matrix([[1.5, 2.5], [0.5, 1.0]]).exp()) == (
        'Matrix([[7.11901824023183, 10.7638891811406], '
        '[2.15277783622812, 4.96624040400371]])')


def test_float_exp_keeps_precision():
    E = _F().exp()
    assert str(E) == (
        'Matrix([[1.14209354891957264311250966407, '
        '0.260350747540692697829791492894], '
        '[0.390526121311039046744687239342, '
        '1.53261967023061168985719690341]])')
    assert {e._prec for e in E} == {103}


def test_float_log():
    assert str(Matrix([[2.0, 1.0], [1.0, 2.0]]).log()) == (
        'Matrix([[0.549306144334055, 0.549306144334055], '
        '[0.549306144334055, 0.549306144334055]])')


def test_float_diagonal_power_exact_zeros():
    G = Matrix([[1.5, 0.0], [0.0, 2.5]])
    assert str(G**0.5) == (
        'Matrix([[1.22474487139159, 0], [0, 1.58113883008419]])')
    assert (G**0.5)[0, 1] is S.Zero


def test_inexact_inverse_matches_rr_and_cc_domains():
    assert str(Matrix([[1.5, 2.5], [0.5, 1.0]]).inv()) == (
        'Matrix([[4.00000000000000, -9.99999999999999], '
        '[-2.00000000000000, 6.00000000000000]])')
    assert str(Matrix([[1.0, 2], [3, 4]]).inv()) == (
        'Matrix([[-2.00000000000000, 1.00000000000000], '
        '[1.50000000000000, -0.500000000000000]])')
    assert str(_F().inv()) == (
        'Matrix([[-20.0000000000000000000000000000, '
        '10.0000000000000000000000000000], '
        '[15.0000000000000000000000000000, '
        '-5.00000000000000000000000000000]])')
    assert str(Matrix([[1.0 + I, 0.5], [0.25, 2.0]]).inv()) == (
        'Matrix([[0.498960498960499 - 0.532224532224532*I, '
        '-0.124740124740125 + 0.133056133056133*I], '
        '[-0.0623700623700624 + 0.0665280665280665*I, '
        '0.515592515592516 - 0.0166320166320166*I]])')
    with pytest.raises(NonInvertibleMatrixError,
                       match='Matrix det == 0; not invertible.'):
        Matrix([[1.0, 2.0], [2.0, 4.0]]).inv()


def test_float_eigenvects_signs():
    assert str(_F().eigenvects()) == (
        '[(-0.0372281323269014, 1, [Matrix([\n[-0.824564840132394],\n'
        '[ 0.565767464968992]])]), (0.537228132326901, 1, [Matrix([\n'
        '[ -0.42222915041526],\n[-0.923052314250193]])])]')


def test_float_is_diagonalizable():
    assert Matrix([[1.5, 2.5], [0.5, 1.0]]).is_diagonalizable() is True
    assert Matrix([[1.0, 2.0], [3.0, 4.0]]).is_diagonalizable() is True
    assert Matrix([[0.5, 1], [0, 0.5]]).is_diagonalizable() is True
    assert Matrix([[1, 1], [0, 1]]).is_diagonalizable() is False


def test_pivot_recognises_function_relations():
    Z = Matrix([[cos(2*x) - cos(x)**2 + sin(x)**2, 1], [0, 0]])
    assert Z.rref() == (Matrix([[0, 1], [0, 0]]), (1,))
    assert Z.nullspace() == [Matrix([[1], [0]])]
    assert Z.rank() == 1
    assert Z.echelon_form() == Matrix([[0, 1], [0, 0]])
    T = Matrix([[sin(x)**2 + cos(x)**2 - 1, x], [x, 1]])
    assert T.rank() == 2
    assert T.rref() == (eye(2), (0, 1))
    assert T.nullspace() == []
    assert str(T.det()) == '-x**2 + sin(x)**2 + cos(x)**2 - 1'
