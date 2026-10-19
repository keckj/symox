"""Numeric skew exponentials compared with an independent matrix power series."""
import mpmath as mp
import numpy as np
import pytest
import sympy as sp
import symox as sx


def _oracle(matrix, digits):
    # Preserve the exact stored floating inputs rather than their short strings.
    with mp.workdps(digits):
        values = []
        for row in matrix.to_sympy().tolist():
            converted = []
            for entry in row:
                numerator, denominator = sp.Rational(entry).as_numer_denom()
                converted.append(mp.mpf(str(numerator))/mp.mpf(str(denominator)))
            values.append(converted)
        return mp.expm(mp.matrix(values))


@pytest.mark.parametrize('size', [2, 3])
def test_skew_exponential_matches_reference_and_numeric_generator(size):
    matrix = sx.Matrix([[0.0, -0.1], [0.1, 0.0]]) if size == 2 else sx.Matrix(
        [[0.0, -0.3, -0.1], [0.3, 0.0, -0.2], [0.1, 0.2, 0.0]])
    result = matrix.exp()
    expected = matrix.to_sympy().exp()
    np.testing.assert_allclose(np.array(result.tolist(), dtype=float),
                               np.array(expected.tolist(), dtype=float), rtol=3e-15, atol=3e-15)
    assert result.shape == matrix.shape
    # The reference generator formats entries directly for C++ test data.
    assert all(format(float(value), '15g') for value in result)


@pytest.mark.parametrize('size', [2, 3])
@pytest.mark.parametrize('angle', ['0.123456789012345678901234567890123456789', '1e-300', '1e40',
    '3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628'])
def test_multiprecision_skew_exponential_keeps_input_accuracy(size, angle):
    x = sx.Float(angle, 80)
    matrix = sx.Matrix([[0, -x], [x, 0]]) if size == 2 else sx.Matrix(
        [[0, -x, -2*x], [x, 0, -3*x], [2*x, 3*x, 0]])
    result = matrix.exp()
    with mp.workdps(110):
        expected = _oracle(matrix, 180)
        for i in range(size):
            for j in range(size):
                value = mp.mpf(str(result[i,j].evalf(90)))
                target = expected[i,j]
                if target:
                    assert abs((value-target)/target) < mp.mpf('1e-75')
                else:
                    assert value == 0


def test_subnormal_angle_is_not_lost_by_squaring_in_machine_precision():
    x = sx.Float(1e-300)
    matrix = sx.Matrix([[0, -x, 0], [x, 0, 0], [0, 0, 0]])
    result = matrix.exp()
    assert result[1,0] != 0
    assert float(result[1,0]) == float(x)
    assert float(result[0,1]) == -float(x)


def test_exact_and_error_paths_remain_available():
    matrix = sx.Matrix([[1,2,0], [0,1,3], [0,0,1]])
    assert matrix.exp().to_sympy() == matrix.to_sympy().exp()
    for entry in [sx.oo, sx.nan, sx.zoo]:
        with pytest.raises(ValueError):
            sx.Matrix([[0, -entry], [entry, 0]]).exp()
    with pytest.raises(sx.NonSquareMatrixError):
        sx.ones(2,3).exp()


def test_matrix_negation_normalizes_floating_zeros_like_reference():
    matrix = sx.Matrix([[0.0, 0.25], [-0.25, 0.0]])
    assert (-matrix).to_sympy() == -matrix.to_sympy()
    assert sx.exp(-sx.Matrix([0.0])).to_sympy() == sp.exp(-sp.Matrix([0.0]))


def test_diagonal_zero_retains_its_existing_large_precision_domain():
    precision = 600_000
    actual = sx.diag(sx.Float(0, precision=precision), sx.Float(0, precision=precision)).exp()
    expected = sp.diag(sp.Float(0, precision=precision), sp.Float(0, precision=precision)).exp()
    for i in range(2):
        assert actual[i,i]-1 == 0
        assert actual[i,i]._prec == expected[i,i]._prec
    assert actual[0,1] == actual[1,0] == 0
