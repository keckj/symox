"""API surface added for SymPy 1.14's array and matrix-expression tests.

Expected values are SymPy 1.14 outputs.
"""
from symox import Array, Matrix, MatrixSymbol, det as top_det
from symox.tensor.array import ImmutableSparseNDimArray
from symox.tensor.array.arrayop import Flatten


def test_flatten():
    assert str(Flatten(Array(range(6)).reshape(2, 3))) == 'Flatten([[0, 1, 2], [3, 4, 5]])'
    assert str(Flatten(Matrix([[1, 2]]))) == 'Flatten(Matrix([[1, 2]]))'
    assert list(Flatten(Array(range(4)).reshape(2, 2))) == [0, 1, 2, 3]
    assert list(Flatten(ImmutableSparseNDimArray([0, 1, 0, 3], (2, 2)))) == [0, 1, 0, 3]
    assert list(Flatten(Matrix([[1, 2], [3, 4]]))) == [1, 2, 3, 4]
    assert list(Flatten([[1, 2], [3, 4]])) == [1, 2, 3, 4]
    for i, v in enumerate(Flatten(Array(range(24)).reshape(4, 6))):
        assert i == v


def test_determinant_module_names():
    from symox.matrices.expressions.determinant import Permanent, det, per
    assert det is top_det
    assert str(det(MatrixSymbol('A', 2, 2))) == 'Determinant(A)'
    assert det(Matrix([[1, 2], [3, 4]])) == -2
    assert per(Matrix([[1, 2], [3, 4]])) == 10
    assert str(Permanent(MatrixSymbol('A', 2, 2))) == 'Permanent(A)'
