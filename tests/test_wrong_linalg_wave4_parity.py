"""SymPy 1.14 parity of non-symbol solve targets, solve_linear and the
matrix structure helpers diag, extract, is_symmetric, is_diagonal and kind.

Expected values are the str() output of SymPy 1.14 for the same calls."""
import pytest

# The calls below are evaluated as strings, hence the unused imports.
from symox import (Function, Identity, ImmutableMatrix, Integral, Inverse,  # noqa: F401
                   Matrix, MatrixSymbol, Rational, SparseMatrix, ZeroMatrix, cos, diag,
                   nonlinsolve, sin, solve, solveset, sqrt)
from symox.abc import t, x, y  # noqa: F401
from symox.solvers.solvers import solve_linear

f = Function('f')
M = MatrixSymbol('M', 2, 2)
N = MatrixSymbol('N', 2, 2)
m43 = Matrix(4, 3, list(range(12)))

CASES = [
    ('solve((x + 2 + sqrt(3), x + 4 + y), y, x + 2)', '{y: -2 + sqrt(3), x + 2: -sqrt(3)}'),
    ('solve(sqrt(2) - 1, 1)', '[sqrt(2)]'),
    ('solve(3 - 2, 1)', '[0]'),
    ('solve(x - 2, 2)', '[x]'),
    ('solve((x + 2 + sqrt(3), x + 4 + y), y, x + 2, dict=True)', '[{y: -2 + sqrt(3), x + 2: -sqrt(3)}]'),
    ('solve((x + 2 + sqrt(3), x + 4 + y), [y, x + 2])', '{y: -2 + sqrt(3), x + 2: -sqrt(3)}'),
    ('solve((x + 2 + sqrt(3), x + 4 + y), {y, x + 2})', '{y: -2 + sqrt(3), x + 2: -sqrt(3)}'),
    ('solve(x + 2*sqrt(3) - 1, sqrt(3))', '[1/2 - x/2]'),
    ('solveset(f(x) - 1, f(x))', '{1}'),
    ('solveset(x + 2 - y, x + 2)', '{y}'),
    ('nonlinsolve([f(x)**2 - 1], [f(x)])', '{(-1,), (1,)}'),
    ('solve_linear(y*cos(x)**2 + y*sin(x)**2 - y)', '(0, 1)'),
    ('solve_linear(x*(cos(y)**2 + sin(y)**2 - 1) + 2)', '(0, 1)'),
    ('solve_linear(2)', 'EXC AttributeError'),
    ('solve_linear(0)', 'EXC AttributeError'),
    ('solve_linear(x + Integral(1, (t, 0, 1)))', '(x, -1)'),
    ('Matrix.diag([1, 2, 3]).shape', '(3, 3)'),
    ('Matrix.diag([1, 2, 3])', 'Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])'),
    ('Matrix.diag([1, 2, 3], unpack=False)', 'Matrix([[1], [2], [3]])'),
    ('Matrix.diag([1, 2, 3], unpack=False).shape', '(3, 1)'),
    ('Matrix.diag([[1, 2, 3]])', 'Matrix([[1], [2], [3]])'),
    ('Matrix.diag([[1, 2, 3]], unpack=False)', 'Matrix([[1, 2, 3]])'),
    ('Matrix.diag(1, 2, rows=3)', 'Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 0]])'),
    ('Matrix.diag(1, 2, rows=3, cols=4)', 'Matrix([[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 0, 0]])'),
    ('Matrix.diag(1, 2, cols=3)', 'Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 0]])'),
    ('Matrix.diag(1, 2, 3, rows=2)', 'EXC ValueError'),
    ('Matrix.diag([1, 2], [3, 4], strict=True)', 'Matrix([[1, 0], [2, 0], [0, 3], [0, 4]])'),
    ('Matrix.diag([1, 2], [3, 4], strict=False)', 'Matrix([[1, 0], [2, 0], [0, 3], [0, 4]])'),
    ('Matrix.diag([[1, 2], [3]], strict=False)', 'Matrix([[1, 0], [2, 0], [0, 3]])'),
    ('Matrix.diag([[1, 2], [3]])', 'Matrix([[1, 0], [2, 0], [0, 3]])'),
    ('SparseMatrix.diag([1, 2, 3]).shape', '(3, 3)'),
    ('ImmutableMatrix.diag([1, 2, 3]).shape', '(3, 3)'),
    ('diag([1, 2, 3], unpack=False)', 'Matrix([[1], [2], [3]])'),
    ('Matrix.diag(rows=2)', 'Matrix([[0, 0], [0, 0]])'),
    ('m43.extract([0, 1, 2, 3], [True, False, True])', 'Matrix([[0, 2], [3, 5], [6, 8], [9, 11]])'),
    ('m43.extract([True, False, True, False], [0, 2])', 'Matrix([[0, 2], [6, 8]])'),
    ('m43.extract(range(4), range(3))', 'Matrix([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])'),
    ('m43.extract([True, False], [0])', 'Matrix([[0]])'),
    ('ImmutableMatrix(m43).extract([0, 1, 2, 3], [True, False, True])', 'Matrix([[0, 2], [3, 5], [6, 8], [9, 11]])'),
    ('SparseMatrix(m43).extract([0, 1, 2, 3], [True, False, True])', 'Matrix([[0, 2], [3, 5], [6, 8], [9, 11]])'),
    ('Matrix(3, 3, [1, x**2 + 2*x + 1, y, (x + 1)**2, 2, 0, y, 0, 3]).is_symmetric()', 'True'),
    ('Matrix(3, 3, [1, x**2 + 2*x + 1, y, (x + 1)**2, 2, 0, y, 0, 3]).is_symmetric(simplify=False)', 'False'),
    ('Matrix(2, 2, [1, sin(x)**2 + cos(x)**2, 1, 2]).is_symmetric()', 'True'),
    ('Matrix(2, 2, [1, sin(x)**2 + cos(x)**2, 1, 2]).is_symmetric(simplify=False)', 'False'),
    ('ImmutableMatrix(2, 2, [1, (x + 1)**2, x**2 + 2*x + 1, 2]).is_symmetric()', 'True'),
    ('SparseMatrix(2, 2, [1, (x + 1)**2, x**2 + 2*x + 1, 2]).is_symmetric()', 'True'),
    ('Matrix(2, 2, [1, (x + 1)**2 - x**2 - 2*x - 1, 0, 2]).is_diagonal()', 'None'),
    ('Matrix(2, 2, [1, 0, sin(x)**2 + cos(x)**2 - 1, 1]).is_diagonal()', 'None'),
    ('Matrix(2, 2, [x, y, y, x]).is_symmetric(simplify=False)', 'True'),
    ('Matrix(2, 2, [1, x*(x + 1), x**2 + x, 2]).is_symmetric()', 'True'),
    ('Matrix(2, 2, [1, x*(x + 1), x**2 + x, 2]).is_symmetric(simplify=False)', 'False'),
    ('Matrix(2, 2, [1, x*(x + 1), x**2 + x, 2]).is_symmetric(simplify=lambda e: e.expand())', 'True'),
    ("MatrixSymbol('A', 2, 2).kind", 'MatrixKind(NumberKind)'),
    ("repr(MatrixSymbol('A', 2, 2).kind)", 'MatrixKind(NumberKind)'),
    ('SparseMatrix([[1, 2], [3, 4]]).kind', 'MatrixKind(NumberKind)'),
    ('(M*N).kind', 'MatrixKind(NumberKind)'),
    ('(M + N).kind', 'MatrixKind(NumberKind)'),
    ('Identity(2).kind', 'MatrixKind(NumberKind)'),
    ('ZeroMatrix(2, 2).kind', 'MatrixKind(NumberKind)'),
    ('M.T.kind', 'MatrixKind(NumberKind)'),
    ('(M**2).kind', 'MatrixKind(NumberKind)'),
    ("Matrix([[M]]).kind if False else MatrixSymbol('A', 2, 2).kind.element_kind", 'NumberKind'),
    ("type(MatrixSymbol('A', 2, 2).kind).__name__", 'MatrixKind'),
    ("MatrixSymbol('A', 2, 2).kind == Matrix([[1]]).kind", 'True'),
    ('(2*M).kind', 'MatrixKind(NumberKind)'),
    ('Inverse(M).kind', 'MatrixKind(NumberKind)'),
]


@pytest.mark.parametrize("call,expected", CASES)
def test_reference_output(call, expected):
    if expected.startswith('EXC '):
        with pytest.raises(Exception) as info:
            eval(call)
        assert type(info.value).__name__ == expected[4:]
    else:
        assert str(eval(call)) == expected


def test_solve_number_target_replaces_its_negative():
    # Number._eval_subs also replaces -2, so x - 2 is linear in the target.
    assert solve(x - 2, 2) == [x]
    assert solve(x + 2, 2, dict=True) == [{2: -x}]


def test_solve_expression_target_keys_are_restored():
    solution = solve((x + 2 + sqrt(3), x + 4 + y), y, x + 2)
    assert list(solution) == [y, x + 2]
    assert solution[x + 2] == -sqrt(3)


def test_solve_linear_error_message():
    with pytest.raises(ValueError, match="solve_linear only handles symbols"):
        solve_linear(x + y, symbols=[x + 1])


def test_matrix_diag_class_and_rows():
    assert isinstance(ImmutableMatrix.diag(1, 2), ImmutableMatrix)
    assert isinstance(Matrix.diag(1, 2, cls=ImmutableMatrix), ImmutableMatrix)
    with pytest.raises(ValueError, match='is not an integer'):
        Matrix.diag(1, rows=Rational(3, 2))


def test_extract_index_errors():
    with pytest.raises(IndexError):
        m43.extract([4], [0])
    assert m43.extract([-1], [-1]) == Matrix([[11]])
    with pytest.raises(TypeError):
        m43.extract(1, [0])
