"""Exact generalized inverses and structured systems, checked independently."""
import pytest
import sympy as sp
import symox as sx


def equal(actual, expected):
    result = actual.to_sympy()
    assert result.shape == expected.shape
    assert all(sp.simplify(a - b) == 0 for a, b in zip(result, expected))


@pytest.mark.parametrize("data", [
    [[1, 2, 3], [2, 4, 6]],
    [[1, 0], [0, 1], [1, 1]],
    [[1, 2, 3]],
    [[0, 0], [0, 0]],
    [[1, sp.I], [sp.I, -1]],
    [[1, sp.I], [0, 2], [2, 1]],
    [[sp.Symbol('x', positive=True), 2 * sp.Symbol('x', positive=True)], [1, 2]],
    [[sp.Symbol('z')], [1]],
])
def test_rank_decomposition_and_penrose(data):
    a = sx.Matrix(data)
    ref = sp.Matrix(data)
    c, f = a.rank_decomposition()
    rc, rf = ref.rank_decomposition()
    equal(c, rc)
    equal(f, rf)
    equal(c * f, ref)
    p = a.pinv()
    equal(p, ref.pinv())
    # All four defining identities; no reference inversion is needed for these.
    equal(a * p * a, ref)
    equal(p * a * p, p.to_sympy())
    equal((a * p).H, (a * p).to_sympy())
    equal((p * a).H, (p * a).to_sympy())


@pytest.mark.parametrize("shape", [(0, 0), (0, 3), (4, 0), (2, 3)])
def test_zero_and_empty_pseudoinverse(shape):
    a = sx.zeros(*shape)
    c, f = a.rank_decomposition()
    assert c.shape == (shape[0], 0)
    assert f.shape == (0, shape[1])
    assert a.pinv().shape == shape[::-1]
    equal(c * f, sp.zeros(*shape))
    b = sx.zeros(shape[0], 2)
    equal(a.pinv_solve(b, sx.zeros(shape[1], 2)), sp.zeros(shape[1], 2))


def test_pinv_solve_arbitrary_nullspace_and_fresh_parameters():
    a = sx.Matrix([[1, 2, 3], [2, 4, 6]])
    b = sx.Matrix([[1, 2], [3, 4]])
    w = sx.Matrix([[1, 2], [3, 4], [5, 6]])
    ref, rb, rw = a.to_sympy(), b.to_sympy(), w.to_sympy()
    equal(a.pinv_solve(b, arbitrary_matrix=w), ref.pinv_solve(rb, arbitrary_matrix=rw))
    minnorm = a.pinv_solve(b, sx.zeros(3, 2))
    equal(minnorm, ref.pinv() * rb)
    first, second = a.pinv_solve(b), a.pinv_solve(b)
    assert first.free_symbols
    assert first.free_symbols.isdisjoint(second.free_symbols)
    equal(a.H * (a * first - b), sp.zeros(3, 2))
    equal(a * (first - minnorm), sp.zeros(2, 2))
    with pytest.raises(sx.ShapeError):
        a.pinv_solve(b, sx.zeros(2, 2))
    with pytest.raises(sx.ShapeError):
        a.pinv_solve(sx.zeros(3, 2))


@pytest.mark.parametrize("method", ["CH", "LDL", "QR", "PINV", "GJ", "GE", "LU", "ADJ", "DM", "DMNC"])
@pytest.mark.parametrize("data", [
    [[1, 0], [1, 1], [1, 2]],
    [[1, sp.I], [0, 2], [2, 1]],
])
def test_least_squares_multiple_rhs(method, data):
    a = sx.Matrix(data)
    b = sx.Matrix([[1, 2], [2, 1], [3, -1]])
    result = a.solve_least_squares(b, method=method)
    ref = a.to_sympy()
    expected = (ref.H * ref).inv() * ref.H * b.to_sympy()
    equal(result, expected)
    equal(a.H * (a * result - b), sp.zeros(2, 2))


@pytest.mark.parametrize("method", ["cholesky_solve", "LDLsolve", "QRsolve"])
@pytest.mark.parametrize("data", [
    [[2, 1], [1, 3]],
    [[0, 1], [1, 0]],
    [[2, sp.I], [-sp.I, 3]],
    [[1, 2], [3, 4]],
    [[sp.I, 1], [0, 2]],
])
def test_structured_square_solves(method, data):
    a = sx.Matrix(data)
    b = sx.Matrix([[1, 2], [3, 4]])
    equal(getattr(a, method)(b), a.to_sympy().inv() * b.to_sympy())


@pytest.mark.parametrize("method", ["cholesky_solve", "LDLsolve", "QRsolve"])
def test_structured_solve_errors_and_empty(method):
    with pytest.raises(NotImplementedError):
        getattr(sx.Matrix([[1, 2]]), method)(sx.Matrix([1]))
    with pytest.raises(sx.NonInvertibleMatrixError):
        getattr(sx.Matrix([[1, 2], [2, 4]]), method)(sx.Matrix([1, 2]))
    with pytest.raises(sx.ShapeError):
        getattr(sx.eye(2), method)(sx.ones(3, 1))
    assert getattr(sx.zeros(0, 0), method)(sx.zeros(0, 3)).shape == (0, 3)
    assert getattr(sx.zeros(2, 0), method)(sx.ones(2, 3)).shape == (0, 3)


@pytest.mark.parametrize("method,data", [
    ('lower_triangular_solve', [[2, 0], [3, 4]]),
    ('upper_triangular_solve', [[2, 3], [0, 4]]),
    ('diagonal_solve', [[2, 0], [0, 4]]),
])
def test_triangular_multiple_rhs(method, data):
    a = sx.Matrix(data)
    b = sx.Matrix([[1, 2], [3, 4]])
    equal(getattr(a, method)(b), a.to_sympy().inv() * b.to_sympy())


def test_structured_error_classes():
    with pytest.raises(ValueError):
        sx.Matrix([[1, 1], [0, 1]]).lower_triangular_solve(sx.ones(2, 1))
    with pytest.raises(TypeError):
        sx.Matrix([[1, 0], [1, 1]]).upper_triangular_solve(sx.ones(2, 1))
    with pytest.raises(TypeError):
        sx.Matrix([[0]]).lower_triangular_solve(sx.ones(1, 1))
    with pytest.raises(ValueError):
        sx.Matrix([[0]]).upper_triangular_solve(sx.ones(1, 1))
    with pytest.raises(sx.NonSquareMatrixError):
        sx.ones(2, 1).lower_triangular_solve(sx.ones(2, 1))
    equal(sx.eye(2).pinv(method='ED'), sp.eye(2))
    with pytest.raises(ValueError):
        sx.eye(2).pinv(method='invalid')
    with pytest.raises(NotImplementedError):
        sx.eye(2).solve_least_squares(sx.ones(2, 1), method='invalid')
    with pytest.raises(NotImplementedError):
        sx.eye(2).rank_decomposition(iszerofunc=lambda x: x == 0)


@pytest.mark.parametrize("operation", [
    lambda a: a.pinv(),
    lambda a: a.rank_decomposition()[0],
    lambda a: a.rank_decomposition()[1],
    lambda a: a.pinv_solve(sx.ones(2, 1)),
    lambda a: a.QRsolve(sx.ones(2, 1)),
    lambda a: a.cholesky_solve(sx.ones(2, 1)),
    lambda a: a.LDLsolve(sx.ones(2, 1)),
    lambda a: a.lower_triangular_solve(sx.ones(2, 1)),
    lambda a: a.upper_triangular_solve(sx.ones(2, 1)),
    lambda a: a.diagonal_solve(sx.ones(2, 1)),
    lambda a: a.solve_least_squares(sx.ones(2, 1)),
])
def test_immutable_results(operation):
    assert isinstance(operation(sx.ImmutableMatrix.eye(2)), sx.ImmutableMatrix)


@pytest.mark.parametrize('method', ['CH', 'LDL', 'QR'])
def test_structured_inverses_and_dispatch(method):
    a = sx.ImmutableMatrix([[1, 2], [3, 4]])
    expected = a.to_sympy().inv()
    equal(a.inv(method=method), expected)
    equal(getattr(a, 'inverse_' + method)(), expected)
    assert isinstance(a.inv(method=method), sx.ImmutableMatrix)
    equal(a.solve(sx.eye(2), method=method), expected)
    with pytest.raises(sx.NonSquareMatrixError):
        sx.ones(3, 2).inv(method=method)


@pytest.mark.parametrize('method', ['QRsolve', 'LDLsolve', 'cholesky_solve'])
def test_symbolic_complex_solve(method):
    z = sx.Symbol('z')
    a = sx.Matrix([[z, 1], [1, 0], [0, 1]])
    b = sx.Matrix([1, 2, 3])
    solution = getattr(a, method)(b)
    ref = a.to_sympy()
    equal(solution, (ref.H * ref).inv() * ref.H * b.to_sympy())
    equal(a.H * (a * solution - b), sp.zeros(2, 1))
