"""Regression coverage for completed matrix APIs and algebraic invariants."""

import pickle
import random

import pytest
import sympy as sp
import symox as sx


def assert_matrix_equal(actual, expected):
    actual = actual.to_sympy()
    assert actual.shape == expected.shape
    assert all(sp.simplify(a - b) == 0 for a, b in zip(actual, expected))


@pytest.mark.parametrize("mutate", [
    lambda m: m.__setitem__((0, 0), 9),
    lambda m: m.__setitem__((slice(None), 0), sx.Matrix([9, 8])),
    lambda m: m.row_del(0),
    lambda m: m.col_del(0),
    lambda m: m.row_swap(0, 1),
    lambda m: m.col_swap(0, 1),
])
def test_immutable_rejects_mutation(mutate):
    matrix = sx.ImmutableMatrix([[1, 2], [3, 4]])
    before = hash(matrix)
    with pytest.raises(TypeError):
        mutate(matrix)
    assert hash(matrix) == before
    assert_matrix_equal(matrix, sp.Matrix([[1, 2], [3, 4]]))


@pytest.mark.parametrize("operation", [
    lambda m: m.T,
    lambda m: m.H,
    lambda m: m[0, :],
    lambda m: m.reshape(1, 4),
    lambda m: m.copy(),
    lambda m: m.inv(),
    lambda m: m**2,
    lambda m: m + sx.eye(2),
    lambda m: sx.eye(2) + m,
    lambda m: m * sx.eye(2),
    lambda m: sx.eye(2) * m,
    lambda m: 3 * m,
    lambda m: sx.Symbol("x") * m,
    lambda m: m / 2,
    lambda m: m.multiply_elementwise(sx.ones(2)),
    lambda m: sx.ones(2).multiply_elementwise(m),
    lambda m: sx.Matrix.hstack(m, m),
    lambda m: sx.Matrix.vstack(m, m),
    lambda m: m.applyfunc(lambda x: x**2),
    lambda m: m.simplify(),
    lambda m: m.expand(),
    lambda m: m.rref()[0],
    lambda m: m.LUdecomposition()[0],
    lambda m: m.QRdecomposition()[0],
    lambda m: m.vec(),
    lambda m: pickle.loads(pickle.dumps(m)),
])
def test_immutable_propagates_to_results(operation):
    result = operation(sx.ImmutableMatrix([[1, 2], [3, 4]]))
    assert isinstance(result, sx.ImmutableMatrix)
    assert isinstance(hash(result), int)
    with pytest.raises(TypeError):
        result[0] = 0


def test_immutable_constructors_and_mutable_copy():
    for matrix in [sx.ImmutableMatrix.eye(2), sx.ImmutableMatrix.zeros(2),
                   sx.ImmutableMatrix.ones(2), sx.ImmutableMatrix.diag(1, 2)]:
        assert isinstance(matrix, sx.ImmutableMatrix)
    original = sx.ImmutableMatrix([[1]])
    assert hash(original) == hash(sx.ImmutableMatrix([[1]]))
    assert {original: "value"}[sx.ImmutableMatrix([[1]])] == "value"
    for mutable in [original.as_mutable(), sx.Matrix(original)]:
        with pytest.raises(TypeError):
            hash(mutable)
        mutable[0] = 2
        assert original[0] == 1
    assert isinstance(original.to_sympy(), sp.ImmutableMatrix)


@pytest.mark.parametrize("immutable", [False, True])
def test_matrix_pickle_preserves_entry_semantics(immutable):
    x, y = sx.symbols("x y", real=True, nonzero=True)
    held = sx.Mul(2, x + y, evaluate=False)
    matrix = sx.Matrix([[x, held, sx.UnevaluatedExpr(x + y)]])
    if immutable:
        matrix = matrix.as_immutable()
    restored = pickle.loads(pickle.dumps(matrix))
    assert restored == matrix
    assert restored[0, 0].is_real is True
    assert restored[0, 0].is_nonzero is True
    assert sx.srepr(restored[0, 1]) == sx.srepr(held)
    assert isinstance(restored, sx.ImmutableMatrix) is immutable


@pytest.mark.parametrize("convention", [None, "math", "maths", "left", "physics", "right"])
def test_hermitian_dot(convention):
    native_left, native_right = sx.Matrix([sx.I, 2]), sx.Matrix([3, 4 * sx.I])
    ref_left, ref_right = sp.Matrix([sp.I, 2]), sp.Matrix([3, 4 * sp.I])
    options = {"hermitian": True, "conjugate_convention": convention}
    assert native_left.dot(native_right, **options).to_sympy() == ref_left.dot(ref_right, **options)
    if convention:
        options.pop("hermitian")
        assert native_left.dot(native_right, **options).to_sympy() == ref_left.dot(ref_right, **options)


def test_dot_validates_convention():
    vector = sx.Matrix([sx.I, 1])
    with pytest.raises(ValueError):
        vector.dot(vector, hermitian=True, conjugate_convention="unknown")
    assert vector.dot(vector, hermitian=False, conjugate_convention="unknown") == 0


@pytest.mark.parametrize("entries", [
    [[-1, 0], [0, 1]],
    [[1, 2], [2, 1]],
    [[4, 2], [2, 3]],
])
def test_symmetric_decompositions_allow_indefinite(entries):
    matrix = sx.Matrix(entries)
    l = matrix.cholesky(hermitian=False)
    assert_matrix_equal(l * l.T, sp.Matrix(entries))
    l, d = matrix.LDLdecomposition(hermitian=False)
    assert_matrix_equal(l * d * l.T, sp.Matrix(entries))
    if sp.Matrix(entries).is_positive_definite is False:
        with pytest.raises(sx.NonPositiveDefiniteMatrixError):
            matrix.cholesky()


def test_diagonalize_options():
    matrix = sx.Matrix([[2, 1], [1, 2]])
    p, d = matrix.diagonalize(normalize=True, sort=True)
    assert_matrix_equal(matrix * p, (p * d).to_sympy())
    assert_matrix_equal(p.H * p, sp.eye(2))
    rotation = sx.Matrix([[0, 1], [-1, 0]])
    assert rotation.is_diagonalizable()
    assert rotation.is_diagonalizable(reals_only=True) is False
    with pytest.raises(sx.MatrixError):
        rotation.diagonalize(reals_only=True)
    assert sx.Matrix([[1, 2]]).is_diagonalizable() is False
    assert sx.Matrix([[1, 1], [0, 1]]).is_diagonalizable() is False


def test_unsupported_eigenproblem_does_not_return_false():
    # Companion matrix for t**5 - t - 1, with no supported radical solution.
    matrix = sx.Matrix([[0, 0, 0, 0, 1], [1, 0, 0, 0, 1],
                        [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]])
    with pytest.raises(NotImplementedError):
        matrix.is_diagonalizable()


def test_gauss_jordan_free_variables_and_rankcheck():
    matrix = sx.Matrix([[0, 1, 2], [0, 2, 4]])
    rhs = sx.Matrix([3, 6])
    solution, parameters, free_columns = matrix.gauss_jordan_solve(rhs, freevar=True)
    assert free_columns == [0, 2]
    assert parameters.shape == (2, 1)
    assert_matrix_equal(matrix * solution, sp.Matrix([3, 6]))
    with pytest.raises(ValueError):
        matrix.LUdecomposition(rankcheck=True)
    l, u, swaps = matrix.LUdecomposition(rankcheck=False)
    assert swaps == []
    assert_matrix_equal(l * u, matrix.to_sympy())


def test_multiple_eigenvalues_and_polynomial_equality():
    assert sx.diag(1, 1, 2).eigenvals(multiple=True) == [1, 1, 2]
    x, y = sx.symbols("x y")
    p = sx.eye(2).charpoly(x)
    q = sx.eye(2).charpoly(y)
    assert p == q and hash(p) == hash(q)
    reference=sp.eye(2).charpoly()
    assert (p == p.as_expr()) == (reference == reference.as_expr())


def test_norm_preserves_symbol_assumptions():
    x = sx.Symbol("x", real=True)
    p = sx.Symbol("p", positive=True)
    result = sx.Matrix([x, 2]).norm(p)
    assert p in result.free_symbols
    assert next(s for s in result.free_symbols if s.name == "p").is_positive is True
    assert sx.Matrix([[1, 2], [3, 4]]).norm("FRO") == sx.sqrt(30)
    with pytest.raises(ValueError):
        sx.Matrix([1, 2]).norm("fro")


@pytest.mark.parametrize("operation, exception", [
    (lambda m: m.eigenvals(chop=True), NotImplementedError),
    (lambda m: m.solve(sx.eye(2), method="unknown"), ValueError),
])
def test_unsupported_options_are_explicit(operation, exception):
    with pytest.raises(exception):
        operation(sx.Matrix([[1, 2], [3, 4]]))
    if exception is ValueError:
        with pytest.raises(ValueError, match="Inversion method unrecognized"):
            sp.Matrix([[1, 2], [3, 4]]).solve(sp.eye(2), method="unknown")


def test_matrix_evalf_accepts_maxn_and_strict():
    source = [[1, 2], [3, 4]]
    actual = sx.Matrix(source).evalf(50, maxn=200, strict=True)
    expected = sp.Matrix(source).evalf(50, maxn=200, strict=True)
    assert [str(entry) for entry in actual] == [str(entry) for entry in expected]
    cancelling = sx.Add(sx.sin(1) ** 2, sx.cos(1) ** 2, -1, evaluate=False)
    with pytest.raises(sx.PrecisionExhausted):
        sx.Matrix([[cancelling, 1]]).evalf(15, strict=True)


def test_random_rational_linalg_matches_sympy():
    rng = random.Random(2749)
    for rows, cols in [(2, 3), (3, 2), (3, 3), (4, 4)]:
        for _ in range(5):
            entries = [[sp.Rational(rng.randrange(-5, 6), rng.randrange(1, 5))
                        for _ in range(cols)] for _ in range(rows)]
            expected = sp.Matrix(entries)
            matrix = sx.Matrix(entries)
            actual_rref, actual_pivots = matrix.rref()
            ref_rref, ref_pivots = expected.rref()
            assert_matrix_equal(actual_rref, ref_rref)
            assert actual_pivots == ref_pivots
            if rows == cols:
                assert matrix.det().to_sympy() == expected.det()
                if expected.det():
                    assert_matrix_equal(matrix.inv(), expected.inv())


@pytest.mark.parametrize("entries", [
    [[1, 2], [3, 4]],
    [[1, 1, 1], [2, 2, 2]],
    [[0, 1], [0, 2], [0, 3]],
    [[sp.I, 1], [1, sp.I]],
    [[0, 0], [0, 0]],
])
def test_qr_decomposition_invariants(entries):
    matrix = sx.Matrix(entries)
    q, r = matrix.QRdecomposition()
    assert_matrix_equal(q * r, sp.Matrix(entries))
    assert_matrix_equal(q.H * q, sp.eye(q.cols))


@pytest.mark.parametrize("build", [
    lambda m: [m.Matrix([1, 2]), m.Matrix([3, 4])],
    lambda m: [[m.Matrix([1, 2]), m.Matrix([3, 4])]],
    lambda m: [m.Matrix([1, 2]), 3],
    lambda m: [[m.eye(2), m.ones(2)], [m.zeros(2), m.eye(2)]],
])
def test_constructor_stacks_flat_matrices_and_joins_block_rows(build):
    assert_matrix_equal(sx.Matrix(build(sx)), sp.Matrix(build(sp)))


def test_matrix_power_rejects_float_integer_overflow():
    with pytest.raises(ValueError, match="exponent too large"):
        sx.Matrix([[-1]]) ** 1e100
    assert sx.Matrix([[2]]) ** 2.0 == sx.Matrix([[4]])


def test_zero_matrix_predicate_preserves_unknown_assumptions():
    assert sx.Matrix([[sx.Symbol("x")]]).is_zero_matrix is None
    assert sx.Matrix([[sx.Symbol("x", nonzero=True)]]).is_zero_matrix is False
    assert sx.Matrix([[sx.Symbol("x", zero=True)]]).is_zero_matrix is True
    assert sx.zeros(2).is_zero_matrix is True
    assert sx.Matrix([[sx.Symbol("x"), 1]]).is_zero_matrix is False
