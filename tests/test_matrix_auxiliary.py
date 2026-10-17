"""Auxiliary MatrixBase methods and top-level matrix helpers compared with SymPy."""
import pickle

import pytest
import sympy as sp
import symox as sx
from symox import _sparse

KINDS = ['dense', 'immutable', 'sparse', 'isparse']
SX = {'dense': sx.Matrix, 'immutable': sx.ImmutableMatrix,
      'sparse': sx.SparseMatrix, 'isparse': sx.ImmutableSparseMatrix}
SP = {'dense': sp.Matrix, 'immutable': sp.ImmutableMatrix,
      'sparse': sp.SparseMatrix, 'isparse': sp.ImmutableSparseMatrix}


def kind_of(value):
    if isinstance(value, sp.MatrixBase):
        sparse = isinstance(value, sp.SparseMatrix.__mro__[1])
        immutable = isinstance(value, (sp.ImmutableDenseMatrix, sp.ImmutableSparseMatrix))
    else:
        sparse = isinstance(value, _sparse._SparseBase)
        immutable = bool(getattr(value, '_immutable', False))
    return ('isparse' if immutable else 'sparse') if sparse else ('immutable' if immutable else 'dense')


def is_matrix(value):
    return isinstance(value, (sp.MatrixBase, sx.MatrixBase))


def norm(value):
    """A library independent, exact rendering of a result."""
    if is_matrix(value):
        return ('M', kind_of(value), value.shape,
                tuple(str(value[i, j]) for i in range(value.rows) for j in range(value.cols)))
    if isinstance(value, (tuple, list)):
        return (type(value).__name__, tuple(norm(v) for v in value))
    if isinstance(value, dict):
        return ('dict', tuple(sorted((str(k), norm(v)) for k, v in value.items())))
    if hasattr(value, 'as_expr') and hasattr(value, 'gens'):
        return ('poly', str(value.as_expr()))
    if isinstance(value, (set, frozenset)):
        return ('set', tuple(sorted(str(v) for v in value)))
    if value is None or isinstance(value, bool):
        return value
    if isinstance(value, int):
        return ('int', value)
    return str(value)


def outcome(fn, lib):
    try:
        return ('ok', norm(fn(lib)))
    except Exception as exc:  # noqa: BLE001 - the exception class is the result
        return ('err', type(exc).__name__)


def same(fn):
    expected = outcome(fn, sp)
    actual = outcome(fn, sx)
    assert actual == expected
    return actual


def build(lib, kind, rows):
    cls = (SX if lib is sx else SP)[kind]
    return cls(rows)


def mk(kind, rows):
    return lambda lib: build(lib, kind, rows)


A3 = [[1, 2, 0], [3, 4, 5], [0, 6, 7]]
SING = [[1, 2, 3], [2, 4, 6], [1, 0, 1]]
RECT = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 12, 11]]


def sym(lib):
    a, b, c, d = lib.symbols('a b c d')
    return [[a, b], [c, d]]


# ---------------------------------------------------------------------------
# Plain kernels on a family of inputs
# ---------------------------------------------------------------------------

UNARY = [
    lambda M: M.berkowitz(),
    lambda M: M.berkowitz_minors(),
    lambda M: M.berkowitz_det(),
    lambda M: M.berkowitz_charpoly(),
    lambda M: M.det_bareis(),
    lambda M: M.det_LU_decomposition(),
    lambda M: M.per(),
    lambda M: M.lower_triangular(),
    lambda M: M.upper_triangular(),
    lambda M: M.lower_triangular(1),
    lambda M: M.upper_triangular(-1),
    lambda M: M.rot90(),
    lambda M: M.rot90(2),
    lambda M: M.rot90(-1),
    lambda M: M.rot90(7),
    lambda M: M.is_nilpotent(),
    lambda M: M.is_symbolic(),
    lambda M: M.get_diag_blocks(),
    lambda M: M.connected_components(),
    lambda M: M.strongly_connected_components(),
    lambda M: M.LUdecomposition_Simple(),
    lambda M: M.LUdecompositionFF(),
    lambda M: M.inverse_BLOCK(),
    lambda M: M.vech(),
    lambda M: M.vech(diagonal=False, check_symmetry=False),
    lambda M: M.cofactorMatrix(),
    lambda M: M.minorEntry(0, 1),
    lambda M: M.minorMatrix(0, 0),
    lambda M: M.is_anti_symmetric(),
    lambda M: M.is_anti_symmetric(simplify=False),
    lambda M: M.dual(),
    lambda M: M.flat(),
    lambda M: M.diagonal(),
    lambda M: M.diagonal(1),
    lambda M: M.diagonal(-1),
    lambda M: M.has(2),
    lambda M: M.todok(),
    lambda M: M.todod(),
    lambda M: list(M.iter_values()),
    lambda M: list(M.iter_items()),
    lambda M: M.atoms(),
    lambda M: M.as_real_imag(),
    lambda M: M.doit(),
]

INPUTS = {
    'A3': lambda lib: A3,
    'singular': lambda lib: SING,
    'rect': lambda lib: RECT,
    'symbolic': sym,
    'one': lambda lib: [[lib.Integer(5)]],
    'nilpotent': lambda lib: [[0, 1, 2], [0, 0, 3], [0, 0, 0]],
    'blocks': lambda lib: [[1, 2, 0, 0], [3, 4, 0, 0], [0, 0, 5, 0], [0, 0, 0, 6]],
    'antisym': lambda lib: [[0, 1, -2], [-1, 0, 3], [2, -3, 0]],
    'symmetric': lambda lib: [[1, 2, 3], [2, 4, 5], [3, 5, 6]],
    'rational': lambda lib: [[lib.Rational(1, 2), 3], [lib.Rational(-2, 3), 7]],
}


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('name', sorted(INPUTS))
@pytest.mark.parametrize('index', range(len(UNARY)))
def test_unary_kernels_match_sympy(kind, name, index):
    op = UNARY[index]
    same(lambda lib: op(build(lib, kind, INPUTS[name](lib))))


@pytest.mark.parametrize('kind', KINDS)
def test_hat_vee_and_dual(kind):
    same(lambda lib: build(lib, kind, [[1], [2], [3]]).hat())
    same(lambda lib: build(lib, kind, [[1, 2, 3]]).hat())
    same(lambda lib: build(lib, kind, [[1], [2]]).hat())
    same(lambda lib: build(lib, kind, [[0, -3, 2], [3, 0, -1], [-2, 1, 0]]).vee())
    same(lambda lib: build(lib, kind, [[0, 3, 2], [3, 0, -1], [-2, 1, 0]]).vee())
    same(lambda lib: build(lib, kind, [[1, 2], [3, 4]]).vee())
    same(lambda lib: build(lib, kind, sym(lib)).dual())
    x = lambda lib: lib.Symbol('x')  # noqa: E731
    same(lambda lib: build(lib, kind, [[0, x(lib)], [-x(lib), 0]]).is_anti_symmetric())
    same(lambda lib: build(lib, kind, [[0, x(lib)**2 - 1], [1 - (x(lib) - 1)*(x(lib) + 1) - 1, 0]])
         .is_anti_symmetric())


@pytest.mark.parametrize('kind', KINDS)
def test_lu_decompositions(kind):
    for rows in (A3, SING, RECT, [[0, 1], [1, 0]], [[2, 1], [4, 2]]):
        same(lambda lib: build(lib, kind, rows).LUdecomposition_Simple())
        same(lambda lib: build(lib, kind, rows).LUdecomposition_Simple(rankcheck=True))
        same(lambda lib: build(lib, kind, rows).LUdecompositionFF())
    P, L, DD, U = sx.SparseMatrix(A3).LUdecompositionFF()
    assert P * build(sx, kind, A3) == L * DD.inv() * U


@pytest.mark.parametrize('kind', KINDS)
def test_inverse_block_and_cramer(kind):
    big = [[(i + 1)*(j + 2) % 7 + (i == j)*9 for j in range(6)] for i in range(6)]
    same(lambda lib: build(lib, kind, big).inverse_BLOCK())
    same(lambda lib: build(lib, kind, big).inverse_BLOCK(iszerofunc=lambda e: e == 0))
    same(lambda lib: build(lib, kind, SING).inverse_BLOCK())
    rhs = lambda lib: lib.Matrix([1, 2, 3])  # noqa: E731
    same(lambda lib: build(lib, kind, A3).cramer_solve(rhs(lib)))
    same(lambda lib: build(lib, kind, A3).cramer_solve(rhs(lib), det_method='bareiss'))
    same(lambda lib: build(lib, kind, A3).cramer_solve(rhs(lib), det_method='berkowitz'))
    same(lambda lib: build(lib, kind, A3).cramer_solve(rhs(lib), det_method=lambda m: m.det()))
    same(lambda lib: build(lib, kind, A3).cramer_solve(lib.Matrix([[1, 2], [3, 4], [5, 6]])))
    same(lambda lib: build(lib, kind, sym(lib)).cramer_solve(lib.Matrix(lib.symbols('u v'))))
    same(lambda lib: build(lib, kind, A3).rref_rhs(rhs(lib)))
    same(lambda lib: build(lib, kind, SING).rref_rhs(rhs(lib)))


@pytest.mark.parametrize('kind', KINDS)
def test_normalized_project_and_equals(kind):
    x = sx.Symbol('x')
    same(lambda lib: build(lib, kind, [[3], [4]]).normalized())
    same(lambda lib: build(lib, kind, [[0], [0]]).normalized())
    same(lambda lib: build(lib, kind, [[1, 2]]).project(lib.Matrix([[3, 4]])))
    same(lambda lib: build(lib, kind, [[1, 2]]).equals(lib.Matrix([[1, 2]])))
    same(lambda lib: build(lib, kind, [[1, 2]]).equals(lib.Matrix([[1, 3]])))
    same(lambda lib: build(lib, kind, [[1, 2]]).equals(lib.Matrix([[1, 2, 3]])))
    same(lambda lib: build(lib, kind, [[1, 2]]).equals(5))
    sxx = lambda lib: lib.Symbol('x')  # noqa: E731
    same(lambda lib: build(lib, kind, [[(sxx(lib) + 1)**2]]).equals(
        lib.Matrix([[sxx(lib)**2 + 2*sxx(lib) + 1]])))
    same(lambda lib: build(lib, kind, [[lib.sin(sxx(lib))**2 + lib.cos(sxx(lib))**2]]).equals(
        lib.Matrix([[1]])))
    same(lambda lib: build(lib, kind, [[sxx(lib)]]).equals(lib.Matrix([[sxx(lib) + 1]])))
    assert build(sx, kind, [[x]]).equals(sx.Matrix([[x]])) is True


@pytest.mark.parametrize('kind', KINDS)
def test_hessenberg_and_bidiagonal_reconstruct(kind):
    def numeric(matrix):
        return matrix.evalf(40).to_sympy() if hasattr(matrix, '_immutable') else matrix.evalf(40)
    for rows in (A3, [[1, 2, 3, 4], [5, 6, 7, 8], [9, 1, 2, 3], [4, 5, 6, 0]]):
        M = build(sx, kind, rows)
        H, P = M.upper_hessenberg_decomposition()
        # Evaluate the independent oracle before decomposition: exact nested
        # radicals trigger expensive adaptive cancellation in SymPy's evalf.
        Hs, Ps = sp.Matrix(rows).evalf(50).upper_hessenberg_decomposition()
        assert H.shape == Hs.shape and P.shape == Ps.shape
        assert (numeric(H) - numeric(Hs)).norm() < 1e-20
        assert (numeric(P) - numeric(Ps)).norm() < 1e-20
        assert (numeric(P)*numeric(H)*numeric(P).H - sp.Matrix(rows)).norm() < 1e-20
        assert all(H[i, j] == 0 for i in range(H.rows) for j in range(i - 1))
    for rows in (A3, RECT):
        for upper in (True, False):
            M = build(sx, kind, rows)
            U, B, V = M.bidiagonal_decomposition(upper=upper)
            Us, Bs, Vs = sp.Matrix(rows).evalf(50).bidiagonal_decomposition(upper=upper)
            assert B.shape == Bs.shape and U.shape == Us.shape and V.shape == Vs.shape
            assert (numeric(B).applyfunc(abs) - numeric(Bs).applyfunc(abs)).norm() < 1e-20
            assert (numeric(U)*numeric(B)*numeric(V) - sp.Matrix(rows)).norm() < 1e-20
            assert (numeric(M.bidiagonalize(upper=upper)) - numeric(B)).norm() < 1e-20
    same(lambda lib: build(lib, kind, A3).bidiagonalize(upper=1))
    same(lambda lib: build(lib, kind, RECT).upper_hessenberg_decomposition())


# ---------------------------------------------------------------------------
# Container protocol and entrywise helpers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('kind', KINDS)
def test_key_protocol(kind):
    for key in [(0, 0), (-1, -1), (slice(None), 1), (0, slice(1, 3)), 5, -1, (3, 0), (0, 4), 12, 'a',
                (slice(0, 2), slice(1, None))]:
        same(lambda lib: build(lib, kind, RECT).key2ij(key))
    for key in [(0, 0), (slice(None), 1), (0, slice(1, 3)), (slice(0, 2), slice(1, None)), (4, 0),
                (0, [1, 2]), ([0, 2], [1, 3])]:
        same(lambda lib: build(lib, kind, RECT).key2bounds(key))
    same(lambda lib: build(lib, kind, RECT).table(lib.StrPrinter()))
    same(lambda lib: build(lib, kind, RECT).table(lib.StrPrinter(), rowstart='<', rowend='>',
                                                  rowsep=';', colsep='|', align='left'))
    same(lambda lib: build(lib, kind, RECT).table(lib.StrPrinter(), align='center'))
    same(lambda lib: build(lib, kind, RECT).table(lib.StrPrinter(), align='bogus'))
    same(lambda lib: build(lib, kind, []).table(lib.StrPrinter()))


@pytest.mark.parametrize('kind', KINDS)
def test_print_nonzero(kind, capsys):
    expected = []
    for lib in (sp, sx):
        build(lib, kind, [[1, 0], [0, 3]]).print_nonzero()
        build(lib, kind, [[1, 0], [0, 3]]).print_nonzero('#')
        expected.append(capsys.readouterr().out)
    assert expected[0] == expected[1]


@pytest.mark.parametrize('kind', ['dense', 'sparse'])
def test_mutators(kind):
    def run(op):
        def fn(lib):
            M = build(lib, kind, [[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            op(M, lib)
            return M
        same(fn)
    run(lambda M, lib: M.copyin_matrix((slice(0, 2), slice(1, 3)), lib.Matrix([[10, 11], [12, 13]])))
    run(lambda M, lib: M.copyin_matrix((slice(0, 2), slice(1, 3)), lib.Matrix([[10, 11]])))
    run(lambda M, lib: M.copyin_list((slice(0, 1), slice(None)), [7, 7, 7]))
    run(lambda M, lib: M.copyin_list((slice(0, 1), slice(None)), 7))
    run(lambda M, lib: M.fill(lib.Symbol('z')))
    run(lambda M, lib: M.fill(0))
    run(lambda M, lib: M.row_add(0, 2, 3))
    run(lambda M, lib: M.row_mult(1, lib.Rational(1, 2)))
    run(lambda M, lib: M.row_mult(1, 0))


@pytest.mark.parametrize('kind', ['immutable', 'isparse'])
def test_immutable_has_no_mutators(kind):
    M = build(sx, kind, [[1]])
    for name in ('fill', 'row_add', 'row_mult', 'copyin_list', 'copyin_matrix'):
        assert hasattr(M, name) == hasattr(build(sp, kind, [[1]]), name)


@pytest.mark.parametrize('kind', KINDS)
def test_entrywise_helpers(kind):
    def m(lib):
        x, y = lib.symbols('x y')
        return build(lib, kind, [[x**2, lib.sin(x)], [lib.Integral(x, x), y + lib.I*x]])
    same(lambda lib: m(lib).atoms())
    same(lambda lib: m(lib).atoms(lib.Symbol))
    same(lambda lib: m(lib).atoms(lib.Number))
    same(lambda lib: build(lib, kind, [[1, 0]]).atoms(lib.Integer))
    same(lambda lib: m(lib).doit())
    same(lambda lib: m(lib).integrate(lib.Symbol('x')))
    same(lambda lib: m(lib).integrate((lib.Symbol('x'), 0, 1)))
    same(lambda lib: m(lib).limit(lib.Symbol('x'), 0))
    same(lambda lib: m(lib).replace(lib.sin, lib.cos))
    actual_parts = m(sx).as_real_imag()
    expected_parts = m(sp).as_real_imag()
    for actual, expected in zip(actual_parts, expected_parts):
        assert kind_of(actual) == kind_of(expected)
        assert actual.shape == expected.shape
        assert all(sp.simplify(sp.expand_complex(actual[i,j].to_sympy() - expected[i,j])) == 0
                   for i in range(actual.rows) for j in range(actual.cols))
    same(lambda lib: m(lib).subs(lib.Symbol('y'), 2).add(m(lib)))
    same(lambda lib: m(lib).rmultiply(build(lib, kind, [[1, 2]])))
    same(lambda lib: m(lib).add(build(lib, kind, [[1, 2]])))


@pytest.mark.parametrize('kind', KINDS)
def test_pow_methods(kind):
    for method in (None, 'multiply', 'mulsimp', 'jordan', 'cayley'):
        same(lambda lib: build(lib, kind, [[1, 1], [0, 1]]).pow(5, method=method))
        same(lambda lib: build(lib, kind, [[2, 0], [0, 3]]).pow(-2, method=method))
    same(lambda lib: build(lib, kind, [[1, 1], [0, 1]]).pow(2, method='bogus'))
    same(lambda lib: build(lib, kind, [[4, 0], [0, 9]]).pow(lib.Rational(1, 2)))
    same(lambda lib: build(lib, kind, [[1, 2]]).pow(2))
    same(lambda lib: build(lib, kind, [[1, 1], [0, 1]]).pow(lib.Symbol('n')))


# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('kind', KINDS)
def test_class_constructors(kind):
    cls = lambda lib: (SX if lib is sx else SP)[kind]  # noqa: E731
    same(lambda lib: cls(lib).jordan_block(3, 2))
    same(lambda lib: cls(lib).jordan_block(size=3, eigenvalue=lib.Symbol('l'), band='lower'))
    same(lambda lib: cls(lib).jordan_block(eigenvalue=2, rows=2, cols=3))
    same(lambda lib: cls(lib).jordan_block(3))
    same(lambda lib: cls(lib).jordan_block(3, 2, band='middle'))
    same(lambda lib: cls(lib).wilkinson(1))
    same(lambda lib: cls(lib).wilkinson(2))
    same(lambda lib: cls(lib).companion(lib.Poly(lib.Symbol('x')**3 - 2*lib.Symbol('x') + 5)))
    same(lambda lib: cls(lib).companion(lib.Poly(2*lib.Symbol('x')**2 + 1)))
    same(lambda lib: cls(lib).companion(lib.Poly(lib.Symbol('x') + 3)))
    same(lambda lib: cls(lib).companion(lib.Poly(lib.Symbol('x')*lib.Symbol('y'))))
    same(lambda lib: cls(lib).irregular(3, lib.ones(2, 1), lib.Matrix([[2, 3]]), lib.eye(2) * 4,
                                        lib.Matrix([[5, 6, 7]])))
    same(lambda lib: cls(lib).irregular(2, lib.ones(2, 2), lib.ones(1, 3)))
    same(lambda lib: cls(lib).orthogonalize(lib.Matrix([1, 2]), lib.Matrix([3, 4])))
    same(lambda lib: cls(lib).orthogonalize(lib.Matrix([1, 2]), lib.Matrix([3, 4]), normalize=True))
    same(lambda lib: cls(lib).orthogonalize(lib.Matrix([[1, 2]]), lib.Matrix([[2, 4]])))
    same(lambda lib: cls(lib).orthogonalize(lib.Matrix([1, 2]), lib.Matrix([2, 4]), rankcheck=True))
    same(lambda lib: cls(lib).orthogonalize())
    same(lambda lib: cls(lib).from_dok(2, 3, {(0, 1): 5, (1, 2): lib.Symbol('t')}))


def test_top_level_functions():
    x, n = sp.symbols('x n')
    same(lambda lib: lib.GramSchmidt([lib.Matrix([1, 1, 0]), lib.Matrix([1, 0, 1]),
                                      lib.Matrix([0, 1, 1])]))
    same(lambda lib: lib.GramSchmidt([lib.Matrix([1, 1, 0]), lib.Matrix([1, 0, 1])], True))
    same(lambda lib: lib.GramSchmidt([lib.Matrix([1, 1]), lib.Matrix([2, 2])]))
    same(lambda lib: lib.casoratian([1, 2**lib.Symbol('n'), lib.factorial(lib.Symbol('n'))],
                                    lib.Symbol('n')))
    same(lambda lib: lib.casoratian([1, 2**lib.Symbol('n')], lib.Symbol('n'), zero=False))
    same(lambda lib: lib.casoratian([3], lib.Symbol('n')))
    same(lambda lib: lib.jordan_cell(4, 3))
    same(lambda lib: lib.Matrix([[1, 2]]).jordan_cell(4, 2))
    same(lambda lib: lib.matrix_multiply_elementwise(lib.Matrix([[1, 2]]), lib.Matrix([[3, 4]])))
    same(lambda lib: lib.matrix_multiply_elementwise(lib.Matrix([[1, 2]]), lib.Matrix([[3]])))
    for fn in ('rot_axis1', 'rot_axis2', 'rot_axis3', 'rot_ccw_axis1', 'rot_ccw_axis2', 'rot_ccw_axis3'):
        same(lambda lib: getattr(lib, fn)(lib.Symbol('t')))
        same(lambda lib: getattr(lib, fn)(lib.pi / 3))
        same(lambda lib: getattr(lib, fn)(0))
    same(lambda lib: lib.rot_givens(1, 3, lib.Symbol('t'), dim=4))
    same(lambda lib: lib.rot_givens(0, 1, lib.pi / 2))
    same(lambda lib: lib.rot_givens(1, 1, 1))
    same(lambda lib: lib.rot_givens(0, 3, 1))
    same(lambda lib: lib.MutableMatrix([[1, 2]]))
    assert sx.MutableMatrix is sx.Matrix


def test_numpy_helpers():
    np = pytest.importorskip('numpy')
    a = sx.list2numpy([1, sx.Symbol('x')])
    b = sp.list2numpy([1, sp.Symbol('x')])
    assert a.dtype == b.dtype and [str(v) for v in a] == [str(v) for v in b]
    a = sx.matrix2numpy(sx.Matrix([[1, 2], [3, 4]]))
    b = sp.matrix2numpy(sp.Matrix([[1, 2], [3, 4]]))
    assert a.shape == b.shape and a.dtype == b.dtype and str(a) == str(b)
    assert sx.matrix2numpy(sx.Matrix([[1, 2]]), dtype=float).dtype == np.float64
    for args, kwargs in [(('a', 3), {}), (('a', (2, 3)), {}), (('x', (2, 2)), {'integer': True}),
                         (('b', 0), {})]:
        a = sx.symarray(*args, **kwargs)
        b = sp.symarray(*args, **kwargs)
        assert a.shape == b.shape and str(a) == str(b)
    assert sx.symarray('x', 2, integer=True)[0].is_integer is True


def test_randmatrix_matches_sympy_seeded():
    for args, kwargs in [((3,), {'seed': 4}), ((2, 3), {'min': -5, 'max': 5, 'seed': 1}),
                         ((4,), {'symmetric': True, 'seed': 9}),
                         ((4,), {'percent': 50, 'seed': 3}),
                         ((5,), {'symmetric': True, 'percent': 30, 'seed': 2}),
                         ((3, 4), {'symmetric': True, 'seed': 1}),
                         ((3,), {'percent': 0, 'seed': 1})]:
        same(lambda lib: lib.randMatrix(*args, **kwargs))


def test_banded_matches_sympy():
    for args, kwargs in [((), {}), (({0: 1, 1: 2},), {}), ((4, {0: 1, -1: 2}), {}),
                         ((3, 5, {1: [1, 2, 3]}), {}), (({0: [1, 2, 3], 2: 4},), {}),
                         (({1: lambda i: i + 1},), {'rows': 3}),
                         (({0: sp.Matrix([[1, 2], [3, 4]])},), {}),
                         ((2, {0: [1, 2, 3]}), {}), ((2, {5: 1}), {}), (({0: 1},), {'rows': 3})]:
        def fn(lib):
            converted = tuple({k: (lib.Matrix(v.tolist()) if isinstance(v, sp.MatrixBase) else v)
                               for k, v in a.items()} if isinstance(a, dict) else a for a in args)
            return lib.banded(*converted, **kwargs)
        same(fn)


def test_matrixbase_isinstance():
    for kind in KINDS:
        assert isinstance(build(sx, kind, [[1]]), sx.MatrixBase)
    assert not isinstance(sx.Symbol('x'), sx.MatrixBase)
    assert not isinstance(1, sx.MatrixBase)
    assert issubclass(sx.Matrix, sx.MatrixBase)


# ---------------------------------------------------------------------------
# Errors, special values and persistence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('kind', KINDS)
def test_errors_match_sympy(kind):
    same(lambda lib: build(lib, kind, RECT).berkowitz())
    same(lambda lib: build(lib, kind, RECT).per())
    same(lambda lib: build(lib, kind, RECT).inverse_BLOCK())
    same(lambda lib: build(lib, kind, RECT).is_nilpotent())
    same(lambda lib: build(lib, kind, A3).vech())
    same(lambda lib: build(lib, kind, A3).dual())
    same(lambda lib: build(lib, kind, A3).lower_triangular(1.5))
    same(lambda lib: build(lib, kind, A3).diagonal(3))
    same(lambda lib: build(lib, kind, A3).diagonal(-3))
    same(lambda lib: build(lib, kind, A3).diagonal(1.0))
    same(lambda lib: build(lib, kind, A3).diagonal('1'))
    same(lambda lib: build(lib, kind, SING).LUdecompositionFF())
    same(lambda lib: build(lib, kind, A3).cramer_solve(lib.Matrix([1, 2, 3]), det_method='bogus'))
    same(lambda lib: build(lib, kind, A3).bidiagonal_decomposition(upper=1))
    same(lambda lib: build(lib, kind, A3).key2ij((5, 5)))
    same(lambda lib: build(lib, kind, A3).key2bounds((5, 0)))
    same(lambda lib: build(lib, kind, A3).equals(lib.Matrix([[1]])))


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('shape', [(0, 0), (0, 3), (3, 0)])
def test_empty_matrices(kind, shape):
    def m(lib):
        return (SX if lib is sx else SP)[kind](*shape, [])
    for op in [lambda M: M.berkowitz(), lambda M: M.berkowitz_det(), lambda M: M.per(),
               lambda M: M.rot90(), lambda M: M.lower_triangular(), lambda M: M.flat(),
               lambda M: M.diagonal(), lambda M: M.todok(), lambda M: M.is_symbolic(),
               lambda M: M.get_diag_blocks(), lambda M: M.connected_components(),
               lambda M: M.LUdecomposition_Simple(), lambda M: M.inverse_BLOCK(),
               lambda M: M.is_nilpotent(), lambda M: M.is_anti_symmetric(), lambda M: M.atoms()]:
        same(lambda lib: op(m(lib)))


@pytest.mark.parametrize('kind', KINDS)
def test_special_values(kind):
    for value in ('oo', 'nan', 'zoo'):
        rows = lambda lib: [[getattr(lib, value), 1], [0, 2]]  # noqa: E731
        same(lambda lib: build(lib, kind, rows(lib)).berkowitz_det())
        same(lambda lib: build(lib, kind, rows(lib)).per())
        same(lambda lib: build(lib, kind, rows(lib)).rot90())
        same(lambda lib: build(lib, kind, rows(lib)).is_symbolic())
    floats = lambda lib: [[lib.Float('1.5', 30), 2], [3, lib.Float('0.25', 10)]]  # noqa: E731
    same(lambda lib: build(lib, kind, floats(lib)).berkowitz())
    same(lambda lib: build(lib, kind, floats(lib)).per())
    same(lambda lib: build(lib, kind, floats(lib)).LUdecompositionFF())
    same(lambda lib: build(lib, kind, floats(lib)).n(5))


@pytest.mark.parametrize('kind', KINDS)
def test_results_pickle_and_print(kind):
    M = build(sx, kind, A3)
    E = build(sp, kind, A3)
    for name in ('rot90', 'lower_triangular', 'upper_triangular', 'inverse_BLOCK'):
        result = getattr(M, name)()
        clone = pickle.loads(pickle.dumps(result))
        assert clone == result and kind_of(clone) == kind_of(result)
        assert str(result) == str(getattr(E, name)())
    assert sx.srepr(M.per()) == sp.srepr(E.per())
    assert sx.srepr(M.berkowitz_det()) == sp.srepr(E.berkowitz_det())


def test_sparse_only_helpers():
    rows = [[1, 0, 3, 2], [0, 0, 1, 0], [4, 0, 0, 5], [0, 6, 7, 0]]
    for kind in ('sparse', 'isparse'):
        same(lambda lib: build(lib, kind, rows).liupc())
        same(lambda lib: build(lib, kind, rows).row_structure_symbolic_cholesky())
        same(lambda lib: build(lib, kind, rows).scalar_multiply(lib.Symbol('k')))
        same(lambda lib: build(lib, kind, rows).multiply(build(lib, kind, rows)))
        same(lambda lib: build(lib, kind, rows).col_insert(1, build(lib, kind, [[9]] * 4)))
        same(lambda lib: build(lib, kind, rows).row_insert(-1, build(lib, kind, [[9] * 4])))
    for kind in KINDS:
        same(lambda lib: build(lib, kind, rows).diagonal(2))


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('offset',[-2**63,2**63-1])
def test_extreme_diagonal_offsets_do_not_overflow(kind,offset):
    original=build(sx,kind,A3)
    reference=build(sp,kind,A3)
    check=lambda actual,expected: norm(actual)==norm(expected)
    assert check(original.lower_triangular(offset),reference.lower_triangular(offset))
    assert check(original.upper_triangular(offset),reference.upper_triangular(offset))
    with pytest.raises(ValueError): original.diagonal(offset)


def test_irregular_zero_top_count_and_constructor_limits():
    code='''
import symox as s
for cls in (s.Matrix,s.ImmutableMatrix,s.SparseMatrix,s.ImmutableSparseMatrix):
    try: cls.irregular(0,s.Matrix([[1]]))
    except ValueError: pass
    else: raise AssertionError('expected a nonempty top row')
    assert cls.irregular(0).shape==(0,0)
    for build in (lambda:cls.wilkinson(2**63),lambda:cls.jordan_block(2**63,1)):
        try: build()
        except NotImplementedError: pass
        else: raise AssertionError('expected bounded allocation')
'''
    import os,subprocess,sys
    subprocess.run([sys.executable,'-c',code],env=os.environ.copy(),check=True,timeout=15)


def test_shared_symbol_free_dags_are_inspected_without_expansion():
    import os,subprocess,sys
    code='''
import symox as s
one=s.Integer(1)
e=s.Function('f')(one)
for _ in range(80):
    e=s.Add(e,e,evaluate=False)
assert s.Matrix([[e]]).is_symbolic() is False
x=s.Symbol('x')
assert s.Matrix([[s.Integral(x,x)]]).is_symbolic() is True
'''
    subprocess.run([sys.executable,'-c',code],env=os.environ.copy(),check=True,timeout=15)


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('rational',[False,True])
def test_rref_rhs_result_classes(kind,rational):
    same(lambda lib: build(lib,kind,[[2 if rational else lib.Symbol('x'),1],[0,3]]).rref_rhs(lib.Matrix([1,2])))


@pytest.mark.parametrize('cls',[sx.SparseMatrix,sx.ImmutableSparseMatrix])
def test_sparse_from_dok_preserves_large_shapes(cls,monkeypatch):
    def fail(*args): raise AssertionError('must not densify')
    monkeypatch.setattr(cls,'as_dense',fail)
    n=10**50
    matrix=cls.from_dok(n,n,{(n-1,0):3})
    assert matrix.shape==(n,n) and matrix.todok()=={(n-1,0):3}


@pytest.mark.parametrize('complex_entries',[False,True])
def test_bidiagonal_zero_pivots_and_complex_reflectors(complex_entries):
    matrix=sx.Matrix([[1,sx.I,1],[sx.I,2,3]]) if complex_entries else sx.Matrix([[0,1],[1,0]])
    u,b,right=matrix.bidiagonal_decomposition()
    un,bn,rn=[m.evalf(40).to_sympy().evalf(40) for m in (u,b,right)]
    assert (un*bn*rn-matrix.to_sympy()).norm() < sp.Float('1e-30')
    assert (un.H*un-sp.eye(u.rows)).norm() < sp.Float('1e-30')
    assert (rn*rn.H-sp.eye(right.rows)).norm() < sp.Float('1e-30')
    assert all(b[i,j]==0 for i in range(b.rows) for j in range(b.cols) if j not in (i,i+1))


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('rows',[2,4])
def test_cramer_rejects_mismatched_rhs_shape(kind,rows):
    with pytest.raises(sx.ShapeError):
        build(sx,kind,A3).cramer_solve(sx.ones(rows,1))


def test_auxiliary_matrix_operations_without_sympy():
    import os,subprocess,sys
    code='''
import sys
sys.modules['sympy']=None
import symox as s
x,n=s.symbols('x n')
for cls in (s.Matrix,s.ImmutableMatrix,s.SparseMatrix,s.ImmutableSparseMatrix):
    a=cls([[1,1],[0,1]])
    assert a.pow(n)[0,1]==n
    assert cls([[s.Abs(x)]]).refine(s.Q.positive(x))[0,0]==x
    assert cls([[x,x**2]]).integrate(x)[0,0]==x**2/2
    assert a.cramer_solve(s.Matrix([3,2]))[0,0]==1
    assert len(cls.orthogonalize(cls([1,0]),cls([1,1])))==2
    assert cls.jordan_block(2,1)==a
'''
    subprocess.run([sys.executable,'-c',code],env=os.environ.copy(),check=True,timeout=20)


@pytest.mark.parametrize('dummy',[False,True])
@pytest.mark.parametrize('collision',[False,True])
def test_charpoly_preserves_symbols_without_capturing_coefficients(dummy,collision):
    x=(sx.Dummy if dummy else sx.Symbol)('x',real=True)
    y,t=sx.symbols('y t')
    a=sx.Matrix([[x if collision else y,1],[2,x if collision else y]])
    p=a.charpoly(x)
    if collision:
        assert p.gens[0] != x
    else:
        assert p.gens[0] == x
    expression=p.as_expr().subs(p.gens[0],t)
    assert sx.expand(expression-(t*sx.eye(2)-a).det())==0
