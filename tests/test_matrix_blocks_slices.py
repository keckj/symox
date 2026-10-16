"""Exact block diagonals and rectangular views retain native matrix-valued DAGs."""
import pickle

import pytest
import sympy as sp
import symox as sx


def explicit(value):
    value = value.to_sympy()
    return value if isinstance(value, sp.MatrixBase) else value.as_explicit()


def same(actual, expected):
    expected = expected if isinstance(expected, sp.MatrixBase) else expected.as_explicit()
    assert explicit(actual).shape == expected.shape
    assert all(sp.simplify(a-b) == 0 for a, b in zip(explicit(actual), expected))


@pytest.mark.parametrize('shapes', [[], [(2, 3)], [(2, 3), (3, 2)], [(2, 2), (1, 1), (3, 2)]])
def test_block_diagonal_shapes_elements_and_interchange(shapes):
    ours = [sx.MatrixSymbol(f'A{i}', *shape) for i, shape in enumerate(shapes)]
    theirs = [sp.MatrixSymbol(f'A{i}', *shape) for i, shape in enumerate(shapes)]
    value, expected = sx.BlockDiagMatrix(*ours), sp.BlockDiagMatrix(*theirs)
    assert value.shape == expected.shape
    same(value, expected)
    assert isinstance(value, sx.BlockDiagMatrix)
    assert value.func(*value.args) == value
    assert sx.sympify(sx.srepr(value)) == value
    assert pickle.loads(pickle.dumps(value)) == value
    assert hash(pickle.loads(pickle.dumps(value))) == hash(value)
    assert value.diag == tuple(ours)
    assert value.blockshape == (len(shapes), len(shapes))
    assert value.rowblocksizes == [shape[0] for shape in shapes]
    assert value.colblocksizes == [shape[1] for shape in shapes]
    assert value.blocks.shape == (len(shapes), len(shapes))


def test_typed_block_grid_preserves_noncommutative_cells():
    A, B, C, D = [sx.MatrixSymbol(name, 2, 2) for name in 'ABCD']
    grid = sx.BlockDiagMatrix(A, B).blocks
    other = sx.BlockDiagMatrix(C, D).blocks
    assert grid[0, 0] == A and grid[1, 1] == B
    assert grid[0, 1] == sx.ZeroMatrix(2, 2)
    assert grid[-1, -1] == B
    assert (grid * other)[0, 0] == A*C
    assert (other * grid)[0, 0] == C*A
    assert (grid * other)[0, 0] != (other * grid)[0, 0]
    assert (grid + other)[1, 1] == B+D
    assert (3 * grid)[0, 0] == 3*A
    assert grid.T[0, 0] == A
    assert grid.free_symbols == {A, B}
    assert pickle.loads(pickle.dumps(grid)) == grid
    assert hash(pickle.loads(pickle.dumps(grid))) == hash(grid)
    assert grid.to_sympy()[0, 0] == sp.MatrixSymbol('A', 2, 2)
    for convert in [sx.Matrix, sx.MatrixExpr]:
        with pytest.raises(TypeError, match='block grid|Block grid'):
            convert(grid)
    with pytest.raises(TypeError):
        sx.Matrix([[1]]) * grid


def test_dense_materialization_intent_and_scalar_entry_guard():
    A, B = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    a, b = sp.MatrixSymbol('A', 2, 2), sp.MatrixSymbol('B', 2, 2)
    same(sx.Matrix(A), sp.Matrix(a))
    same(sx.Matrix([[A, B]]), sp.Matrix([[a, b]]))
    # Native dense matrices contain scalar expressions; matrix-valued scalar cells
    # are represented by the typed block grid instead of commutative proxies.
    for entries in [[A], lambda i, j: A]:
        with pytest.raises(TypeError, match='scalar entries'):
            sx.Matrix(1, 1, entries)
    n = sx.Symbol('n', integer=True, positive=True)
    with pytest.raises(NotImplementedError):
        sx.Matrix(sx.MatrixSymbol('N', n, n))


def test_block_diagonal_algebra_calculus_and_rectangular_singularity():
    A, B = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 1, 1)
    x = sx.Symbol('x', real=True)
    value = sx.BlockDiagMatrix(x*A, B)
    same(value.T, value.to_sympy().T)
    same(value.H, value.to_sympy().adjoint())
    same(value.diff(x), sp.BlockDiagMatrix(sp.MatrixSymbol('A', 2, 2), sp.ZeroMatrix(1, 1)))
    assert value.inverse() == sx.BlockDiagMatrix((x*A).inverse(), B.inverse())
    assert sx.det(value) == sx.det(x*A)*sx.det(B)
    assert sx.trace(value) == sx.trace(x*A)+sx.trace(B)
    singular = sx.BlockDiagMatrix(sx.MatrixSymbol('R', 2, 3), sx.MatrixSymbol('S', 3, 2))
    assert sx.det(singular) == 0
    with pytest.raises(Exception) as error:
        singular.inverse()
    assert type(error.value).__name__ == 'NonInvertibleMatrixError'
    n, m = sx.symbols('n m', integer=True, positive=True)
    unknown = sx.BlockDiagMatrix(sx.MatrixSymbol('U', n, m), sx.MatrixSymbol('V', m, n))
    assert isinstance(sx.det(unknown), sx.Determinant)
    assert isinstance(unknown.inverse(), sx.Inverse)
    assert sx.det(sx.BlockDiagMatrix()) == 1


@pytest.mark.parametrize('rows,cols', [(slice(None), slice(None)), (slice(1, 4), slice(2, 6)), (slice(-3, -1), slice(-5, None)), (2, slice(1, 5)), (slice(2, 4), -1)])
def test_matrix_slice_against_sympy(rows, cols):
    A, a = sx.MatrixSymbol('A', 5, 6), sp.MatrixSymbol('A', 5, 6)
    value, expected = A[rows, cols], a[rows, cols]
    assert isinstance(value, sx.MatrixSlice)
    assert value.parent == A
    assert value.shape == expected.shape
    same(value, expected)
    assert value.func(*value.args) == value
    assert pickle.loads(pickle.dumps(value)) == value
    assert sx.sympify(sx.srepr(value)) == value
    same(value.T, expected.T)
    same(value.H, expected.adjoint())


def test_slice_nested_negative_element_and_boundaries():
    A = sx.MatrixSymbol('A', 5, 6)
    view = A[2:4, 1:5]
    assert view[-1, -1] == A[3, 4]
    assert view[1:, 1:].parent == A
    assert view[1:, 1:][0, 0] == A[3, 2]
    assert A[1:3].shape == (2, 6)
    assert A[-1] == A[4, 5]
    assert A[7] == A[1, 1]
    for row in [slice(0, 0), slice(-7, 2), slice(1, 9), slice(None, None, -2)]:
        selected = list(range(5))[row]
        assert A[row, :].shape == (len(selected), 6)
        assert [A[row, :][i, 0] for i in range(len(selected))] == [A[i, 0] for i in selected]
    with pytest.raises(ValueError):
        A[::0, :]
    assert A[::2, :].shape == (3, 6)
    with pytest.raises(IndexError):
        view[2, 0]


def test_symbolic_slice_substitution_and_derivative():
    n, k = sx.symbols('n k', integer=True, positive=True)
    x = sx.Symbol('x', real=True)
    A = sx.MatrixSymbol('A', n, n)
    view = A[1:k, :]
    assert view.rows.is_integer is True
    assert view.rows.is_nonnegative is True
    assert view.xreplace({n:10, k:4}).shape == (3, 10)
    assert view.free_symbols == {A, k, n}
    replacement = sx.MatrixSymbol('B', 5, 5)
    assert view.xreplace({A: replacement, k: 4}).cols == n
    concrete = view.xreplace({A: replacement, k: 4, n: 5})
    assert concrete.shape == (3, 5)
    assert concrete[0, 0] == sx.MatrixSymbol('B', 5, 5)[1, 0]
    assert ((x*A)[1:k, :]).diff(x) == A[1:k, :]
    assert A[1:k, :].diff(x) == sx.ZeroMatrix(*view.shape)
    with pytest.raises(NotImplementedError, match='moving|start|bound'):
        A[k:n, :].diff(k)


def test_import_paths():
    from symox.matrices.expressions.blockmatrix import BlockDiagMatrix
    from symox.matrices.expressions.slice import MatrixSlice
    assert BlockDiagMatrix is sx.BlockDiagMatrix
    assert MatrixSlice is sx.MatrixSlice


@pytest.mark.parametrize('mode', ['blockdiag', 'slice'])
def test_numpy_source_and_tensor_source_lowering(mode):
    import ast
    import numpy as np
    A, B = sx.MatrixSymbol('A', 3, 3), sx.MatrixSymbol('B', 2, 2)
    values = {'A': np.arange(9.).reshape(3, 3), 'B': np.arange(4.).reshape(2, 2)}
    if mode == 'blockdiag':
        expression = sx.BlockDiagMatrix(A, B)
        expected = np.zeros((5, 5)); expected[:3, :3] = values['A']; expected[3:, 3:] = values['B']
    else:
        expression = A[1:3, :2]; expected = values['A'][1:3, :2]
    entries = expression.as_explicit()
    rows, cols = entries.shape
    source = 'numpy.array([' + ','.join('[' + ','.join(sx.numpycode(entries[i, j]) for j in range(cols)) + ']' for i in range(rows)) + '])'
    np.testing.assert_allclose(eval(source, {'numpy': np, **values}), expected)
    for generate in [sx.torchcode_function, sx.jaxcode_function]:
        source = generate('matrix_view', (A, B), expression)
        ast.parse(source)
        assert 'BlockDiagMatrix' not in source and 'MatrixSlice' not in source
        assert 'symox' not in source and 'sympy' not in source


def test_native_block_and_slice_runtime_without_sympy():
    import subprocess
    import sys
    source = '''
import sys
sys.modules['sympy'] = None
import symox as s
A = s.MatrixSymbol('A', 3, 3)
B = s.MatrixSymbol('B', 2, 2)
D = s.BlockDiagMatrix(A, B)
assert D.shape == (5, 5)
assert D[3, 4] == B[0, 1]
assert D[0, 4] == 0
assert D[1:4, 2:5].as_explicit().shape == (3, 3)
assert (D.blocks * D.blocks)[0, 0] == A**2
assert s.sympify(s.srepr(D)) == D
'''
    result = subprocess.run([sys.executable, '-c', source], capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr


def test_matrix_metadata_dependency_mask_is_cache_independent():
    A = sx.MatrixSymbol('MaskCache', 2, 2)
    entry = A[0, 1]
    hash(A)
    _ = A.free_symbols
    assert entry.has(A)
    assert entry.xreplace({A: sx.ImmutableMatrix([[1, 2], [3, 4]])}) == 2
    independently_built = sx.MatrixSymbol('MaskCache', 2, 2)
    assert independently_built == A
    assert independently_built[0, 1].has(A)
    assert independently_built[0, 1].xreplace({A: sx.ImmutableMatrix([[1, 2], [3, 4]])}) == 2


def test_sympy_matrixexpr_dense_conversion_preserves_shape_and_scalar_entries():
    a, b = sp.MatrixSymbol('A', 2, 2), sp.MatrixSymbol('B', 2, 2)
    same(sx.Matrix(a), sp.Matrix(a))
    same(sx.Matrix([[a, b]]), sp.Matrix([[a, b]]))
    for source in [[a], lambda i, j: a]:
        with pytest.raises(TypeError, match='scalar entries'):
            sx.Matrix(1, 1, source)


@pytest.mark.parametrize('foreign', [False, True])
def test_dense_entry_mutations_never_store_commutative_matrix_proxies(foreign):
    A = sp.MatrixSymbol('A', 2, 2) if foreign else sx.MatrixSymbol('A', 2, 2)
    dense = sx.Matrix([[1]])
    with pytest.raises(TypeError):
        dense[0, 0] = A
    with pytest.raises(TypeError):
        dense.applyfunc(lambda entry: A)
    same(sx.diag(A), sp.Matrix(sp.MatrixSymbol('A', 2, 2)))
