"""Ordered block elimination, including rectangular and zero-dimensional factors."""
import pickle
import subprocess
import sys

import pytest
import sympy as sm
import symox as sx


def symbolic_blocks(module, p, q):
    return (module.MatrixSymbol('A', p, p), module.MatrixSymbol('B', p, q),
            module.MatrixSymbol('C', q, p), module.MatrixSymbol('D', q, q))


def grid(module, blocks):
    a, b, c, d = blocks
    return module.BlockMatrix([[a, b], [c, d]])


def dense_equal(actual, expected):
    actual = actual.to_sympy() if hasattr(actual, 'to_sympy') else actual
    assert actual.shape == expected.shape
    assert all(sm.simplify(a-b) == 0 for a, b in zip(actual, expected))


def explicit_blocks(module, p, r, q, pivot):
    make = lambda rows, cols, f: module.ImmutableMatrix(rows, cols, [f(i, j) for i in range(rows) for j in range(cols)])
    if pivot == 'A':
        a = make(p, p, lambda i, j: (2+module.I*(i+1)) if i == j else 0)
        b = make(p, q, lambda i, j: i+2*j+1)
        c = make(r, p, lambda i, j: 2*i-j+module.I)
        d = make(r, q, lambda i, j: i+3*j-2)
    else:
        a = make(p, r, lambda i, j: i+3*j-2)
        b = make(p, q, lambda i, j: i+2*j+1)
        c = make(q, r, lambda i, j: 2*i-j+module.I)
        d = make(q, q, lambda i, j: (2+module.I*(i+1)) if i == j else 0)
    return a, b, c, d


def assemble_dense(blocks):
    a, b, c, d = blocks
    return sm.ImmutableMatrix(a.rows+c.rows, a.cols+b.cols,
        lambda i, j: (a[i, j] if j < a.cols else b[i, j-a.cols]) if i < a.rows
        else (c[i-a.rows, j] if j < a.cols else d[i-a.rows, j-a.cols]))


def test_square_symbolic_oracle_and_reconstruction():
    n, m = sx.symbols('n m', integer=True, positive=True)
    ns, ms = sm.symbols('n m', integer=True, positive=True)
    value = grid(sx, symbolic_blocks(sx, n, m))
    reference = grid(sm, symbolic_blocks(sm, ns, ms))
    assert value.schur() == value.schur(mat='A')
    for pivot in ['A', 'D']:
        assert value.schur(pivot).to_sympy() == reference.schur(pivot)
    for method in ['LDUdecomposition', 'UDLdecomposition']:
        factors = getattr(value, method)()
        expected = getattr(reference, method)()
        assert isinstance(factors, tuple) and len(factors) == 3
        assert all(isinstance(f, sx.BlockMatrix) for f in factors)
        assert isinstance(factors[1], sx.BlockDiagMatrix)
        assert [f.to_sympy() for f in factors] == list(expected)
        assert sx.block_collapse(factors[0]*factors[1]*factors[2]) == value
        for factor in factors:
            assert sx.sympify(sx.srepr(factor)) == factor
            assert pickle.loads(pickle.dumps(factor)) == factor


@pytest.mark.parametrize('pivot', ['A', 'D'])
@pytest.mark.parametrize('p,r,q', [(2, 3, 4), (0, 3, 2), (2, 0, 3), (2, 3, 0), (0, 0, 0), (1, 1, 1)])
def test_rectangular_complex_dense_certificate_and_empty_shapes(pivot, p, r, q):
    blocks = explicit_blocks(sx, p, r, q, pivot)
    reference = explicit_blocks(sm, p, r, q, pivot)
    value = grid(sx, blocks)
    a, b, c, d = reference
    # Independent scalar-matrix formula; no reliance on SymPy's block LDU
    # identities, which use the wrong identity size for rectangular matrices.
    expected_schur = d-c*a.inv()*b if pivot == 'A' else a-b*d.inv()*c
    dense_equal(value.schur(pivot).as_explicit(), expected_schur)
    method = 'LDUdecomposition' if pivot == 'A' else 'UDLdecomposition'
    factors = getattr(value, method)()
    rows, cols = value.shape
    assert factors[0].shape == (rows, rows)
    assert factors[1].shape == value.shape
    assert factors[2].shape == (cols, cols)
    dense_factors = [f.as_explicit() for f in factors]
    dense_equal(dense_factors[0]*dense_factors[1]*dense_factors[2], assemble_dense(reference))
    collapsed = sx.block_collapse(factors[0]*factors[1]*factors[2])
    dense_equal(collapsed.as_explicit(), assemble_dense(reference))
    assert collapsed.shape == value.shape


@pytest.mark.parametrize('method', ['LDUdecomposition', 'UDLdecomposition'])
def test_singular_complement_and_whole_matrix_are_allowed(method):
    a = sx.Identity(2)
    value = sx.BlockMatrix([[a, a], [a, a]])
    factors = getattr(value, method)()
    assert factors[1].blocks[1, 1] == sx.ZeroMatrix(2, 2) if method == 'LDUdecomposition' else factors[1].blocks[0, 0] == sx.ZeroMatrix(2, 2)
    assert sx.block_collapse(factors[0]*factors[1]*factors[2]) == value


@pytest.mark.parametrize('p,q', [(0, 0), (0, 2), (2, 0), (2, 3)])
def test_two_block_diagonal_and_zero_dimensional_pivots(p, q):
    a, d = sx.Identity(p), sx.Identity(q)
    value = sx.BlockDiagMatrix(a, d)
    assert value.schur() == d and value.schur('D') == a
    for method in ['LDUdecomposition', 'UDLdecomposition']:
        factors = getattr(value, method)()
        dense_equal((factors[0]*factors[1]*factors[2]).as_explicit(), value.as_explicit().to_sympy())
    # Unique empty linear map is invertible, regardless of its spelling.
    empty_zero = sx.BlockDiagMatrix(sx.ZeroMatrix(0, 0), sx.Identity(2))
    assert empty_zero.schur() == sx.Identity(2)
    assert empty_zero.LDUdecomposition()[1].shape == (2, 2)


def test_pivot_singularity_and_square_validation_are_local():
    a, b, c, d = symbolic_blocks(sx, 2, 3)
    non_square_a = sx.BlockMatrix([[sx.MatrixSymbol('R', 2, 1), b], [sx.MatrixSymbol('S', 3, 1), d]])
    for operation in [lambda: non_square_a.schur(), non_square_a.LDUdecomposition]:
        with pytest.raises(sx.NonSquareMatrixError):
            operation()
    non_square_d = sx.BlockMatrix([[a, sx.MatrixSymbol('R', 2, 4)], [c, sx.MatrixSymbol('S', 3, 4)]])
    for operation in [lambda: non_square_d.schur('D'), non_square_d.UDLdecomposition]:
        with pytest.raises(sx.NonSquareMatrixError):
            operation()
    singular = sx.BlockDiagMatrix(sx.ImmutableMatrix([[0]]), sx.Identity(2))
    for operation in [lambda: singular.schur(), singular.LDUdecomposition]:
        with pytest.raises(sx.NonInvertibleMatrixError):
            operation()
    # The opposite pivot need not be invertible.
    singular.UDLdecomposition()
    swapped = sx.BlockDiagMatrix(sx.Identity(2), sx.ImmutableMatrix([[0]]))
    swapped.LDUdecomposition()
    for operation in [lambda: swapped.schur('D'), swapped.UDLdecomposition]:
        with pytest.raises(sx.NonInvertibleMatrixError):
            operation()


def test_unknown_pivot_squareness_is_not_assumed_for_factorization():
    p, r = sx.symbols('p r', integer=True, positive=True)
    value = sx.BlockMatrix([[sx.MatrixSymbol('A', p, r), sx.MatrixSymbol('B', p, 2)],
                            [sx.MatrixSymbol('C', 2, r), sx.Identity(2)]])
    with pytest.raises(NotImplementedError, match='certified'):
        value.LDUdecomposition()
    # A symbolic Schur expression retains the same formal inverse policy as
    # MatrixExpr.inverse; factorization must additionally certify partitions.
    assert value.schur().shape == (2, 2)


def test_grid_and_options_fail_explicitly():
    a = sx.MatrixSymbol('A', 2, 2)
    for value in [sx.BlockMatrix([[a]]), sx.BlockDiagMatrix(), sx.BlockDiagMatrix(a, a, a)]:
        for operation in [value.schur, value.LDUdecomposition, value.UDLdecomposition]:
            with pytest.raises(sx.ShapeError, match='2x2'):
                operation()
    value = sx.BlockDiagMatrix(a, a)
    with pytest.raises(NotImplementedError, match='Generalized'):
        value.schur(generalized=True)
    for pivot in ['B', 'C']:
        with pytest.raises(NotImplementedError, match='off-diagonal'):
            value.schur(pivot)
    with pytest.raises(ValueError, match='pivot'):
        value.schur('unknown')


def test_native_runtime_without_sympy():
    code = '''
import sys
sys.modules['sympy']=None
import symox as s
A=s.MatrixSymbol('A',2,2)
B=s.MatrixSymbol('B',2,4)
C=s.MatrixSymbol('C',3,2)
D=s.MatrixSymbol('D',3,4)
x=s.BlockMatrix([[A,B],[C,D]])
l,d,u=x.LDUdecomposition()
assert s.block_collapse(l*d*u)==x
assert x.schur()==D-C*A.inv()*B
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=15)
