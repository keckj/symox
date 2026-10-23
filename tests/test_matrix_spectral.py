"""Exact spectral matrix algorithms: oracle checks and defining identities."""
import pytest
import sympy as sp
import symox as sx


def equal(actual, expected):
    actual = actual.to_sympy() if hasattr(actual, 'to_sympy') else actual
    if hasattr(expected, 'shape'):
        assert actual.shape == expected.shape
        assert all(sp.simplify(a-b) == 0 for a,b in zip(actual,expected))
    else:
        assert actual == expected or sp.simplify(actual-expected) == 0


@pytest.mark.parametrize('data', [
    [[3,0],[0,4]],
    [[1,2],[3,4]],
    [[1,2,3],[2,4,6]],
    [[1,0],[0,1],[1,1]],
    [[1,sp.I],[sp.I,-1]],
    [[0,0],[0,0]],
    [[1,0,0]],
    [[1],[0],[0]],
    [[sp.Symbol('a',positive=True),0],[0,2]],
])
def test_singular_values_and_condition_number(data):
    a, reference = sx.Matrix(data), sp.Matrix(data)
    values = a.singular_values()
    expected = reference.singular_values()
    assert len(values) == len(expected) == a.cols
    for actual, ref in zip(values, expected):
        equal(actual, ref)
    equal(a.condition_number(), reference.condition_number())
    if a.rows != 1 and a.cols != 1:
        equal(a.norm(2), reference.norm(2))
        equal(a.norm(-2), reference.norm(-2))


@pytest.mark.parametrize('shape', [(0,0),(0,3),(4,0)])
def test_empty_spectral_quantities(shape):
    a = sx.zeros(*shape)
    assert len(a.singular_values()) == shape[1]
    assert all(value == 0 for value in a.singular_values())
    assert a.condition_number() == 0
    assert a.pinv(method='ED').shape == shape[::-1]


@pytest.mark.parametrize('data', [
    [[1,2],[3,4]],
    [[1,0],[0,2],[0,0]],
    [[1,2,3],[2,4,6]],
    [[1,sp.I],[sp.I,-1]],
    [[1,sp.I],[0,2],[2,1]],
    [[2,0,0],[0,2,0]],
    [[0,0],[0,0]],
    [[sp.Symbol('a',positive=True),0],[0,2]],
])
def test_ed_pseudoinverse(data):
    a, reference = sx.Matrix(data), sp.Matrix(data)
    p = a.pinv(method='ED')
    equal(p, reference.pinv())
    equal(p, a.pinv().to_sympy())
    equal(a*p*a, reference)
    equal(p*a*p, p.to_sympy())
    equal((a*p).H, (a*p).to_sympy())
    equal((p*a).H, (p*a).to_sympy())


def jordan_block(size, value):
    result = sp.eye(size)*value
    for i in range(size-1):
        result[i,i+1] = 1
    return result


@pytest.mark.parametrize('blocks', [
    [(1,0)],
    [(3,0)],
    [(3,2),(2,2),(1,2),(2,-1)],
    [(2,sp.I),(1,sp.I),(2,-sp.I)],
    [(2,sp.Symbol('a')),(1,sp.Symbol('a'))],
    [(1,2),(1,3),(1,3)],
])
def test_jordan_blocks_and_change_of_basis(blocks):
    reference_j = sp.diag(*(jordan_block(size,value) for size,value in blocks))
    n = reference_j.rows
    # Nontrivial exact triangular similarity, deterministic and invertible.
    transform = sp.eye(n)
    for i in range(n-1):
        transform[i,i+1] = i+1
    reference_a = transform*reference_j*transform.inv()
    a = sx.ImmutableMatrix(reference_a)
    p,j = a.jordan_form()
    expected = reference_a.jordan_form(calc_transform=False)
    equal(j, expected)
    equal(a*p, p.to_sympy()*j.to_sympy())
    assert p.det() != 0
    assert isinstance(p, sx.ImmutableMatrix)
    assert isinstance(j, sx.ImmutableMatrix)
    equal(a.jordan_form(calc_transform=False), expected)
    q,cells = a.jordan_cells()
    assert all(isinstance(cell,sx.ImmutableMatrix) for cell in cells)
    equal(q* sx.diag(*cells), a.to_sympy()*q.to_sympy())
    equal(sx.diag(*cells), expected)


def test_jordan_distinct_algebraic_eigenvalues():
    a = sx.Matrix([[1,1],[1,0]])
    p,j = a.jordan_form()
    equal(j, a.to_sympy().jordan_form(calc_transform=False))
    equal(a*p, p.to_sympy()*j.to_sympy())


def test_empty_jordan_and_invalid_arguments():
    a = sx.zeros(0,0)
    p,j = a.jordan_form()
    assert p.shape == j.shape == (0,0)
    assert a.jordan_cells()[1] == []
    with pytest.raises(sx.NonSquareMatrixError):
        sx.ones(2,3).jordan_form()
    # SymPy 1.14 accepts chop (used only for Float matrices).
    assert sx.eye(2).jordan_form(chop=True) == (sx.eye(2), sx.eye(2))
    with pytest.raises(ValueError):
        sx.eye(2).pinv(method='bad')


def test_unresolved_eigenproblem_is_explicit():
    # No numerical approximation should masquerade as an exact Jordan form.
    a = sx.Matrix([[0,0,0,0,-1],[1,0,0,0,1],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0]])
    with pytest.raises((NotImplementedError,sx.MatrixError)):
        a.jordan_form()


@pytest.mark.parametrize('value',[sx.oo,-sx.oo,sx.zoo,sx.nan])
def test_dense_nonfinite_products_do_not_discard_zero_terms(value):
    zero=sx.Matrix([[0]])
    special=sx.Matrix([[value]])
    assert (zero*special)[0] == sx.nan
    assert (special*zero)[0] == sx.nan
    a=sx.Matrix([[0,2],[1,0]])
    b=sx.Matrix([[value,0,1],[0,3,0]])
    expected=a.to_sympy()*b.to_sympy()
    assert a*b == sx.Matrix(expected)
    assert b.T*a.T == sx.Matrix(b.T.to_sympy()*a.T.to_sympy())
    # The same scalar semantics apply to the sparse/dense reflected routes.
    assert sx.SparseMatrix(a)*b == sx.Matrix(expected)
    assert a*sx.SparseMatrix(b) == sx.Matrix(expected)
