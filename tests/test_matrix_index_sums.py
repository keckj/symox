"""Recover typed matrix products from full, explicit index contractions."""
import itertools
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('directions', list(itertools.product([False, True], repeat=4)))
@pytest.mark.parametrize('reverse', [False, True])
def test_rectangular_chains(directions, reverse):
    indices = sx.symbols('i0:5', integer=True)
    sizes = [2, 3, 4, 2, 5]
    factors, matrices = [], []
    for j, transposed in enumerate(directions):
        rows, columns = (sizes[j+1], sizes[j]) if transposed else (sizes[j], sizes[j+1])
        a = sx.MatrixSymbol('A'+str(j), rows, columns)
        factors.append(a[indices[j+1], indices[j]] if transposed else a[indices[j], indices[j+1]])
        matrices.append(a.T if transposed else a)
    expression = sx.Sum(sx.Mul(*factors),
                        *((indices[j], 0, sizes[j]-1) for j in range(1, 4)))
    first, last = (indices[-1], indices[0]) if reverse else (indices[0], indices[-1])
    actual = sx.MatrixExpr.from_index_summation(expression, first_index=first, last_index=last)
    expected = sx.MatMul(*matrices).doit()
    if reverse:
        expected = expected.T.doit()
    assert actual == expected
    reference = sp.MatrixExpr.from_index_summation(expression.to_sympy(), first_index=first.to_sympy(),
                                                    last_index=last.to_sympy())
    assert actual.to_sympy().doit() == reference.doit()


def test_default_orientation_and_official_examples():
    n = sx.Symbol('n', positive=True, integer=True)
    i, j, k, l = sx.symbols('i j k l', integer=True)
    a, b = sx.MatrixSymbol('A', n, n), sx.MatrixSymbol('B', n, n)
    cases = [
        (sx.Sum(a[i, j]*b[j, k], (j, 0, n-1)), a*b),
        (sx.Sum(a[j, i]*b[j, k], (j, 0, n-1)), a.T*b),
        (sx.Sum(a[i, i], (i, 0, n-1)), sx.Trace(a)),
        (sx.Sum(a[i, j]*b[k, j]*a[l, k], (j, 0, n-1), (k, 0, n-1)), a*b.T*a.T),
    ]
    for expression, expected in cases:
        assert sx.MatrixExpr.from_index_summation(expression) == expected


def test_hadamard_products_and_powers():
    i, j = sx.symbols('i j', integer=True)
    a, b = sx.MatrixSymbol('A', 2, 3), sx.MatrixSymbol('B', 3, 2)
    assert sx.MatrixExpr.from_index_summation(a[i, j]*b[j, i]) == sx.HadamardProduct(a, b.T)
    assert sx.MatrixExpr.from_index_summation(a[i, j]**2) == sx.HadamardPower(a, 2)
    assert sx.MatrixExpr.from_index_summation(a[i, j], first_index=j) == a.T


def test_traces_and_disconnected_scalar_components():
    i, j, k, l = sx.symbols('i j k l', integer=True)
    a, b, c = [sx.MatrixSymbol(name, 2, 2) for name in 'ABC']
    expr = sx.Sum(a[i, j]*b[j, i], (i, 0, 1), (j, 0, 1))
    result = sx.MatrixExpr.from_index_summation(expr)
    assert result == sx.Trace(a*b)
    expression = sx.Sum(a[i, j]*b[j, k]*c[l, l], (j, 0, 1), (l, 0, 1))
    assert sx.MatrixExpr.from_index_summation(expression) == sx.Trace(c)*a*b
    expression = sx.Sum(a[i, i]+b[i, i], (i, 0, 1))
    assert sx.MatrixExpr.from_index_summation(expression) == sx.Trace(a)+sx.Trace(b)


def test_linear_combinations_and_unused_sums():
    i, j, k, l = sx.symbols('i j k l', integer=True)
    a, b, c = [sx.MatrixSymbol(name, 2, 2) for name in 'ABC']
    expression = sx.Sum(2*a[i, j]*b[j, k]+3*a[i, j]*c[j, k], (j, 0, 1))
    assert sx.MatrixExpr.from_index_summation(expression) == 2*a*b+3*a*c
    expression = sx.Sum(a[i, j]*b[j, k], (j, 0, 1), (l, 2, 4))
    assert sx.MatrixExpr.from_index_summation(expression) == 3*a*b
    assert sx.MatrixExpr.from_index_summation(sx.Sum(3, (l, 2, 4))) == 9


def test_bounds_dimensions_and_unsupported_arrays_are_checked():
    i, j, k = sx.symbols('i j k', integer=True)
    a, b = sx.MatrixSymbol('A', 2, 3), sx.MatrixSymbol('B', 3, 4)
    for bounds in [(j, 1, 2), (j, 0, 1), (j, 0, 3)]:
        with pytest.raises(ValueError):
            sx.MatrixExpr.from_index_summation(sx.Sum(a[i, j]*b[j, k], bounds))
    with pytest.raises(ValueError):
        sx.MatrixExpr.from_index_summation(a[i, j], dimensions={i: 3})
    with pytest.raises(ValueError):
        sx.MatrixExpr.from_index_summation(a[i, j], first_index=k)
    with pytest.raises(NotImplementedError):
        sx.MatrixExpr.from_index_summation(a[i, j]*b[j, k])
    with pytest.raises(NotImplementedError):
        sx.MatrixExpr.from_index_summation(sx.Sum(j*a[i, j]*b[j, k], (j, 0, 2)))


def test_fromiter_uses_concrete_constructor():
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    assert sx.MatAdd.fromiter(iter([a, b])) == sx.MatAdd(a, b)
    assert (a*b).fromiter([b, a]) == b*a
    assert sx.MatrixSymbol.fromiter(['C', 2, 3]) == sx.MatrixSymbol('C', 2, 3)


def test_contraction_without_sympy():
    code = '''
import sys
sys.modules['sympy']=None
import symox as s
i,j,k=s.symbols('i j k',integer=True)
a=s.MatrixSymbol('A',2,3); b=s.MatrixSymbol('B',3,4)
expr=s.Sum(a[i,j]*b[j,k],(j,0,2))
assert s.MatrixExpr.from_index_summation(expr)==a*b
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
