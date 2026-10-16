"""Held calculus operators have presentation arguments distinct from storage."""
import pytest
import sympy as sp
import symox as sx


def normalized(key):
    if isinstance(key, tuple):
        return tuple(normalized(item) for item in key)
    return key.to_sympy() if isinstance(key, sx.Expr) else key


@pytest.mark.parametrize('operation', [
    lambda s, x, y, f: s.Derivative(f(x), x),
    lambda s, x, y, f: s.Derivative(f(x), (x, 3)),
    lambda s, x, y, f: s.Derivative(f(x, y), x, (y, 2)),
    lambda s, x, y, f: 2*s.Derivative(f(x), (x, 2)),
    lambda s, x, y, f: s.Derivative(f(x), (x, 2))**3,
    lambda s, x, y, f: s.Subs(f(x), x, 2),
    lambda s, x, y, f: s.Subs(f(x, y), (x, y), (y, 3)),
    lambda s, x, y, f: 2*s.Subs(f(x), x, 2),
    lambda s, x, y, f: s.Function('sin')(x),
    lambda s, x, y, f: s.Function('Tuple')(x),
    lambda s, x, y, f: s.UnevaluatedExpr(x+y),
    lambda s, x, y, f: s.Limit(f(x)+y, x, 0),
    lambda s, x, y, f: s.Limit(f(x), x, s.oo),
    lambda s, x, y, f: 2*s.Limit(f(x), x, 0, dir='-'),
    lambda s, x, y, f: s.Function('Limit')(x),
])
def test_held_calculus_keys_match_sympy(operation):
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    expression = operation(sx, x, y, sx.Function('f'))
    expected = operation(sp, X, Y, sp.Function('f'))
    assert normalized(sx.default_sort_key(expression)) == sp.default_sort_key(expected)


def test_matrix_key_preserves_names_of_alpha_equivalent_substitutions():
    x, y = sx.symbols('x y')
    f = sx.Function('f')
    left, right = sx.Subs(f(x), x, 0), sx.Subs(f(y), y, 0)
    assert left == right
    matrix = sx.MatrixExpr(sx.ImmutableMatrix([[left, right]]))
    key = sx.default_sort_key(matrix)
    cell_keys = key[1][1][2][1][1]
    assert cell_keys == (sx.default_sort_key(left), sx.default_sort_key(right))
    assert cell_keys[0] != cell_keys[1]
    assert 'f(x)' in str(matrix) and 'f(y)' in str(matrix)


@pytest.mark.parametrize('operation', [
    lambda s, a, b, d: s.BlockMatrix([[a, b], [b.T, d]]),
    lambda s, a, b, d: s.Transpose(s.BlockMatrix([[a, b], [b.T, d]])),
    lambda s, a, b, d: 2*s.BlockMatrix([[a, b], [b.T, d]]),
    lambda s, a, b, d: s.BlockMatrix([[a]]),
])
def test_block_matrix_presentation_keys_match_sympy(operation):
    actual = operation(sx, sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 1),
                       sx.MatrixSymbol('D', 1, 1))
    expected = operation(sp, sp.MatrixSymbol('A', 2, 2), sp.MatrixSymbol('B', 2, 1),
                         sp.MatrixSymbol('D', 1, 1))
    assert normalized(sx.default_sort_key(actual)) == sp.default_sort_key(expected)


@pytest.mark.parametrize('name', ['DiagMatrix', 'DiagonalOf', 'DiagonalMatrix'])
@pytest.mark.parametrize('shape', [(3,1), (1,3), (1,1)])
def test_diagonal_matrix_presentation_keys_match_sympy(name, shape):
    actual = getattr(sx, name)(sx.MatrixSymbol('A', *shape))
    expected = getattr(sp, name)(sp.MatrixSymbol('A', *shape))
    assert normalized(sx.default_sort_key(actual)) == sp.default_sort_key(expected)


@pytest.mark.parametrize('dummy',[False,True])
def test_single_lambda_key_uses_tuple_of_variables_and_preserves_identity(dummy):
    if dummy:
        x=sx.Dummy('x',dummy_index=37,integer=True)
        X=sp.Dummy('x',dummy_index=37,integer=True)
    else:
        x=sx.Symbol('x');X=sp.Symbol('x')
    actual=sx.Lambda(x,2*x+1)._e
    expected=sp.Lambda(X,2*X+1)
    for expression,oracle in [(actual,expected),(2*actual,2*expected),(actual**2,expected**2)]:
        assert normalized(sx.default_sort_key(expression))==sp.default_sort_key(oracle)
    if dummy:
        other=sx.Dummy('x',dummy_index=38,integer=True)
        assert sx.default_sort_key(sx.Lambda(other,2*other+1)._e)!=sx.default_sort_key(actual)


@pytest.mark.parametrize('operation',[
    lambda d,x:d,
    lambda d,x:-d,
    lambda d,x:2*d,
    lambda d,x:d**2,
    lambda d,x:2*d**2,
    lambda d,x:d**x,
    lambda d,x:2*d**x,
    lambda d,x:d**-1,
])
def test_dummy_atom_and_coefficient_power_keys_have_distinct_argument_shapes(operation):
    from sympy.core.cache import clear_cache
    clear_cache()
    d=sx.Dummy('key_dummy',dummy_index=113,real=True)
    D=sp.Dummy('key_dummy',dummy_index=113,real=True)
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=operation(d,x);expected=operation(D,X)
    assert normalized(sx.default_sort_key(actual))==sp.default_sort_key(expected)
