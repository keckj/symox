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
