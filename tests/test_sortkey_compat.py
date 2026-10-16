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
])
def test_held_calculus_keys_match_sympy(operation):
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    expression = operation(sx, x, y, sx.Function('f'))
    expected = operation(sp, X, Y, sp.Function('f'))
    assert normalized(sx.default_sort_key(expression)) == sp.default_sort_key(expected)
