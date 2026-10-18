"""Held derivative free symbols come from the body, not count metadata."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source',[
    'Derivative(f(x),x,evaluate=False)',
    'Derivative(f(y),x,3,evaluate=False)',
    'Derivative(f(x),y,evaluate=False)',
    'Derivative(1,x,evaluate=False)',
    'Derivative(x,x,2,evaluate=False)',
    'Derivative(f(x),f(y),evaluate=False)',
    'Lambda(x,Derivative(f(y),x,evaluate=False))',
])
def test_held_derivative_free_symbols_match_reference(source):
    x,y=sx.symbols('x y');f=sx.Function('f')
    X,Y=sp.symbols('x y');F=sp.Function('f')
    actual=eval(source,vars(sx)|{'x':x,'y':y,'f':f})
    expected=eval(source,vars(sp)|{'x':X,'y':Y,'f':F})
    assert {value.to_sympy() for value in actual.free_symbols}==expected.free_symbols


def test_billion_order_symbol_queries_do_not_enumerate_variables():
    x,y=sx.symbols('x y');f=sx.Function('f')
    derivative=sx.Derivative(f(y),x,2**32-1,evaluate=False)
    assert derivative.free_symbols=={y}
    assert derivative.has(x)
    assert derivative.has(sx.Integer(2**32-1))
    assert derivative.has(sx.Tuple(x,2**32-1))
    assert sx.Derivative(1,x,2**32-1,evaluate=False).free_symbols==set()
    assert sx.Subs(derivative,x,0).free_symbols=={y}


def test_undefined_function_head_queries_use_compact_derivative_metadata():
    x=sx.Symbol('x');f=sx.Function('f');g=sx.Function('g')
    held=sx.Derivative(f(x),(x,10**9),evaluate=False)
    assert held.has(f)
    assert not held.has(g)
    expression=x
    for _ in range(80):expression=f(expression,expression)
    assert expression.has(f) and not expression.has(g)
