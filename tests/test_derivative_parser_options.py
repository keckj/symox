"""Textual derivative specifications use the same validation as Python calls."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source',[
    'Derivative(sin(x),x,2)','Derivative(sin(x),[x,2])',
    'Derivative(sin(x),(x,2))','Derivative(sin(x),x,0)',
    'Derivative(sin(x),x,2,evaluate=True)','Derivative(sin(x))',
    'diff(sin(x),x,2)','diff(x**3)','diff(sin(x),x,2,evaluate=False)',
    'Derivative(exp(x*y),(x,2),y)','diff(exp(x*y),x,y)',
    'Derivative(2)','diff(2)','Derivative(f(x),x,Integer(2))',
])
def test_textual_derivative_forms_match_reference(source):
    actual=sx.sympify(source).to_sympy()
    expected=sp.sympify(source,locals={'f':sp.Function('f')})
    assert actual==expected or actual.expand()==expected.expand()


@pytest.mark.parametrize('source',[
    'Derivative(sin(x),x,0,0)','Derivative(sin(x),x,2,3)',
    'Derivative(sin(x),x,-1)','Derivative(sin(x),x,2.0)',
    'Derivative(sin(x),sin(x))','diff(x*y)',
    'Min(I,1)','Max(nan,1)','Min(zoo,1)',
])
def test_invalid_textual_derivatives_and_nonreal_extrema_raise(source):
    with pytest.raises(ValueError):sx.sympify(source)
