"""Simplify controls change candidate selection and evaluation behavior."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('rational',[False,True,None])
def test_rationalization_and_restoration_match_reference(rational):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=sx.Float('.5')*x+sx.Float('.25')*x
    expected=sp.simplify(expression.to_sympy(),rational=rational)
    result=sx.simplify(expr=expression,rational=rational)
    assert result.to_sympy()==expected
    assert expression.simplify(rational=rational)==result
    assert sx.simplify(sx.Float('.5'),rational=rational).is_Float


def test_inverse_combination_is_an_explicit_option():
    x=sx.Symbol('x')
    expression=sx.asin(sx.sin(x))
    assert sx.simplify(expression)==expression
    assert sx.simplify(expression,inverse=True)==x
    assert expression.simplify(inverse=True)==x


def test_held_calculus_is_evaluated_only_with_doit():
    x=sx.Symbol('x')
    held=sx.Integral(x,(x,0,1))
    assert sx.simplify(held,doit=False).has(sx.Integral)
    assert sx.simplify(held,doit=True)==sx.Rational(1,2)
    assert sx.simplify(held)==sx.Rational(1,2)


def test_ratio_limits_complexity_and_custom_measure_is_used():
    x=sx.Symbol('x')
    expression=(x+1)**2
    expanded=sx.expand(expression)
    calls=[]
    def prefer_expanded(candidate):
        calls.append(candidate)
        return 0 if candidate==expanded else 10
    result=sx.simplify(expression,measure=prefer_expanded,ratio=sx.oo)
    assert result==expanded and expanded in calls
    assert sx.simplify(expression,ratio=0)==expression


def test_matrix_controls_preserve_exact_entries():
    x=sx.Symbol('x')
    matrix=sx.Matrix([[sx.Float('.5')*x,sx.asin(sx.sin(x))]])
    result=sx.simplify(matrix,rational=True,inverse=True)
    assert result==sx.ImmutableMatrix([[x/2,x]])
