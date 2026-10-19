"""Derivative construction and constructor-level methods use actual native nodes."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('options',[{}, {'evaluate':False}, {'evaluate':True},
    {'simplify':False}, {'unused_option':True}])
def test_derivative_reference_keyword_and_type(options):
    x=sx.Symbol('x');X=sp.Symbol('x')
    value=sx.Derivative(expr=x*x,**options)
    expected=sp.Derivative(expr=X*X,**options)
    assert value.to_sympy()==expected
    assert isinstance(value,sx.Derivative)==isinstance(expected,sp.Derivative)
    assert sx.Derivative.doit(value)==value.doit()


@pytest.mark.parametrize('name',['Integral','Sum','Product'])
def test_calculus_constructors_expose_unbound_methods(name):
    x=sx.Symbol('x');X=sp.Symbol('x')
    constructor=getattr(sx,name)
    value=constructor(x,(x,1,3))
    expected=getattr(sp,name)(X,(X,1,3))
    assert constructor.doit(value).to_sympy()==expected.doit()
    assert constructor.has(value,x)
    assert constructor.class_key()==getattr(sp,name).class_key()


@pytest.mark.parametrize('options',[{'unused_option':True},{'evaluate':False},{'commutative':False}])
def test_subs_accepts_reference_constructor_assumptions(options):
    x=sx.Symbol('x');X=sp.Symbol('x')
    value=sx.Subs(expr=x*x,variables=x,point=2,**options)
    expected=sp.Subs(expr=X*X,variables=X,point=2,**options)
    assert value.to_sympy()==expected
    assert sx.Subs.doit(value).to_sympy()==expected.doit()
