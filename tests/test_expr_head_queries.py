"""Structural has queries retain native literal and undefined-head matching."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('name',['Expr','Symbol','Integer','Rational','Add','Mul',
    'Pow','sin','cos','log','Derivative','Integral','Sum','Piecewise','Function',
    'Eq','Ne','Lt','Le','Gt','Ge'])
def test_class_and_function_queries_match_reference(name):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y');f=sx.Function('f');F=sp.Function('f')
    expressions=[x,sx.Integer(2),x+y,x*y,x**3,sx.sin(x)+sx.cos(y),sx.log(x),f(x),
        sx.Derivative(f(x),x,2),sx.Integral(sx.sin(x),x),sx.Sum(x,(x,1,3)),
        sx.Piecewise((x,x>0),(0,True)),sx.Eq(x,y),sx.Ne(x,y),x<y,x<=y,x>y,x>=y]
    references=[X,sp.Integer(2),X+Y,X*Y,X**3,sp.sin(X)+sp.cos(Y),sp.log(X),F(X),
        sp.Derivative(F(X),X,2),sp.Integral(sp.sin(X),X),sp.Sum(X,(X,1,3)),
        sp.Piecewise((X,X>0),(0,True)),sp.Eq(X,Y),sp.Ne(X,Y),X<Y,X<=Y,X>Y,X>=Y]
    for native,reference in zip(expressions,references):
        assert native.has(getattr(sx,name))==reference.has(getattr(sp,name))


def test_undefined_functions_do_not_alias_builtin_heads():
    x=sx.Symbol('x');unknown=sx.Function('sin');other=sx.Function('sin',real=True)
    assert not unknown(x).has(sx.sin)
    assert not sx.sin(x).has(unknown)
    assert unknown(x).has(unknown)
    assert not unknown(x).has(other)
    assert unknown(x).has(unknown(x).func)
    assert not sx.sin(x).has(unknown(x).func)
    assert sx.sin(x).has(sx.sin(x).func)
    assert not unknown(x).has(sx.sin(x).func)


def test_mixed_queries_literals_and_algebraic_subexpressions():
    x,y,z=sx.symbols('x y z');f=sx.Function('f')
    expression=x+y+z+sx.sin(x)
    assert expression.has(x+y)
    assert expression.has(f,sx.cos,x)
    assert expression.has(sx.sin,z)
    assert not expression.has(f,sx.cos,sx.Derivative)
    # A string pattern is skipped by the reference (strict ``_sympify``).
    assert not expression.has('x')
    assert not expression.has()
    assert f(x).has(sx.Function)
    assert sx.Derivative(f(x),x,2).has(sx.Integer)


@pytest.mark.parametrize('source',['3*a*x+4*x+5','a*x*y+b*y-x','a*sin(x)+b*cos(x)+3',
    '1/(x+a)','x+y','0','3.5*x+2','a*exp(x)+b-exp(x)-y'])
@pytest.mark.parametrize('generators',[(),('x',),('y',),('x','y')])
def test_coefficients_dict_generator_arguments_and_missing_entries(source,generators):
    actual=sx.sympify(source).as_coefficients_dict(*(sx.Symbol(name) for name in generators))
    expected=sp.sympify(source).as_coefficients_dict(*(sp.Symbol(name) for name in generators))
    assert {key.to_sympy():value.to_sympy() for key,value in actual.items()}==expected
    assert actual[sx.Symbol('missing_generator')]==0


def test_held_arithmetic_coefficient_dictionary_uses_public_arguments():
    x=sx.Symbol('x');X=sp.Symbol('x')
    for native,reference in [(sx.Add(3*x,4,evaluate=False),sp.Add(3*X,4,evaluate=False)),
        (sx.Mul(2,3*x,evaluate=False),sp.Mul(2,3*X,evaluate=False))]:
        for args in [(),(x,)]:
            actual=native.as_coefficients_dict(*args)
            expected=reference.as_coefficients_dict(*(value.to_sympy() for value in args))
            assert {key.to_sympy():value.to_sympy() for key,value in actual.items()}==expected
