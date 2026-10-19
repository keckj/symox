"""Documented circle constructors and curve-parameter substitution domains."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source',[
    'x**2+y**2-25','2*x**2+2*y**2-4*x+8*y-40',
    '(x-3)**2+(y+2)**2-7','x**2+y**2',
    'a*(x**2+y**2-1)',
])
def test_circle_from_polynomial_and_equality(source):
    actual=sx.sympify(source);expected=sp.sympify(source)
    for a,b in [(actual,expected),(sx.Eq(actual,0),sp.Eq(expected,0))]:
        result=sx.Circle(a)
        assert result.to_sympy()==sp.Circle(b)


def test_circle_named_variables_and_degenerate_three_points():
    a,b=sx.symbols('a b');A,B=sp.symbols('a b')
    assert sx.Circle(sx.Eq(a*a+b*b,25),x='a',y=b).to_sympy()==sp.Circle(sp.Eq(A*A+B*B,25),x='a',y=B)
    for points in [((0,0),(1,0),(2,0)),((0,0),(0,0),(0,0))]:
        assert sx.Circle(*points).to_sympy()==sp.Circle(*points)
    for source in ['x*x+2*y*y-1','x*y-1','x**3+y*y-1','x*x+y*y+1']:
        with pytest.raises(sx.GeometryError):sx.Circle(sx.sympify(source))


@pytest.mark.parametrize('value',[4,'s','2*s+1'])
def test_curve_parameter_substitution_evaluates_to_a_point(value):
    t=sx.Symbol('t');T=sp.Symbol('t')
    a=sx.Curve((t,t*t),(t,0,1));b=sp.Curve((T,T*T),(T,0,1))
    value=sx.sympify(value);reference=value.to_sympy()
    assert a.subs(t,value).to_sympy()==b.subs(T,reference)
    assert a.subs({t:value}).to_sympy()==b.subs({T:reference})


def test_curve_other_substitutions_preserve_parameter_binding():
    t,x=sx.symbols('t x');T,X=sp.symbols('t x')
    a=sx.Curve((x*t,t*t),(t,0,x));b=sp.Curve((X*T,T*T),(T,0,X))
    assert a.subs(x,3).to_sympy()==b.subs(X,3)
    assert a.subs([(x,3),(t,2)]).to_sympy()==b.subs([(X,3),(T,2)])
