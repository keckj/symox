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


def test_integral_class_hook_preserves_held_calculus():
    x=sx.Symbol('x')
    X=sp.Symbol('x')
    held=sx.Integral(x,(x,0,1))
    reference=sp.Integral(X,(X,0,1))
    for doit in (False,True):
        assert sx.to_sympy(sx.simplify(held,doit=doit))==sp.simplify(reference,doit=doit)
    assert held.doit()==sx.Rational(1,2)


@pytest.mark.parametrize('doit',[False,True])
@pytest.mark.parametrize('source',[
    'Integral(Heaviside(t-x),(x,0,t))',
    'Integral(sin(x)**2+cos(x)**2,x)',
    'Integral(2*a*x,(x,0,1))',
    'Integral(0,x)',
    '2*Integral(x,(x,0,1))',
])
def test_integral_simplification_matches_reference_class_protocol(source,doit):
    namespaces=[]
    for module in (sx,sp):
        namespace={name:getattr(module,name) for name in ('Integral','Heaviside','sin','cos')}
        namespace.update(zip(('x','t','a'),module.symbols('x t a',real=True)))
        namespaces.append(namespace)
    actual=sx.simplify(eval(source,namespaces[0]),doit=doit)
    reference=sp.simplify(eval(source,namespaces[1]),doit=doit)
    assert sx.to_sympy(actual)==reference


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
