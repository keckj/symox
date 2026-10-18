import pytest
import symox as s
import sympy as sp


@pytest.mark.parametrize('value',[0,1,-2,1.0,1.5,'2',s.Rational(3,2),s.Symbol('n'),s.oo,s.I])
def test_intlike(value):
    from sympy.functions.special.gamma_functions import intlike
    reference=value.to_sympy() if hasattr(value,'to_sympy') else value
    try:
        expected=intlike(reference)
    except TypeError:
        with pytest.raises(TypeError):s.intlike(value)
    else:
        assert s.intlike(value)==expected


@pytest.mark.parametrize('p',range(1,8))
@pytest.mark.parametrize('x',[s.Rational(7,2),s.Rational(3,4),s.Symbol('x')])
def test_multigamma_product_and_derivative(p,x):
    actual=s.multigamma(x,p)
    expected=sp.multigamma(x.to_sympy(),p)
    assert sp.simplify(actual.to_sympy()-expected)==0
    if x.is_Symbol:
        assert sp.simplify(actual.diff(x).to_sympy()-expected.diff(x.to_sympy()))==0


def test_multigamma_symbolic_order_and_invalid_dimensions():
    x,p=s.symbols('x p',positive=True)
    n=s.Symbol('n',positive=True,integer=True)
    f=s.multigamma(x,n)
    assert f.has(s.Product)
    assert s.simplify(f.subs(n,2).doit()-s.multigamma(x,2))==0
    for order in (0,-1,s.Rational(1,2)):
        with pytest.raises(ValueError):s.multigamma(x,order)
