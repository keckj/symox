"""Ordered term data and leading coefficients on real and complex rays."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def converted(value):
    if isinstance(value,(list,tuple)):return type(value)(converted(item) for item in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


@pytest.mark.parametrize('source',['x*x+2*x*y+1','x-y','sin(x)+cos(y)+2',
    '3*x*x+Rational(1,2)*x','x*x+y**3+x*y','0','2*I*x+x*x'])
def test_ordered_terms_with_coefficient_data_match_reference(source):
    expression=sx.sympify(source);reference=sp.sympify(source)
    assert converted(expression.as_ordered_terms(order='lex'))==reference.as_ordered_terms(order='lex')
    assert converted(expression.as_ordered_terms(data=True))==reference.as_ordered_terms(data=True)
    assert converted(expression.as_ordered_factors(order='lex'))==reference.as_ordered_factors(order='lex')


def test_ordering_preserves_alpha_equal_constructor_presentations():
    x,y=sx.symbols('x y');f=sx.Function('f')
    first=sx.Subs(f(x),x,0);second=sx.Subs(f(y),y,0)
    expression=sx.Add(first,second,evaluate=False)
    terms,generators=expression.as_ordered_terms(data=True)
    assert len(terms)==2
    assert {str(term) for term,_ in terms}=={str(first),str(second)}
    with pytest.raises(ValueError):expression.as_ordered_terms(order='invalid')
    for order in ['grlex','grevlex']:
        assert expression.as_ordered_factors(order=order)==[expression]


@pytest.mark.parametrize('order',['lex','grlex','grevlex','ilex','igrlex','igrevlex',
    'rev-lex','rev-grlex','rev-grevlex','rev-ilex','rev-igrlex','rev-igrevlex'])
@pytest.mark.parametrize('source',['x*x+y**3+x*y+1','x**4+x*y*y+y**3+2',
    'sin(x)**2*cos(x)+sin(x)**2+1','2*I*x+x*x-y','1-x','x**(-2)+x+y**(-3)'])
def test_monomial_term_orders_and_data_match_reference(order,source):
    expression=sx.sympify(source);reference=sp.sympify(source)
    assert converted(expression.as_ordered_terms(order=order))==reference.as_ordered_terms(order=order)
    assert converted(expression.as_ordered_terms(order=order,data=True))==reference.as_ordered_terms(order=order,data=True)


def test_callable_and_symbol_monomial_orders():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    expression=x*x+y**3+x*y;reference=X*X+Y**3+X*Y
    order=lambda monomial:(sum(monomial),tuple(reversed(monomial)))
    assert converted(expression.as_ordered_terms(order=order,data=True))==reference.as_ordered_terms(order=order,data=True)
    assert converted(expression.as_ordered_terms(order=sx.Symbol('grlex')))==reference.as_ordered_terms(order='grlex')
    with pytest.raises(ValueError):expression.as_ordered_terms(order=17)


@pytest.mark.parametrize('order',[None,'lex','grlex','grevlex','rev-lex'])
def test_asymptotic_order_terms_remain_after_ordinary_terms(order):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=1+x+x*x+sx.Order(x**3,x)
    reference=1+X+X*X+sp.Order(X**3,X)
    assert converted(expression.as_terms())==reference.as_terms()
    assert converted(expression.as_ordered_terms(order=order,data=True))==reference.as_ordered_terms(order=order,data=True)


@pytest.mark.parametrize('source',['x*x','exp(x)-1','sin(x)','1/x','sqrt(x)','x*x*log(x)',
    'log(x)','log(-x)','(1-cos(x))/x'])
@pytest.mark.parametrize('direction',[0,1,-1,sx.I,-sx.I])
def test_leading_coefficients_match_reference_along_rays(source,direction):
    x=sx.Symbol('x');X=sp.Symbol('x');L=sx.Symbol('L');LL=sp.Symbol('L')
    actual=sx.sympify(source).leadterm(x,logx=L,cdir=direction)
    expected=sp.sympify(source).leadterm(X,logx=LL,cdir=sx.sympify(direction).to_sympy())
    assert actual[1].to_sympy()==expected[1]
    assert sp.simplify(actual[0].to_sympy()-expected[0])==0


@pytest.mark.parametrize('direction',[1,-1,sx.I,-sx.I,1+sx.I])
def test_fractional_negative_base_reconstructs_the_specified_ray(direction):
    x=sx.Symbol('x');value=sx.sqrt(-x)
    leading=value.as_leading_term(x,cdir=direction)
    point=direction*sx.Rational(1,10**6)
    assert abs(complex((leading.subs(x,point)/value.subs(x,point)).evalf(40))-1)<1e-12


def test_leading_terms_multiple_variables_centers_and_invalid_options():
    x,y,L=sx.symbols('x y L')
    value=(x+x*x)*(y+y*y)
    assert value.as_leading_term() is value
    assert value.as_leading_term(x,y)==x*y
    assert (sx.exp(x)-1).leadterm(x,dir='-')==(1,1)
    assert ((x-2)**3).leadterm(x,x0=2,dir='-')==(1,3)
    assert sx.log(x-2).leadterm(x,logx=L,x0=2)==(L,0)
    with pytest.raises(ValueError):value.leadterm(x,cdir=x)
    with pytest.raises(ValueError):value.leadterm(x,cdir=sx.oo)
    with pytest.raises(TypeError):value.as_leading_term(x,unknown=True)
    with pytest.raises(ValueError):value.leadterm(x,dir='both')


def test_matrix_leading_terms_use_original_variable_coefficients():
    x,L=sx.symbols('x L');A=sx.MatrixSymbol('A',2,2)
    for direction in [-1,sx.I,-sx.I]:
        assert ((sx.exp(x)-1)*A).leadterm(x,cdir=direction)==(A,1)
        coefficient,power=(sx.log(-x)*A).leadterm(x,logx=L,cdir=direction)
        scalar,_=sx.log(-x).leadterm(x,logx=L,cdir=direction)
        assert coefficient==scalar*A and power==0


def test_refine_accepts_the_sympy_method_keyword():
    x=sx.Symbol('x')
    expression=sx.sqrt(x*x)
    assert expression.refine(assumption=sx.Q.positive(x))==x
    assert expression.refine(assumptions=sx.Q.positive(x))==x
    assert expression.refine(assumption=sx.Q.negative(x))==-x
    with pytest.raises(TypeError):expression.refine(unknown=True)


def test_method_options_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');L=Symbol('L')
assert (exp(x)-1).leadterm(x,cdir=-1)==(1,1)
assert log(x).leadterm(x,logx=L,cdir=I)==(L,0)
assert sqrt(x*x).refine(assumption=Q.positive(x))==x
terms,gens=(x*x+1).as_ordered_terms(data=True)
assert terms[0][0]==x*x and gens==[x]
'''
    subprocess.run([sys.executable,'-c',source],check=True)
