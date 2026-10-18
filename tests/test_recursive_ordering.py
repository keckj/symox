"""Explicit monomial orders reach nested factors and structural metadata."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def converted(value):
    if isinstance(value,tuple):return tuple(converted(child) for child in value)
    if isinstance(value,list):return [converted(child) for child in value]
    return value.to_sympy() if hasattr(value,'to_sympy') else value


def expression(module,source):
    x,y=module.symbols('x y')
    return eval(source,vars(module)|{'x':x,'y':y})


@pytest.mark.parametrize('order',[None,'lex','grlex','grevlex','ilex','igrlex','igrevlex',
    'rev-lex','rev-grlex','rev-grevlex'])
@pytest.mark.parametrize('source',[
    '1-x','x*x+y**3+x*y','(x*x+y**3)*(x**3+y)',
    'sin(x*x+y**3)','2*(x*x+y**3)**(x+y)','exp(x*x+y**3)',
    'Derivative(sin(x*x+y**3),x,evaluate=False)',
    'Piecewise((x*x+y**3,x>0),(0,True))',
    'Tuple(x*x+y**3,x)','FiniteSet(x*x+y**3,0)',
    'Interval(x*x+y**3,oo)','Lambda(x,x*x+y**3)',
    'Subs(x*x+y**3,x,y)','Order(x*x,x)','I','-oo','Rational(1,3)',
    'MatrixSymbol("A",x*x+y**3,2)',
    '(x*x+y**3)*MatrixSymbol("A",2,2)+x*MatrixSymbol("B",2,2)',
    'Trace((x*x+y**3)*MatrixSymbol("A",2,2))',
    'AlgebraicNumber(sqrt(2),[1,2])',
    'conjugate(MatrixSymbol("A",2,2))',
    '(conjugate(x)**2+conjugate(y)**3)*conjugate(MatrixSymbol("A",2,2))',
    'BlockMatrix([[MatrixSymbol("A",2,2),MatrixSymbol("B",2,2)]])',
])
def test_recursive_sort_keys_match_reference(source,order):
    actual=expression(sx,source);expected=expression(sp,source)
    assert converted(sx.default_sort_key(actual,order))==sp.default_sort_key(expected,order)


@pytest.mark.parametrize('order',['lex','grlex','grevlex','rev-lex','igrlex'])
@pytest.mark.parametrize('source',[
    '(x*x+y**3)*(x**3+y)*sin(x*x+y**3)',
    '2*(x*x+y**3)**2/(x+y)',
    '(x*x+y**3)*MatrixSymbol("A",2,2)*MatrixSymbol("B",2,2)',
])
def test_ordered_factors_recurse_without_reordering_matrix_products(order,source):
    actual=expression(sx,source);expected=expression(sp,source)
    assert converted(actual.as_ordered_factors(order))==expected.as_ordered_factors(order)


@pytest.mark.parametrize('order',['lex','grlex','grevlex','rev-lex','igrlex'])
def test_matrix_terms_share_scalar_monomial_order_rules(order):
    source='x*x*MatrixSymbol("A",2,2)+y**3*MatrixSymbol("B",2,2)+x*y*MatrixSymbol("A",2,2)'
    actual=expression(sx,source);expected=expression(sp,source)
    assert converted(actual.as_ordered_terms(order=order,data=True))==expected.as_ordered_terms(order=order,data=True)


def test_callable_orders_and_container_keys_propagate_to_scalar_children():
    native=expression(sx,'[x*x+y**3, sin(x*x+y**3)]')
    foreign=expression(sp,'[x*x+y**3, sin(x*x+y**3)]')
    order=lambda powers:(sum(powers),tuple(reversed(powers)))
    assert converted(sx.default_sort_key(native,order))==sp.default_sort_key(foreign,order)
    assert converted(sx.default_sort_key({'a':native},'grlex'))==sp.default_sort_key({'a':foreign},'grlex')


def test_nondefault_keys_retain_shared_subgraphs_and_bound_names():
    x,y=sx.symbols('x y');pair=sx.Function('pair')
    value=x
    for _ in range(80):value=pair(value,value)
    key=value.sort_key('grlex')
    for _ in range(80):
        arguments=key[1][1]
        assert arguments[0] is arguments[1]
        key=arguments[0]
    f=sx.Function('f');first=sx.Subs(f(x),x,0);second=sx.Subs(f(y),y,0)
    assert first==second
    assert first.sort_key('grlex')!=second.sort_key('grlex')


def test_invalid_orders_raise_when_a_monomial_order_is_needed():
    x,y=sx.symbols('x y')
    with pytest.raises(ValueError):(x+y).sort_key('invalid')
    with pytest.raises(ValueError):((x+y)*(x*x+y)).as_ordered_factors('invalid')


def test_recursive_ordering_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y')
assert (x*x+y**3).as_ordered_terms(order='grlex')==[y**3,x*x]
assert (x*x+y**3).sort_key('grlex')[1][1][0]==(y**3).sort_key('grlex')
A=MatrixSymbol('A',2,2);B=MatrixSymbol('B',2,2)
assert (x*x*A+y**3*B).as_ordered_terms(order='grlex')==[y**3*B,x*x*A]
'''
    subprocess.run([sys.executable,'-c',source],check=True)
