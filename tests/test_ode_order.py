"""Differential order uses structural derivative counts without unfolding."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def expression(module,source):
    x,y=module.symbols('x y');f,g=module.Function('f'),module.Function('g')
    return module.sympify(eval(source,vars(module)|{'x':x,'y':y,'f':f,'g':g}))


@pytest.mark.parametrize('function',['f(x)','g(x)'])
@pytest.mark.parametrize('source',[
    '0','f(x)','g(x)','sin(f(x))',
    'Derivative(f(x),x,2,evaluate=False)+Derivative(f(x),x,evaluate=False)**2',
    'Derivative(f(x),x,2,evaluate=False)+Derivative(g(x),x,3,evaluate=False)',
    'Derivative(f(x)**2,x,3,evaluate=False)',
    'Derivative(sin(f(x)),x,2,evaluate=False)',
    'Derivative(Derivative(f(x),x,2,evaluate=False)**2,x,3,evaluate=False)',
    'Eq(Derivative(f(x),x,4,evaluate=False),g(x))',
    'Derivative(g(f(x)),x,2,evaluate=False)',
    'Derivative(f(x)+g(x),x,2,evaluate=False)',
    'Derivative(f(x,y),x,2,y,3,evaluate=False)',
    'exp(Derivative(f(x),x,2,evaluate=False))',
    'Piecewise((Derivative(f(x),x,3,evaluate=False),x>0),(0,True))',
])
def test_held_differential_orders_match_reference(source,function):
    native=expression(sx,source);foreign=expression(sp,source)
    assert sx.ode_order(native,expression(sx,function))==sp.ode_order(foreign,expression(sp,function))


def test_total_mixed_order_and_independent_function_orders():
    x,y=sx.symbols('x y');f,g=sx.Function('f'),sx.Function('g')
    equation=sx.Derivative(f(x,y),x,2,y,3,evaluate=False)+sx.Derivative(g(x,y),x,8,evaluate=False)
    assert sx.ode_order(equation,f(x,y))==5
    assert sx.ode_order(equation,g(x,y))==8


def test_huge_held_counts_are_not_enumerated():
    x=sx.Symbol('x');f=sx.Function('f')
    derivative=sx.Derivative(f(x),x,2**32-1,evaluate=False).diff(x,evaluate=False)
    assert sx.ode_order(derivative,f(x))==2**32


def test_repeated_subgraphs_are_not_unfolded():
    x=sx.Symbol('x');f,pair=sx.Function('f'),sx.Function('pair')
    expression=sx.Derivative(f(x),x,3,evaluate=False)
    for _ in range(80):expression=pair(expression,expression)
    assert sx.ode_order(expression,f(x))==3


def test_invalid_nonscalar_equations_and_dependent_functions_raise():
    x=sx.Symbol('x')
    with pytest.raises(TypeError):sx.ode_order(sx.Matrix([[x]]),x)
    with pytest.raises(TypeError):sx.ode_order(sx.FiniteSet(x),x)
    with pytest.raises(TypeError):sx.ode_order(x,sx.Matrix([[x]]))


def test_differential_order_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');f=Function('f');g=Function('g')
equation=Derivative(f(x),x,2,evaluate=False)+Derivative(g(x),x,3,evaluate=False)
assert ode_order(equation,f(x))==2
assert ode_order(equation,g(x))==3
'''
    subprocess.run([sys.executable,'-c',source],check=True)
