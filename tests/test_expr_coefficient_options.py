"""Structural coefficients preserve dependent residual factors."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source',[
    'x*(x+1)','x*(x+y)','x*(y+1)','(x+1)**2','x+z*(x+x*y)',
    '3+2*x+4*x**2','x+sqrt(2)*x','-x-2*y','3*x*y+2*x*z+x*y*z',
    'x**3+x**2+1/x','sin(x)+x','x*sin(x)+x*cos(x)','x+y+2',
])
@pytest.mark.parametrize('target',['x','x*x','x*y','sin(x)','x+y','1','-1','2','sqrt(2)'])
@pytest.mark.parametrize('power',[-1,0,1,2])
def test_structural_coefficients_match_reference(source,target,power):
    expression=sx.sympify(source);reference=sp.sympify(source)
    factor=sx.sympify(target);F=sp.sympify(target)
    for right in (False,True):
        actual=expression.coeff(factor,power,right=right)
        assert actual.to_sympy()==reference.coeff(F,power,right=right)


def test_explicit_factor_selection_keeps_a_dependent_residual():
    x,y=sx.symbols('x y')
    assert (x*(x+1)).coeff(x)==x+1
    assert (x*(x+y)).coeff(x)==x+y
    assert ((x+1)**2).coeff(x)==0
    assert (x*(x+y)).coeff(x+y)==x
    assert (x*x).coeff(x)==0
    assert (x*x).coeff(x*x)==1


def test_numeric_coefficients_select_explicit_factors():
    x,y=sx.symbols('x y')
    assert (-x-2*y).coeff(2)==-y
    assert (x-2*y).coeff(-1)==2*y
    assert (x+2*y).coeff(1)==x
    assert (3+2*x+4*x*x).coeff(1)==0
    assert (sx.Rational(2,3)*x).coeff(sx.Rational(2,3))==x


def test_held_argument_and_repeated_factor_views_follow_reference():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    for expression,reference in [
        (sx.Add(3*x,4,evaluate=False),sp.Add(3*X,4,evaluate=False)),
        (sx.Mul(x,x,y,evaluate=False),sp.Mul(X,X,Y,evaluate=False)),
    ]:
        for target,T in [(x,X),(x*x,X*X),(sx.Integer(3),sp.Integer(3))]:
            assert expression.coeff(target).to_sympy()==reference.coeff(T)


def test_invalid_powers_and_resource_limits():
    x=sx.Symbol('x')
    for power in [True,False,None,.5,2.0,sx.Symbol('n')]:
        with pytest.raises((ValueError,TypeError)):x.coeff(x,power)
    with pytest.raises(NotImplementedError):x.coeff(2,10**9)
    assert x.coeff([x])==0
    assert x.coeff(0)==0


def test_matrix_scalar_coefficients_preserve_parent_nodes():
    x=sx.Symbol('x');A=sx.MatrixSymbol('A',2,2)
    assert (x*sx.Trace(x*A)).coeff(x)==sx.Trace(x*A)
    assert (x*A[0,1]+2*x*A[1,0]).coeff(x)==A[0,1]+2*A[1,0]


def test_coefficient_options_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y')
assert (x*(x+1)).coeff(x,right=True)==x+1
assert (-x-2*y).coeff(2)==-y
assert (x+2*y).coeff(1)==x
'''
    subprocess.run([sys.executable,'-c',source],check=True)
