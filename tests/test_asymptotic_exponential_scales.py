"""Comparable exponential classes retain slower scales as coefficient functions."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('n',[2,3,5])
@pytest.mark.parametrize('source',[
    '1/(1+exp(-x)+exp(-x*x))',
    '(1+exp(-x*x))/(1+exp(-x))',
    '(1+exp(-x**3))/(1+exp(-x))',
])
def test_slower_coefficients_expand_before_faster_scale_terms(n,source):
    x=sx.Symbol('x',positive=True);X=sp.Symbol('x',positive=True)
    expression=eval(source,dict(vars(sx),x=x));foreign=eval(source,dict(vars(sp),x=X))
    actual=expression.aseries(x,n=n)
    expected=foreign.aseries(X,n=n)
    assert actual.to_sympy()==expected


@pytest.mark.parametrize('hir',[False,True])
def test_a_finite_sum_of_incomparable_powers_can_still_be_exact(hir):
    x=sx.Symbol('x',positive=True)
    expression=sx.exp(-x)+sx.exp(-x*x)+sx.exp(-x**3)
    assert expression.aseries(x,n=4,hir=hir)==expression


@pytest.mark.parametrize('n',[2,3,5])
def test_hierarchical_geometric_expansion_has_an_exact_remainder_certificate(n):
    x=sx.Symbol('x',positive=True);w,z=sx.symbols('w z',positive=True)
    expression=1/(1+sx.exp(-x)+sx.exp(-x*x))
    result=expression.aseries(x,n=n,hir=True)
    # Replace the distinct exponential monomials by independent formal
    # variables. The remaining rational identity certifies all coefficients.
    rules={sx.exp(-k*x*x):w**k for k in range(1,n+1)}
    rules.update({sx.exp(-k*x):z**k for k in range(1,n+1)})
    polynomial=result.removeO().xreplace(rules)
    expected=sum((-w)**k/(1+z)**(k+1) for k in range(n))
    assert sx.cancel(polynomial-expected)==0
    remainder=(-w)**n/((1+z)**n*(1+z+w))
    assert sx.cancel(1/(1+z+w)-polynomial-remainder)==0
    assert result.getO()==sx.O(sx.exp(-n*x*x),(x,sx.oo))


@pytest.mark.parametrize('hir',[False,True])
def test_same_leading_scale_keeps_the_exact_exponent_residual(hir):
    x=sx.Symbol('x',positive=True);w,z=sx.symbols('w z',positive=True)
    expression=1/(1+sx.exp(-x*x)+sx.exp(-x*x-x))
    result=expression.aseries(x,n=3,hir=hir)
    expected=1-(1+sx.exp(-x))*sx.exp(-x*x)+(1+sx.exp(-x))**2*sx.exp(-2*x*x)
    assert sx.expand(result.removeO()-expected)==0
    assert result.getO()==sx.O(sx.exp(-3*x*x),(x,sx.oo))


def test_different_leading_coefficients_and_subleading_exponents():
    x=sx.Symbol('x',positive=True)
    expression=1/(1+sx.exp(-2*x*x)+sx.exp(-x*x-x))
    result=expression.aseries(x,n=3)
    expected=1-sx.exp(-x)*sx.exp(-x*x)+(-1+sx.exp(-2*x))*sx.exp(-2*x*x)
    assert result.removeO()==expected
    assert result.getO()==sx.O(sx.exp(-3*x*x),(x,sx.oo))


def test_hierarchical_finite_rational_coefficient_is_not_reexpanded():
    x=sx.Symbol('x',positive=True)
    expression=(1+sx.exp(-x*x))/(1+sx.exp(-x))
    result=expression.aseries(x,n=3,hir=True)
    assert sx.cancel(result-expression)==0
    assert result.getO() is None


def test_scale_expansion_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x',positive=True)
result=(1/(1+exp(-x)+exp(-x*x))).aseries(x,n=3)
assert result==1-exp(-x)+exp(-2*x)+O(exp(-3*x),(x,oo))
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
