"""Weak sign bounds and lattice absorption preserve symbolic extrema."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('facts',[{'nonnegative':True},{'nonpositive':True},
    {'integer':True,'nonnegative':True},{'zero':True}])
@pytest.mark.parametrize('source',['Max(0,n)','Min(0,n)','Max(x,x+n)','Min(x,x+n)',
    'Max(-n,n)','Min(-n,n)','Max(x,x+n,x+2*n)','Min(x,x+n,x+2*n)'])
def test_weak_symbolic_bounds_match_reference(facts,source):
    n=sx.Symbol('n',**facts);N=sp.Symbol('n',**facts);x=sx.Symbol('x');X=sp.Symbol('x')
    environment=dict(vars(sx),x=x,n=n);reference=dict(vars(sp),x=X,n=N)
    actual=eval(source,environment);expected=eval(source,reference)
    if facts.get('zero'):
        assert (actual.to_sympy()-expected).is_zero is True
    else:
        assert actual.to_sympy()==expected
    assert actual.subs(n,0).to_sympy()==expected.subs(N,0)


@pytest.mark.parametrize('source',[
    'Min(x,Max(x,y))','Max(x,Min(x,y))','Min(Max(x,y),Max(x,y,z))',
    'Max(Min(x,y),Min(x,y,z))','Min(x,y,Max(x,y,z))',
    'Max(x,y,Min(x,y,z))','Min(x,Max(x,y),Max(x,z))',
    'Max(x,Min(x,y),Min(x,z))',
])
def test_lattice_absorption_matches_reference(source):
    x,y,z=sx.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    expression=eval(source,dict(vars(sx),x=x,y=y,z=z))
    expected=eval(source,dict(vars(sp),x=X,y=Y,z=Z))
    assert expression.to_sympy()==expected
    for values in [(-1,2,3),(0,0,0),(4,-3,1)]:
        assert expression.subs(dict(zip((x,y,z),values))).to_sympy()==expected.subs(dict(zip((X,Y,Z),values)))


def test_equal_candidates_do_not_all_disappear():
    x=sx.Symbol('x')
    held=sx.Add(x,0,evaluate=False)
    assert sx.Min(x,held).subs(x,3)==3
    assert sx.Max(x,held).subs(x,3)==3
    for value in [sx.Min(x,sx.Max(x,x)),sx.Max(x,sx.Min(x,x))]:
        assert value==x


def test_nonnegative_matrix_dimensions_simplify_clipped_shapes():
    n=sx.Symbol('n',integer=True,nonnegative=True)
    assert sx.Max(0,n)==n
    A=sx.MatrixSymbol('A',n,2)
    assert A[:,:].shape==(n,2)
    assert sx.Min(0,n)==0


def test_symbolic_extrema_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y');n=Symbol('n',nonnegative=True)
assert Max(0,n)==n
assert Min(0,n)==0
assert Min(x,Max(x,y))==x
'''
    subprocess.run([sys.executable,'-c',source],check=True)
