"""Function images retain endpoint attainment and gaps between components."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.calculus.util import function_range
from sympy.calculus.util import function_range as oracle


@pytest.mark.parametrize('source,domain',[
    ('x*x','Reals'),('x','Reals'),('-x','Interval(0,1)'),
    ('1/x','Reals'),('1/(x*x)','Reals'),('x+1/x','Interval.open(0,oo)'),
    ('exp(x)','Reals'),('log(x)','Reals'),('sqrt(x)','Interval(-5,9)'),
    ('x*x','Interval(-2,3)'),('x*x','Interval.open(0,1)'),
    ('x*x','Interval.Lopen(0,1)'),('x*x','Interval.Ropen(0,1)'),
    ('x**3-3*x','Interval(-2,2)'),
    ('sin(x)','Reals'),('sin(x)*cos(x)','Reals'),('sin(x)**2','Reals'),
    ('sin(x)','Interval.open(0,pi/2)'),('sin(x)','Interval(2,3)'),
    ('cos(x)','Interval.open(0,2*pi)'),('tan(x)','Reals'),('cot(x)','Reals'),
    ('exp(sin(x))','Reals'),('Abs(x)','Reals'),
    ('sinh(x)','Reals'),('cosh(x)','Reals'),('tanh(x)','Reals'),('coth(x)','Reals'),
    ('Abs(x*x-1)','Interval(-2,2)'),
    ('x*x','Union(Interval(-3,-2),Interval(1,4))'),
    ('x*x','FiniteSet(-2,1,3)'),('1/x','FiniteSet(0,1,2)'),('2','Reals'),
])
def test_function_ranges_match_real_reference(source,domain):
    x=sx.Symbol('x');X=sp.Symbol('x',real=True)
    actual=function_range(sx.sympify(source),x,eval(domain,vars(sx)))
    expected=oracle(sp.sympify(source,locals={'x':X}),X,eval(domain,vars(sp)))
    assert actual.to_sympy()==expected


def test_one_open_period_can_omit_an_extreme_value():
    x=sx.Symbol('x')
    assert function_range(sx.cos(x),x,sx.Interval.open(0,2*sx.pi))==sx.Interval.Ropen(-1,1)
    assert function_range(sx.sin(x),x,sx.Interval.open(0,2*sx.pi))==sx.Interval(-1,1)
    assert function_range(sx.cos(x),x,sx.Interval.Lopen(0,2*sx.pi))==sx.Interval(-1,1)
    assert function_range(sx.cos(x),x,sx.Interval.Ropen(0,2*sx.pi))==sx.Interval(-1,1)
    assert sx.Interval.Lopen(0,2*sx.pi).contains(2*sx.pi)==sx.true
    assert sx.cos(2*sx.pi)==1
    assert function_range(sx.cos(x),x,sx.Interval.open(0,4*sx.pi))==sx.Interval(-1,1)


def test_noncompact_images_keep_unattained_bounds_and_disconnected_gaps():
    x=sx.Symbol('x')
    assert function_range(sx.exp(x),x,sx.S.Reals)==sx.Interval.open(0,sx.oo)
    assert function_range(x*x,x,sx.S.Reals)==sx.Interval(0,sx.oo)
    assert function_range(1/x,x,sx.S.Reals)==sx.Union(sx.Interval.open(-sx.oo,0),sx.Interval.open(0,sx.oo))
    assert function_range(sx.tan(x),x,sx.S.Reals)==sx.S.Reals


def test_empty_images_constants_and_real_parameters():
    x=sx.Symbol('x');a=sx.Symbol('a',real=True);p=sx.Symbol('p',positive=True)
    assert function_range(x,x,sx.S.EmptySet)==sx.S.EmptySet
    assert function_range(sx.sqrt(x),x,sx.Interval(-2,-1))==sx.S.EmptySet
    assert function_range(2,x,sx.S.Integers)==sx.FiniteSet(2)
    assert function_range(x*x+a*x+a,x,sx.S.Reals)==sx.Interval(a-a*a/4,sx.oo)
    assert function_range(p*x,x,sx.Interval.open(0,1))==sx.Interval.open(0,p)


def test_large_periodic_domains_are_reduced_without_enumerating_all_roots():
    x=sx.Symbol('x')
    assert function_range(sx.sin(x),x,sx.Interval(0,10**100*sx.pi))==sx.Interval(-1,1)


def test_unresolved_proofs_and_invalid_domains_raise():
    x,y=sx.symbols('x y')
    with pytest.raises(ValueError):function_range(x,x*x,sx.S.Reals)
    with pytest.raises(ValueError):function_range(x,x,sx.Interval(0,x))
    with pytest.raises(NotImplementedError):function_range(sx.Function('f')(x),x,sx.S.Reals)
    with pytest.raises(NotImplementedError):function_range(x+y,x,sx.S.Reals)
    with pytest.raises(NotImplementedError):function_range(x,x,sx.S.Complexes)
    with pytest.raises(NotImplementedError):function_range(sx.sin(x)/x,x,sx.Interval(1,2))


def test_real_coefficient_proofs_reject_complex_or_undecided_offsets():
    x=sx.Symbol('x');a=sx.Symbol('a',real=True);p=sx.Symbol('p',positive=True)
    for expression in [x+sx.log(a),x+sx.sqrt(a)]:
        with pytest.raises(NotImplementedError):function_range(expression,x,sx.S.Reals)
        with pytest.raises(NotImplementedError):sx.minimum(expression,x)
        with pytest.raises(NotImplementedError):sx.is_increasing(expression,symbol=x)
    assert function_range(x+sx.log(p),x,sx.S.Reals)==sx.S.Reals
    for expression in [sx.true,sx.Eq(x,0),sx.FiniteSet(1,2)]:
        with pytest.raises(TypeError):function_range(expression,x,sx.S.Reals)
        with pytest.raises(TypeError):sx.minimum(expression,x)
        with pytest.raises(TypeError):sx.is_increasing(expression,symbol=x)


def test_function_range_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
from symox.calculus.util import function_range
x=Symbol('x')
assert function_range(exp(x),x,S.Reals)==Interval.open(0,oo)
assert function_range(cos(x),x,Interval.open(0,2*pi))==Interval.Ropen(-1,1)
assert function_range(tan(x),x,S.Reals)==S.Reals
assert function_range(tanh(x),x,S.Reals)==Interval.open(-1,1)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
