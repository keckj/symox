"""Elementary extrema use complete critical sets and one-sided endpoint limits."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('function',['minimum','maximum'])
@pytest.mark.parametrize('source,domain',[
    ('x*x+2*x+5','Reals'),('-x*x+2*x+5','Reals'),
    ('x','Interval(0,1)'),('x*x','Interval(-2,3)'),('x*x','Interval.open(0,1)'),
    ('x**3-3*x','Interval(-2,2)'),('exp(x)','Reals'),
    ('exp(x)','Interval(0,1)'),('log(x)','Interval.open(0,1)'),
    ('sqrt(x)','Interval(0,4)'),('1/x','Interval(1,2)'),
    ('1/x','Reals'),('1/(x*x)','Reals'),('x+1/x','Interval.open(0,oo)'),
    ('(x*x-1)/(x*x+1)','Reals'),
    ('sin(x)','Reals'),('sin(x)*cos(x)','Reals'),('sin(x)**2','Reals'),
    ('sin(x)','Interval(2,3)'),('sin(x)','Interval(-pi,pi/4)'),
    ('exp(sin(x))','Reals'),('tan(x)','Reals'),('cot(x)','Reals'),
    ('sinh(x)','Reals'),('cosh(x)','Reals'),('tanh(x)','Reals'),('coth(x)','Reals'),
    ('Abs(x)','Reals'),('Abs(x*x-1)','Interval(-2,2)'),
    ('x*x','Union(Interval(-3,-2),Interval(1,4))'),
    ('x*x','FiniteSet(-2,1,3)'),('1/x','FiniteSet(0,1,2)'),('2','Reals'),
])
def test_supported_extrema_match_reference(function,source,domain):
    x=sx.Symbol('x');X=sp.Symbol('x',real=True)
    result=getattr(sx,function)(sx.sympify(source),x,eval(domain,vars(sx)))
    expected=getattr(sp,function)(sp.sympify(source,locals={'x':X}),X,eval(domain,vars(sp)))
    assert result.to_sympy()==expected


def test_endpoint_values_are_infima_and_suprema_even_when_not_attained():
    x=sx.Symbol('x')
    assert sx.minimum(sx.exp(x),x)==0
    assert sx.maximum(x,x,sx.Interval.open(0,1))==1
    assert sx.minimum(x*x,x,sx.Interval.open(0,1))==0
    assert sx.maximum(1/x,x)==sx.oo
    assert sx.minimum(1/x,x)==-sx.oo


def test_periodic_ranges_reduce_before_enumerating_a_huge_domain():
    x=sx.Symbol('x');domain=sx.Interval(0,10**100*sx.pi)
    assert sx.minimum(sx.sin(x),x,domain)==-1
    assert sx.maximum(sx.sin(x),x,domain)==1
    assert sx.minimum(sx.cos(x),x,sx.Interval.open(123,sx.oo))==-1
    assert sx.maximum(sx.cos(x),x,sx.Interval.open(-sx.oo,456))==1


def test_unbounded_limits_can_prove_one_extremum_without_all_critical_points():
    x=sx.Symbol('x')
    assert sx.maximum(x**4+sx.sin(x),x)==sx.oo
    assert sx.minimum(-x**4+sx.sin(x),x)==-sx.oo


def test_exact_tiny_and_large_finite_bounds_do_not_round_together():
    x=sx.Symbol('x');tiny=sx.Rational(1,10**1000)
    assert sx.minimum(x,x,sx.Interval(tiny,2*tiny))==tiny
    assert sx.maximum(x,x,sx.Interval(tiny,2*tiny))==2*tiny
    large=sx.Integer(10**100)
    assert sx.minimum(x,x,sx.Interval(large,large+1))==large
    assert sx.maximum(x,x,sx.Interval(large,large+1))==large+1


def test_parameters_and_standard_aliases():
    from symox.calculus.util import minimum,maximum
    assert minimum is sx.minimum and maximum is sx.maximum
    x=sx.Symbol('x');a=sx.Symbol('a',real=True);p=sx.Symbol('p',positive=True)
    assert sx.minimum(x*x+a*x+a,x)==a-a*a/4
    assert sx.minimum(p*x*x,x)==0
    assert sx.maximum(p*x*x,x)==sx.oo
    assert sx.minimum(2,x,sx.S.Integers)==2


def test_invalid_domains_and_unresolved_critical_sets_raise():
    x,y=sx.symbols('x y')
    for function in [sx.minimum,sx.maximum]:
        with pytest.raises(ValueError):function(x,x,sx.S.EmptySet)
        with pytest.raises(ValueError):function(x,x*x)
        with pytest.raises(ValueError):function(x,x,sx.Interval(0,x))
        with pytest.raises(ValueError):function(sx.sqrt(x),x,sx.Interval(-2,-1))
        with pytest.raises(NotImplementedError):function(sx.Function('f')(x),x)
        with pytest.raises(NotImplementedError):function(x+y,x)
        with pytest.raises(NotImplementedError):function(x,x,sx.S.Complexes)
    with pytest.raises(NotImplementedError):sx.minimum(x**4+sx.sin(x),x)


def test_function_extrema_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x')
assert minimum(x*x+2*x+5,x)==4
assert maximum(sin(x),x)==1
assert minimum(1/x,x)==-oo
assert maximum(x,x,Interval.open(0,1))==1
'''
    subprocess.run([sys.executable,'-c',source],check=True)
