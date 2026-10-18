"""Finite-set membership is an existential condition with real domains."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('members,target',[
    ('[x/2]','Interval(0,1)'),
    ('[x,x*x]','Interval(1,2)'),
    ('[x*x/(x+2)]','Interval(1,oo)'),
    ('[1/x]','Interval(1,2)'),
    ('[x*x]','Interval.open(0,1)'),
    ('[x*x]','Interval.Lopen(0,1)'),
    ('[x*x]','Interval.Ropen(0,1)'),
    ('[x+1,2*x]','Union(Interval(-3,-1),Interval.open(2,4))'),
    ('[x/(x-1)]','Interval(-1,1)'),
    ('[x**3]','Interval(-8,27)'),
    ('[2*x+1]','Interval.open(-oo,3)'),
    ('[x*x+1]','Interval(0,1)'),
    ('[x*x+1]','Interval.open(0,1)'),
])
def test_rational_membership_matches_reference(members,target):
    x=sx.Symbol('x');X=sp.Symbol('x')
    native=vars(sx)|{'x':x};reference=vars(sp)|{'x':X}
    specification=sx.Intersection(sx.FiniteSet(*eval(members,native)),eval(target,native))
    foreign=sp.Intersection(sp.FiniteSet(*eval(members,reference)),eval(target,reference))
    assert sx.not_empty_in(specification,x).to_sympy()==sp.not_empty_in(foreign,X)


@pytest.mark.parametrize('member,expected',[
    ('x','Reals'),('1/x','Union(Interval.open(-oo,0),Interval.open(0,oo))'),
    ('log(x)','Interval.open(0,oo)'),('sqrt(x)','Interval(0,oo)'),
    ('asin(x)','Interval(-1,1)'),('atanh(x)','Interval.open(-1,1)'),
    ('floor(x)','Reals'),('ceiling(x)','Reals'),('sign(x)','Reals'),
    ('acot(x)','Reals'),('Heaviside(x)','Reals'),
    ('sign(1/x)','Union(Interval.open(-oo,0),Interval.open(0,oo))'),
    ('floor(sqrt(x))','Interval(0,oo)'),
])
def test_members_need_to_be_defined_rather_than_continuous(member,expected):
    x=sx.Symbol('x');namespace=vars(sx)|{'x':x}
    assert sx.not_empty_in(sx.FiniteSet(eval(member,namespace)),x)==eval(expected,namespace)


def test_finite_membership_targets_and_complements():
    x=sx.Symbol('x')
    specification=sx.Intersection(sx.FiniteSet(x,x*x),sx.FiniteSet(1,4))
    assert sx.not_empty_in(specification,x)==sx.FiniteSet(-2,-1,1,2,4)
    target=sx.Complement(sx.S.Reals,sx.FiniteSet(0))
    assert sx.not_empty_in(sx.Intersection(sx.FiniteSet(x),target),x)==target
    assert sx.not_empty_in(sx.Intersection(sx.FiniteSet(x),sx.FiniteSet(1/x)),x)==sx.FiniteSet(-1,1)


def test_one_defined_constant_member_is_a_sufficient_witness():
    x=sx.Symbol('x');f=sx.Function('f')
    assert sx.not_empty_in(sx.Intersection(sx.FiniteSet(0,x),sx.Interval(0,1)),x)==sx.S.Reals
    assert sx.not_empty_in(sx.Intersection(sx.FiniteSet(0,f(x)),sx.Interval(0,1)),x)==sx.S.Reals
    assert sx.not_empty_in(sx.Union(sx.FiniteSet(0),sx.Intersection(sx.FiniteSet(f(x)),sx.Interval(0,1))),x)==sx.S.Reals
    assert sx.not_empty_in(sx.S.EmptySet,x)==sx.S.EmptySet


def test_real_parameters_and_independent_intersection_order():
    x=sx.Symbol('x');a=sx.Symbol('a',real=True);p=sx.Symbol('p',positive=True)
    assert sx.not_empty_in(sx.FiniteSet(a),x)==sx.S.Reals
    assert sx.not_empty_in(sx.Intersection(sx.Interval(0,1),sx.FiniteSet(p*x)),x)==sx.Interval(0,1/p)
    assert sx.not_empty_in(sx.Intersection(sx.FiniteSet(x),sx.Interval(x-1,x+1)),x)==sx.S.Reals


def test_unsupported_proofs_and_invalid_arguments_raise():
    x,y=sx.symbols('x y')
    with pytest.raises(ValueError):sx.not_empty_in(sx.FiniteSet(x))
    with pytest.raises(ValueError):sx.not_empty_in(sx.FiniteSet(x),x*x)
    with pytest.raises(ValueError):sx.not_empty_in(sx.Interval(0,1),x)
    with pytest.raises(ValueError):sx.not_empty_in(x,x)
    with pytest.raises(NotImplementedError):sx.not_empty_in(sx.FiniteSet(x+y),x,y)
    with pytest.raises(NotImplementedError):sx.not_empty_in(sx.FiniteSet(x+y),x)
    with pytest.raises(NotImplementedError):sx.not_empty_in(sx.FiniteSet(sx.Function('f')(x)),x)
    with pytest.raises(NotImplementedError):sx.not_empty_in(sx.FiniteSet(sx.I*x),x)


def test_native_nonempty_domains_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x')
assert not_empty_in(FiniteSet(x/2).intersect(Interval(0,1)),x)==Interval(0,2)
assert not_empty_in(FiniteSet(0,x).intersect(Interval(0,1)),x)==S.Reals
assert not_empty_in(FiniteSet(log(x)),x)==Interval.open(0,oo)
assert not_empty_in(FiniteSet(floor(x)),x)==S.Reals
'''
    subprocess.run([sys.executable,'-c',source],check=True)
