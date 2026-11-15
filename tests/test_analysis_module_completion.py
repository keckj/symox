"""Integral common multiples and derivative-predicate interval certificates."""
import importlib
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.calculus.util import lcim
from symox.calculus.singularities import monotonicity_helper
from sympy.calculus.util import lcim as reference_lcim
from sympy.calculus.singularities import monotonicity_helper as reference_helper


@pytest.mark.parametrize('source',[
    '[Rational(1,2),Rational(3,4),Rational(5,6)]',
    '[2*pi,3*pi,pi,pi/2]','[sqrt(2),sqrt(8),3*sqrt(2)/2]',
    '[2,3,5]','[-2,-3]','[Rational(-1,2),Rational(3,4)]',
    '[1,2*pi]','[sqrt(2),sqrt(3)]','[0,2,3]','[0]',
])
def test_lcim_matches_reference_and_has_integral_quotients(source):
    values=eval(source,vars(sx));reference=eval(source,vars(sp))
    actual=lcim(values);expected=reference_lcim(list(map(sp.sympify,reference)))
    if expected is None:assert actual is None
    else:
        assert sp.simplify(actual.to_sympy()-expected)==0
        if actual!=0:
            assert all(sx.simplify(actual/value).is_integer is True for value in values)
    assert lcim([])==1


@pytest.mark.parametrize('source,domain',[
    ('x*x','Reals'),('x*x','Interval(0,2)'),('x*x','Interval.open(0,2)'),
    ('x**3','Reals'),('-x','Reals'),('exp(x)','Reals'),
    ('log(x)','Interval.open(0,oo)'),('sqrt(x)','Interval.open(0,oo)'),
    ('7','Reals'),
])
@pytest.mark.parametrize('predicate',['lambda d:d>=0','lambda d:d>0','lambda d:d<=0'])
def test_derivative_predicates_match_reference_on_supported_domains(source,domain,predicate):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=monotonicity_helper(sx.sympify(source),eval(predicate),eval(domain,vars(sx)),x)
    expected=reference_helper(sp.sympify(source),eval(predicate),eval(domain,vars(sp)),X)
    if source=='x**3' and predicate=='lambda d:d>0' and expected is None:
        # SymPy's set inclusion can remain undecided despite the witness x=0.
        assert actual is False
        assert sp.diff(X**3,X).subs(X,0)==0
    else:assert actual is expected


@pytest.mark.parametrize('source,domain',[
    ('x*x','Union(Interval(-2,-1),Interval(1,2))'),
    ('x','Union(Interval(-2,-1),Interval(1,2))'),
    ('1/x','Union(Interval(-2,-1),Interval(1,2))'),
    ('exp(x)','Union(Interval(-2,-1),Interval(1,2))'),
    ('x**3-3*x','Union(Interval(-2,-1),Interval(1,2))'),
    ('x','FiniteSet(1,2)'),('x*x','Complement(Reals,FiniteSet(0))'),
    ('x','Complexes'),('x','Naturals'),('x','Interval(0,x)'),
    ('x+y','Reals'),('Abs(x)','Interval(0,1)'),
])
@pytest.mark.parametrize('predicate',['lambda d:d>=0','lambda d:d>0','lambda d:d<0'])
def test_derivative_predicates_fall_back_to_reference_off_native_domains(source,domain,predicate):
    # Disconnected or non-interval domains, bounds depending on the variable,
    # parameters not known real and unsupported functions give SymPy's answer
    # instead of an internal exception.
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=monotonicity_helper(sx.sympify(source),eval(predicate),eval(domain,{**vars(sx),'x':x}),x)
    expected=reference_helper(sp.sympify(source),eval(predicate),eval(domain,{**vars(sp),'x':X}),X)
    assert actual is expected


def test_multivariate_without_symbol_raises_reference_error():
    x,y=sx.symbols('x y')
    with pytest.raises(NotImplementedError,match='multivariate'):
        monotonicity_helper(x+y,lambda d:d>=0)
    assert monotonicity_helper(sx.Integer(3),lambda d:d>0) is False
    # SymPy overflows the recursion limit on undefined functions; the native
    # NotImplementedError is kept rather than deferring to it.
    with pytest.raises(NotImplementedError):
        monotonicity_helper(sx.Function('f')(x),lambda d:d>=0,symbol=x)


def test_singletons_empty_domains_uncertainty_and_validation():
    x,a=sx.symbols('x a',real=True)
    assert monotonicity_helper(x*x,lambda d:d>0,sx.FiniteSet(0),x) is False
    assert monotonicity_helper(x*x,lambda d:d>0,sx.S.EmptySet,x) is True
    assert monotonicity_helper(a*x,lambda d:d>0,symbol=x) is None
    assert lcim([sx.I,1]) is None
    with pytest.raises(NotImplementedError):monotonicity_helper(x+a,lambda d:d>0)
    with pytest.raises(ValueError):monotonicity_helper(x,lambda d:d>0,symbol=x*x)


def test_compound_derivative_predicates_use_exact_set_operations():
    x=sx.Symbol('x')
    assert monotonicity_helper(x*x,lambda d:sx.And(d>=0,d<=4),sx.Interval(0,2),x) is True
    assert monotonicity_helper(x*x,lambda d:sx.Or(d<0,d>0),sx.Interval(0,2),x) is False
    assert monotonicity_helper(x*x,lambda d:sx.Or(d<0,d>0),sx.Interval.open(0,2),x) is True
    assert monotonicity_helper(x*x,lambda d:sx.Eq(d,0),sx.S.Reals,x) is False


def test_analysis_completion_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
from symox.calculus.util import lcim
from symox.calculus.singularities import monotonicity_helper
x=Symbol('x')
assert lcim([2*pi,3*pi,pi/2])==6*pi
assert monotonicity_helper(x*x,lambda d:d>=0,Interval(0,2),x) is True
'''
    subprocess.run([sys.executable,'-c',source],check=True)
