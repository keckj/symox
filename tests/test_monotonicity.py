"""Native monotonicity proofs respect domains and isolated stationary points."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


FUNCTIONS=['is_increasing','is_strictly_increasing','is_decreasing',
           'is_strictly_decreasing','is_monotonic']


@pytest.mark.parametrize('name',FUNCTIONS)
@pytest.mark.parametrize('source,domain',[
    ('x','Reals'),('-x','Reals'),('x*x','Reals'),('x**4','Reals'),
    ('x**3+3*x','Reals'),('x**3-3*x*x+4*x','Reals'),
    ('exp(x)','Reals'),('exp(-x)','Reals'),
    ('x*x','Interval.open(0,oo)'),('x*x','Interval.open(-oo,0)'),
    ('1/x','Interval.open(0,oo)'),('1/x','Interval.open(-oo,0)'),
    ('log(x)','Interval.open(0,oo)'),('sqrt(x)','Interval.open(0,oo)'),
])
def test_supported_monotonicity_cases_match_reference(name,source,domain):
    x=sx.Symbol('x');X=sp.Symbol('x')
    D=eval(domain,vars(sp));d=eval(domain,vars(sx))
    actual=getattr(sx,name)(sx.sympify(source),d,x)
    expected=getattr(sp,name)(sp.sympify(source),D,X)
    assert actual is expected


def test_mathematical_strictness_allows_isolated_zero_derivatives():
    x=sx.Symbol('x')
    assert sx.is_strictly_increasing(x**3,symbol=x) is True
    assert sx.is_strictly_decreasing(-x**3,symbol=x) is True
    assert sx.is_monotonic(x**3,symbol=x) is True
    # This derivative vanishes at x=1/3 without creating a plateau.
    assert sx.is_strictly_increasing((3*x-1)**3,symbol=x) is True
    assert sx.is_monotonic(2,symbol=x) is True
    assert sx.is_increasing(2,symbol=x) is True
    assert sx.is_decreasing(2,symbol=x) is True
    assert sx.is_strictly_increasing(2,symbol=x) is False
    assert sx.is_strictly_decreasing(2,symbol=x) is False


def test_closed_endpoints_need_no_finite_derivative():
    x=sx.Symbol('x');positive=sx.Interval(0,sx.oo)
    assert sx.is_strictly_increasing(x*x,positive,x) is True
    assert sx.is_strictly_increasing(sx.sqrt(x),positive,x) is True
    assert sx.is_strictly_decreasing(-sx.sqrt(x),positive,x) is True
    assert sx.is_strictly_decreasing(x*x,sx.Interval(-1,0),x) is True


def test_function_must_be_defined_and_real_throughout_the_requested_interval():
    x=sx.Symbol('x')
    for function in FUNCTIONS:
        assert getattr(sx,function)(1/x,sx.S.Reals,x) is False
        assert getattr(sx,function)(sx.log(x),sx.Interval(0,1),x) is False
        assert getattr(sx,function)(sx.sqrt(x),sx.Interval(-1,1),x) is False
        assert getattr(sx,function)(1/x,sx.FiniteSet(0),x) is False
    assert sx.is_monotonic(1/x,sx.Interval.open(0,sx.oo),x) is True


def test_empty_and_single_point_domains_are_vacuously_monotonic():
    x=sx.Symbol('x')
    for name in FUNCTIONS:
        function=getattr(sx,name)
        assert function(x*x,sx.S.EmptySet,x) is True
        assert function(x*x,sx.FiniteSet(1),x) is True


def test_parameter_signs_and_inference_retain_uncertainty():
    x=sx.Symbol('x');a=sx.Symbol('a',real=True);p=sx.Symbol('p',positive=True)
    assert sx.is_strictly_increasing(p*x,symbol=x) is True
    assert sx.is_strictly_decreasing(-p*x,symbol=x) is True
    assert sx.is_increasing(a*x,symbol=x) is None
    assert sx.is_decreasing(a*x,symbol=x) is None
    # Reference: no turning point of a*x in the reals, hence monotonic.
    assert sx.is_monotonic(a*x,symbol=x) is True
    assert sx.is_strictly_increasing(x+a,symbol=x) is True
    assert sx.is_strictly_increasing(sx.exp(x)) is True
    positive_x=sx.Symbol('positive_x',positive=True)
    assert sx.is_monotonic(positive_x*positive_x,sx.S.Reals,positive_x) is False


def test_standard_module_aliases_and_unsupported_proofs_are_explicit():
    from symox.calculus.singularities import is_monotonic
    from symox.calculus import is_strictly_increasing
    assert is_monotonic is sx.is_monotonic
    assert is_strictly_increasing is sx.is_strictly_increasing
    x,y=sx.symbols('x y')
    with pytest.raises(NotImplementedError):sx.is_increasing(x+y)
    # Values below are those of the reference monotonicity_helper:
    # interval.is_subset(solveset(predicate(derivative), symbol, S.Reals)).
    assert sx.is_increasing(x+y,symbol=x) is True
    # The reference overflows the recursion limit here; an explicit
    # NotImplementedError is kept instead of a crash.
    with pytest.raises(NotImplementedError):sx.is_increasing(sx.Function('f')(x),symbol=x)
    assert sx.is_increasing(sx.sin(x),symbol=x) is False
    assert sx.is_increasing(x,sx.S.Complexes,x) is False
    assert sx.is_increasing(x,sx.FiniteSet(0,1),x) is True
    assert sx.is_increasing(x,sx.Union(sx.Interval(0,1),sx.Interval(2,3)),x) is True
    assert sx.is_increasing(x,sx.Interval(0,x),x) is True
    with pytest.raises(ValueError):sx.is_increasing(x,symbol=x*x)


@pytest.mark.parametrize('source,domain',[
    ('x*x','Reals'),('-x*x','Reals'),('x**4','Reals'),('x**3','Reals'),
    ('exp(x)','Reals'),('log(x)','Interval.open(0,oo)'),
    ('1/x','Interval.open(0,oo)'),('1/x','Interval.open(-oo,0)'),
    ('sqrt(x)','Interval.open(0,oo)'),('x','Reals'),('2','Reals'),
])
def test_univariate_convexity_matches_reference(source,domain):
    x=sx.Symbol('x');X=sp.Symbol('x')
    assert sx.is_convex(sx.sympify(source),x,domain=eval(domain,vars(sx))) is \
        sp.is_convex(sp.sympify(source),X,domain=eval(domain,vars(sp)))


def test_convexity_uses_the_interval_interior_and_validates_scope():
    from symox.calculus.util import is_convex
    assert is_convex is sx.is_convex
    x,y=sx.symbols('x y')
    assert is_convex(-sx.sqrt(x),x,domain=sx.Interval(0,sx.oo)) is True
    assert is_convex(x**3,x,domain=sx.Interval(0,sx.oo)) is True
    assert is_convex(1/x,x) is False
    with pytest.raises(ValueError):is_convex(x*x)
    with pytest.raises(NotImplementedError):is_convex(x*x+y*y,x,y)
    with pytest.raises(NotImplementedError):is_convex(sx.Abs(x),x)


def test_monotonicity_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x')
assert is_strictly_increasing(x**3,symbol=x) is True
assert is_strictly_increasing(sqrt(x),Interval(0,oo),x) is True
assert is_decreasing(1/x,Interval.open(0,oo),x) is True
assert is_monotonic(1/x,S.Reals,x) is False
assert is_convex(-sqrt(x),x,domain=Interval(0,oo)) is True
'''
    subprocess.run([sys.executable,'-c',source],check=True)
