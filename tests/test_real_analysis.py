"""Native singularities, continuity domains, periods and stationary points."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox import singularities,periodicity,stationary_points
from symox.calculus.util import continuous_domain


def canonical_lattices(value):
    """Normalize affine integer families by exact period-offset arithmetic."""
    index=sp.Symbol('_lattice_index',integer=True)
    def convert(image):
        if len(image.lamda.variables)!=1 or image.base_set!=sp.S.Integers:return image
        old=image.lamda.variables[0]
        period=sp.expand(image.lamda.expr).coeff(old)
        offset=sp.expand(image.lamda.expr).subs(old,0)
        if period.could_extract_minus_sign():period=-period
        ratio=sp.cancel(offset/period)
        if ratio.is_Rational:offset=period*(ratio%1)
        return sp.ImageSet(sp.Lambda(index,period*index+offset),sp.S.Integers)
    return value.replace(lambda item:isinstance(item,sp.ImageSet),convert)


@pytest.mark.parametrize('source',[
    'x**2+3*x','1/x','1/(x*x-1)','exp(1/x)','sin(1/x)','log(x)',
    'log(x*x-1)','atanh(x)','tan(x)','cot(x)','sec(x)','csc(x)',
    'tanh(x)','coth(x)','sech(x)','csch(x)','1/sqrt(x)',
])
@pytest.mark.parametrize('domain',['Reals','Complexes'])
def test_singularity_sets_match_reference(source,domain):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=sx.sympify(source);reference=sp.sympify(source)
    actual=singularities(expression,x,getattr(sx.S,domain)).to_sympy()
    expected=sp.singularities(reference,X,getattr(sp.S,domain))
    assert canonical_lattices(actual)==canonical_lattices(expected)


@pytest.mark.parametrize('source,expected',[
    ('x*x','Reals'),('1/x','Union(Interval.open(-oo,0),Interval.open(0,oo))'),
    ('sqrt(x)','Interval(0,oo)'),('1/sqrt(x)','Interval.open(0,oo)'),
    ('log(x)','Interval.open(0,oo)'),('log(x*x-1)','Union(Interval.open(-oo,-1),Interval.open(1,oo))'),
    ('asin(x)','Interval(-1,1)'),('acos(x)','Interval(-1,1)'),('acosh(x)','Interval(1,oo)'),
    ('atanh(x)','Interval.open(-1,1)'),('exp(x)','Reals'),
    ('log(x)+sqrt(1-x)','Interval.Lopen(0,1)'),
    ('x**pi','Interval(0,oo)'),
])
def test_real_continuity_domains(source,expected):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=continuous_domain(sx.sympify(source),x,sx.S.Reals)
    reference=sp.sympify(expected)
    assert actual.to_sympy()==reference
    if source!='x**pi':
        from sympy.calculus.util import continuous_domain as oracle
        assert actual.to_sympy()==oracle(sp.sympify(source),X,sp.S.Reals)


def test_poles_on_bounded_domains_and_principal_real_powers():
    x=sx.Symbol('x')
    assert singularities(sx.tan(x),x,sx.Interval(0,sx.pi))==sx.FiniteSet(sx.pi/2)
    assert continuous_domain(sx.tan(x),x,sx.Interval(0,sx.pi))==sx.Union(
        sx.Interval.Ropen(0,sx.pi/2),sx.Interval.Lopen(sx.pi/2,sx.pi))
    # Scalar Pow uses the principal complex branch for negative bases.
    assert continuous_domain(x**sx.Rational(1,3),x,sx.S.Reals)==sx.Interval(0,sx.oo)


@pytest.mark.parametrize('source,expected',[
    ('x*x','FiniteSet(0)'),('x**3-3*x','FiniteSet(-1,1)'),
    ('1/x','EmptySet'),('log(x)','EmptySet'),('exp(x)','EmptySet'),
    ('x+1/x','FiniteSet(-1,1)'),('sqrt(x)','EmptySet'),
])
def test_stationary_points_match_reference(source,expected):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=stationary_points(sx.sympify(source),x)
    assert actual.to_sympy()==sp.sympify(expected)
    assert actual.to_sympy()==sp.stationary_points(sp.sympify(source),X)


def test_stationary_points_restrict_real_branches_and_open_endpoints():
    x=sx.Symbol('x')
    assert stationary_points(sx.sin(x),x,sx.Interval(0,4*sx.pi))==sx.FiniteSet(
        sx.pi/2,3*sx.pi/2,5*sx.pi/2,7*sx.pi/2)
    assert stationary_points(x*x,x,sx.Interval.open(0,2))==sx.S.EmptySet
    assert stationary_points(sx.log(x*x-1),x)==sx.S.EmptySet
    assert stationary_points(x,x,sx.S.EmptySet)==sx.S.EmptySet
    with pytest.raises(NotImplementedError):stationary_points(sx.Abs(x),x)


@pytest.mark.parametrize('source,expected',[
    ('sin(x)','2*pi'),('cos(3*x)','2*pi/3'),('tan(2*x)','pi/2'),
    ('sin(x)**2','pi'),('cos(x)**4','pi'),('Abs(sin(x))','pi'),
    ('sin(x)+cos(2*x)','2*pi'),('sin(x/2)+sin(x/3)','12*pi'),
    ('sin(cos(x))','2*pi'),('exp(sin(x))','2*pi'),('exp(2*I*x)','pi'),
    ('sin(x)+sin(sqrt(2)*x)',None),('x*x',None),('sinh(x)',None),
    ('2','0'),('(sin(x)+cos(x))**2','pi'),('Mod(2*x,3)','3/2'),
])
def test_structural_periods_and_symbolic_verification(source,expected):
    x=sx.Symbol('x')
    result=periodicity(sx.sympify(source),x,check=True)
    if expected is None:
        assert result is None
    else:
        assert sx.simplify(result-sx.sympify(expected))==0


def test_unknown_functions_domains_and_invalid_variables_are_explicit():
    x,a=sx.symbols('x a');f=sx.Function('f')
    with pytest.raises(NotImplementedError):singularities(f(x),x)
    with pytest.raises(NotImplementedError):singularities(x**a,x)
    with pytest.raises(NotImplementedError):singularities(sx.gamma(x),x)
    with pytest.raises(NotImplementedError):continuous_domain(f(x),x,sx.S.Reals)
    with pytest.raises(NotImplementedError):continuous_domain(sx.Abs(sx.sqrt(x)),x,sx.S.Reals)
    assert continuous_domain(f(x),x,sx.S.EmptySet)==sx.S.EmptySet
    assert singularities(f(x),x,sx.S.EmptySet)==sx.S.EmptySet
    with pytest.raises(NotImplementedError):continuous_domain(x,x,sx.S.Complexes)
    with pytest.raises(ValueError):continuous_domain(x,x,sx.Interval(0,x))
    for function in [singularities,periodicity,stationary_points]:
        with pytest.raises(ValueError):function(x,x*x)
    assert periodicity(f(sx.sin(x)),x)==2*sx.pi
    assert periodicity(f(x),x) is None


def test_no_sympy_import():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols,sin,pi,S,FiniteSet,Interval,log
from symox.calculus.singularities import singularities
from symox.calculus.util import continuous_domain,periodicity,stationary_points
x=symbols('x')
assert singularities(1/(x*x-1),x,S.Reals)==FiniteSet(-1,1)
assert continuous_domain(log(x),x,S.Reals)==Interval.open(0,S.Infinity)
assert stationary_points(x*x,x)==FiniteSet(0)
assert periodicity(sin(2*x),x)==pi
'''
    subprocess.run([sys.executable,'-c',source],check=True)
