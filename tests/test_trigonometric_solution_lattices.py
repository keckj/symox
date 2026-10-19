"""Compare entire periodic solution families, including complex targets."""
import math
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


def lattice_cosets(value):
    """Convert rational-pi lattices to cosets of one common period."""
    if value is sp.S.EmptySet:return []
    def inverses(expr):
        return expr.replace(sp.asec,lambda z:sp.acos(1/z)).replace(
            sp.acsc,lambda z:sp.asin(1/z)).replace(sp.acot,lambda z:sp.atan(1/z))
    families=list(value.args) if isinstance(value,sp.Union) else [value]
    lattices=[]
    for family in families:
        assert isinstance(family,sp.ImageSet),family
        index=family.lamda.variables[0]
        expression=family.lamda.expr
        period=sp.simplify(sp.diff(expression,index)/sp.pi)
        center=sp.simplify(inverses(expression.subs(index,0)))
        assert period.is_Rational and period!=0
        lattices.append((abs(period),center))
    return lattices


def assert_same_lattices(actual,expected):
    native,reference=lattice_cosets(actual),lattice_cosets(expected)
    if not native or not reference:
        assert native==reference;return
    periods=[p for p,_ in native+reference]
    common=sp.Rational(math.lcm(*(int(p.p) for p in periods)),
                       math.gcd(*(int(p.q) for p in periods)))
    def expand(lattices):
        return [sp.expand(center+i*period*sp.pi) for period,center in lattices
                for i in range(int(common/period))]
    ncosets,rcosets=expand(native),expand(reference)
    def same(a,b):return sp.simplify((a-b)/(common*sp.pi)).is_integer is True
    assert all(any(same(a,b) for b in rcosets) for a in ncosets),(actual,expected)
    assert all(any(same(a,b) for a in ncosets) for b in rcosets),(actual,expected)


@pytest.mark.parametrize('head',['sin','cos','tan','cot','sec','csc'])
@pytest.mark.parametrize('target',[-2,-1,0,1,2,sx.I,1+sx.I,sx.Rational(1,3)])
def test_complex_trigonometric_families_match_reference_modulo_integer_shifts(head,target):
    if (head in ['sec','csc'] and target==0) or (head in ['tan','cot'] and target==sx.I):
        pytest.skip('reference returns nan/infinite pseudo-roots; independently checked below')
    x=sx.Symbol('x');X=sp.Symbol('x');value=sx.sympify(target)
    actual=sx.solveset(getattr(sx,head)(2*x+1)-value,x,sx.S.Complexes).to_sympy()
    expected=sp.solveset(getattr(sp,head)(2*X+1)-value.to_sympy(),X,sp.S.Complexes)
    assert_same_lattices(actual,expected)


@pytest.mark.parametrize('head',['sin','cos'])
@pytest.mark.parametrize('digits',[30,100,1000])
@pytest.mark.parametrize('sign',[-1,1])
def test_exact_real_boundaries_have_no_machine_epsilon_band(head,digits,sign):
    x=sx.Symbol('x');denominator=10**digits
    outside=sign*(sx.Integer(1)+sx.Rational(1,denominator))
    assert sx.solveset(getattr(sx,head)(x)-outside,x,sx.S.Reals)==sx.S.EmptySet
    inside=sign*(sx.Integer(1)-sx.Rational(1,denominator))
    roots=sx.solveset(getattr(sx,head)(x)-inside,x,sx.S.Reals)
    assert roots!=sx.S.EmptySet and not isinstance(roots,sx.ConditionSet)


@pytest.mark.parametrize('head',['sec','csc'])
@pytest.mark.parametrize('domain',['Reals','Complexes'])
def test_reciprocal_trigonometric_functions_never_equal_zero(head,domain):
    # 1/cos(z), 1/sin(z) have poles but no zeros at finite z.
    x=sx.Symbol('x')
    assert sx.solveset(getattr(sx,head)(x),x,getattr(sx.S,domain))==sx.S.EmptySet


@pytest.mark.parametrize('head',['tan','cot'])
@pytest.mark.parametrize('sign',[-1,1])
def test_tangent_omitted_values_are_not_infinite_roots(head,sign):
    # tan(z)=-i*(exp(2iz)-1)/(exp(2iz)+1) omits ±i since exp never vanishes.
    x=sx.Symbol('x')
    assert sx.solveset(getattr(sx,head)(x)-sign*sx.I,x,sx.S.Complexes)==sx.S.EmptySet


@pytest.mark.parametrize('head',['sin','cos','tan','cot','sec','csc'])
@pytest.mark.parametrize('slope',['I','1+I','-2'])
def test_complex_affine_argument_has_complete_transformed_lattices(head,slope):
    x=sx.Symbol('x')
    a=sx.sympify(slope);b=sx.Rational(1,3)+sx.I
    roots=sx.solveset(getattr(sx,head)(a*x+b)-2,x,sx.S.Complexes).to_sympy()
    expected=sp.solveset(getattr(sp,head)(sp.Symbol('x'))-2,sp.Symbol('x'),sp.S.Complexes)
    transformed=roots.replace(lambda e:isinstance(e,sp.ImageSet),
        lambda e:sp.ImageSet(sp.Lambda(e.lamda.variables[0],a.to_sympy()*e.lamda.expr+b.to_sympy()),e.base_set))
    assert_same_lattices(transformed,expected)


def test_real_domain_does_not_imply_a_complex_affine_argument_is_real():
    x=sx.Symbol('x')
    roots=sx.solveset(sx.sin(2*x+sx.asin(2))-2,x,sx.S.Reals)
    assert roots!=sx.S.EmptySet
    for index in [-100,-2,0,3,100]:assert roots.contains(index*sx.pi)==sx.true
    periodic=sx.solveset(sx.exp(sx.I*x)-1,x,sx.S.Reals).to_sympy()
    expected=sp.ImageSet(sp.Lambda(sp.Symbol('n',integer=True),2*sp.pi*sp.Symbol('n',integer=True)),sp.S.Integers)
    assert_same_lattices(periodic,expected)
    assert sx.solveset(sx.exp((1+sx.I)*x)-1,x,sx.S.Reals)==sx.FiniteSet(0)
    hyperbolic=sx.solveset(sx.sinh(sx.I*x+1)-sx.sinh(1),x,sx.S.Reals).to_sympy()
    assert_same_lattices(hyperbolic,expected)


def test_real_exponential_does_not_accept_a_tiny_imaginary_target():
    x=sx.Symbol('x')
    for digits in [30,1000]:
        target=1+sx.I*sx.Rational(1,10**digits)
        assert sx.solveset(sx.exp(x)-target,x,sx.S.Reals)==sx.S.EmptySet


def test_trigonometric_solving_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x')
assert not isinstance(solveset(sin(x)-2,x,S.Complexes),ConditionSet)
assert solveset(sec(x),x,S.Complexes)==S.EmptySet
assert solveset(tan(x)-I,x,S.Complexes)==S.EmptySet
assert solveset(sin(x)-1-Rational(1,10**1000),x,S.Reals)==S.EmptySet
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=20)
