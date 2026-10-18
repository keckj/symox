import itertools
import pickle
import pytest
import symox as s
import sympy as sp
from symox.sets.ordinals import omega as w
from sympy.sets.ordinals import omega as sw


@pytest.mark.parametrize('n',range(7))
def test_power_set_elements_membership_and_rewrite(n):
    values=s.FiniteSet(*range(n));p=s.PowerSet(values)
    expected=sp.PowerSet(sp.FiniteSet(*range(n)))
    assert len(p)==2**n
    assert {v.to_sympy() for v in p}==set(expected)
    assert p.rewrite(s.FiniteSet).to_sympy()==expected.rewrite(sp.FiniteSet)
    assert p.contains(s.FiniteSet(0))==(s.true if n else s.false)
    assert s.EmptySet in p
    assert pickle.loads(pickle.dumps(p))==p


def test_infinite_powers_and_disjoint_union():
    assert s.PowerSet(s.S.Integers).contains(s.FiniteSet(2,3))==s.true
    assert s.PowerSet(s.Interval(0,1)).contains(s.FiniteSet(2))==s.false
    a=s.DisjointUnion(s.FiniteSet(1,2),s.FiniteSet(2),s.EmptySet)
    b=sp.DisjointUnion(sp.FiniteSet(1,2),sp.FiniteSet(2),sp.EmptySet)
    assert len(a)==3
    assert {(int(x),int(i)) for x,i in a}==set(b)
    assert a.contains((2,0))==s.true and a.contains((2,1))==s.true
    assert a.contains((2,2))==s.false and a.contains(2)==s.false
    assert a.rewrite(s.Union).to_sympy()==b.rewrite(sp.Union)
    i=s.Symbol('i',integer=True)
    assert a.contains((1,i)).subs(i,0)==s.true


@pytest.mark.parametrize('point',[0,1,s.I,1+s.I,2+s.I,s.Rational(1,2)+s.I,-s.I])
def test_complex_rectangles(point):
    bases=s.ProductSet(s.Interval(0,1),s.Interval(0,2))
    a=s.ComplexRegion(bases);b=sp.ComplexRegion(bases.to_sympy())
    assert a.contains(point)==s.sympify(b.contains(s.sympify(point).to_sympy()))
    assert a.a_interval==s.Interval(0,1) and a.b_interval==s.Interval(0,2)
    assert a.measure==2
    assert a.to_sympy()==b


@pytest.mark.parametrize('point',[0,1,s.I,-1,-s.I,2,1+s.I])
def test_complex_polar_regions(point):
    bases=s.ProductSet(s.Interval(0,1),s.Interval(0,s.pi))
    a=s.ComplexRegion(bases,polar=True);b=sp.ComplexRegion(bases.to_sympy(),polar=True)
    assert a.contains(point)==s.sympify(b.contains(s.sympify(point).to_sympy()))
    assert a.to_sympy()==b


def test_complex_finite_regions_and_images():
    r=s.ComplexRegion(s.ProductSet(s.FiniteSet(0,1),s.FiniteSet(1,2)))
    assert set(r)=={s.I,2*s.I,1+s.I,1+2*s.I}
    x,y=s.symbols('x y')
    base=s.FiniteSet(-2,-1,1,2)
    assert s.imageset(x,x*x,base)==s.FiniteSet(1,4)
    assert s.imageset(lambda x:x*x,base)==s.FiniteSet(1,4)
    assert s.imageset(s.Lambda((x,y),x+y),base,s.FiniteSet(0,1))==s.FiniteSet(-2,-1,0,1,2,3)
    assert isinstance(s.imageset(x,2*x,s.S.Integers),s.ImageSet)
    with pytest.raises(ValueError):s.imageset(s.Lambda((x,y),x+y),base)


@pytest.mark.parametrize('a',range(7))
@pytest.mark.parametrize('b',range(7))
def test_cantor_ordinal_arithmetic(a,b):
    x=w**2*a+w*b+3;y=w*(b+1)+a
    sx=sw**2*a+sw*b+3;sy=sw*(b+1)+a
    assert x.to_sympy()==sx and y.to_sympy()==sy
    assert (x+y).to_sympy()==sx+sy
    assert (x*y).to_sympy()==sx*sy
    assert (y*x).to_sympy()==sy*sx
    assert (x<y)==(sx<sy)
    assert s.Ordinal.convert(a).to_sympy()==sp.Ordinal.convert(a)


def test_transfinite_exponents_and_validation():
    a=w**(w+1)+w**3*2+1
    b=sw**(sw+1)+sw**3*2+1
    assert a.to_sympy()==b
    assert (a+w**w).to_sympy()==b+sw**sw
    assert (a*(w+1)).to_sympy()==b*(sw+1)
    assert s.ord0+5==5 and 5+s.ord0==5
    assert w.is_limit_ordinal and (w+1).is_successor_ordinal
    with pytest.raises(ValueError):s.Ordinal(s.OmegaPower(0,1),s.OmegaPower(1,1))
    with pytest.raises(TypeError):s.OmegaPower(0,0)


def test_region_union_area_and_ordinal_hashes():
    a=s.ProductSet(s.Interval(0,2),s.Interval(0,2))
    b=s.ProductSet(s.Interval(1,3),s.Interval(1,3))
    assert s.ComplexRegion(s.Union(a,b)).measure==7
    assert s.ComplexRegion(s.ProductSet(s.Interval(1,2),s.Interval(0,s.pi)),polar=True).measure==3*s.pi/2
    assert hash(s.Ordinal.convert(3))==hash(3)
    assert hash(s.ord0)==hash(0)
