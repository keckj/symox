import math
import mpmath as mp
import pytest
import symox as s
import sympy as sp
from symox.functions.special.bessel import (BesselBase, SphericalBesselBase,
 SphericalHankelBase, AiryBase, assume_integer_order)


@pytest.mark.parametrize('kind',[1,2])
@pytest.mark.parametrize('n',range(5))
def test_spherical_hankel_rules_and_values(kind,n):
    z=s.Symbol('z',positive=True);f=getattr(s,'hn%d'%kind)(n,z)
    assert isinstance(f,BesselBase) and isinstance(f,SphericalBesselBase) and isinstance(f,SphericalHankelBase)
    assert f.order==n and f.argument==z
    sign=1 if kind==1 else -1
    expanded=s.expand_func(f)
    assert sp.simplify(expanded.to_sympy()-sp.expand_func(getattr(sp,'hn%d'%kind)(n,z.to_sympy())))==0
    assert f.diff(z)==getattr(s,'hn%d'%kind)(n-1,z)-(n+1)*f/z
    assert f.fdiff()==f.diff(z)
    assert s.conjugate(f)==getattr(s,'hn%d'%(3-kind))(n,z)
    assert sp.simplify(f.rewrite(s.besselj).to_sympy()-getattr(sp,'hn%d'%kind)(n,z.to_sympy()).rewrite(sp.besselj))==0
    for x in (0.5,2.5,7.0):
        actual=s._core._numeric('hn%d'%kind,float(n),x)
        want=complex(mp.sqrt(mp.pi/(2*x))*(mp.hankel1 if kind==1 else mp.hankel2)(n+0.5,x))
        assert abs(actual-want)<2e-12*max(1,abs(want))


@pytest.mark.parametrize('n',[0,1,2,5,20,64,128])
def test_native_spherical_bessel_zeros(n):
    roots=s.jn_zeros(n,4)
    expected=[float(mp.besseljzero(n+0.5,k)) for k in range(1,5)]
    assert all(float(b)>float(a) for a,b in zip(roots,roots[1:]))
    assert max(abs(float(a)-b) for a,b in zip(roots,expected))<2e-11
    for zero in roots:
        assert abs(s._core._numeric('jn',float(n),float(zero)))<2e-12


@pytest.mark.parametrize('m,a,b',[(0,0.5,1.2),(1,0.5,1.2),(2,2.0,3.0),(5,5.0,1.0),(2.5,1.3,2.2),(1,20.0,20.1)])
def test_marcum_numeric_probability(m,a,b):
    actual=s._core._numeric('marcumq',m,a,b)
    with mp.workdps(40):
        aa,bb,mm=map(mp.mpf,(a,b,m))
        integrand=lambda t:aa**(1-mm)*t**mm*mp.exp(-(aa*aa+t*t)/2)*mp.besseli(mm-1,aa*t)
        expected=mp.quad(integrand,[bb,max(aa,bb)+10,mp.inf])
    assert actual is not None and abs(actual.real-float(expected))<3e-12


def test_marcum_symbolic_rules_and_base_protocols():
    m,a,b,z=s.symbols('m a b z')
    f=s.marcumq(m,a,b)
    assert f.diff(a)==a*(s.marcumq(m+1,a,b)-f)
    assert f.diff(b)==-a**(1-m)*b**m*s.exp(-(a*a+b*b)/2)*s.besseli(m-1,a*b)
    for args in [(0,0,0),(0,a,0),(1,a,a),(2,a,a),(m,0,b)]:
        actual=s.marcumq(*args).to_sympy()
        expected=sp.marcumq(*(v.to_sympy() if hasattr(v,'to_sympy') else v for v in args))
        assert sp.simplify(actual-expected)==0
    assert isinstance(s.besselj(m,z),BesselBase)
    assert isinstance(s.airyai(z),AiryBase)
    class Helper:
        @assume_integer_order
        def action(self,nu,z):return nu+z
    assert Helper().action(s.Integer(2),z)==2+z
    assert Helper().action(s.Rational(1,2),z) is None


def test_marcum_integral_series_and_finite_bessel_rewrites():
    m,a,b,t,k=s.symbols('m a b t k')
    f=s.marcumq(m,a,b)
    expected=a**(1-m)*s.Integral(t**m*s.exp(-(t*t+a*a)/2)*s.besseli(m-1,a*t),(t,b,s.oo))
    assert f.rewrite(s.Integral,x=t)==expected
    expected=s.exp(-(a*a+b*b)/2)*s.Sum((a/b)**k*s.besseli(k,a*b),(k,1-m,s.oo))
    assert f.rewrite(s.Sum,k=k)==expected
    assert sp.simplify(s.marcumq(4,a,a).rewrite(s.besseli).to_sympy()-sp.marcumq(4,a.to_sympy(),a.to_sympy()).rewrite(sp.besseli))==0
    n,z=s.symbols('n z')
    # The defining relation, correcting SymPy 1.14's direct-rewrite order.
    rewritten=s.hn1(n,z).rewrite(s.hankel1)
    assert rewritten==s.sqrt(s.pi/(2*z))*s.hankel1(n+s.Rational(1,2),z)
