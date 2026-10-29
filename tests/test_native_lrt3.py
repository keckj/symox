"""Native general log_to_real for auxiliary polynomials of degree >= 3."""
import importlib
import re

import pytest
import sympy as sp
import symox as sx

rt=importlib.import_module('symox.integrals.rationaltools')
ref=importlib.import_module('sympy.integrals.rationaltools')

# Integrands whose logpart has a degree >= 3 auxiliary polynomial with
# rational resultant roots and rational or pure quadratic pairs.
ADMITTED=['1/(x**4+1)','1/(x**4+4)','1/(x**4+1/4)','1/(x**4-x**2+1)',
          '1/((x**2+1)*(x**2+4)*(x-1))','x**4/(x**6+1)','1/(x**6+1)','x**2/(x**4+1)',
          '(x**3+2*x+5)/((x**2+3)*(x**2+2*x+2))']
# Cube roots, Cardano roots, nested fields and incomplete real root lists
# (where the exact RootSum is kept) stay on the general path.
DECLINED=['1/(x**3-2)','1/(x**3+x+1)','1/(x**5+1)','1/(x**6-2)','1/(x**8+1)',
          'x**2/(x**4+x+1)','1/(x**3-3*x+1)']


def groups(source):
    x,t=sx.symbols('x t');f=sx.sympify(source,locals={'x':x})
    n,d=sx.fraction(sx.cancel(f));_,r=sx.Poly(n,x).div(sx.Poly(d,x))
    pairs=rt.ratint_logpart(r,sx.Poly(d,x),x,t)
    return x,t,[(h,q) for h,q in pairs if q.degree()>=3]


def complete(h,q,x,t,monkeypatch):
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_log_to_real',lambda *args:None)
        return rt.log_to_real(h,q,x,t)


def normalize(e):
    return sp.expand(e.replace(lambda f:f.func in (sp.log,sp.atan),
        lambda f:f.func(sp.cancel(f.args[0],extension=True))))


def trace_certificate(value,h,q):
    # d/dx log_to_real(h, q) is the trace of t*h_x/h modulo q: no roots.
    X,T=sp.symbols('x t');K=sp.QQ.frac_field(X)
    hs,qs=h.as_expr().to_sympy(),q.as_expr().to_sympy()
    qpoly=sp.Poly(qs,T,domain=K);hpoly=sp.Poly(hs,T,domain=K)
    r=(sp.Poly(T*sp.diff(hs,X),T,domain=K)*sp.invert(hpoly,qpoly)).rem(qpoly)
    trace=sum((r*sp.Poly(T**j,T,domain=K)).rem(qpoly).nth(j) for j in range(qpoly.degree()))
    return sp.cancel(sp.diff(value.to_sympy(),X)-trace,extension=True)==0


@pytest.mark.parametrize('source',ADMITTED)
def test_admitted_groups_match_complete_path_and_reference(source,monkeypatch):
    x,t,pairs=groups(source);X,T=sp.symbols('x t')
    assert pairs
    for h,q in pairs:
        direct=sx._core._log_to_real(h,q,x,t)
        assert direct is not None
        assert rt.log_to_real(h,q,x,t)==direct
        assert complete(h,q,x,t,monkeypatch)==direct
        expected=ref.log_to_real(sp.Poly(h.as_expr().to_sympy(),X,domain='QQ[t]'),
                                 sp.Poly(q.as_expr().to_sympy(),T,domain='QQ'),X,T)
        assert normalize(direct.to_sympy())==normalize(expected)
        assert trace_certificate(direct,h,q)


@pytest.mark.parametrize('source',DECLINED)
def test_undecided_groups_decline_to_complete_path(source):
    x,t,pairs=groups(source)
    assert pairs
    for h,q in pairs:
        assert sx._core._log_to_real(h,q,x,t) is None


@pytest.mark.parametrize('source,rootsum',[(s,False) for s in ADMITTED]+[
    ('1/(x**3-2)',False),('1/(x**5+1)',False),('1/(x**6-2)',True),('1/(x**8+1)',True)])
def test_public_ratint_matches_complete_path_and_is_an_antiderivative(source,rootsum,monkeypatch):
    x=sx.Symbol('x');f=sx.sympify(source,locals={'x':x})
    value=rt.ratint(f,x)
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_log_to_real',lambda *args:None)
        patch.setattr(sx._core,'_ratint_driver',lambda *args:None)
        old=rt.ratint(f,x)
    strip=lambda e:re.sub(r'dummy_index=\d+','',sx.srepr(e))
    assert strip(value)==strip(old)
    assert value.has(sx.RootSum)==rootsum
    derivative=sx.diff(value,x)-f
    for point in (sx.Rational(3,7),sx.Rational(-5,3),sx.Rational(11,4)):
        assert abs(complex(derivative.subs(x,point).evalf(30)))<1e-20
