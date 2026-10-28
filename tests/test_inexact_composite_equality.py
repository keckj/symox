"""Typed inexact parameter elements retain equality and immutable hash laws."""
import copy
import pickle
from fractions import Fraction as F

import pytest
import symox as sx
import sympy as sp

P=(1,53,139,140,211)
t,u=sx.symbols('t u')
rt,ru=sp.symbols('t u')


def dyad(value):
    if value.is_Float:
        sign,man,exp,_=value._mpf_
        return F((-1 if sign else 1)*man)*F(2)**exp
    return F(int(value.p),int(value.q))


def domains(p,cc,field,variables=(t,),ref_variables=(rt,)):
    g=sx.ComplexField(p) if cc else sx.RealField(p)
    rg=sp.ComplexField(p) if cc else sp.RealField(p)
    return (g.frac_field(*variables) if field else g.poly_ring(*variables),
            rg.frac_field(*ref_variables) if field else rg.poly_ring(*ref_variables),g,rg)


@pytest.mark.parametrize('p',P)
@pytest.mark.parametrize('cc',(False,True))
@pytest.mark.parametrize('field,variant',[(False,'uni'),(False,'multi'),(True,'uni'),(True,'multi'),(True,'fraction')])
def test_same_binary_inputs_parent_values_and_reference_equalities(p,cc,field,variant):
    variables=(t,u) if variant=='multi' else (t,)
    ref_variables=(rt,ru) if variant=='multi' else (rt,)
    k,r,g,rg=domains(p,cc,field,variables,ref_variables)
    c=g.convert(sx.Rational(9,8)+(sx.I*sx.Rational(3,8) if cc else 0))
    re=dyad(sx.re(c.as_expr()));im=dyad(sx.im(c.as_expr()))
    rc=sp.Rational(re.numerator,re.denominator)+sp.I*sp.Rational(im.numerator,im.denominator)
    # The reference receives precisely the selected native stored dyades.
    reference_coefficient=rg.convert(rc);coefficient_expr=rg.to_sympy(reference_coefficient)
    assert dyad(sp.re(coefficient_expr))==re and dyad(sp.im(coefficient_expr))==im
    source=c.as_expr()*t+2;reference_source=rc*rt+2
    if variant=='multi':source+=u*(t+1);reference_source+=ru*(rt+1)
    if variant=='fraction':source/=(t+1);reference_source/=(rt+1)
    a=k.convert(source);ar=r.convert(reference_source)
    assert a.parent()==k and ar.parent()==r
    values=[a,k.convert(a),k.convert(a.as_expr()),k.from_sympy(a.as_expr()),a/k.one,a+k.zero,
            copy.copy(a),copy.deepcopy(a),pickle.loads(pickle.dumps(a))]
    references=[ar,r.convert(ar),r.convert(ar.as_expr()),r.from_sympy(ar.as_expr()),ar/r.one,ar+r.zero,
                copy.copy(ar),ar,ar]
    for actual,reference in zip(values,references):
        assert actual.parent()==k and reference.parent()==r
        assert (actual==a)==(reference==ar) is True
        assert (a==actual)==(ar==reference) is True
        assert (actual!=a)==(reference!=ar) is False
        assert actual.as_expr()==a.as_expr()
        assert hash(actual)==hash(a)
    assert len(set(values))==1
    assert {a:'first',values[-1]:'last'}=={a:'last'}
    assert {a:19}[values[4]]==19


@pytest.mark.parametrize('p',P)
@pytest.mark.parametrize('cc',(False,True))
@pytest.mark.parametrize('field',(False,True))
def test_zero_one_and_cross_parent_constants_keep_scalar_rules(p,cc,field):
    k,r,_,_=domains(p,cc,field)
    for value in (0,1,-1,2):
        a=k.convert(value);reference=r.convert(value)
        assert (a==value)==(reference==value) is True
        assert (value==a)==(value==reference) is True
        assert hash(a)==hash(value)
        assert bool(a)==bool(reference)==bool(value)
        assert (a!=value) is False
    assert k.zero+k.zero==k.zero
    assert k.one*k.one==k.one
    assert k.one/k.one==k.one
    assert hash(k.zero)==hash(0) and hash(k.one)==hash(1)
    other,ref_other,_,_=domains(53,not cc,not field,(u,),(ru,))
    assert (k.one==other.one)==(r.one==ref_other.one) is True
    assert (k.zero==other.zero)==(r.zero==ref_other.zero) is True


@pytest.mark.parametrize('p',P)
@pytest.mark.parametrize('cc',(False,True))
@pytest.mark.parametrize('field',(False,True))
def test_distinct_coefficient_dyades_remain_unequal_at_same_domain(p,cc,field):
    k,r,g,rg=domains(p,cc,field)
    next_value=sx.Rational(2**(p-1)+1,2**(p-1))
    a=k.convert(t+1);b=k.convert(next_value*t+1)
    ar=r.convert(rt+1);br=r.convert(sp.Rational(2**(p-1)+1,2**(p-1))*rt+1)
    assert (a==b)==(ar==br) is False
    assert (a!=b)==(ar!=br) is True
    assert len({a,b})==2
    assert {f._prec for f in a.as_expr().atoms(sx.Float)}=={f._prec for f in ar.as_expr().atoms(sp.Float)}
    assert {f._prec for f in b.as_expr().atoms(sx.Float)}=={f._prec for f in br.as_expr().atoms(sp.Float)}
    assert all(f._mpf_[3]<=p for f in a.as_expr().atoms(sx.Float))
    assert all(f._mpf_[3]<=p for f in b.as_expr().atoms(sx.Float))


@pytest.mark.parametrize('field',(False,True))
@pytest.mark.parametrize('difference',('precision','generators','ground','order'))
def test_nonconstant_other_parents_remain_distinct(field,difference):
    k,r,_,_=domains(139,True,field,(t,u),(rt,ru))
    g=sx.ComplexField(140 if difference=='precision' else 139)
    rg=sp.ComplexField(140 if difference=='precision' else 139)
    if difference=='ground':g=sx.RealField(139);rg=sp.RealField(139)
    vs,rs=((u,t),(ru,rt)) if difference=='generators' else ((t,u),(rt,ru))
    order='grlex' if difference=='order' else 'lex'
    b=(g.frac_field(*vs,order=order) if field else g.poly_ring(*vs,order=order)).convert(t+u+1)
    br=(rg.frac_field(*rs,order=order) if field else rg.poly_ring(*rs,order=order)).convert(rt+ru+1)
    a=k.convert(t+u+1);ar=r.convert(rt+ru+1)
    assert a.parent()!=b.parent() and ar.parent()!=br.parent()
    assert (a==b)==(ar==br) is False
    assert (b==a)==(br==ar) is False


@pytest.mark.parametrize('p',P)
def test_ex_zero_payload_precision_and_hash_are_untouched(p):
    z=sx.EX.convert(sx.Float(0,precision=p))
    e=z.as_expr()
    assert e.is_Float and e._prec==p and e._mpf_[1]==0
    assert z==z and z==sx.EX.convert(e)
    assert hash(z)==hash(sx.EX.convert(e))
    assert z.as_expr()._prec==p


def test_typed_equality_does_not_call_public_parent_or_expression_overrides(monkeypatch):
    k=sx.ComplexField(140).poly_ring(t);a=k.convert((1+sx.I)*t+2);b=a/k.one
    calls=[]
    def forbidden(*args):
        calls.append(1);raise AssertionError('typed equality must not call public metadata methods')
    monkeypatch.setattr(sx._core.DomainElement,'as_expr',forbidden)
    monkeypatch.setattr(sx._core.DomainElement,'parent',forbidden)
    assert a==b and not a!=b and hash(a)==hash(b)
    assert calls==[]


def test_custom_equality_and_foreign_sympify_callbacks_are_not_bypassed(monkeypatch):
    cls=sx._core.DomainElement;native_eq=cls.__eq__;calls=[]
    k=sx.ComplexField(140).poly_ring(t);a=k.convert(t+1)
    def equality(self,other):
        calls.append('eq');return native_eq(self,other)
    monkeypatch.setattr(cls,'__eq__',equality)
    assert a==a and calls==['eq']
    class Foreign:
        def _sympy_(self):
            calls.append('sympify');return t+1
    assert (a==Foreign()) is False and calls[-2:]==['eq','sympify']
    class Broken:
        def _sympy_(self):
            calls.append('broken');raise ValueError('foreign conversion')
    assert (a==Broken()) is False and calls[-2:]==['eq','broken']


@pytest.mark.parametrize('field',(False,True))
def test_ordering_and_failed_conversions_keep_errors(field):
    k,r,_,_=domains(140,True,field)
    a=k.convert(t+1);ar=r.convert(rt+1)
    with pytest.raises(TypeError):a<a
    assert (ar<ar) is False  # Existing native ordering gap is outside Eq/Ne.
    with pytest.raises(sx.CoercionFailed):k.convert(sx.nan)
    with pytest.raises(ZeroDivisionError):a/k.zero
    assert a!=k.convert(t+2)
