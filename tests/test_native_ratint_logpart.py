"""LRT pairs fix logarithm branches; a derivative certificate is additional."""
import importlib
import random

import pytest
import sympy as sp
import symox as sx

rt=importlib.import_module('symox.integrals.rationaltools')
reference=importlib.import_module('sympy.integrals.rationaltools')
x,t=sx.symbols('x t')
xp,tp=sp.symbols('x t')


def pairs(f,g,domain='ZZ'):
    return sx.Poly(sx.sympify(str(f),locals={'x':x}),x,domain=domain),sx.Poly(sx.sympify(str(g),locals={'x':x}),x,domain=domain)


def payload(value):
    return [(sp.sympify(str(h.as_expr())),str(h.domain),sp.sympify(str(q.as_expr())),str(q.domain)) for h,q in value]


def certificate(value,f,g):
    # Trace of multiplication by t*h_x/h in QQ(x)[t]/q is precisely the
    # sum over all roots of q. No floating roots or branch samples are used.
    field=sp.QQ.frac_field(xp)
    total=0
    for h,_,q,_ in value:
        q=sp.Poly(q,tp,domain=field)
        hpoly=sp.Poly(h,tp,domain=field)
        inv=sp.invert(hpoly,q)
        r=(sp.Poly(tp*sp.diff(h,xp),tp,domain=field)*inv).rem(q)
        total+=sum((r*sp.Poly(tp**j,tp,domain=field)).rem(q).nth(j) for j in range(q.degree()))
    assert sp.cancel(total-f/g)==0


@pytest.mark.parametrize('domain',['ZZ','QQ'])
@pytest.mark.parametrize('seed',range(48))
def test_exact_pairs_domains_and_trace_certificate(domain,seed,monkeypatch):
    rng=random.Random(714231+seed)
    degree=seed%7+1
    # Squarefree denominators include irreducible and reducible factors,
    # repeated residues, missing PRS degrees and nonmonic leading terms.
    g=sp.Poly(sum(rng.randrange(-4,5)*xp**j for j in range(degree))+rng.choice([-3,-1,1,2])*xp**degree,xp)
    g=g.sqf_part().as_expr()*rng.choice([-3,-1,1,2])
    f=sum(sp.Rational(rng.randrange(-5,6),rng.randrange(1,6) if domain=='QQ' else 1)*xp**j for j in range(sp.degree(g,xp)))
    f=f or 1
    if sp.gcd(sp.Poly(f,xp),sp.Poly(g,xp)).degree():f=1
    a,b=pairs(f,g,domain)
    native=sx._core._ratint_logpart(a,b,x,t)
    assert native is not None
    expected=reference.ratint_logpart(sp.Poly(f,xp,domain=domain),sp.Poly(g,xp,domain=domain),xp,tp)
    actual=rt.ratint_logpart(a,b,x,t)
    assert payload(native)==payload(expected)==payload(actual)
    certificate(payload(native),f,g)
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert payload(rt.ratint_logpart(a,b,x,t))==payload(native)


@pytest.mark.parametrize('f,g',[(1,xp*xp+xp+1),(12,xp*xp-xp-2),(1,2*xp+2),(xp,xp*xp+1),(2*xp,xp*xp-1),(3*xp*xp,xp**3+1),(1,(xp-1)*(xp+1)*(xp*xp+1)),(xp**3+1,xp**4+xp+1)])
def test_signed_squarefree_groups_and_unchanged_whole_denominator(f,g):
    a,b=pairs(f,g)
    result=rt.ratint_logpart(a,b,x,t)
    expected=reference.ratint_logpart(sp.Poly(f,xp),sp.Poly(g,xp),xp,tp)
    assert payload(result)==payload(expected)
    certificate(payload(result),f,g)
    if len(result)==1 and result[0][0].degree()==b.degree():
        assert result[0][0]==b


@pytest.mark.parametrize('f,g',[(0,1),(0,xp),(1,1),(1,0),(xp,xp*xp),(1,xp*xp),(xp*xp+1,xp*xp+1)])
def test_nonpreconditions_keep_full_existing_result_or_error(f,g,monkeypatch):
    a,b=pairs(f,g)
    assert sx._core._ratint_logpart(a,b,x,t) is None
    def run():
        try:return ('value',payload(rt.ratint_logpart(a,b,x,t)))
        except Exception as error:return ('error',type(error),str(error))
    actual=run()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert run()==actual


@pytest.mark.parametrize('domain',['EX','RR','QQ<sqrt(2)>','GF(5)'])
def test_other_domains_retain_complete_algorithm(domain,monkeypatch):
    if domain=='QQ<sqrt(2)>':a,b=sx.Poly(1,x,extension=sx.sqrt(2)),sx.Poly(x*x+1,x,extension=sx.sqrt(2))
    elif domain=='GF(5)':a,b=sx.Poly(1,x,modulus=5),sx.Poly(x*x+1,x,modulus=5)
    else:a,b=pairs(1,xp*xp+xp+1,domain)
    assert sx._core._ratint_logpart(a,b,x,t) is None
    def run():
        try:return ('value',payload(rt.ratint_logpart(a,b,x,t)))
        except Exception as error:return ('error',type(error),str(error))
    actual=run();monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert run()==actual


def test_t_alias_false_and_custom_parameters_keep_existing_normalization(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1)
    for parameter in (x,t+1,sx.Integer(2)):
        assert sx._core._ratint_logpart(a,b,x,parameter) is None
    # Dummy creation remains in the public adapter. A falsy t has the same
    # observable constructor protocol as the omitted argument.
    calls=[];old=rt.Dummy
    def dummy(*args,**kwargs):calls.append((args,kwargs));return old(*args,**kwargs)
    monkeypatch.setattr(rt,'Dummy',dummy)
    assert len(rt.ratint_logpart(a,b,x,0))==1 and len(calls)==1


def test_protocol_refusal_does_not_touch_custom_objects():
    class Protocol:
        def __getattribute__(self,name):raise AssertionError('eligibility called a protocol')
    obj=Protocol()
    assert sx._core._ratint_logpart(obj,obj,obj,obj) is None
    class Derived(sx.Poly):pass
    a,b=pairs(1,xp*xp+1)
    assert sx._core._ratint_logpart(Derived(a),b,x,t) is None


@pytest.mark.parametrize('method',['sqf_list','LC','primitive','coeffs','monoms','invert','diff'])
def test_actual_poly_overrides_are_observed(method,monkeypatch):
    a,b=pairs(1,xp*xp+xp+1)
    old=getattr(sx.Poly,method);calls=[]
    def override(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    monkeypatch.setattr(sx.Poly,method,override)
    assert sx._core._ratint_logpart(a,b,x,t) is None
    rt.ratint_logpart(a,b,x,t)
    assert calls


def test_resultant_override_is_observed_and_native_skipped(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1);calls=[];old=rt.resultant
    def result(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    monkeypatch.setattr(rt,'resultant',result)
    rt.ratint_logpart(a,b,x,t)
    assert calls==[1]


def test_pending_interrupt_is_never_a_refusal():
    import _thread
    a,b=pairs(1,xp*xp+xp+1)
    with pytest.raises(KeyboardInterrupt):
        _thread.interrupt_main()
        sx._core._ratint_logpart(a,b,x,t)


def test_held_evaluation_preserves_complete_path(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1)
    # The held host route is allowed to retain its existing errors.
    def run():
        with sx.evaluate(False):
            try:return ('value',payload(rt.ratint_logpart(a,b,x,t)))
            except Exception as error:return ('error',type(error),str(error))
    actual=run();monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert run()==actual


@pytest.mark.parametrize('f,g',[(1,xp*xp+xp+1),(12,xp*xp-xp-2),(1,xp**4+1),(xp**3+1,xp**4+xp+1)])
def test_downstream_integral_retains_identical_logarithm_branches(f,g,monkeypatch):
    a,b=pairs(f,g)
    actual=rt.ratint((a.as_expr(),b.as_expr()),x)
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    original=rt.ratint((a.as_expr(),b.as_expr()),x)
    canonical=lambda e:e.xreplace({d:t for d in e.atoms(sx.Dummy)})
    assert canonical(actual)==canonical(original)
    assert sx.cancel(actual.diff(x)-a.as_expr()/b.as_expr())==0
    for value in (-3,sx.Rational(-1,2),0,sx.Rational(1,2),3):
        if b.eval(value)==0:continue
        assert canonical(actual.subs(x,value))==canonical(original.subs(x,value))


def test_native_stock_route_avoids_the_host_lrt_recurrence():
    import sys
    a,b=pairs(1,xp*xp+xp+1)
    calls=[]
    def profile(frame,event,arg):
        if event=='call' and frame.f_code.co_name=='_include_sign':calls.append(1)
    old=sys.getprofile()
    try:
        sys.setprofile(profile)
        rt.ratint_logpart(a,b,x,t)
    finally:sys.setprofile(old)
    assert calls==[]


def test_poly_constructor_override_retains_all_intermediate_calls(monkeypatch):
    import symox._poly_constructor as construction
    a,b=pairs(1,xp*xp+xp+1);calls=[];old=construction.construct
    def construct(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    monkeypatch.setattr(construction,'construct',construct)
    actual=payload(rt.ratint_logpart(a,b,x,t));count=len(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert payload(rt.ratint_logpart(a,b,x,t))==actual
    assert len(calls)==count and count>2


def test_degree_bound_selects_the_original_complete_algorithm(monkeypatch):
    a,b=pairs(33*xp**32,xp**33+1)
    assert sx._core._ratint_logpart(a,b,x,t) is None
    actual=payload(rt.ratint_logpart(a,b,x,t))
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert payload(rt.ratint_logpart(a,b,x,t))==actual


def test_constructor_replaced_before_module_reload_is_not_blessed(monkeypatch):
    import symox._poly_constructor as construction
    a,b=pairs(1,xp*xp+xp+1);old=construction.construct;calls=[]
    def replacement(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    with monkeypatch.context() as patch:
        patch.setattr(construction,'construct',replacement)
        importlib.reload(rt)
        rt.ratint_logpart(a,b,x,t)
        assert len(calls)>2
    importlib.reload(rt)


@pytest.mark.parametrize('residues',[(2,2,3),(1,-1,-1),(3,3,3),(1,2,3),(sp.Rational(1,2),sp.Rational(-2,3),sp.Rational(1,2))])
@pytest.mark.parametrize('irreducible',[False,True])
def test_distinct_multiplicity_groups_preserve_order_content_and_domains(residues,irreducible):
    factors=(xp*xp+1,xp*xp+xp+1,xp-2) if irreducible else (xp-1,xp+1,xp-2)
    g=sp.prod(factors)
    f=sp.cancel(sum(r*sp.diff(h,xp)*g/h for r,h in zip(residues,factors)))
    a,b=pairs(f,g,'QQ')
    expected=reference.ratint_logpart(sp.Poly(f,xp,domain='QQ'),sp.Poly(g,xp,domain='QQ'),xp,tp)
    actual=sx._core._ratint_logpart(a,b,x,t)
    assert actual is not None
    assert payload(actual)==payload(expected)
    certificate(payload(actual),f,g)


@pytest.mark.parametrize('bits',[80,200])
@pytest.mark.parametrize('rational',[False,True])
def test_large_exact_coefficients_do_not_use_float_decisions(bits,rational):
    c=sp.Integer(2)**bits+7
    f=sp.Rational(1,c) if rational else c
    g=3*xp*xp+xp+1
    a,b=pairs(f,g,'QQ' if rational else 'ZZ')
    native=sx._core._ratint_logpart(a,b,x,t)
    assert native is not None
    expected=reference.ratint_logpart(sp.Poly(f,xp),sp.Poly(g,xp),xp,tp)
    assert payload(native)==payload(expected)
    certificate(payload(native),f,g)


def test_noncommutative_parameter_preserves_original_error(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1);nc=sx.Symbol('n',commutative=False)
    assert sx._core._ratint_logpart(a,b,x,nc) is None
    def run():
        try:return ('value',payload(rt.ratint_logpart(a,b,x,nc)))
        except Exception as error:return ('error',type(error),str(error))
    actual=run();monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert run()==actual and actual[0]=='error'


def test_expression_as_poly_override_keeps_intermediate_callbacks(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1);old=sx.Expr.as_poly;calls=[]
    def as_poly(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    monkeypatch.setattr(sx.Expr,'as_poly',as_poly)
    actual=payload(rt.ratint_logpart(a,b,x,t));count=len(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert payload(rt.ratint_logpart(a,b,x,t))==actual
    assert count==len(calls)>0


def test_constructor_method_override_before_reload_is_not_blessed(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1);old=sx.Poly._from_normalized;calls=[]
    def construct(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    with monkeypatch.context() as patch:
        patch.setattr(sx.Poly,'_from_normalized',staticmethod(construct))
        importlib.reload(rt)
        assert sx._core._ratint_logpart(a,b,x,t) is None
        rt.ratint_logpart(a,b,x,t)
        assert len(calls)>2
    importlib.reload(rt)


def test_expression_override_before_reload_is_not_blessed(monkeypatch):
    from functools import wraps
    a,b=pairs(1,xp*xp+xp+1);old=sx.Expr.as_poly;calls=[]
    @wraps(old)
    def as_poly(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    with monkeypatch.context() as patch:
        patch.setattr(sx.Expr,'as_poly',as_poly)
        importlib.reload(rt)
        rt.ratint_logpart(a,b,x,t)
        assert calls
    importlib.reload(rt)


def test_resultant_override_before_reload_is_not_blessed(monkeypatch):
    import symox.polys as polys
    a,b=pairs(1,xp*xp+xp+1);old=polys.resultant;calls=[]
    def result(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    with monkeypatch.context() as patch:
        patch.setattr(polys,'resultant',result)
        importlib.reload(rt)
        rt.ratint_logpart(a,b,x,t)
        assert calls==[1]
    importlib.reload(rt)


def test_poly_init_keeps_every_intermediate_callback(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1);calls=[];old=sx.Poly.__init__
    def init(self,*args,**kwargs):calls.append(1);return old(self)
    monkeypatch.setattr(sx.Poly,'__init__',init)
    actual=payload(rt.ratint_logpart(a,b,x,t));count=len(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert payload(rt.ratint_logpart(a,b,x,t))==actual
    assert len(calls)==count and count>2


def test_public_poly_alias_keeps_as_poly_constructor_callbacks(monkeypatch):
    a,b=pairs(1,xp*xp+xp+1);calls=[];old=sx.Poly
    def poly(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    monkeypatch.setattr(sx,'Poly',poly)
    actual=payload(rt.ratint_logpart(a,b,x,t));count=len(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert payload(rt.ratint_logpart(a,b,x,t))==actual
    assert count==len(calls)>0


def test_constructor_poly_alias_keeps_native_factory_callbacks(monkeypatch):
    import symox._poly_constructor as construction
    a,b=pairs(1,xp*xp+xp+1);calls=[];old=construction.Poly
    class Meta(type):
        def __instancecheck__(cls,value):return isinstance(value,old)
        def __getattr__(cls,name):
            if name.startswith('_from_'):calls.append(name)
            return getattr(old,name)
    class Proxy(metaclass=Meta):pass
    monkeypatch.setattr(construction,'Poly',Proxy)
    actual=payload(rt.ratint_logpart(a,b,x,t));count=len(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    assert payload(rt.ratint_logpart(a,b,x,t))==actual
    assert count==len(calls)>2
