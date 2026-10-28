"""Certified domains and endpoint-log branches for retained rational integrals."""
import importlib
import pytest
import sympy as sp
import symox as sx

x,a=sx.symbols('x a')
interval=sx.Interval(0,1)


def held(f):return sx.Integral(f,(x,0,1))
def pieces(e):
    assert e.is_Piecewise
    return [(pair[0],pair[1]) for pair in e.args]


@pytest.mark.parametrize('sign',[1,-1])
def test_exact_pole_domain_and_endpoint_formula_are_structural_certificates(sign):
    c=sign*a;f=1/(x+c);result=sx.integrate(f,(x,0,1));branches=pieces(result)
    value,condition=branches[0]
    expected_condition=sx.And(sx.Contains(c,sx.S.Complexes),sx.Not(sx.Contains(-c,interval)))
    assert value==sx.log(1+c)-sx.log(c) and condition==expected_condition
    assert branches[-1]==(held(f),sx.true)
    # The excluded point is the unique zero of the denominator, exactly;
    # no numerical root location or quadrature is used for this certificate.
    xp,cp=sp.symbols('x c')
    polynomial=sp.Poly(xp+cp,xp)
    assert polynomial.degree()==1 and polynomial.eval(-cp)==0
    assert sp.cancel(sp.diff(sp.log(xp+cp),xp)-1/(xp+cp))==0
    assert value==sx.log((x+c).subs(x,1))-sx.log((x+c).subs(x,0))


@pytest.mark.parametrize('sign',[1,-1])
@pytest.mark.parametrize('c',[sx.Integer(-1),sx.Rational(-3,4),sx.Rational(-1,2),sx.Integer(0)])
def test_every_selected_real_pole_keeps_the_original_integral_until_doit(sign,c):
    f=1/(x+sign*a);result=sx.integrate(f,(x,0,1))
    at=result.subs(a,sign*c)
    assert at==held(1/(x+c))
    want=sx.oo if c==0 else (-sx.oo if c==-1 else sx.nan)
    assert at.doit()==want
    assert sx.integrate(1/(x+c),(x,0,1))==want


@pytest.mark.parametrize('c',[sx.Integer(-3),sx.Rational(-5,4),sx.Rational(1,7),sx.Integer(2),sx.I,-2+sx.I,2-3*sx.I])
def test_finite_regular_substitutions_select_the_same_principal_endpoint_logs(c):
    result=sx.integrate(1/(x+a),(x,0,1)).subs(a,c)
    expected=sx.log(1+c)-sx.log(c)
    assert result==expected and not result.has(sx.Integral)
    cp=sp.sympify(str(c));reference=sp.log(1+cp)-sp.log(cp)
    assert sp.simplify(result.to_sympy()-reference)==0
    if c in (-3,sx.Rational(-5,4)):
        assert result.is_real is True


def test_symbolic_log_branch_and_linear_factorization_certificates():
    # The branch proof is the horizontal path described in the module and
    # evidence file. These exact identities certify the relevant primitive
    # and unique pole without replacing the proof by sampled values.
    xp,ap=sp.symbols('x a',nonzero=True)
    assert sp.cancel((1+ap*xp)-ap*(xp+1/ap))==0
    assert sp.cancel(sp.diff(sp.log(xp+1/ap)/ap,xp)-1/(1+ap*xp))==0
    assert sp.Poly(1+ap*xp,xp).degree()==1
    assert sp.cancel((1+ap*xp).subs(xp,-1/ap))==0
    result=sx.integrate(1/(1+a*x),(x,0,1));branches=pieces(result)
    assert branches[0]==(sx.Integer(1),sx.Eq(a,0))
    assert branches[1][0]==(sx.log(1+1/a)-sx.log(1/a))/a
    assert branches[-1]==(held(1/(1+a*x)),sx.true)
    # The nonzero condition may be removed as redundant after Eq(a,0);
    # the effective branch still explicitly excludes that preceding case.
    actual=sx.And(branches[1][1],sx.Ne(a,0))
    expected=sx.And(sx.Contains(a,sx.S.Complexes),sx.Ne(a,0),sx.Not(sx.Contains(-1/a,interval)))
    assert actual==expected


@pytest.mark.parametrize('value',[0,1,2,sx.Rational(-1,2),sx.I,-2+sx.I])
def test_affine_slope_specializations_preserve_zero_and_separate_log_branches(value):
    value=sx.sympify(value)
    result=sx.integrate(1/(1+a*x),(x,0,1)).subs(a,value)
    expected=sx.Integer(1) if value==0 else (sx.log(1+1/value)-sx.log(1/value))/value
    assert result==expected and not result.has(sx.Integral)


@pytest.mark.parametrize('value',[-1,-2,-3])
def test_affine_slope_poles_are_never_interpreted_as_principal_values(value):
    result=sx.integrate(1/(1+a*x),(x,0,1)).subs(a,value)
    assert result==held(1/(1+value*x))
    expected=sx.oo if value==-1 else sx.nan
    assert result.doit()==expected


@pytest.mark.parametrize('assumptions',[{}, {'real':True},{'complex':True},{'nonzero':True}])
@pytest.mark.parametrize('value',[sx.oo,-sx.oo,sx.zoo,sx.nan])
def test_nonfinite_substitutions_never_select_a_finite_continued_value(assumptions,value):
    c=sx.Symbol('c',**assumptions)
    for f in (1/(x+c),1/(x-c),1/(1+c*x)):
        result=sx.integrate(f,(x,0,1))
        at=result.subs(c,value)
        # The fallback retains the original substituted integrand, including
        # any NaN or zero produced by that substitution.
        assert at==held(f).subs(c,value)
        assert at.doit()==sx.integrate(f.subs(c,value),(x,0,1))


@pytest.mark.parametrize('deep',[False,True])
def test_repeated_doit_is_idempotent_and_does_not_embed_its_own_piecewise(deep):
    for f in (1/(x+a),1/(x-a),1/(1+a*x)):
        result=sx.integrate(f,(x,0,1));want=result
        for _ in range(5):
            result=result.doit(deep=deep)
            assert result==want and sx.count_ops(result)==sx.count_ops(want)
        assert result.has(sx.Integral)


def test_piecewise_fixed_point_guard_still_evaluates_other_branches():
    p=sx.Symbol('p',real=True)
    conditional=sx.integrate(1/(x+a),(x,0,1))
    expression=sx.Piecewise((sx.Integral(x,(x,0,1)),p>0),(conditional,True))
    result=expression.doit()
    assert result.args[0][0]==sx.Rational(1,2)
    assert result==result.doit()


def test_piecewise_doit_preserves_fatal_host_errors(monkeypatch):
    module=importlib.import_module('symox.integrals._methods');calls=[]
    def interrupted(*args):calls.append(args);raise KeyboardInterrupt('host integral')
    monkeypatch.setattr(module,'_native_host',interrupted)
    f=sx.Function('f')
    expression=sx.Piecewise((sx.Integral(f(x),x),a>0),(0,True))
    with pytest.raises(KeyboardInterrupt,match='host integral'):expression.doit()
    assert calls


def test_unchanged_raw_piecewise_preserves_its_unevaluated_presentation():
    value=sx.Piecewise((1,a>0),(2,a>0),(3,True),evaluate=False)
    assert len(value.args)==3
    assert value.doit().args==value.args


def test_conditional_rules_do_not_require_runtime_sympy():
    import os,subprocess,sys
    script="""import sys
sys.modules['sympy']=None
import symox as s
x,a=s.symbols('x a')
for f in (1/(x+a),1/(x-a),1/(1+a*x)):
 result=s.integrate(f,(x,0,1))
 assert result.is_Piecewise and result.has(s.Integral)
 assert result.doit()==result==result.doit().doit()
 assert s.integrate(f,(x,0,1),conds='none').has(s.Integral)
assert s.integrate(1/(1+a*x),(x,0,1)).subs(a,0)==1
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


@pytest.mark.parametrize('assumptions',[{}, {'real':True},{'complex':True},{'nonzero':True},{'negative':True}])
def test_slope_zero_substitution_is_safe_even_when_it_contradicts_declared_assumptions(assumptions):
    c=sx.Symbol('c',**assumptions)
    result=sx.integrate(1/(1+c*x),(x,0,1))
    for _ in range(3):
        assert result.subs(c,0).doit()==1
        result=result.doit()
