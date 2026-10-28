"""Principal E1 preserves the analytic Laplace value on the whole right half-plane.

The reference Laplace table has a known prefactor error for the reciprocal.
Expected values here follow the exact E1 substitution/parts proof, not that table.
"""
import pytest
import sympy as sp
import symox as sx

t,s=sx.symbols('t s')
tp,spv=sp.symbols('t s')


def data(a,c,logarithm=False,frequency=0):
    a,c,frequency=map(sx.sympify,(a,c,frequency))
    q=c*(s+frequency)/a
    tail=sx.exp(q)*sx.expint(1,q)
    kernel=sx.log(a*t+c) if logarithm else 1/(a*t+c)
    expected=(sx.log(c)+tail)/(s+frequency) if logarithm else tail/a
    return sx.exp(-frequency*t)*kernel,expected


@pytest.mark.parametrize('a,c',[(1,1),(2,3),(sx.Rational(3,2),sx.Rational(2,5))])
@pytest.mark.parametrize('logarithm',[False,True])
@pytest.mark.parametrize('simplify',[False,True])
def test_exact_e1_value_and_roc_for_every_existing_positive_affine_rule(a,c,logarithm,simplify):
    kernel,expected=data(a,c,logarithm)
    actual=sx.laplace_transform(kernel,t,s,simplify=simplify)
    assert sx.cancel(actual[0]-expected)==0
    assert actual[1:]==(0,sx.true)
    assert not actual[0].has(sx.Ei)
    assert sx.laplace_transform(kernel,t,s,noconds=True,simplify=simplify)==actual[0]


@pytest.mark.parametrize('facts',[{}, {'real':True},{'positive':True},{'complex':True},{'nonzero':True}])
@pytest.mark.parametrize('logarithm',[False,True])
def test_frequency_assumptions_never_discard_e1_before_complex_substitution(facts,logarithm):
    z=sx.Symbol('z',**facts)
    kernel,expected=data(1,1,logarithm)
    value=sx.laplace_transform(kernel,t,z,noconds=True)
    expected=expected.subs(s,z)
    assert value==expected and not value.has(sx.Ei)
    # Even a specialization contrary to initial real/positive assumptions
    # retains the same E1 expression, rather than crossing a hidden Ei cut.
    for point in (1+sx.I,1-sx.I,sx.Rational(1,3)+2*sx.I):
        assert value.subs(z,point)==expected.subs(z,point)


@pytest.mark.parametrize('a,c',[(1,1),(2,3)])
@pytest.mark.parametrize('logarithm',[False,True])
@pytest.mark.parametrize('frequency',[sx.Integer(2),1+sx.I,-sx.Rational(1,2)+sx.I])
def test_complex_frequency_shift_preserves_value_and_abscissa(a,c,logarithm,frequency):
    kernel,expected=data(a,c,logarithm,frequency)
    value,plane,condition=sx.laplace_transform(kernel,t,s)
    assert sx.cancel(value-expected)==0
    assert plane==-sx.re(frequency) and condition==sx.true
    assert not value.has(sx.Ei)


@pytest.mark.parametrize('a,c',[(1,1),(2,3)])
@pytest.mark.parametrize('logarithm',[False,True])
@pytest.mark.parametrize('degree',[1,2])
def test_time_moments_use_exact_frequency_differentiation(a,c,logarithm,degree):
    kernel,expected=data(a,c,logarithm)
    value,plane,condition=sx.laplace_transform(t**degree*kernel,t,s)
    assert sx.cancel(value-(-1)**degree*expected.diff(s,degree))==0
    assert plane==0 and condition==sx.true


def test_symbolic_affine_scaling_has_exact_integral_and_differential_certificates():
    a,c=sx.symbols('a c',positive=True)
    f,expected=data(a,c)
    value=sx.laplace_transform(f,t,s,noconds=True)
    assert value==expected
    # The substitution u=t+c/a yields exp(beta*s)*E1(beta*s)/a;
    # its ROC and the absence of branch jumps are proved in the artifact.
    assert sx.cancel(a*value.diff(s)-c*value+1/s)==0
    g,logvalue=data(a,c,True)
    got=sx.laplace_transform(g,t,s,noconds=True)
    assert sx.cancel(s*got-sx.log(c)-a*value)==0
    # Independent exact symbolic identities, not numeric agreement alone.
    ap,cp=sp.symbols('a c',positive=True)
    vp=sp.exp(cp*spv/ap)*sp.expint(1,cp*spv/ap)/ap
    assert sp.cancel(ap*sp.diff(vp,spv)-cp*vp+1/spv)==0
    assert sp.cancel(sp.diff(sp.log(ap*tp+cp),tp)-ap/(ap*tp+cp))==0


@pytest.mark.parametrize('imaginary',[1,-1])
@pytest.mark.parametrize('digits',[40,80])
@pytest.mark.parametrize('logarithm',[False,True])
def test_upper_and_lower_complex_half_planes_keep_requested_precision(imaginary,digits,logarithm):
    kernel,expected=data(2,3,logarithm)
    value=sx.laplace_transform(kernel,t,s,noconds=True)
    point=sx.Rational(2,3)+imaginary*sx.I/7
    got=value.subs(s,point).evalf(digits,strict=True,maxn=digits*3)
    zp=sp.Rational(2,3)+imaginary*sp.I/7
    tail=sp.exp(3*zp/2)*sp.expint(1,3*zp/2)
    want=(sp.log(3)+tail)/zp if logarithm else tail/2
    # Precision regression only; the analytic proof fixes the branch.
    error=abs(sp.N(sp.sympify(str(got))-want,digits+15))
    assert error<sp.Rational(1,10**(digits-5))


def test_unrelated_and_unsupported_rules_keep_established_paths():
    f=sx.Function('f')
    assert sx.laplace_transform(f(t),t,s,noconds=True)==sx.LaplaceTransform(f(t),t,s)
    assert sx.laplace_transform(1,t,s)==(1/s,0,sx.true)


def test_e1_laplace_rules_do_not_import_runtime_sympy():
    import os,subprocess,sys
    script="""import sys
class NoReference:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoReference())
import symox as S
t,s=S.symbols('t s')
for kernel,want in [(1/(t+1),S.exp(s)*S.expint(1,s)),(S.log(t+1),S.exp(s)*S.expint(1,s)/s)]:
 value,plane,condition=S.laplace_transform(kernel,t,s)
 assert value==want and plane==0 and condition==S.true
 assert not value.has(S.Ei)
assert 'sympy' not in sys.modules
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


@pytest.mark.parametrize('a,c',[(-1,-2),(sx.I,2*sx.I),(-2*sx.I,-3*sx.I)])
def test_reciprocal_rule_keeps_complex_scalars_when_the_offset_ratio_is_positive(a,c):
    kernel,expected=data(a,c)
    value,plane,condition=sx.laplace_transform(kernel,t,s)
    assert sx.cancel(value-expected)==0
    assert plane==0 and condition==sx.true and not value.has(sx.Ei)


def test_held_transform_doit_uses_the_same_principal_branch():
    for logarithm in (False,True):
        kernel,expected=data(2,3,logarithm)
        held=sx.LaplaceTransform(kernel,t,s)
        assert sx.cancel(held.doit(noconds=True)-expected)==0
