"""Gamma moments use the sign of a real exponential power, without degree caps."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('power',['1/2','3/2','-1','-2','-1/2'])
@pytest.mark.parametrize('rate',[1,2])
@pytest.mark.parametrize('conds',['piecewise','none','separate'])
def test_rational_signed_powers_match_reference_formula_and_domain(power,rate,conds):
    x,a=sx.symbols('x a');q=sx.sympify(power)
    expression=x**a*sx.exp(-rate*x**q)
    actual=sx.integrate(expression,(x,0,sx.oo),conds=conds)
    expected=sp.integrate(expression.to_sympy(),(x.to_sympy(),0,sp.oo),conds=conds)
    if conds=='piecewise':
        converted=actual.to_sympy()
        assert isinstance(converted,sp.Piecewise)
        assert converted.args[0][1]==expected.args[0][1]
        assert sp.simplify(converted.args[0][0]-expected.args[0][0])==0
        assert converted.args[1]==expected.args[1]
    elif conds=='separate':
        assert actual[1].to_sympy()==expected[1]
        assert sp.simplify(actual[0].to_sympy()-expected[0])==0
    else:assert sp.simplify(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('positive',[True,False])
def test_symbolic_signed_real_power_uses_absolute_jacobian(positive):
    x,a=sx.symbols('x a');q=sx.Symbol('q',**({'positive':True} if positive else {'negative':True}))
    beta=sx.Symbol('b',positive=True);s=(a+1)/q
    value,condition=sx.integrate(x**a*sx.exp(-beta*x**q),(x,0,sx.oo),conds='separate')
    assert value==sx.gamma(s)/((q if positive else -q)*beta**s)
    assert condition==((sx.re(a)>-1) if positive else (sx.re(a)<-1))


def test_negative_even_power_has_convergent_whole_line_moments():
    x=sx.Symbol('x')
    assert sx.integrate(x**(-2)*sx.exp(-x**(-2)),(x,-sx.oo,sx.oo))==sx.sqrt(sx.pi)
    assert sx.integrate(x**(-3)*sx.exp(-x**(-2)),(x,-sx.oo,sx.oo))==0


def test_zero_input_keeps_reference_scalar_shape_with_separate():
    x=sx.Symbol('x')
    assert sx.integrate(0,(x,0,sx.oo),conds='separate')==0
    assert sx.integrate(x*sx.exp(-x*x),(x,-sx.oo,sx.oo),conds='separate')==(0,sx.true)


def test_large_or_symbolic_even_power_does_not_require_an_i64_degree():
    script='''
import sys
class NoReference:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoReference())
import symox as sx
x,a=sx.symbols('x a');q=sx.Integer(2)**100
value,condition=sx.integrate(x**a*sx.exp(-x**q),(x,0,sx.oo),conds='separate')
assert value==sx.gamma((a+1)/q)/q
assert condition==(sx.re(a)>-1)
assert sx.integrate(sx.exp(-x**q),(x,-sx.oo,sx.oo))==2*sx.gamma(1/q)/q
q=sx.Symbol('q',positive=True,even=True)
assert sx.integrate(sx.exp(-x**q),(x,-sx.oo,sx.oo))==2*sx.gamma(1/q)/q
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=20)
