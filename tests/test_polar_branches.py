"""Polar arithmetic must preserve turns and principal branch boundaries."""
import builtins
import pytest
import symox as s
import sympy as sp

x=s.Symbol('x')
p=s.Symbol('p',positive=True)
z=s.Symbol('z',polar=True)


@pytest.mark.parametrize('coefficient',[-9,-5,-3,-2,-1,0,1,2,3,5,9])
@pytest.mark.parametrize('period',[s.pi,2*s.pi,4*s.pi,s.oo])
def test_periodic_arguments_keep_winding(coefficient,period):
    value=s.exp_polar(coefficient*s.I*s.pi)
    actual=s.periodic_argument(value,period)
    expected=sp.periodic_argument(sp.exp_polar(coefficient*sp.I*sp.pi),period.to_sympy())
    assert sp.simplify(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('value',[1,-1,4,-4,s.I,-s.I,2*s.I,-2*s.I,1+s.I])
def test_polar_lift(value):
    value=s.sympify(value)
    actual=s.polar_lift(value)
    expected=sp.polar_lift(value.to_sympy())
    assert actual.to_sympy()==expected
    assert s.Abs(actual)==s.Abs(value)
    assert s.unpolarify(actual)==value


def test_polar_structure_and_substitution():
    lifted=s.polar_lift(x)
    assert lifted.class_name=='polar_lift'
    assert not lifted.is_UndefinedFunction
    assert lifted.is_polar is True
    assert s.polar_lift(z)==z
    assert s.polar_lift(p)==p*s.exp_polar(0)
    assert s.polar_lift(p*x)==p*s.polar_lift(x)
    assert s.polar_lift(x).subs(x,-2)==2*s.exp_polar(s.I*s.pi)
    assert (-s.exp_polar(x)).is_polar is not True
    assert s.periodic_argument(p*s.exp_polar(3*s.I*s.pi),2*s.pi)==s.pi
    assert s.periodic_argument(p*x,2*s.pi)==s.periodic_argument(x,2*s.pi)


def test_polar_exp_log_power_and_conjugate():
    ep=s.exp_polar(x)
    assert s.log(ep)==x
    assert ep**s.Rational(1,2)==s.exp_polar(x/2)
    assert ep**0==1
    assert s.diff(ep,x)==ep
    assert s.conjugate(ep)==s.exp_polar(s.conjugate(x))
    assert s.Abs(ep)==s.exp(s.re(x))
    assert s.unpolarify(s.exp_polar(3*s.I*s.pi))==-1
    assert s.unbranched_argument(s.exp_polar(3*s.I*s.pi))==3*s.pi


@pytest.mark.parametrize('angle',[-5,-3,-1,0,1,3,5])
def test_principal_branch(angle):
    ep=2*s.exp_polar(angle*s.I*s.pi)
    actual=s.principal_branch(ep,2*s.pi)
    target=2*s.exp_polar((0 if angle%2==0 else 1)*s.I*s.pi)
    assert actual==target
    assert actual.is_polar is True
    assert s.principal_branch(ep,s.oo)==ep
    assert s.principal_branch(s.polar_lift(x),2*s.pi)==s.principal_branch(x,2*s.pi)


def test_principal_branch_positive_and_unknown():
    assert s.principal_branch(p,2*s.pi)==p*s.exp_polar(0)
    unknown=s.principal_branch(x,2*s.pi)
    assert unknown.class_name=='principal_branch'
    assert unknown.subs(x,-2)==2*s.exp_polar(s.I*s.pi)
    assert s.periodic_argument(s.exp_polar(s.I*x),2*s.pi).class_name=='periodic_argument'


def test_polarify_roundtrip_and_explicit_lift():
    expr=(x+1)**p+s.sin(x)
    result,reverse=s.polarify(expr)
    assert reverse and all(v.is_polar for v in reverse)
    assert s.unpolarify(result,subs=reverse)==expr
    assert s.polarify(x,lift=True)==s.polar_lift(x)
    assert s.polarify(x+1,lift=True)==s.polar_lift(x+1)
    assert not isinstance(s.polarify(x,subs=False),tuple)
    assert s.polarify(s.exp(x),subs=False)==s.exp(x)


def test_unpolarify_protects_branch_sensitive_arguments():
    lifted=s.polar_lift(x)
    root=s.sqrt(lifted)
    assert s.unpolarify(root)==root
    assert s.unpolarify(s.log(lifted))==s.log(lifted)
    assert s.unpolarify(lifted**2)==x**2
    assert s.unpolarify(s.sin(lifted))==s.sin(x)
    assert s.unpolarify(s.principal_branch(x,2*s.pi))==x
    assert s.unpolarify(s.exp_polar(0),exponents_only=True)==1
    assert s.unpolarify(True) is True


def test_native_computation_has_no_sympy_import(monkeypatch):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('native polar arithmetic imported SymPy')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    lifted=s.polar_lift(-3)
    assert s.Abs(lifted)==3
    assert s.unpolarify(lifted)==-3
    assert s.periodic_argument(s.exp_polar(5*s.I*s.pi),2*s.pi)==s.pi


def test_floating_exponent_stays_polar():
    value=s.exp_polar(s.Float('0.25'))
    assert value.class_name=='exp_polar'
    assert value.is_polar is True
    assert not value.is_UndefinedFunction
