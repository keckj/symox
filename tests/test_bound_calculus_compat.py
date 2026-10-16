"""Native held limits, safe alpha identity, and exact sign exclusions."""
import pickle
import subprocess
import sys
import pytest
import sympy as sp
import symox as sm


@pytest.mark.parametrize('flags', [
    {'nonnegative': False}, {'nonpositive': False},
    {'nonnegative': False, 'nonpositive': False},
    {'nonnegative': False, 'real': True}, {'nonpositive': False, 'real': True},
    {'nonnegative': False, 'negative': False}, {'nonpositive': False, 'positive': False},
    {'real': False}, {'imaginary': True},
])
def test_sign_exclusions_are_not_opposite_real_signs(flags):
    x=sm.Symbol('x')
    for actual,expected in [(sm.Symbol('q',**flags),sp.Symbol('q',**flags)),
                            (sm.Function('q',**flags)(x),sp.Function('q',**flags)(sp.Symbol('x')))]:
        for fact in ['real','positive','negative','nonnegative','nonpositive','zero','nonzero']:
            assert getattr(actual,'is_'+fact)==getattr(expected,'is_'+fact),(flags,fact)
        assert sm.from_sympy(expected)==actual
        assert sm.from_sympy(sm.to_sympy(actual))==actual
        assert pickle.loads(pickle.dumps(actual))==actual
        assert sm.sympify(actual.srepr())==actual
    with pytest.raises(ValueError):sm.Symbol('bad',real=True,nonnegative=False,nonpositive=False)


def test_limit_is_native_arithmetic_with_exact_public_args_and_scope():
    x,y,z=sm.symbols('x y z'); f=sm.Function('f')
    held=sm.Limit(f(x)+y,x,0)
    assert isinstance(held,sm.Expr) and isinstance(held,sm.Limit)
    assert held.is_Limit and not held.is_Function
    assert held.args==(f(x)+y,x,0,sm.Symbol('+'))
    assert held.func(*held.args)==held
    assert held.expr==f(x)+y and held.var==x and held.point==0 and held.dir=='+'
    assert held.free_symbols=={y} and held.has(x) and held.has(y)
    assert (held+1).has(held)
    assert held.subs(y,z)==sm.Limit(f(x)+z,x,0)
    assert held.subs(x,z)==sm.Limit(f(z)+y,z,0)
    assert held.xreplace({y:z})==sm.Limit(f(x)+z,x,0)
    assert held.doit()==held and (held+1).doit()==held+1
    assert held.is_commutative is None
    assert sm.Function('Limit')(x).is_Limit is False
    assert not isinstance(sm.Function('Limit')(x),sm.Limit)
    with pytest.raises(NotImplementedError):held.subs(x,2)
    with pytest.raises(NotImplementedError):sm.Limit(f(x),x,x)
    with pytest.raises(ValueError):sm.Limit(f(x),x,0,dir='invalid')


@pytest.mark.parametrize('point,direction',[(0,'+'),(0,'-'),(0,'+-'),(sm.oo,'-'),(-sm.oo,'+')])
def test_limit_roundtrips_and_unknown_continuity(point,direction):
    x=sm.Symbol('x');f=sm.Function('f')
    held=sm.Limit(f(x),x,point,dir=direction)
    assert held.doit()==held
    assert str(held)==str(sp.Limit(sp.Function('f')(sp.Symbol('x')),sp.Symbol('x'),sm.to_sympy(point),dir=direction))
    assert held.args[-1]==sm.Symbol(direction)
    assert pickle.loads(pickle.dumps(held))==held
    assert sm.sympify(held.srepr())==held
    assert sm.sympify(str(held))==held
    assert sm.from_sympy(sm.to_sympy(held))==held
    assert '\\lim' in sm.latex(held)
    with pytest.raises((ValueError,TypeError,NotImplementedError)):sm.lambdify(x,held,'numpy')


def test_native_limit_doit_preserves_existing_math_and_errors():
    x=sm.Symbol('x')
    assert sm.Limit(sm.sin(x)/x,x,0,dir='+-').doit()==1
    assert (1+sm.Limit(sm.sin(x)/x,x,0)).doit()==2
    assert float(sm.Limit(sm.sin(x)/x,x,0).evalf(50))==1.0
    assert sm.Limit(1/x,x,0,dir='+-').doit()==sm.zoo
    with pytest.raises(ValueError):sm.Limit(sm.sign(x),x,0,dir='+-').doit()
    with pytest.raises(ValueError):(1+sm.Limit(sm.sign(x),x,0,dir='+-')).doit()
    assert sm.Limit(x,x,sm.oo,dir='+').args[-1]==sm.Symbol('-')


def test_safe_alpha_subs_hash_equality_and_structural_rewrites():
    x,y,z=sm.symbols('x y z');f=sm.Function('f'); pair=sm.Function('pair')
    a,b=sm.Subs(f(x),x,0),sm.Subs(f(y),y,0)
    assert a==b and hash(a)==hash(b)
    assert len({a,b})==1 and {a:1}[b]==1
    assert a+b==2*a
    assert a.variables==(x,) and b.variables==(y,)
    assert (a+z).has(b) and (a+z).xreplace({b:1})==z+1
    assert (a+z).subs(b,1)==z+1
    assert pair(a,b).xreplace({x:z}).args[0].variables==(z,)
    assert pair(a,b).xreplace({x:z}).args[1].variables==(y,)
    assert sm.Subs(f(x),x,y)==sm.Subs(f(z),z,y)
    assert sm.Subs(f(x),x,x)!=sm.Subs(f(y),y,y)
    assert sm.Subs(1,x,0)!=sm.Subs(1,y,0)
    assert sm.Subs(x+y,(x,y),(y,z))!=sm.Subs(x+y,(x,y),(z,y))
    for held in [a,b]:
        restored=pickle.loads(pickle.dumps(held))
        assert restored==held and restored.args==held.args


def test_nested_binders_and_unsupported_binders_keep_safe_identity():
    x,y,z=sm.symbols('x y z');f=sm.Function('f');g=sm.Function('g')
    a=sm.Subs(g(sm.Subs(f(x)+y,y,0),x),x,1)
    b=sm.Subs(g(sm.Subs(f(z)+y,y,0),z),z,1)
    assert a==b and hash(a)==hash(b)
    # A Limit is a different binder; alpha identity deliberately stays structural.
    a=sm.Subs(sm.Limit(f(x),x,0),x,1)
    b=sm.Subs(sm.Limit(f(y),y,0),y,1)
    assert a!=b


def test_bound_operations_work_without_sympy_imports():
    subprocess.run([sys.executable,'-c',"""
import sys
sys.modules['sympy']=None
import symox as s,pickle
x,y=s.symbols('x y');f=s.Function('f')
a,b=s.Subs(f(x),x,0),s.Subs(f(y),y,0)
assert a==b and hash(a)==hash(b)
q=s.Symbol('q',nonnegative=False)
assert q.is_real is None and q.is_negative is None and q.is_nonnegative is False
held=s.Limit(f(x),x,0)
assert held.free_symbols==set() and (held+1).has(x)
assert pickle.loads(pickle.dumps(held))==held
assert s.Limit(s.sin(x)/x,x,0).doit()==1
"""],check=True,timeout=30)
