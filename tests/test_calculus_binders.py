"""Native held calculus: ordered limits are binding operators, not functions."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm


@pytest.mark.parametrize('operation',[
    lambda s,x,y,z:s.Integral(x,x),
    lambda s,x,y,z:s.Integral(1,x),
    lambda s,x,y,z:s.Integral(x,(x,2)),
    lambda s,x,y,z:s.Integral(x,(x,0,y)),
    lambda s,x,y,z:s.Integral(x,(x,0,x)),
    lambda s,x,y,z:s.Integral(x+y,(x,0,y),(y,0,z)),
    lambda s,x,y,z:s.Integral(s.Function('f')(x),x,x),
    lambda s,x,y,z:s.Integral(s.Function('f')(x),(x,0,1)),
    lambda s,x,y,z:s.Sum(x*y,(x,1,z)),
    lambda s,x,y,z:s.Sum(x+y,(x,1,y),(y,1,z)),
    lambda s,x,y,z:s.Sum(s.Function('f')(x),(x,0,y)),
])
def test_held_calculus_has_exact_arguments_scope_properties_and_roundtrips(operation):
    x,y,z=sm.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    actual=operation(sm,x,y,z);expected=operation(sp,X,Y,Z)
    constructor=sm.Integral if actual.is_Integral else sm.Sum
    assert isinstance(actual,sm.Expr) and isinstance(actual,constructor)
    assert not actual.is_Function
    assert actual.func is constructor and actual.func(*actual.args)==actual
    assert tuple(sm.to_sympy(arg) for arg in actual.args)==expected.args
    assert sm.to_sympy(actual.function)==expected.function
    assert [tuple(sm.to_sympy(value) for value in limit) for limit in actual.limits]==[tuple(limit) for limit in expected.limits]
    assert [sm.to_sympy(value) for value in actual.variables]==expected.variables
    assert [sm.to_sympy(value) for value in actual.bound_symbols]==expected.bound_symbols
    assert {sm.to_sympy(value) for value in actual.free_symbols}==expected.free_symbols
    assert actual.is_number==expected.is_number
    assert actual.is_commutative==expected.is_commutative
    assert str(actual)==str(expected)
    assert actual.srepr()==sp.srepr(expected)
    assert sm.to_sympy(actual)==expected
    assert sm.from_sympy(expected)==actual
    assert sm.sympify(actual.srepr())==actual
    assert sm.sympify(str(actual))==actual
    restored=pickle.loads(pickle.dumps(actual))
    assert restored==actual and restored.args==actual.args
    assert sm.default_sort_key(actual)==sp.default_sort_key(expected)
    assert sm.default_sort_key(2*actual)==sp.default_sort_key(2*expected)


@pytest.mark.parametrize('operation',[
    lambda s,x,y:s.Integral(x**2,x),
    lambda s,x,y:s.Integral(x,(x,2)),
    lambda s,x,y:s.Integral(x,(x,0,y)),
    lambda s,x,y:s.Integral(x+y,(x,0,y),(y,0,2)),
    lambda s,x,y:s.Sum(x,(x,2,-1)),
    lambda s,x,y:s.Sum(x,(x,3,0)),
    lambda s,x,y:s.Sum(s.sin(x),(x,3,0)),
    lambda s,x,y:s.Sum(s.Function('f')(x),(x,1,0)),
    lambda s,x,y:s.Sum(s.Function('f')(x),(x,y,y+2)),
    lambda s,x,y:s.Sum(s.sin(x),(x,y+2,y-1)),
    lambda s,x,y:s.Sum(x*y,(x,1,2),(y,1,3)),
])
def test_doit_uses_native_engines_and_exact_karr_endpoint_convention(operation):
    x,y=sm.symbols('x y');X,Y=sp.symbols('x y')
    actual=operation(sm,x,y).doit();expected=operation(sp,X,Y).doit()
    assert sp.simplify(sm.to_sympy(actual)-expected)==0


@pytest.mark.parametrize('operation,old,new',[
    (lambda s,x,y,z:s.Integral(s.Function('f')(x),(x,0,y)),'x','z'),
    (lambda s,x,y,z:s.Integral(s.Function('f')(x),(x,0,x)),'x','z'),
    (lambda s,x,y,z:s.Integral(s.Function('f')(x),x),'x',2),
    (lambda s,x,y,z:s.Integral(s.Function('f')(x),x,x),'x',2),
    (lambda s,x,y,z:s.Integral(s.Function('f')(x),x,x),'x','z'),
    (lambda s,x,y,z:s.Integral(x+y,(x,0,y),(y,0,z)),'y',2),
    (lambda s,x,y,z:s.Integral(x+y,(x,0,y),(y,0,z)),'z','y'),
    (lambda s,x,y,z:s.Sum(s.Function('f')(x)+y,(x,1,z)),'x','z'),
    (lambda s,x,y,z:s.Sum(s.Function('f')(x)+y,(x,1,z)),'y','x'),
])
def test_substitution_observes_definite_and_indefinite_binding_scopes(operation,old,new):
    x,y,z=sm.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    native={'x':x,'y':y,'z':z};oracle={'x':X,'y':Y,'z':Z}
    actual=operation(sm,x,y,z).subs(native[old],native.get(new,new))
    expected=operation(sp,X,Y,Z).subs(oracle[old],oracle.get(new,new))
    assert sm.to_sympy(actual)==expected


@pytest.mark.parametrize('kind',['Integral','Sum'])
def test_simultaneous_and_exact_replacements_have_distinct_scope_rules(kind):
    x,y=sm.symbols('x y');X,Y=sp.symbols('x y')
    expression=getattr(sm,kind)(x+y,(x,0,y))
    expected=getattr(sp,kind)(X+Y,(X,0,Y))
    assert sm.to_sympy(expression.subs({x:y,y:x},simultaneous=True))==expected.subs({X:Y,Y:X},simultaneous=True)
    assert sm.to_sympy(expression.xreplace({x:y}))==expected.xreplace({X:Y})
    with pytest.raises(ValueError):expression.xreplace({x:2})


@pytest.mark.parametrize('operation,variable',[
    (lambda s,x,y:s.Integral(s.Function('f')(x),x),'x'),
    (lambda s,x,y:s.Integral(s.Function('f')(x),(x,0,x)),'x'),
    (lambda s,x,y:s.Integral(s.Function('f')(x),(x,0,y)),'y'),
    (lambda s,x,y:s.Integral(y*s.Function('f')(x),(x,0,y)),'y'),
    (lambda s,x,y:s.Integral(s.Function('f')(x),(x,y,2*y)),'y'),
    (lambda s,x,y:s.Integral(x**2,(x,0,y)),'y'),
    (lambda s,x,y:s.Integral(s.Function('f')(x),(x,0,1)),'x'),
    (lambda s,x,y:s.Sum(x*y,(x,1,4)),'y'),
    (lambda s,x,y:s.Sum(x,(x,1,y)),'y'),
])
def test_ftc_leibniz_and_held_discrete_endpoint_derivatives(operation,variable):
    x,y=sm.symbols('x y');X,Y=sp.symbols('x y')
    actual=sm.diff(operation(sm,x,y),{'x':x,'y':y}[variable])
    expected=sp.diff(operation(sp,X,Y),{'x':X,'y':Y}[variable])
    assert sm.to_sympy(actual)==expected


def test_unknown_operations_remain_typed_and_do_not_collide_with_user_functions():
    x=sm.Symbol('x');f=sm.Function('f')
    held=sm.Integral(f(x),x)
    assert sm.integrate(f(x),x)==held and held.doit()==held
    unknown=sm.Sum(f(x),(x,0,sm.oo))
    assert sm.summation(f(x),(x,0,sm.oo))==unknown and unknown.doit()==unknown
    for name in ['Integral','Sum']:
        user=sm.Function(name)(x)
        assert user.is_Function and user.is_UndefinedFunction
        assert not isinstance(user,getattr(sm,name))
        assert sm.from_sympy(sm.to_sympy(user))==user
        assert pickle.loads(pickle.dumps(user))==user
        for attr in ['function','limits','bound_symbols']:
            with pytest.raises(AttributeError):getattr(user,attr)
    assert '\\int' in sm.latex(held)
    assert '\\sum' in sm.latex(unknown)
    with pytest.raises((ValueError,TypeError,NotImplementedError)):sm.lambdify(x,held,'numpy')


def test_limits_keep_symbol_assumptions_and_dummy_identity():
    bound=sm.Dummy('x',positive=True);free=sm.Symbol('x')
    held=sm.Integral(bound+free,(bound,0,1))
    assert held.free_symbols=={free}
    assert sm.from_sympy(sm.to_sympy(held))==held
    assert pickle.loads(pickle.dumps(held)).limits[0][0]==bound
    plain=sm.Symbol('x');positive=sm.Symbol('x',positive=True)
    expression=sm.Sum(positive,(plain,1,2))
    assert expression.free_symbols=={positive}


def test_constructor_inference_none_limits_and_explicit_unsupported_inputs():
    x,y=sm.symbols('x y');X,Y=sp.symbols('x y')
    assert sm.Integral(x)==sm.Integral(x,x)
    assert sm.Integral(sm.Integral(x,x),y).limits==((x,),(y,))
    for bounds in [(None,y),(y,None),(None,None)]:
        assert sm.to_sympy(sm.Integral(x,(x,*bounds)))==sp.Integral(X,(X,*(Y if b==y else b for b in bounds)))
    for call in [lambda:sm.Integral(1),lambda:sm.Integral(x+y),lambda:sm.Integral(x,2),lambda:sm.Sum(x,x),lambda:sm.Sum(x,(x,1)),lambda:sm.Sum(x,(x,1,2,3))]:
        with pytest.raises((TypeError,ValueError)):call()
    with pytest.raises(TypeError,match=r"^Basic.__new__\(\) got an unexpected keyword argument 'risch'$"):sm.Integral(x,x,risch=True)


def test_native_binders_work_without_importing_sympy():
    program=r'''
import sys,pickle
class BlockSympy:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('native calculus imported SymPy')
sys.meta_path.insert(0,BlockSympy())
import symox as s
x,y=s.symbols('x y')
f=s.Function('f')
e=s.Integral(f(x),(x,0,y))
assert e.free_symbols=={y}
assert s.diff(e,y)==f(y)
assert e.subs(x,2)==e
assert s.sympify(e.srepr())==e
assert pickle.loads(pickle.dumps(e))==e
assert s.Sum(x,(x,2,-1)).doit()==-1
assert s.Integral(x,(x,0,2)).doit()==2
'''
    subprocess.run([sys.executable,'-c',program],check=True,capture_output=True,text=True,timeout=15)


def test_composite_substitutions_preserve_binding_scopes_and_reject_new_dummy_dependencies():
    x,y,z=sm.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    f=sm.Function('f');F=sp.Function('f')
    for old,OLD in [(f(x),F(X)),(f(x)+y,F(X)+Y),(x*y,X*Y),(sm.Integer(1),sp.S.One)]:
        actual=sm.Integral(f(x)+x*y+y,(x,0,y)).subs(old,z)
        expected=sp.Integral(F(X)+X*Y+Y,(X,0,Y)).subs(OLD,Z)
        assert sm.to_sympy(actual)==expected
    with pytest.raises(ValueError):sm.Integral(f(y),(x,0,1)).subs(f(y),x)
    # Like the reference, a variable replaced by a function application
    # becomes the new integration variable.
    actual=sm.Integral(f(x),x).subs(x,f(y))
    assert sm.to_sympy(actual)==sp.Integral(F(X),X).subs(X,F(Y))


def test_functional_endpoint_integral_and_large_integer_karr_sum():
    x=sm.Symbol('x')
    assert sm.integrate(x,(x,2))==sm.Integral(x,(x,2)).doit()==2
    with pytest.raises(ValueError):sm.summation(x)
    offset=10**100
    assert sm.Sum(x,(x,offset+2,offset-1)).doit()==-2*offset-1
    held=sm.Sum(1,(x,sm.oo,sm.oo))
    assert held.doit()==held  # unsupported doubly infinite convention stays explicit


def test_internal_apply_and_legacy_restore_validate_and_normalize_limits():
    from symox import _core
    from symox._compat import _restore_expr
    x=sm.Symbol('x');f=sm.Function('f')
    assert _core._apply('Integral',f(x),x)==sm.Integral(f(x),x)
    assert _restore_expr('Integral',(f(x),x))==sm.Integral(f(x),x)
    assert _restore_expr('Sum',(f(x),(x,0,sm.oo)))==sm.Sum(f(x),(x,0,sm.oo))
    for name,args in [('Integral',()),('Sum',()),('Integral',(x,2)),('Sum',(x,x)),('Sum',(x,(x,1)) )]:
        with pytest.raises((TypeError,ValueError)):_core._apply(name,*args)
