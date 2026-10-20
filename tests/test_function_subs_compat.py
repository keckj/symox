"""Undefined function identities, held calculus, and exact substitutions."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm


@pytest.mark.parametrize("name",["f","sin","exp","Tuple","Subs","UnevaluatedExpr","Order"])
def test_function_heads_keep_identity_including_builtin_collisions(name):
    x=sm.Symbol("x")
    f=sm.Function(name)
    applied=f(x)
    assert applied.func==f and hash(applied.func)==hash(f)
    assert applied.func(*applied.args)==applied
    assert applied.has(f)
    assert applied.is_UndefinedFunction
    assert not isinstance(applied,sm.Subs)
    assert sm.diff(applied,x).is_Derivative
    assert pickle.loads(pickle.dumps(f))==f
    assert pickle.loads(pickle.dumps(applied))==applied
    assert sm.sympify(applied.srepr())==applied
    assert sm.from_sympy(sm.to_sympy(applied))==applied
    assert sm.from_sympy(sm.to_sympy(f))==f


@pytest.mark.parametrize("flags",[{"real":True},{"real":False},{"positive":True}])
def test_function_metadata_nargs_and_head_hash_roundtrip(flags):
    x,y=sm.symbols("x y")
    f=sm.Function("f",nargs=(2,1,2),**flags)
    assert f==sm.Function("f",nargs=(1,2),**flags)
    assert f.nargs==sm.FiniteSet(1,2)
    assert f(x).func==f and hash(f(x).func)==hash(f)
    expected=sp.Function("f",nargs=(1,2),**flags)(sp.Symbol("x"))
    assert f(x).is_real==expected.is_real
    assert sm.from_sympy(expected)==f(x)
    assert sm.from_sympy(sm.to_sympy(f))==f
    for value in [f(x),f(x,y)]:
        assert pickle.loads(pickle.dumps(value))==value
        assert sm.sympify(value.srepr())==value
        assert sm.to_sympy(value)==sp.Function("f",nargs=(1,2),**flags)(*map(sm.to_sympy,value.args))
    with pytest.raises(TypeError):f()
    with pytest.raises(TypeError):f(x,y,x)
    with pytest.raises(ValueError):sm.Function("f",real=True,imaginary=True)
    with pytest.raises(ValueError):sm.Function("f",nargs=-1)


def test_explicit_undefined_sin_does_not_evaluate_or_print_as_builtin():
    x=sm.Symbol("x")
    f=sm.Function("sin")
    assert f(x)!=sm.sin(x)
    assert f(1).evalf(50).is_UndefinedFunction
    with pytest.raises((ValueError,TypeError,NotImplementedError)):sm.lambdify(x,f(x),"numpy")
    # SymPy's code printers look known functions up by name.
    assert sm.cxxcode(f(x))=="sin(x)"
    from symox.printing.cxx import CXX17CodePrinter
    assert CXX17CodePrinter().doprint(f(x))=="sin(x)"
    assert CXX17CodePrinter({"user_functions":{"sin":"custom_sin"}}).doprint(f(x))=="custom_sin(x)"


def test_typed_subs_binding_properties_evaluation_and_roundtrips():
    x,y,z=sm.symbols("x y z")
    held=sm.Subs(x*x+y,x,z)
    assert isinstance(held,sm.Subs) and held.is_Subs and not held.is_Function
    assert held.free_symbols=={y,z}
    assert held.variables==(x,) and held.point==(z,) and held.expr==x*x+y
    assert held.func(*held.args)==held
    assert held.doit()==z*z+y
    assert sm.Subs(x+y,(x,y),(y,2)).doit()==4  # prescribed Subs evaluation order
    assert sm.Subs(sm.Derivative(x*x,x),x,2).doit(deep=False)==4
    for value in [held,sm.Subs(sm.Derivative(sm.Function("f")(x),x),x,0)]:
        assert pickle.loads(pickle.dumps(value))==value
        assert sm.sympify(value.srepr())==value
        assert sm.from_sympy(sm.to_sympy(value))==value
    with pytest.raises(ValueError):sm.Subs(x,(x,y),(1,))
    with pytest.raises(ValueError):sm.Subs(x,(x,x),(1,2))
    with pytest.raises(NotImplementedError):sm.Subs(x*x,x*x,y)
    assert "right|" in sm.latex(held)


def test_derivative_grouped_args_evaluate_and_substitution():
    x,y=sm.symbols("x y")
    f=sm.Function("f")
    held=sm.Derivative(f(x,y),x,y,x)
    assert held.variable_count==((x,1),(y,1),(x,1))
    assert held.variables==(x,y,x) and held.derivative_count==3
    assert held.args[1].args==(x,1)
    assert held.func(*held.args)==held
    assert sm.Derivative(x*x,x,evaluate=True)==2*x
    assert sm.Derivative(x*x,x,evaluate=False).doit()==2*x
    assert sm.Derivative(x*x,x).subs(x,2).doit()==4
    assert sm.diff(f(x*x),x).subs(x,0).doit()==0
    assert sm.Derivative(f(x),x).subs(x,y)==sm.Derivative(f(y),y)
    assert isinstance(sm.Derivative(f(x,y),x).subs(x,y),sm.Subs)
    assert sm.Derivative(f(x,y),x).subs(y,x)==sm.Subs(sm.Derivative(f(x,y),x),y,x)
    with pytest.raises(ValueError):sm.Derivative(x*x,x*x)


def test_simultaneous_symbol_substitution_preserves_rhs_and_poles():
    x,y=sm.symbols("x y")
    assert (x+y).subs({x:y,y:x},simultaneous=True)==x+y
    assert (x/y).subs({x:0,y:0},simultaneous=True)==sm.nan
    assert (x+y).subs([(x,y),(y,2)],simultaneous=True)==y+2
    assert (x+y).subs([(x,y),(y,2)])==4
    assert (x*x).subs({x*x:y},simultaneous=True) == y
    with pytest.raises(NotImplementedError): x.subs({x:y},unrecognized=True)


def test_simultaneous_substitution_preserves_bound_operators():
    x,y=sm.symbols("x y")
    f=sm.Function("f")
    derivative=sm.Derivative(f(x),x)
    result=derivative.subs({x:2},simultaneous=True)
    assert result==sm.Subs(derivative,x,2)
    assert result.free_symbols==set()
    bound=sm.Subs(f(x),x,y)
    result=bound.subs({y:x,x:y},simultaneous=True)
    assert result==sm.Subs(f(x),x,x)
    assert result.free_symbols=={x}
    sx,sy=sp.symbols("x y")
    sf=sp.Function("f")
    assert sm.to_sympy(result)==sp.Subs(sf(sx),sx,sy).subs({sy:sx,sx:sy},simultaneous=True)
    z=sm.Symbol("z")
    ordered=sm.Subs(x+y,(x,y),(y,z))
    assert ordered.subs({y:x},simultaneous=True)==ordered
    assert sm.Subs(f(x),x,x).subs({x:z},simultaneous=True)==sm.Subs(f(x),x,z)


def test_undefined_order_has_no_order_properties():
    x=sm.Symbol("x")
    value=sm.Function("Order")(x)
    assert not value.is_Order
    for attribute in ["expr","variables","point"]:
        with pytest.raises(AttributeError):getattr(value,attribute)


@pytest.mark.parametrize("name",["Symbol","Add","Mul","Pow","UnevaluatedExpr"])
def test_printer_builtin_handlers_do_not_capture_undefined_names(name):
    from symox.printing.cxx import CXX17CodePrinter
    from symox.printing.str import StrPrinter
    x=sm.Symbol("x")
    applied=sm.Function(name)(x)
    assert StrPrinter().doprint(applied)==name+"(x)"
    with pytest.raises(NotImplementedError,match="Unsupported by"):CXX17CodePrinter().doprint(applied)
    assert CXX17CodePrinter({"user_functions":{name:"custom"}}).doprint(applied)=="custom(x)"


def test_undefined_function_limits_stay_held_without_continuity_assumptions():
    x=sm.Symbol("x")
    for f in [sm.Function("f"),sm.Function("f",real=True),sm.Function("sin")]:
        result=sm.limit(f(x),x,0)
        assert isinstance(result,sm.Limit)
        assert result.doit()==result
    assert sm.limit(sm.sin(x)/x,x,0)==1


def test_new_symbolic_apis_have_no_sympy_runtime_dependency():
    subprocess.run([sys.executable,"-c","""
import sys
sys.modules['sympy']=None
import symox as s
x,y=s.symbols('x y')
f=s.Function('sin',real=True,nargs=1)
assert f(x).func==f and f(x).is_real is True
assert s.diff(f(x),x).is_Derivative
assert s.Subs(s.Derivative(x*x,x),x,2).doit()==4
assert (x+y).subs({x:y,y:x},simultaneous=True)==x+y
assert isinstance(s.limit(f(x),x,0),s.Limit)
"""],check=True,timeout=30)
