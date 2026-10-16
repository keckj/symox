"""Lambda signatures, simultaneous application and native interoperation."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def same(actual, expected):
    assert sx.to_sympy(actual) == expected


@pytest.mark.parametrize("signature", ["x", "(x,y)", "(x,(y,z))", "((x,y),z)", "()", "((),x)"])
def test_signature_structure_free_symbols_and_reconstruction(signature):
    actual = sx.sympify(f"Lambda({signature},x+2*y+3*z+a)")
    expected = sp.sympify(f"Lambda({signature},x+2*y+3*z+a)")
    assert isinstance(actual, sx.Lambda)
    same(actual.signature, expected.signature)
    assert tuple(sx.to_sympy(v) for v in actual.variables) == expected.variables
    assert actual.bound_symbols == actual.variables
    same(actual.nargs, expected.nargs)
    assert {sx.to_sympy(v) for v in actual.free_symbols} == expected.free_symbols
    same(actual.expr, expected.expr)
    assert actual.func(*actual.args) == actual
    assert sx.from_sympy(expected) == actual
    assert sx.sympify(actual.srepr()) == actual
    assert sx.sympify(str(actual)) == actual
    assert sp.sympify(actual.srepr()) == expected
    assert str(actual) == str(expected)
    assert sx.default_sort_key(actual) == sp.default_sort_key(expected)
    assert actual.is_identity == expected.is_identity
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        restored = pickle.loads(pickle.dumps(actual, protocol))
        assert restored == actual
        assert restored.args == actual.args


@pytest.mark.parametrize("signature,arguments", [
    ("x", "(y,)"), ("(x,y)", "(y,x)"), ("(x,y)", "(y+1,x-1)"),
    ("(x,(y,z))", "(1,(2,3))"), ("((x,y),z)", "((1,2),3)"),
    ("()", "()"), ("((),x)", "((),4)"),
])
def test_application_is_exact_simultaneous_and_destructures(signature, arguments):
    x, y, z = sx.symbols("x y z")
    X, Y, Z = sp.symbols("x y z")
    native_args = eval(arguments, {}, dict(x=x, y=y, z=z))
    oracle_args = eval(arguments, {}, dict(x=X, y=Y, z=Z))
    actual = sx.sympify(f"Lambda({signature},x+2*y+3*z)")
    expected = sp.sympify(f"Lambda({signature},x+2*y+3*z)")
    same(actual(*native_args), expected(*oracle_args))


@pytest.mark.parametrize("signature", ["(x,x)", "(x,(y,x))", "(x+1,y)", "(1,y)", "3"])
def test_bad_signature_rejected_by_public_constructor(signature):
    x, y = sx.symbols("x y")
    X, Y = sp.symbols("x y")
    native_signature = eval(signature, {}, dict(x=x, y=y))
    oracle_signature = eval(signature, {}, dict(x=X, y=Y))
    with pytest.raises(sp.core.function.BadSignatureError):
        sp.Lambda(oracle_signature, X + Y)
    with pytest.raises(sx.BadSignatureError):
        sx.Lambda(native_signature, x + y)


@pytest.mark.parametrize("arguments", [(), (1,), (1, 2, 3)])
def test_bad_top_level_arity(arguments):
    x, y = sx.symbols("x y")
    value = sx.Lambda((x, y), x-y)
    with pytest.raises(sx.BadArgumentsError, match="takes exactly 2 arguments"):
        value(*arguments)


@pytest.mark.parametrize("arguments", [(1, 2), (1, (2,)), (1, (2, 3, 4))])
def test_tuple_argument_requires_matching_shape(arguments):
    x, y, z = sx.symbols("x y z")
    with pytest.raises(sx.BadArgumentsError):
        sx.Lambda((x, (y, z)), x+y+z)(*arguments)


def test_public_constructor_keywords_native_tuple_and_empty_signature():
    x, y = sx.symbols("x y")
    value = sx.Lambda(signature=sx.Tuple(x, y), expr=x-y)
    assert value(y, x) == y-x
    assert sx.Lambda((), 7)() == 7
    with pytest.raises(TypeError):
        sx.Lambda(x, x+1, evaluate=False)


@pytest.mark.parametrize("method", ["subs", "xreplace"])
def test_structural_substitution_updates_signature_and_preserves_validation(method):
    x, y, z, a = sx.symbols("x y z a")
    X, Y, Z, A = sp.symbols("x y z a")
    actual = sx.Lambda((x, y), x+y+z)
    expected = sp.Lambda((X, Y), X+Y+Z)
    same(getattr(actual, method)({x: a}), getattr(expected, method)({X: A}))
    same(getattr(actual, method)({z: x}), getattr(expected, method)({Z: X}))
    for replacement in (2, y):
        with pytest.raises(sx.BadSignatureError):
            getattr(actual, method)({x: replacement})
    pair = sx.Function("pair")(actual, actual)
    with pytest.raises(sx.BadSignatureError):
        getattr(pair, method)({x: 2})


def test_call_preserves_sympy_capture_contract_and_nested_free_symbols():
    x, y, z = sx.symbols("x y z")
    X, Y, Z = sp.symbols("x y z")
    actual = sx.Lambda(x, sx.Lambda(y, x+y+z))
    expected = sp.Lambda(X, sp.Lambda(Y, X+Y+Z))
    assert actual.free_symbols == {z}
    same(actual(y), expected(Y))
    f, F = sx.Function("f"), sp.Function("f")
    actual = sx.Lambda(x, sx.Integral(f(x, y), (y, 0, 1)))
    expected = sp.Lambda(X, sp.Integral(F(X, Y), (Y, 0, 1)))
    same(actual(y), expected(Y))
    same(sx.Lambda(x, sx.Subs(f(x), x, y))(z), sp.Lambda(X, sp.Subs(F(X), X, Y))(Z))


def test_has_diff_and_doit_use_lambda_structure_and_free_variables():
    x, y, z = sx.symbols("x y z")
    X, Y, Z = sp.symbols("x y z")
    actual = sx.Lambda((x, y), x*y+z)
    expected = sp.Lambda((X, Y), X*Y+Z)
    assert actual.has(x) and actual.has(y) and actual.has(z)
    assert actual.diff(x) == expected.diff(X) == 0
    same(actual.diff(z), expected.diff(Z))
    actual = sx.Lambda(x, x+sx.Integral(y, (y, 0, 1)))
    expected = sp.Lambda(X, X+sp.Integral(Y, (Y, 0, 1)))
    same(actual.doit(), expected.doit())


def test_lambda_equality_is_structural_with_independent_dummy_identity():
    x, y, z = sx.symbols("x y z")
    assert sx.Lambda(x, x+z) != sx.Lambda(y, y+z)
    d, e = sx.Dummy("x"), sx.Dummy("x")
    assert sx.Lambda(d, d+z) != sx.Lambda(e, e+z)
    assert sx.Lambda(d, d+z).free_symbols == {z}
    assert len({sx.Lambda(x, x+z), sx.Lambda(y, y+z)}) == 2


def test_identity_function_singleton_native_api_interop_and_pickles():
    x, y = sx.symbols("x y")
    actual = sx.Lambda(x, x)
    expected = sp.S.IdentityFunction
    assert actual is sx.S.IdentityFunction
    assert sx.Lambda(y, y) is actual
    assert sx.IdentityFunction() is actual
    assert isinstance(actual, sx.Lambda) and isinstance(actual, sx.IdentityFunction)
    assert actual.func is sx.IdentityFunction
    assert actual.args == ()
    assert actual.free_symbols == set()
    assert actual.nargs == sx.FiniteSet(1)
    assert actual.variables == actual.bound_symbols
    assert len(actual.variables) == 1 and actual.variables[0].is_Dummy
    assert actual.signature == sx.Tuple(*actual.variables)
    assert actual.expr == actual.variables[0]
    assert not actual.is_identity
    assert actual(x+y) == x+y
    assert not actual.has(actual.variables[0])
    assert actual.xreplace({actual.variables[0]: 2}) == actual
    assert sx.to_sympy(actual) is expected
    assert sx.from_sympy(expected) == actual
    assert actual.srepr() == sp.srepr(expected)
    assert str(actual) == str(expected)
    assert sx.default_sort_key(actual) == sp.default_sort_key(expected)
    assert sx.sympify(actual.srepr()) == actual
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        assert pickle.loads(pickle.dumps(actual, protocol)) == actual
    with pytest.raises(sx.BadArgumentsError):
        actual()
    with pytest.raises(TypeError):
        sx.IdentityFunction(1)


@pytest.mark.parametrize("name", ["Lambda", "IdentityFunction", "FunctionMatrix"])
def test_undefined_function_head_collisions_remain_ordinary_functions(name):
    x, y = sx.symbols("x y")
    X, Y = sp.symbols("x y")
    actual = sx.Function(name, nargs=2)(x, y)
    expected = sp.Function(name, nargs=(2,))(X, Y)
    assert not actual.is_Lambda
    assert actual.func.name == name
    assert actual.func.nargs == sx.FiniteSet(2)
    assert not isinstance(actual, sx.Lambda)
    assert actual.free_symbols == {x, y}
    same(actual, expected)
    assert sx.from_sympy(expected) == actual
    assert pickle.loads(pickle.dumps(actual)) == actual
    assert sx.sympify(actual.srepr()) == actual
    for field in ("signature", "variables", "bound_symbols", "expr"):
        with pytest.raises(AttributeError):
            getattr(actual, field)


def test_lambda_native_without_sympy_and_shared_dag_application():
    program = r'''
import builtins, pickle, sys
original = builtins.__import__
def block(name,*args,**kwargs):
    if name == 'sympy' or name.startswith('sympy.'):
        raise AssertionError('SymPy runtime fallback')
    return original(name,*args,**kwargs)
builtins.__import__=block
import symox as s
x,y,z=s.symbols('x y z')
f=s.Function('pair')
body=x+y
for _ in range(40):body=f(body,body)
value=s.Lambda((x,y),body)
result=value(y,x)
for _ in range(40):
    assert result.args[0]==result.args[1]
    result=result.args[0]
assert result==x+y
value=pickle.loads(pickle.dumps(s.Lambda((x,(y,z)),x+2*y+z)))
assert value(1,(2,3))==8
assert s.Lambda(x,x)(7)==7
assert not any(k=='sympy' or k.startswith('sympy.') for k in sys.modules)
'''
    subprocess.run([sys.executable, "-c", program], check=True, timeout=15)


def test_tuple_body_is_retained_and_roundtrips_through_both_printers():
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    value = sx.Lambda((x, y), sx.Tuple(x, y))
    expected = sp.Lambda((X, Y), (X, Y))
    same(value, expected)
    assert value.is_identity == expected.is_identity
    same(value(1, 2), expected(1, 2))
    assert sx.sympify(str(value)) == value
    assert sx.sympify(value.srepr()) == value


def test_application_rejects_invalid_calculus_binding_variables():
    x = sx.Symbol('x')
    X = sp.Symbol('x')
    f, F = sx.Function('f'), sp.Function('f')
    with pytest.raises(ValueError):
        sp.Lambda(X, sp.Integral(F(X), X))(2)
    with pytest.raises(ValueError):
        sx.Lambda(x, sx.Integral(f(x), x))(2)


def test_tuple_factory_preserves_symbolic_type_and_func_reconstruction():
    x, y = sx.symbols('x y')
    value = sx.Tuple(x, (y, 2))
    assert isinstance(value, sx.Tuple)
    assert value.is_Tuple
    assert value.func is sx.Tuple
    assert value.func(*value.args) == value
    same(value, sp.Tuple(sp.Symbol('x'), (sp.Symbol('y'), 2)))
    undefined = sx.Function('Tuple')(x, y)
    assert not isinstance(undefined, sx.Tuple)
    assert not undefined.is_Tuple
    assert pickle.loads(pickle.dumps(value)) == value
