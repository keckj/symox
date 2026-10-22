"""Held elementwise Lambda application on native noncommutative matrix DAGs."""
import ast
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sm
from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction as SymPyApply

import symox as sx


def equal_entries(actual, expected):
    assert actual.shape == expected.shape
    assert all(sx.expand(a-b) == 0 for a,b in zip(actual, expected))


def test_class_arguments_properties_and_native_constructor_method():
    z, t = sx.symbols('z t')
    A = sx.MatrixSymbol('A', 2, 3)
    function = sx.Lambda(z, z*z+t)
    value = sx.ElementwiseApplyFunction(function, A)
    assert isinstance(value, sx.ElementwiseApplyFunction)
    assert isinstance(value, sx.MatrixExpr)
    assert value._kind_tag == 'ElementwiseApplyFunction'
    assert value.is_ElementwiseApplyFunction
    assert value.args == (function, A)
    assert value.function == function and value.expr == A
    assert value.shape == A.shape
    assert value.func is sx.ElementwiseApplyFunction
    assert value.func(*value.args) == value
    assert A.applyfunc(function) == value
    assert value[0, 1] == A[0, 1]**2+t
    assert value[-1, -1] == A[1, 2]**2+t
    assert value.free_symbols == {A, t}
    assert z not in value.free_symbols
    assert sx.sympify(sx.srepr(value)) == value


@pytest.mark.parametrize('rows,cols', [(0, 0), (0, 3), (3, 0), (1, 1), (2, 3), (3, 2)])
def test_dense_evaluation_empty_shapes_and_negative_indices(rows, cols):
    z = sx.Symbol('z')
    parent = sx.MatrixExpr(sx.ImmutableMatrix(rows, cols, [r+10*c for r in range(rows) for c in range(cols)]))
    value = parent.applyfunc(sx.Lambda(z, z*z+1))
    expected = sx.ImmutableMatrix(rows, cols, [(r+10*c)**2+1 for r in range(rows) for c in range(cols)])
    assert value.shape == (rows, cols)
    assert value.as_explicit() == expected
    assert value.doit() == expected
    assert value.doit().shape == (rows, cols)
    assert value.doit(deep=False) == expected
    if rows and cols:
        for r in range(-rows, rows):
            for c in range(-cols, cols):
                assert value[r, c] == expected[r, c]
    for index in [(rows, 0), (0, cols), (-rows-1, 0), (0, -cols-1)]:
        with pytest.raises(IndexError):
            value[index]


def test_sympy_interchange_with_named_lambda_and_held_one_by_one_parent():
    z, t = sx.symbols('z t')
    u, v = sm.symbols('z t')
    for rows, cols in [(2, 3), (0, 3), (1, 1)]:
        A, a = sx.MatrixSymbol('A', rows, cols), sm.MatrixSymbol('A', rows, cols)
        value = A.applyfunc(sx.Lambda(z, z*z+t))
        # SymPy's public constructor can replace a 1x1 elementwise expression
        # by a matrix power. Preserve the held class explicitly for the oracle.
        reference = sm.Basic.__new__(SymPyApply, sm.Lambda(u, u*u+v), a)
        assert value.to_sympy() == reference
        assert sx.sympify(reference) == value
        assert sx.MatrixExpr(reference) == value
        assert sx.default_sort_key(value) == sm.default_sort_key(reference)
        assert sx.sympify(sx.srepr(value)) == value


def test_symbolic_dimensions_substitutions_and_bound_parameter_scope():
    n, m = sx.symbols('n m', integer=True, positive=True)
    z, t, u = sx.symbols('z t u')
    A = sx.MatrixSymbol('A', n, m)
    value = A.applyfunc(sx.Lambda(z, z*z+t))
    assert value.shape == (n, m)
    assert value.free_symbols == {A, t}
    B = sx.MatrixSymbol('B', 2, 3)
    transformed = value.xreplace({A: B, t: u})
    assert transformed.shape == (2, 3)
    assert transformed[1, 2] == B[1, 2]**2+u
    assert transformed.free_symbols == {B, u}
    assert value.subs({A: B, t: u}) == transformed
    # The same-named variable is bound in the function but free in its parent.
    dependent = (t*B).applyfunc(sx.Lambda(t, t*t))
    equal_entries(dependent.diff(t).as_explicit(), dependent.as_explicit().applyfunc(lambda entry: entry.diff(t)))
    assert dependent.free_symbols == {B, t}


def test_chain_rule_includes_explicit_function_parameter_dependence():
    z = sx.Symbol('z')
    t = sx.Symbol('t', real=True)
    i, j = sx.symbols('i j', integer=True)
    parent = sx.FunctionMatrix(2, 3, sx.Lambda((i, j), t*(i+2*j+1)))
    value = parent.applyfunc(sx.Lambda(z, z*z+t*z+t*t))
    for order in [1, 2]:
        actual = value.diff(t, order).as_explicit()
        expected = value.as_explicit().applyfunc(lambda entry: entry.diff(t, order))
        equal_entries(actual, expected)
    A = sx.MatrixSymbol('A', 2, 3)
    # SymPy currently returns ZeroMatrix for this parameter-only derivative.
    # Its scalar dense entries give an independent mathematical certificate.
    parameter_only = A.applyfunc(sx.Lambda(z, z+t))
    equal_entries(parameter_only.diff(t).as_explicit(), sx.ImmutableMatrix.ones(2, 3))
    assert parameter_only.diff(sx.Symbol('missing')) == sx.ZeroMatrix(2, 3)
    constant = parent.applyfunc(sx.Lambda(z, t))
    equal_entries(constant.diff(t).as_explicit(), sx.ImmutableMatrix.ones(2, 3))


def test_identity_and_nested_composition_evaluate_exactly_without_alpha_equality():
    z, t, u = sx.symbols('z t u')
    A = sx.MatrixSymbol('A', 2, 3)
    identity = A.applyfunc(sx.Lambda(z, z))
    assert isinstance(identity, sx.ElementwiseApplyFunction)
    assert identity.doit() == A
    assert identity.doit(deep=False) == A
    assert identity.diff(t) == A.diff(t)
    inner = A.applyfunc(sx.Lambda(z, z+u))
    nested = inner.applyfunc(sx.Lambda(z, z*z+t))
    result = nested.doit()
    assert isinstance(result, sx.ElementwiseApplyFunction)
    assert result.expr == A
    assert result.free_symbols == {A, t, u}
    equal_entries(result.as_explicit(), nested.as_explicit())
    equal_entries(result.as_explicit(), A.as_explicit().applyfunc(lambda entry: (entry+u)**2+t))
    assert sx.sympify(sx.srepr(result)) == result
    assert pickle.loads(pickle.dumps(result)) == result


def test_transpose_adjoint_conjugation_and_reverse_slice_keep_exact_entries():
    z, t = sx.symbols('z t')
    parent = sx.MatrixExpr(sx.ImmutableMatrix([[1+sx.I, 2-sx.I], [3*sx.I, 4]]))
    value = parent.applyfunc(sx.Lambda(z, sx.conjugate(z)+z*z+t))
    expected = value.as_explicit()
    assert value.T.as_explicit() == expected.T
    assert value.H.as_explicit() == expected.H
    assert value.conjugate().as_explicit() == expected.conjugate()
    assert value[::-1, ::-1].as_explicit() == sx.ImmutableMatrix([[expected[1, 1], expected[1, 0]], [expected[0, 1], expected[0, 0]]])
    assert value.T.function == value.function
    assert value.T.expr == parent.T


def test_expansion_visits_function_and_parent_without_evaluating_the_wrapper():
    z, x, y = sx.symbols('z x y')
    A, B, C = [sx.MatrixSymbol(name, 2, 2) for name in 'ABC']
    value = ((A+B)*C).applyfunc(sx.Lambda(z, (x+y)**2+z))
    expanded = value.expand()
    assert isinstance(expanded, sx.ElementwiseApplyFunction)
    assert expanded.expr == A*C+B*C
    assert expanded.function(0) == x*x+2*x*y+y*y
    assert value.expand(deep=False) == value
    assert value.expand(mul=False).expr == value.expr
    assert value.expand(multinomial=False).function == value.function
    assert sx.expand(value) == expanded
    equal_entries(expanded.as_explicit(), value.as_explicit())


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL+1))
def test_flat_codec_keeps_held_application_and_mixed_scalar_metadata(protocol):
    z, t = sx.symbols('z t')
    A = sx.MatrixSymbol('A', 2, 3)
    value = A.applyfunc(sx.Lambda(z, z*z+t))
    for expression in [value, value.T, value[::-1, ::-1], sx.Trace(value.T*value), sx.Determinant(value.T*value)]:
        restored = pickle.loads(pickle.dumps(expression, protocol=protocol))
        assert restored == expression
        assert restored.args == expression.args
        assert sx.sympify(sx.srepr(restored)) == restored


@pytest.mark.parametrize('name', ['sin', 'exp', 'Abs', 'conjugate', 'log'])
def test_known_symbolic_heads_are_normalized_to_native_unary_lambdas(name):
    head = getattr(sx, name)
    A = sx.MatrixSymbol('A', 2, 3)
    value = A.applyfunc(head)
    assert isinstance(value, sx.ElementwiseApplyFunction)
    assert value.function(2) == head(2)
    assert value[1, 2] == head(A[1, 2])
    assert value.free_symbols == {A}
    f = sx.Function('f', nargs=1, real=True)
    assert A.applyfunc(f)[1, 2] == f(A[1, 2])


def test_invalid_functions_parent_and_arities_fail_explicitly():
    z, u = sx.symbols('z u')
    A = sx.MatrixSymbol('A', 2, 3)
    for function in [sx.Lambda((z, u), z+u), sx.Integer(1)]:
        with pytest.raises(ValueError, match='Lambda|one|1'):
            A.applyfunc(function)
    with pytest.raises(ValueError, match='1'):
        A.applyfunc(sx.Function('f', nargs=2))
    with pytest.raises(NotImplementedError, match='destructur|scalar'):
        A.applyfunc(sx.Lambda(((z, u),), z+u))
    with pytest.raises(TypeError):
        sx.ElementwiseApplyFunction(sx.Lambda(z, z*z), 1)
    with pytest.raises(NotImplementedError, match='string'):
        A.applyfunc('lambda z:z*z')
    called = []
    def callback(z):
        called.append(z)
        return z*z
    with pytest.raises(TypeError, match='Lambda|symbolic'):
        A.applyfunc(callback)
    assert called == []


def test_numpy_source_and_existing_tensor_codegen_lower_applications_to_entries():
    z, t = sx.symbols('z t', real=True)
    A = sx.MatrixSymbol('A', 2, 3)
    value = (t*A).applyfunc(sx.Lambda(z, z*z+sx.sin(t)))
    entries = value.as_explicit()
    source = 'numpy.array([' + ','.join('['+','.join(sx.numpycode(entries[r, c]) for c in range(3))+']' for r in range(2))+'])'
    data = np.arange(6.).reshape(2, 3)
    np.testing.assert_allclose(eval(source, {'numpy': np, 't': .25, 'A': data}), (.25*data)**2+np.sin(.25))
    for generate in [sx.extra.torchcode_function, sx.extra.jaxcode_function]:
        source = generate('applied_matrix', (A, t), value)
        ast.parse(source)
        assert 'ElementwiseApplyFunction' not in source and 'Lambda' not in source
        assert 'symox' not in source and 'sympy' not in source


def test_native_runtime_and_shared_codec_without_sympy():
    from symox.matrices.expressions.applyfunc import ElementwiseApplyFunction
    assert ElementwiseApplyFunction is sx.ElementwiseApplyFunction
    code = '''
import sys,pickle
sys.modules['sympy']=None
import symox as s
z,t=s.symbols('z t')
A=s.MatrixSymbol('A',2,3)
f=A.applyfunc(s.Lambda(z,z+t))
assert f[1,2] == A[1,2]+t
assert f.diff(t).as_explicit() == s.ImmutableMatrix.ones(2,3)
for i in range(45):
    f=s.ElementwiseApplyFunction(s.IdentityFunction(),f)
for protocol in range(pickle.HIGHEST_PROTOCOL+1):
    payload=pickle.dumps(f,protocol=protocol)
    assert len(payload)<50000
    r=pickle.loads(payload)
    assert r==f
assert f.doit().expr==A
assert s.sympify(s.srepr(A.applyfunc(s.Lambda(z,z*z)))) == A.applyfunc(s.Lambda(z,z*z))
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=15)
