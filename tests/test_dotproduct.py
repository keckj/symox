"""Native scalar bilinear contractions of typed row and column vectors."""
import ast
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx


def orient(values, row, module):
    return module.ImmutableMatrix(1, len(values), values) if row else module.ImmutableMatrix(len(values), 1, values)


def scalar_equal(actual, expected):
    expected = sx.sympify(expected)
    assert sx.simplify(sx.expand(actual-expected)) == 0


@pytest.mark.parametrize('arow,brow', [(True, True), (True, False), (False, True), (False, False)])
def test_class_orientation_scope_interchange_and_protocol(arow, brow):
    a = sx.MatrixSymbol('a', *( (1, 3) if arow else (3, 1)))
    b = sx.MatrixSymbol('b', *( (1, 3) if brow else (3, 1)))
    value = sx.DotProduct(a, b)
    assert isinstance(value, sx.DotProduct) and isinstance(value, sx.Expr)
    assert not isinstance(value, sx.MatrixExpr)
    assert value.class_name == 'DotProduct'
    assert value.args == (a, b)
    assert value.func is sx.DotProduct
    assert value.func(*value.args) == value
    assert value.free_symbols == {a, b}
    assert value.has(a) and value.has(b)
    assert value.is_number is False and value.is_Number is False
    assert value.is_Function is False and value.is_commutative is None
    assert value.count_ops() == 0
    assert value != sx.DotProduct(b, a)  # held constructor preserves order
    reference = sp.DotProduct(a.to_sympy(), b.to_sympy())
    assert value.to_sympy() == reference
    assert sx.sympify(reference) == value
    assert sx.sympify(sx.srepr(value)) == value
    assert sx.default_sort_key(value) == sp.default_sort_key(reference)
    assert sx.latex(value) == sp.latex(reference)
    result = sum((a[0, i] if arow else a[i, 0])*(b[0, i] if brow else b[i, 0]) for i in range(3))
    assert value.doit() == result
    assert value.doit(deep=False, expand=False) == result
    assert value.doit(expand=True) == result


@pytest.mark.parametrize('arow,brow', [(True, True), (True, False), (False, True), (False, False)])
@pytest.mark.parametrize('complex_values', [False, True])
def test_exact_and_floating_bilinear_values_without_implicit_conjugation(arow, brow, complex_values):
    if complex_values:
        aa = [1+2*sx.I, -3+sx.I, 2*sx.I]
        bb = [2-sx.I, sx.I, -1+sx.I]
        av = np.asarray([1+2j, -3+1j, 2j])
        bv = np.asarray([2-1j, 1j, -1+1j])
    else:
        aa, bb = [sx.Rational(1, 3), -2, sx.Rational(5, 7)], [4, sx.Rational(2, 5), -3]
        av, bv = np.asarray([1/3, -2., 5/7]), np.asarray([4., 2/5, -3.])
    value = sx.DotProduct(orient(aa, arow, sx), orient(bb, brow, sx))
    expected = sum(a*b for a, b in zip(aa, bb))
    scalar_equal(value.doit(), expected)
    assert value.is_number is False and value.free_symbols == set()
    np.testing.assert_allclose(complex(value.evalf()), np.dot(av, bv), rtol=1e-13)
    if complex_values:
        assert np.dot(av, bv) != np.vdot(av, bv)
        hermitian = sx.DotProduct(value.args[0].conjugate(), value.args[1])
        np.testing.assert_allclose(complex(hermitian.doit().evalf()), np.vdot(av, bv), rtol=1e-13)
    # Binary floating values have a separate independent NumPy oracle.
    afloat = [sx.Float(float(v)) for v in [.125, -.5, 1.25]]
    bfloat = [sx.Float(float(v)) for v in [1.5, .75, -.25]]
    floated = sx.DotProduct(orient(afloat, arow, sx), orient(bfloat, brow, sx))
    np.testing.assert_allclose(float(floated.doit()), np.dot([.125, -.5, 1.25], [1.5, .75, -.25]))


@pytest.mark.parametrize('arow,brow', [(True, True), (True, False), (False, True), (False, False)])
@pytest.mark.parametrize('length', [0, 1])
def test_empty_and_scalar_vectors(length, arow, brow):
    a = orient(list(range(length)), arow, sx)
    b = orient([3]*length, brow, sx)
    result = sx.DotProduct(a, b)
    assert result.args[0].shape == a.shape and result.args[1].shape == b.shape
    assert result.doit() == 0
    assert result.to_sympy().doit() == 0
    assert sx.sympify(sx.srepr(result)) == result


def test_symbolic_length_sum_binding_substitution_and_scope():
    n = sx.Symbol('n', integer=True, positive=True)
    a, b = sx.MatrixSymbol('a', 1, n), sx.MatrixSymbol('b', n, 1)
    held = sx.DotProduct(a, b)
    assert held.free_symbols == {a, b}  # matrix dimensions are metadata
    assert held.has(n)
    total = held.doit()
    assert isinstance(total, sx.Sum)
    assert total.free_symbols == {a, b, n}
    index, lower, upper = total.args[1].args
    assert index.is_Dummy and lower == 0 and upper == n-1
    assert total.args[0] == a[0, index]*b[index, 0]
    concrete = held.xreplace({n: 3})
    assert concrete.args[0].shape == (1, 3) and concrete.args[1].shape == (3, 1)
    assert concrete.doit() == sum(concrete.args[0][0, i]*concrete.args[1][i, 0] for i in range(3))
    assert held.subs(n, 3) == concrete
    assert held.diff(n) == 0


def test_native_scalar_chain_rule_and_complex_real_parameter():
    t = sx.Symbol('t', real=True)
    a = sx.MatrixExpr(sx.ImmutableMatrix([[1+t, 2*sx.I*t, t*t]]))
    b = sx.MatrixExpr(sx.ImmutableMatrix([[t+sx.I], [3-t], [1+2*t]]))
    value = sx.DotProduct(a, b)
    expected = (1+t)*(t+sx.I)+(2*sx.I*t)*(3-t)+t*t*(1+2*t)
    for order in [1, 2, 3]:
        scalar_equal(value.diff(t, order).doit(), expected.diff(t, order))
    scalar_equal(value.doit(), expected)
    scalar_equal(sx.conjugate(value).doit(), sx.conjugate(expected))
    assert value.diff(sx.Symbol('missing')) == 0


def test_substitutions_are_exact_simultaneous_and_validate_nested_lengths():
    n = sx.Symbol('n', integer=True, positive=True)
    a, b = sx.MatrixSymbol('a', 1, n), sx.MatrixSymbol('b', n, 1)
    c, d = sx.MatrixSymbol('c', 1, n), sx.MatrixSymbol('d', n, 1)
    value = sx.DotProduct(a, b)
    replaced = value.xreplace({a: c, b: d, n: 2})
    assert replaced.args == (c, d)  # replacement matrices are not revisited
    assert value.xreplace({value: sx.Integer(7), n: 2}) == 7
    with pytest.raises(TypeError, match='length'):
        value.xreplace({a: sx.MatrixSymbol('bad', 1, 2)})
    with pytest.raises(TypeError):
        value.xreplace({a: sx.MatrixSymbol('bad', 2, 2)})
    outer = sx.DotProduct(sx.MatrixExpr(sx.ImmutableMatrix([[value]])), sx.MatrixExpr(sx.ImmutableMatrix([[1]])))
    with pytest.raises(TypeError, match='length'):
        outer.xreplace({a: sx.MatrixSymbol('bad', 1, 2)})


def test_expansion_preserves_held_contraction_and_applies_hints_to_parents():
    t, u = sx.symbols('t u')
    a, b = sx.MatrixSymbol('a', 1, 3), sx.MatrixSymbol('b', 3, 1)
    value = sx.DotProduct((t+u)**2*a, b)
    expanded = value.expand()
    assert isinstance(expanded, sx.DotProduct)
    assert expanded.args[0] == (t*t+2*t*u+u*u)*a
    assert value.expand(deep=False) == value
    assert value.expand(multinomial=False) == value
    scalar_equal(expanded.doit(), value.doit())
    assert sx.DotProduct(2*a, 3*b).count_ops() == 2
    scalar_equal(sx.block_collapse(value).doit(), value.doit())


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL+1))
def test_mixed_flat_codec_with_bindings_precision_and_matrix_scalar_parents(protocol):
    t = sx.Symbol('t', real=True)
    z = sx.Symbol('z')
    f = sx.Function('f', nargs=1, positive=True)
    bound = sx.Subs(f(z), z, t)
    value = sx.DotProduct(sx.MatrixExpr(sx.ImmutableMatrix([[bound, sx.Float('0.125', 80)]])),
                          sx.MatrixExpr(sx.ImmutableMatrix([[t], [sx.sin(t)]])))
    mixed = sx.Trace(sx.MatrixExpr(sx.ImmutableMatrix([[value]])))+value**2
    for expression in [value, mixed]:
        restored = pickle.loads(pickle.dumps(expression, protocol=protocol))
        assert restored == expression and restored.args == expression.args
        assert restored.free_symbols == expression.free_symbols
        assert sx.sympify(sx.srepr(restored)) == restored
    assert value.free_symbols == {t} and z not in value.free_symbols


@pytest.mark.parametrize('nonfinite', [sx.oo, -sx.oo, sx.zoo, sx.nan])
def test_nonfinite_zero_terms_are_not_skipped(nonfinite):
    a, b = sx.ImmutableMatrix([[0]]), sx.ImmutableMatrix([[nonfinite]])
    assert sx.DotProduct(a, b).doit() == sx.nan
    assert sx.DotProduct(b, a).doit() == sx.nan
    nested = sx.ImmutableMatrix([[sx.Symbol('x')*nonfinite]])
    assert sx.DotProduct(a, nested).doit() == sx.nan


def test_invalid_arguments_lengths_and_options():
    row = sx.MatrixSymbol('a', 1, 3)
    for a, b in [(1, row), (row, 1), (sx.ZeroMatrix(0, 0), row),
                 (sx.MatrixSymbol('bad', 2, 2), row), (row, sx.MatrixSymbol('bad', 4, 1))]:
        with pytest.raises(TypeError):
            sx.DotProduct(a, b)
    n, m = sx.symbols('n m', integer=True, positive=True)
    with pytest.raises(TypeError, match='length'):
        sx.DotProduct(sx.MatrixSymbol('a', 1, n), sx.MatrixSymbol('b', m, 1))
    value = sx.DotProduct(row, row)
    with pytest.raises(NotImplementedError):
        value.doit(unknown=True)
    with pytest.raises(TypeError):
        value.doit(expand='yes')
    undefined = sx.Function('DotProduct')(sx.Symbol('x'), sx.Symbol('y'))
    assert not isinstance(undefined, sx.DotProduct)
    assert value.diff(row)==sx.Array([[2*row[0,j] for j in range(3)]])


def test_zero_dimensional_matrix_inputs_remain_explicitly_unsupported_in_codegen():
    a, b = sx.MatrixSymbol('a', 1, 0), sx.MatrixSymbol('b', 0, 1)
    value = sx.DotProduct(a, b)
    assert value.doit() == 0
    for generate in [sx.torchcode_function, sx.jaxcode_function]:
        with pytest.raises(NotImplementedError, match='nonempty'):
            generate('empty_contract', (a, b), value)


def test_finite_numpy_and_tensor_codegen_sources_lower_typed_parents():
    a, b = sx.MatrixSymbol('a', 1, 3), sx.MatrixSymbol('b', 3, 1)
    t = sx.Symbol('t', real=True)
    value = sx.DotProduct(t*a, b)
    source = sx.numpycode(value)
    av, bv = np.asarray([[1+2j, 3, -1j]]), np.asarray([[2-1j], [1j], [4]])
    actual = eval(source, {'numpy': np, 'a': av, 'b': bv, 't': .3})
    np.testing.assert_allclose(actual, .3*np.dot(av.reshape(-1), bv.reshape(-1)))
    for generate in [sx.torchcode_function, sx.jaxcode_function]:
        source = generate('contract', (a, b, t), value)
        ast.parse(source)
        assert 'DotProduct' not in source and 'symox' not in source and 'sympy' not in source
    assert 'DotProduct' not in sx.ccode(value)
    assert 'DotProduct' not in sx.cxxcode(value)


def test_shared_mixed_dag_codec_and_expansion_without_sympy():
    from symox.matrices.expressions.dotproduct import DotProduct
    assert DotProduct is sx.DotProduct
    code = '''
import sys,pickle
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x',real=True)
y=x
for i in range(35):
    matrix=s.MatrixExpr(s.ImmutableMatrix([[y]]))
    y=s.DotProduct(matrix,matrix)
assert y.free_symbols=={x}
for protocol in range(pickle.HIGHEST_PROTOCOL+1):
    payload=pickle.dumps(y,protocol=protocol)
    assert len(payload)<60000
    restored=pickle.loads(payload)
    assert restored==y and restored.args==y.args
assert y.expand()==y
assert y.xreplace({x:s.Symbol('z',real=True)}).free_symbols=={s.Symbol('z',real=True)}
a=s.MatrixSymbol('a',1,3);b=s.MatrixSymbol('b',3,1)
assert s.DotProduct(a,b).doit()==sum(a[0,i]*b[i,0] for i in range(3))
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=20)
