"""Native FunctionMatrix and its exact scalar Lambda entry definition."""
import ast
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sm
import symox as sx


def context(rows=2, cols=3):
    i, j, t = sx.symbols('i j t', integer=True)
    u, v, z = sm.symbols('i j t', integer=True)
    native = sx.FunctionMatrix(rows, cols, sx.Lambda((i, j), i+10*j+t*t))
    reference = sm.FunctionMatrix(rows, cols, sm.Lambda((u, v), u+10*v+z*z))
    return native, reference, (i, j, t), (u, v, z)


@pytest.mark.parametrize('rows,cols', [(1, 1), (2, 3), (3, 2), (4, 4)])
def test_functionmatrix_construction_entries_and_sympy_interchange(rows, cols):
    value, reference, _, _ = context(rows, cols)
    assert isinstance(value, sx.FunctionMatrix)
    assert isinstance(value, sx.MatrixExpr)
    assert value.is_FunctionMatrix and value._kind_tag == 'FunctionMatrix'
    assert value.args == (rows, cols, value.lamda)
    assert value.func is sx.FunctionMatrix
    assert value.func(*value.args) == value
    assert value.to_sympy() == reference
    assert sx.sympify(reference) == value
    assert sx.MatrixExpr(reference) == value
    assert sx.sympify(sx.srepr(value)) == value
    assert value.as_explicit().to_sympy() == reference.as_explicit()
    for i in range(rows):
        for j in range(cols):
            assert value[i, j].to_sympy() == reference[i, j]


@pytest.mark.parametrize('rows,cols', [(0, 0), (0, 3), (3, 0), (1, 1), (2, 3)])
def test_empty_shapes_and_negative_indices_use_independent_dense_certificate(rows, cols):
    value, _, (_, _, t), _ = context(rows, cols)
    expected = sx.ImmutableMatrix(rows, cols, [i+10*j+t*t for i in range(rows) for j in range(cols)])
    assert value.shape == (rows, cols)
    assert value.as_explicit() == expected
    assert value.as_explicit().shape == (rows, cols)
    assert value.T.as_explicit() == expected.T
    assert value.H.as_explicit() == expected.H
    assert value.conjugate().as_explicit() == expected.conjugate()
    if rows and cols:
        for i in range(-rows, rows):
            for j in range(-cols, cols):
                assert value[i, j] == expected[i, j]
    for index in [(rows, 0), (0, cols), (-rows-1, 0), (0, -cols-1)]:
        with pytest.raises(IndexError):
            value[index]


def test_symbolic_shape_free_variables_and_substitution_follow_lambda_scope():
    n, m = sx.symbols('n m', integer=True, positive=True)
    i, j, t = sx.symbols('i j t', integer=True)
    value = sx.FunctionMatrix(n, m, sx.Lambda((i, j), i+10*j+t))
    assert value.shape == (n, m)
    assert value.free_symbols == {n, m, t}
    assert i not in value.free_symbols and j not in value.free_symbols
    assert value[-1, -1] == n-1+10*(m-1)+t
    assert value[i, j] == i+10*j+t
    for replace in [lambda v: v.subs({n: 2, m: 3, t: 4}),
                    lambda v: v.xreplace({n: 2, m: 3, t: 4})]:
        concrete = replace(value)
        assert concrete.shape == (2, 3)
        assert concrete[1, 2] == 25
        assert list(concrete.as_explicit()) == [4, 14, 24, 5, 15, 25]
    # SymPy 1.14: MatrixExpr.as_explicit rejects a symbolic shape with ValueError.
    with pytest.raises(ValueError, match='symbolic shape'):
        value.as_explicit()
    # Lambda application/substitution is simultaneous and permits capture,
    # matching SymPy rather than applying an alpha-renaming convention.
    captured = value.xreplace({t: i})
    assert captured[1, 2] == 22
    k = sx.Symbol('k', integer=True)
    renamed = value.xreplace({i: k})
    assert renamed != value
    assert renamed[1, 2] == value[1, 2]
    assert renamed.free_symbols == value.free_symbols


def test_free_matrix_symbols_inside_lambda_and_exact_scalar_derivatives():
    i, j = sx.symbols('i j', integer=True)
    t = sx.Symbol('t', real=True)
    A = sx.MatrixSymbol('A', 2, 3)
    value = sx.FunctionMatrix(2, 3, sx.Lambda((i, j), t*A[i, j]+i*j))
    assert value.free_symbols == {A, t}
    assert value[1, 2] == t*A[1, 2]+2
    derivative = value.diff(t)
    assert isinstance(derivative, sx.FunctionMatrix)
    assert derivative.shape == value.shape
    assert derivative.as_explicit() == A.as_explicit()
    assert value.diff(i) == sx.ZeroMatrix(2, 3)
    assert value.diff(j) == sx.ZeroMatrix(2, 3)
    assert value.diff(t, 2) == sx.ZeroMatrix(2, 3)
    assert value.T.diff(t).as_explicit() == A.T.as_explicit()
    assert value.H.diff(t).as_explicit() == A.H.as_explicit()
    # SymPy FunctionMatrix.diff currently raises AttributeError for this case;
    # certify the new native derivative through independent dense calculus.
    assert value.as_explicit().applyfunc(lambda e: e.diff(t)) == derivative.as_explicit()


def test_functionmatrix_scalar_expansion_and_trace_keep_native_nodes():
    i, j, x, y = sx.symbols('i j x y')
    value = sx.FunctionMatrix(2, 2, sx.Lambda((i, j), (x+y)**2+i+j))
    expanded = value.expand()
    assert isinstance(expanded, sx.FunctionMatrix)
    assert expanded[0, 0] == x*x+2*x*y+y*y
    assert value.expand(deep=False) == value
    assert value.expand(multinomial=False) == value
    assert sx.expand(value) == expanded
    assert value.trace() == 2*(x+y)**2+2
    assert value.T.as_explicit() == value.as_explicit().T
    assert (value+value).as_explicit() == value.as_explicit()+value.as_explicit()
    actual_product = (value*value).as_explicit()
    dense_product = value.as_explicit()*value.as_explicit()
    assert actual_product.shape == dense_product.shape
    assert all(sx.expand(a-b) == 0 for a,b in zip(actual_product,dense_product))


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL+1))
def test_functionmatrix_pickle_and_mixed_scalar_metadata(protocol):
    value, _, _, _ = context()
    expressions = [value, value.T, value[::-1, ::-1], sx.Trace(value.T*value),
                   sx.Determinant(value.T*value), value[1, 2]]
    for expression in expressions:
        restored = pickle.loads(pickle.dumps(expression, protocol=protocol))
        assert restored == expression
        assert restored.args == expression.args
        assert sx.sympify(sx.srepr(restored)) == restored


def test_undefined_and_known_native_function_heads_become_real_two_argument_lambdas():
    f = sx.Function('f', nargs=(1, 2), real=True)
    value = sx.FunctionMatrix(2, 3, f)
    assert value.lamda.func.__name__ == 'Lambda'
    assert value.lamda(1, 2) == f(1, 2)
    assert value[1, 2] == f(1, 2)
    assert value.free_symbols == set()
    assert pickle.loads(pickle.dumps(value)) == value
    assert sx.sympify(sx.srepr(value)) == value
    assert sx.FunctionMatrix(2, 3, sx.atan2)[1, 2] == sx.atan2(1, 2)
    choose = sx.FunctionMatrix(3, 3, sx.binomial)
    assert list(choose.as_explicit()) == [sx.binomial(i, j) for i in range(3) for j in range(3)]
    i, j = sx.symbols('i j', integer=True)
    diagonal = sx.FunctionMatrix(3, 3, sx.Lambda((i, j), sx.Piecewise((1, sx.Eq(i, j)), (0, True))))
    assert diagonal.as_explicit() == sx.ImmutableMatrix(3, 3, [int(i == j) for i in range(3) for j in range(3)])
    with pytest.raises(ValueError, match='two|2'):
        sx.FunctionMatrix(2, 3, sx.Function('g', nargs=1))


@pytest.mark.parametrize('rows,cols', [(-1, 2), (2, -1), (sx.Rational(1, 2), 2)])
def test_invalid_matrix_dimensions(rows, cols):
    i, j = sx.symbols('i j')
    with pytest.raises(ValueError):
        sx.FunctionMatrix(rows, cols, sx.Lambda((i, j), i+j))


def test_invalid_lambda_signatures_and_python_callables_fail_explicitly():
    i, j, k = sx.symbols('i j k')
    for lamda in [sx.Lambda(i, i+1), sx.Lambda((i, j, k), i+j), sx.Symbol('L'), sx.Integer(1)]:
        with pytest.raises(ValueError, match='Lambda|two|2'):
            sx.FunctionMatrix(2, 3, lamda)
    with pytest.raises(NotImplementedError, match='destructur|scalar'):
        sx.FunctionMatrix(2, 3, sx.Lambda(((i, j), k), i+j+k))
    # SymPy sympifies strings: a lambda string is accepted, an expression is not.
    assert sx.FunctionMatrix(2, 3, 'lambda i,j:i+j').as_explicit() == sx.Matrix([[0, 1, 2], [1, 2, 3]])
    with pytest.raises(ValueError, match='compatible'):
        sx.FunctionMatrix(2, 3, 'i + j')
    called = []
    def python_function(i, j):
        called.append((i, j))
        return i+j
    with pytest.raises(TypeError, match='Lambda|symbolic'):
        sx.FunctionMatrix(2, 3, python_function)
    assert called == []


def test_functionmatrix_codegen_uses_scalar_entries_and_existing_tensor_output_lowering():
    i, j = sx.symbols('i j', integer=True)
    t = sx.Symbol('t', real=True)
    value = sx.FunctionMatrix(2, 3, sx.Lambda((i, j), sx.sin(t)+i+2*j))
    entries = value.as_explicit()
    for scalar in entries:
        for print_expression in [sx.numpycode, sx.cxxcode, sx.extra.cudacode, sx.rust_code, sx.extra.torchcode, sx.extra.jaxcode]:
            source = print_expression(scalar)
            assert 'Lambda' not in source and 'FunctionMatrix' not in source
    source = 'numpy.array([' + ','.join('['+','.join(sx.numpycode(entries[r, c]) for c in range(3))+']' for r in range(2))+'])'
    expected = np.array([[np.sin(0.25)+r+2*c for c in range(3)] for r in range(2)])
    np.testing.assert_allclose(eval(source, {'numpy': np, 't': 0.25}), expected)
    for generate in [sx.extra.torchcode_function, sx.extra.jaxcode_function]:
        source = generate('function_matrix', (t,), value)
        ast.parse(source)
        assert 'FunctionMatrix' not in source and 'Lambda' not in source
        assert 'symox' not in source and 'sympy' not in source


def test_import_path_and_runtime_are_native_without_sympy():
    from symox.matrices.expressions.funcmatrix import FunctionMatrix
    assert FunctionMatrix is sx.FunctionMatrix
    code = '''
import sys,pickle
sys.modules['sympy']=None
import symox as s
i,j,t=s.symbols('i j t')
F=s.FunctionMatrix(2,3,s.Lambda((i,j),i+10*j+t*t))
assert F[-1,-1] == 21+t*t
assert F.diff(t)[1,2] == 2*t
assert F.as_explicit().shape == (2,3)
assert F.xreplace({t:4})[1,2] == 37
assert pickle.loads(pickle.dumps(F)) == F
assert s.sympify(s.srepr(F)) == F
assert s.FunctionMatrix(0,3,F.lamda).as_explicit().shape == (0,3)
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=15)
