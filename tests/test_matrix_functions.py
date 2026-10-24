"""Exact dense functional calculus, independent polynomial/differential certificates."""
import ast
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx


def equal(actual, expected):
    actual = actual.to_sympy() if hasattr(actual, 'to_sympy') else actual
    expected = expected.to_sympy() if hasattr(expected, 'to_sympy') else expected
    assert actual.shape == expected.shape
    assert all(sp.simplify(sp.expand(a-b)) == 0 for a, b in zip(actual, expected))


@pytest.mark.parametrize('data', [
    [[0]], [[1]], [[0, 1], [0, 0]], [[2, 3], [0, 2]],
    [[0, 1], [-1, 0]], [[1, 2], [3, 4]],
    [[1, 1, 0], [0, 1, 1], [0, 0, 1]],
    [[sp.I, 2], [0, -sp.I]],
    [[sp.Symbol('a', real=True), 1], [0, sp.Symbol('a', real=True)]],
    [[sp.Symbol('z'), 1], [0, sp.Symbol('z')]],
])
@pytest.mark.parametrize('immutable', [False, True])
def test_exponential_matches_sympy_exact_and_preserves_matrix_type(data, immutable):
    reference = sp.ImmutableMatrix(data) if immutable else sp.Matrix(data)
    value = sx.ImmutableMatrix(data) if immutable else sx.Matrix(data)
    result = value.exp()
    equal(result, reference.exp())
    assert isinstance(result, sx.ImmutableMatrix) is immutable
    assert result.shape == value.shape
    assert pickle.loads(pickle.dumps(result)) == result


@pytest.mark.parametrize('name', ['exp', 'sin', 'cos', 'log', 'sqrt', 'resolvent'])
def test_analytic_func_hermite_matches_oracle_with_repeated_eigenvalues(name):
    x, X = sx.Symbol('x'), sp.Symbol('x')
    a = sx.ImmutableMatrix([[2, 3, 0], [0, 2, 1], [0, 0, 3]])
    reference = a.to_sympy()
    if name == 'resolvent':
        f, F = 1/(x-5), 1/(X-5)
    else:
        f, F = getattr(sx, name)(x), getattr(sp, name)(X)
    actual = a.analytic_func(f, x)
    equal(actual, reference.analytic_func(F, X))
    equal(a*actual, actual*a)
    assert isinstance(actual, sx.ImmutableMatrix)
    if name == 'resolvent':
        equal(actual, (a-5*sx.eye(3)).inv())


def test_polynomial_path_does_not_require_solvable_eigenvalues():
    # Companion of an irreducible quintic: analytic nonpolynomial functions
    # remain unsupported, but polynomial substitution needs no eigensystem.
    a = sx.ImmutableMatrix([[0, 0, 0, 0, -1], [1, 0, 0, 0, 1], [0, 1, 0, 0, 0],
                            [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]])
    x, t = sx.symbols('x t')
    f = x**9-2*x**3+t*x+3
    # Independent Horner evaluation of the coefficient sequence.
    result = sx.zeros(5)
    coefficients = [1, 0, 0, 0, 0, 0, -2, 0, t, 3]
    for coefficient in coefficients:
        result = result*a+coefficient*sx.eye(5)
    equal(a.analytic_func(f, x), result)
    with pytest.raises(NotImplementedError, match='exact|spectrum'):
        a.exp()
    with pytest.raises(NotImplementedError, match='exact|spectrum'):
        a.analytic_func(sx.sin(x), x)


def test_similarity_equivariance_and_known_jordan_block_formula():
    p = sx.ImmutableMatrix([[1, 2, -1], [0, 1, 3], [0, 0, 1]])
    j = sx.ImmutableMatrix([[2, 1, 0], [0, 2, 1], [0, 0, 2]])
    a = p*j*p.inv()
    expected_j = sx.exp(2)*sx.ImmutableMatrix([[1, 1, sx.Rational(1, 2)], [0, 1, 1], [0, 0, 1]])
    equal(j.exp(), expected_j)
    equal(a.exp(), p*expected_j*p.inv())
    x = sx.Symbol('x')
    equal(a.analytic_func(sx.sin(x), x), p*j.analytic_func(sx.sin(x), x)*p.inv())


def test_exponential_differential_equation_initial_value_inverse_and_complex_generator():
    t = sx.Symbol('t', real=True)
    generator = sx.ImmutableMatrix([[1, sx.I, 0], [0, 1, 2], [0, 0, 1]])
    exponential = (t*generator).exp()
    equal(exponential.diff(t), generator*exponential)
    equal(exponential.subs(t, 0), sx.eye(3))
    equal(exponential*(-t*generator).exp(), sx.eye(3))
    equal(exponential.H, (t*generator.H).exp())


@pytest.mark.parametrize('shape', [(0, 0), (0, 3), (3, 0), (2, 3)])
def test_empty_functions_and_nonsquare_errors(shape):
    x = sx.Symbol('x')
    a = sx.ImmutableMatrix(shape[0], shape[1], [0]*(shape[0]*shape[1]))
    if shape[0] == shape[1]:
        assert a.exp().shape == a.analytic_func(sx.log(x), x).shape == shape
        assert isinstance(a.exp(), sx.ImmutableMatrix)
    else:
        with pytest.raises(sx.NonSquareMatrixError):
            a.exp()
        with pytest.raises(sx.NonSquareMatrixError):
            a.analytic_func(sx.exp(x), x)


@pytest.mark.parametrize('f_source', ['1/x', 'log(x)', '1/(x-1)', 'sqrt(x)'])
def test_explicit_poles_and_nonanalytic_derivatives_are_rejected(f_source):
    x = sx.Symbol('x')
    function = eval(f_source, dict(vars(sx), x=x))
    matrix = sx.ImmutableMatrix([[0, 1], [0, 0]]) if f_source != '1/(x-1)' else sx.eye(2)
    with pytest.raises(ValueError, match='analytic|pole|nonfinite'):
        matrix.analytic_func(function, x)


@pytest.mark.parametrize('value', [sx.oo, -sx.oo, sx.zoo, sx.nan])
def test_nonfinite_matrix_entries_and_functions_raise_value_error(value):
    x = sx.Symbol('x')
    matrix = sx.ImmutableMatrix([[value, 1], [0, 0]])
    with pytest.raises(ValueError, match='finite'):
        matrix.exp()
    with pytest.raises(ValueError, match='finite'):
        matrix.analytic_func(x*x, x)
    with pytest.raises(ValueError):
        sx.eye(2).analytic_func(value*x, x)


def test_parameter_validation_scope_and_scalar_arguments():
    x, y = sx.symbols('x y')
    a = sx.eye(2)
    for function, parameter in [(sx.exp(x), 1), (sx.Integer(2), x), (sx.exp(y), x)]:
        with pytest.raises(ValueError, match='symbol|parameter'):
            a.analytic_func(function, parameter)
    with pytest.raises(ValueError, match='matrix'):
        sx.ImmutableMatrix([[x]]).analytic_func(sx.exp(x), x)
    with pytest.raises(TypeError):
        a.analytic_func(sx.MatrixExpr(a), x)
    with pytest.raises(TypeError):
        a.analytic_func(sx.exp(x), sx.MatrixSymbol('A', 2, 2))
    # Free symbols govern collisions; a bound integration variable does not
    # turn the scalar function's parameter into a free matrix parameter.
    integral = sx.Integral(x, (x, 0, 1))
    equal(sx.ImmutableMatrix([[integral]]).analytic_func(sx.exp(x), x), sx.ImmutableMatrix([[sx.exp(sx.Rational(1, 2))]]))


def test_diagonal_float_evaluation_preserves_declared_precision_and_polynomial_float_path():
    number = sx.Float('0.125', 80)
    a = sx.ImmutableMatrix([[number, 0], [0, -number]])
    result = a.exp()
    assert result[0, 0]._prec >= number._prec and result[1, 1]._prec >= number._prec
    np.testing.assert_allclose([float(result[0, 0]), float(result[1, 1])], np.exp([.125, -.125]), rtol=1e-14)
    x = sx.Symbol('x')
    floating = sx.ImmutableMatrix([[.125, .25], [-.5, .75]])
    equal(floating.analytic_func(x*x+x, x), floating*floating+floating)
    # SymPy 1.14 goes through the Jordan form of the Float matrix
    assert str(floating.exp()) == (
        'Matrix([[1.04589638195854, 0.385445367798660], '
        '[-0.770890735597320, 2.00950980145519]])')


def test_finite_codegen_uses_existing_scalar_printers_without_matrix_runtime():
    t = sx.Symbol('t', real=True)
    generator = sx.ImmutableMatrix([[1, 1], [0, 1]])
    result = (t*generator).exp()
    source = 'numpy.array(['+','.join('['+','.join(sx.numpycode(result[i, j]) for j in range(2))+']' for i in range(2))+'])'
    for point in [-.3, 0., .8]:
        actual = eval(source, {'numpy': np, 't': point})
        np.testing.assert_allclose(actual, np.exp(point)*np.asarray([[1., point], [0., 1.]]), rtol=1e-13)
    for generate in [sx.extra.torchcode_function, sx.extra.jaxcode_function]:
        source = generate('matrix_exponential', t, result)
        ast.parse(source)
        assert 'symox' not in source and 'sympy' not in source


def test_native_runtime_without_sympy():
    subprocess.run([sys.executable, '-c', '''
import sys
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x')
a=s.ImmutableMatrix([[1,2,0],[0,1,3],[0,0,1]])
expected=s.exp(1)*s.ImmutableMatrix([[1,2,3],[0,1,3],[0,0,1]])
assert a.exp()==expected
assert a.analytic_func(x*x+2*x,x)==a*a+2*a
assert a.analytic_func(1/(x-3),x)==(a-3*s.eye(3)).inv()
assert s.zeros(0).exp().shape==(0,0)
'''], check=True, timeout=20)


def test_resource_limits_reject_before_evaluation_and_keep_easy_large_powers():
    # The subprocess deadline detects accidental evaluation of a billion-bit
    # held number without ever constructing that number in the parent process.
    subprocess.run([sys.executable, '-c', '''
import symox as s
x=s.Symbol('x')
def limited(call):
    try:
        call()
    except NotImplementedError as error:
        assert 'resource limit' in str(error)
    else:
        raise AssertionError('missing exact matrix function resource guard')

exponent=2**100
nil=s.ImmutableMatrix([[0,1],[0,0]])
assert nil.analytic_func(x**exponent,x)==s.zeros(2)
assert s.ImmutableMatrix(s.eye(2)).analytic_func(x**exponent,x)==s.eye(2)
assert s.ImmutableMatrix([[0,1],[1,0]]).analytic_func(x**exponent,x)==s.eye(2)
a=s.ImmutableMatrix([[2,1],[0,2]])
limited(lambda:a.analytic_func(x**(2**50),x))
held=s.Pow(2,10**9,evaluate=False)
raw=s._core._restore_arithmetic('Mul',s.Integer(1),[(s.Integer(2),s.Integer(10**9))],False)
limited(lambda:s.ImmutableMatrix([[raw]]).exp())
# Parent metadata is not present in the scalar .args tuple.
trace=s.Trace(s.MatrixExpr(s.ImmutableMatrix([[held]])))
limited(lambda:s.ImmutableMatrix([[trace]]).exp())
limited(lambda:s.ImmutableMatrix([[held]]).exp())
limited(lambda:a.analytic_func(s.Mul(x,held,evaluate=False),x))
# This power becomes numerical only during scalar spectral sampling.
limited(lambda:s.ImmutableMatrix([[2]]).analytic_func(s.sin(x**(10**9)),x))
limited(lambda:s.ImmutableMatrix([[10**9]]).analytic_func(s.gamma(x),x))
# Raw mixed DAGs retain an exponentially large numeric coefficient even
# though their leaves include a free symbolic parameter.
one=s.Integer(1)
y=s.Symbol('y')
mixed=s._core._restore_arithmetic('Mul',s.Integer(2),[(y,one)],False)
for _ in range(80):
    mixed=s._core._restore_arithmetic('Mul',one,[(mixed,one),(mixed,one)],False)
limited(lambda:s.ImmutableMatrix([[mixed]]).analytic_func(s.exp(x),x))
# Repeated shared products have no single large exponent to screen.
f=x
for _ in range(50):
    f=s.Mul(f,f,evaluate=False)
limited(lambda:a.analytic_func(f,x))
'''], check=True, timeout=20)


@pytest.mark.parametrize('kind',['Symbol','Dummy'])
@pytest.mark.parametrize('flags',[{'finite':False},{'infinite':True},{'extended_real':True,'finite':False}])
def test_nonfinite_declared_symbols_are_rejected_in_matrix_and_function(kind,flags):
    x=sx.Symbol('x')
    z=getattr(sx,kind)('z',**flags)
    assert z.is_finite is False
    a=sx.ImmutableMatrix([[z]])
    with pytest.raises(ValueError,match='finite'):
        a.exp()
    with pytest.raises(ValueError,match='finite'):
        a.analytic_func(sx.exp(x),x)
    # The concrete nonfinite parameter stays present even if future scalar
    # canonicalization can absorb finite summands next to an infinity.
    function=sx.Add(x,z,evaluate=False)
    with pytest.raises(ValueError,match='nonfinite'):
        sx.eye(1).analytic_func(function,x)


@pytest.mark.parametrize('flags',[{}, {'extended_real':True}])
def test_unknown_finiteness_keeps_symbolic_matrix_functions_available(flags):
    z=sx.Symbol('z',**flags)
    x=sx.Symbol('x')
    assert z.is_finite is None
    equal(sx.ImmutableMatrix([[z]]).exp(),sx.ImmutableMatrix([[sx.exp(z)]]))
    equal(sx.eye(1).analytic_func(x+z,x),sx.ImmutableMatrix([[1+z]]))
