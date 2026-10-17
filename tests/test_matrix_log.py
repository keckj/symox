"""Exact primary matrix logarithms and independent reconstruction certificates."""
import ast
import os
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx


def equal(actual, expected):
    a = actual.to_sympy() if hasattr(actual, 'to_sympy') else actual
    b = expected.to_sympy() if hasattr(expected, 'to_sympy') else expected
    assert a.shape == b.shape
    assert all(sp.simplify(sp.expand(x-y)) == 0 for x, y in zip(a, b))


def jordan_log_series(matrix, value):
    """Finite log(lambda*I+N) series, with N**n=0; no spectral solver."""
    n = matrix.rows
    identity = sx.eye(n)
    nil = matrix-value*identity
    power = identity
    result = sx.log(value)*identity
    assert nil**n == sx.zeros(n)
    for degree in range(1, n):
        power = power*nil
        result += sx.Rational((-1)**(degree+1), degree)*power/value**degree
    return result


@pytest.mark.parametrize('immutable', [False, True])
@pytest.mark.parametrize('value', [1, 2, -2, sp.I, -sp.I, 2+sp.I])
@pytest.mark.parametrize('size', [1, 2, 3])
def test_jordan_logarithm_matches_finite_series_sympy_and_reconstructs(value, size, immutable):
    matrix_class = sx.ImmutableMatrix if immutable else sx.Matrix
    value = sx.sympify(value)
    nil = sx.zeros(size)
    for i in range(size-1):
        nil[i, i+1] = i+2
    matrix = matrix_class(value*sx.eye(size)+nil)
    result = matrix.log()
    equal(result, jordan_log_series(matrix, value))
    equal(result, matrix.to_sympy().log())
    equal(result.exp(), matrix)
    equal(matrix*result, result*matrix)
    assert isinstance(result, sx.ImmutableMatrix) is immutable
    assert pickle.loads(pickle.dumps(result)) == result
    # SymPy may evaluate scalar logs or rational complex coefficients while
    # crossing the boundary; certify the mathematical roundtrip instead.
    equal(sx.ImmutableMatrix(result.to_sympy()), result)


@pytest.mark.parametrize('flags', [{}, {'nonzero': True}, {'positive': True}, {'negative': True}, {'real': True}])
def test_generic_parameter_formula_uses_scalar_branch_without_assuming_global_identities(flags):
    z = sx.Symbol('z', **flags)
    matrix = sx.ImmutableMatrix([[z, 2, 0], [0, z, 3], [0, 0, z]])
    result = matrix.log()
    equal(result, jordan_log_series(matrix, z))
    # Generic formulas are valid wherever their scalar denominators are nonzero.
    # Recompute the matrix logarithm when specializing a generic parameter to 0.
    points = [-2, -1] if flags.get('negative') else [1, 2]
    if flags == {}:
        points += [1+2*sx.I]
    for point in points:
        equal(result.subs(z, point).exp(), matrix.subs(z, point))
    if flags.get('positive'):
        equal(result.exp(), matrix)


@pytest.mark.parametrize('data', [
    [[2, 1], [1, 2]],
    [[0, -1], [1, 0]],
    [[sp.Rational(3, 4), sp.Rational(5, 4)], [sp.Rational(5, 4), sp.Rational(3, 4)]],
    [[1, 0], [0, -1]],
])
def test_nontriangular_and_negative_spectra_match_sympy_without_real_projection(data):
    matrix = sx.ImmutableMatrix(data)
    result = matrix.log()
    equal(result, matrix.to_sympy().log())
    # Independent exact reconstruction also covers complex discriminants
    # whose principal square roots the native exp solver can leave opaque.
    equal(result.to_sympy().exp(), matrix)
    reconstructed = result.exp().to_sympy().evalf(70)
    target = np.asarray(matrix.to_sympy(), dtype=np.complex128)
    np.testing.assert_allclose(np.asarray(reconstructed, dtype=np.complex128), target, rtol=1e-13, atol=1e-13)
    equal(matrix*result, result*matrix)


def test_similarity_complex_coefficients_and_nilpotent_series_certificate():
    p = sx.ImmutableMatrix([[1, sx.I, 1], [0, 1, 2-sx.I], [0, 0, 1]])
    j = sx.ImmutableMatrix([[2, 2, 0], [0, 2, 3], [0, 0, 2]])
    matrix = p*j*p.inv()
    expected = p*jordan_log_series(j, sx.Integer(2))*p.inv()
    result = matrix.log()
    equal(result, expected)
    equal(result.exp(), matrix)
    equal(matrix*result, result*matrix)


def test_commuting_family_derivative_is_inverse_with_assumptions_explicit():
    t = sx.Symbol('t', positive=True)
    nil = sx.ImmutableMatrix([[0, 2, 0], [0, 0, 3], [0, 0, 0]])
    matrix = t*sx.eye(3)+nil
    result = matrix.log()
    equal(result.diff(t), matrix.inv())
    equal(result.exp(), matrix)
    # This derivative certificate uses A'=I; A^-1*A' is not a general rule
    # for the logarithm of a family whose matrices do not commute.


def test_scalar_negative_imaginary_unit_has_exact_principal_branch():
    assert sx.log(-sx.I) == -sx.I*sx.pi/2
    assert sx.exp(sx.log(-sx.I)) == -sx.I
    z = sx.Symbol('z')
    assert str(sx.log(-z).func) == 'log'
    assert sx.log(-z).args == (-z,)


def test_branch_cut_does_not_enable_invalid_log_exp_or_product_rewrites():
    matrix = 2*sx.pi*sx.I*sx.eye(2)
    assert matrix.exp().log() == sx.zeros(2)
    assert matrix.exp().log() != matrix
    minus = -sx.eye(2)
    assert (minus*minus).log() == sx.zeros(2)
    equal(minus.log(), sx.pi*sx.I*sx.eye(2))
    assert (minus*minus).log() != 2*minus.log()


@pytest.mark.parametrize('immutable', [False, True])
def test_simplify_options_are_explicit_and_callable_postprocesses_once(immutable):
    matrix_class = sx.ImmutableMatrix if immutable else sx.Matrix
    matrix = matrix_class([[2, 1], [0, 2]])
    default = matrix.log()
    for option in [False, None, sx.cancel, sp.cancel]:
        equal(matrix.log(simplify=option), default)
    calls = []
    def postprocess(value):
        calls.append(value)
        assert isinstance(value, sx.ImmutableMatrix) is immutable
        return value.applyfunc(sx.expand)
    result = matrix.log(postprocess)
    assert len(calls) == 1
    equal(result, default)
    assert isinstance(result, sx.ImmutableMatrix) is immutable
    for option in [True, 2, 'cancel', object()]:
        with pytest.raises(TypeError, match='simplify'):
            matrix.log(simplify=option)
    with pytest.raises(TypeError):
        matrix.log(unknown=True)
    def raises(_):
        raise RuntimeError('caller callback failure')
    with pytest.raises(RuntimeError, match='caller callback failure'):
        matrix.log(raises)


def test_explicit_simplifier_callback_is_not_confused_with_a_rebound_public_cancel(monkeypatch):
    matrix = sx.ImmutableMatrix([[1, 2], [0, 1]])
    calls = []
    def custom(value):
        calls.append(value)
        return value
    monkeypatch.setattr(sx, 'cancel', custom)
    result = matrix.log(simplify=sx.cancel)
    assert calls == [result]
    assert result == sx.ImmutableMatrix([[0, 2], [0, 0]])


@pytest.mark.parametrize('shape', [(0, 0), (0, 3), (3, 0), (2, 3)])
@pytest.mark.parametrize('immutable', [False, True])
def test_empty_shapes_and_square_validation(shape, immutable):
    cls = sx.ImmutableMatrix if immutable else sx.Matrix
    matrix = cls(shape[0], shape[1], [0]*(shape[0]*shape[1]))
    if shape == (0, 0):
        for option in [False, None, sx.cancel]:
            result = matrix.log(simplify=option)
            assert result.shape == shape
            assert isinstance(result, sx.ImmutableMatrix) is immutable
        assert matrix.log().shape == shape
    else:
        with pytest.raises(sx.NonSquareMatrixError):
            matrix.log()


@pytest.mark.parametrize('data', [[[0]], [[0, 1], [0, 0]], [[1, 0], [0, 0]], [[1, 2], [2, 4]]])
def test_zero_eigenvalue_has_no_finite_matrix_logarithm(data):
    with pytest.raises(ValueError, match='pole|nonfinite|analytic'):
        sx.ImmutableMatrix(data).log()


@pytest.mark.parametrize('value', [sx.oo, -sx.oo, sx.zoo, sx.nan])
def test_nonfinite_entries_rejected(value):
    with pytest.raises(ValueError, match='finite'):
        sx.ImmutableMatrix([[value]]).log()


@pytest.mark.parametrize('kind', ['Symbol', 'Dummy'])
@pytest.mark.parametrize('flags', [{'finite': False}, {'infinite': True}, {'extended_real': True, 'finite': False}])
def test_explicit_nonfinite_symbol_facts_rejected(kind, flags):
    value = getattr(sx, kind)('z', **flags)
    with pytest.raises(ValueError, match='finite'):
        sx.ImmutableMatrix([[value]]).log()


def test_float_diagonal_keeps_scalar_precision_and_general_float_spectrum_is_unsupported():
    value = sx.Float('0.125', 80)
    matrix = sx.ImmutableMatrix([[value, 0], [0, -value]])
    result = matrix.log()
    assert result[0, 0]._prec >= value._prec
    np.testing.assert_allclose(complex(result[0, 0]), np.log(.125), rtol=1e-14)
    np.testing.assert_allclose(complex(result[1, 1]), np.log(.125)+1j*np.pi, rtol=1e-14)
    with pytest.raises(NotImplementedError):
        sx.ImmutableMatrix([[.125, .25], [-.5, .75]]).log()


def test_unresolved_exact_spectrum_is_not_approximated():
    matrix = sx.ImmutableMatrix([[0, 0, 0, 0, -1], [1, 0, 0, 0, 1], [0, 1, 0, 0, 0],
                                 [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]])
    with pytest.raises(NotImplementedError, match='exact|spectrum'):
        matrix.log()


@pytest.mark.parametrize('immutable', [False, True])
def test_sparse_delegation_is_bounded_and_preserves_mutability(immutable, monkeypatch):
    cls = sx.ImmutableSparseMatrix if immutable else sx.SparseMatrix
    value = cls([[1, 2], [0, 1]])
    result = value.log()
    assert isinstance(result, cls)
    assert result.todok() == {(0, 1): sx.Integer(2)}
    assert pickle.loads(pickle.dumps(result)) == result
    equal(result.as_dense().exp(), value)
    for shape in [(0, 0), (0, 3), (3, 0)]:
        if shape == (0, 0):
            assert cls(*shape, {}).log().shape == shape
        else:
            with pytest.raises(sx.NonSquareMatrixError):
                cls(*shape, {}).log()
    huge = cls(10**6, 10**6, {})
    monkeypatch.setattr(type(huge), 'as_dense', lambda _: pytest.fail('densified before guard'))
    with pytest.raises(NotImplementedError, match='resource limit.*densification'):
        huge.log()


def test_finite_scalar_printers_produce_independent_numpy_source():
    t = sx.Symbol('t', positive=True)
    matrix = sx.ImmutableMatrix([[t, 2, 0], [0, t, sx.I], [0, 0, t]])
    result = matrix.log()
    sources = [sx.numpycode(result[i, j]) for i in range(3) for j in range(3)]
    for point in [.1, 1., 2.5]:
        actual = np.asarray([eval(source, {'numpy': np, 't': point}) for source in sources]).reshape(3, 3)
        nil = np.asarray([[0, 2, 0], [0, 0, 1j], [0, 0, 0]])
        expected = np.log(point)*np.eye(3)+nil/point-(nil@nil)/(2*point**2)
        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
    for generate in [sx.torchcode_function, sx.jaxcode_function]:
        source = generate('exact_log', t, result)
        ast.parse(source)
        assert 'symox' not in source and 'sympy' not in source


def test_native_runtime_and_preflights_reject_latent_costs_before_work():
    subprocess.run([sys.executable, '-c', '''
import sys
sys.modules['sympy']=None
import symox as s
m=s.ImmutableMatrix([[1,2,0],[0,1,3],[0,0,1]])
l=m.log()
assert l==s.ImmutableMatrix([[0,2,-3],[0,0,3],[0,0,0]])
assert l.exp()==m
assert m.log(simplify=s.cancel)==l

def limited(value):
    try:
        value.log()
    except NotImplementedError as error:
        assert 'resource limit' in str(error)
    else:
        raise AssertionError('matrix log did not enforce preflight')
held=s.Pow(2,10**9,evaluate=False)
limited(s.ImmutableMatrix([[held]]))
trace=s.Trace(s.MatrixExpr(s.ImmutableMatrix([[held]])))
limited(s.ImmutableMatrix([[trace]]))
one=s.Integer(1)
y=s.Symbol('y')
raw=s._core._restore_arithmetic('Mul',s.Integer(2),[(y,one)],False)
for _ in range(80):
    raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
limited(s.ImmutableMatrix([[raw]]))
entries=[s.Add(*(s.Symbol('a%d_%d'%(i,j)) for j in range(6))) for i in range(36)]
limited(s.ImmutableMatrix(6,6,entries))
'''], check=True, timeout=20, env=os.environ.copy())
