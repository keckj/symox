"""Condensed exact SVD: defining certificates, not nonunique singular bases."""
import ast
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx


def equal(a, b):
    a = a.to_sympy() if hasattr(a, 'to_sympy') else a
    b = b.to_sympy() if hasattr(b, 'to_sympy') else b
    assert a.shape == b.shape
    for left, right in zip(a, b):
        delta = sp.simplify(left-right)
        if delta != 0:
            # SymPy leaves adjoint(1/sqrt(z*conjugate(z)+1)) opaque even
            # for complex z. Real/imaginary coordinates independently prove
            # positivity of re(z)**2+im(z)**2+1 and the exact certificate.
            delta = sp.simplify(sp.expand_complex(delta))
        assert delta == 0


def check(reference, immutable=True):
    a = sx.ImmutableMatrix(reference) if immutable else sx.Matrix(reference)
    u, sigma, v = a.singular_value_decomposition()
    rank = reference.rank()
    assert u.shape == (a.rows, rank) and sigma.shape == (rank, rank) and v.shape == (a.cols, rank)
    assert all(isinstance(factor, sx.ImmutableMatrix) == immutable for factor in (u, sigma, v))
    us, ss, vs = (factor.to_sympy() for factor in (u, sigma, v))
    equal(us.H*us, sp.eye(rank))
    equal(vs.H*vs, sp.eye(rank))
    equal(us*ss*vs.H, reference)
    assert all(ss[i, j] == 0 for i in range(rank) for j in range(rank) if i != j)
    # SymPy bases can differ by phases and by arbitrary orthonormal rotations
    # inside repeated eigenspaces. Match only the positive singular values.
    expected = [value for value in reference.singular_values() if value.is_zero is not True]
    assert len(expected) == rank
    for i in range(rank):
        value = ss[i, i]
        assert value.is_positive is not False
        matches = [j for j, target in enumerate(expected) if sp.simplify(value-target) == 0]
        assert matches, (value, expected)
        expected.pop(matches[0])
    return a, u, sigma, v


def projector(vector):
    vector = sp.Matrix(vector)
    return sp.eye(vector.rows)-vector*vector.H/(vector.H*vector)[0]


CASES = [
    sp.Matrix([[3]]), sp.Matrix([[-3]]), sp.Matrix([[1+sp.I]]),
    sp.Matrix([[1, 2], [2, 4]]), sp.Matrix([[1, 2], [3, 4]]),
    sp.Matrix([[1, 2, 3], [2, 4, 6]]), sp.Matrix([[1, 0], [0, 1], [1, 1]]),
    sp.Matrix([[1, sp.I]]), sp.Matrix([[1], [sp.I]]),
    sp.Matrix([[1, sp.I], [sp.I, 1]]),
    sp.Matrix([[1, sp.I, 2], [sp.I, 1, 0]]),
    sp.Matrix([[1, sp.I, 2], [sp.I, 1, 0]]).H,
    projector([1, 1, 1]), projector([1, sp.I, 1]),
    sp.Matrix([[1, 0, 0, 0, 2], [0, 0, 3, 0, 0], [0, 0, 0, 0, 0], [0, 2, 0, 0, 0]]),
]


@pytest.mark.parametrize('reference', CASES)
@pytest.mark.parametrize('immutable', [False, True])
def test_condensed_exact_svd_certificates_and_singular_values(reference, immutable):
    check(reference, immutable)


@pytest.mark.parametrize('shape', [(0, 0), (0, 3), (4, 0), (2, 3), (3, 2)])
@pytest.mark.parametrize('immutable', [False, True])
def test_zero_and_empty_shapes_preserve_condensed_dimensions(shape, immutable):
    constructor = sx.ImmutableMatrix if immutable else sx.Matrix
    a = constructor(sx.zeros(*shape))
    u, sigma, v = a.singular_value_decomposition()
    assert (u.shape, sigma.shape, v.shape) == ((shape[0], 0), (0, 0), (shape[1], 0))
    assert all(isinstance(factor, sx.ImmutableMatrix) == immutable for factor in (u, sigma, v))
    equal(u*sigma*v.H, sp.zeros(*shape))
    equal(u.H*u, sp.eye(0))
    equal(v.H*v, sp.eye(0))
    # SymPy1.14 raises on the empty-input case. These shapes are instead
    # certified independently by the matrix product and orthonormality laws.


def test_exact_parameter_domains_and_generic_complex_known_rank():
    p = sp.Symbol('p', positive=True)
    r = sp.Symbol('r', nonzero=True)
    z = sp.Symbol('z', complex=True)
    for a in [p*sp.Matrix([[1, 0, 2], [0, 3, 0]]), sp.diag(r, 0),
              sp.Matrix([[z, 1]]), sp.Matrix([[z], [1]])]:
        check(a)
    # An unknown sign does not require a branch when the parameter is nonzero.
    a, u, sigma, v = check(sp.diag(r, 0))
    equal(sigma, sp.Matrix([[sp.Abs(r)]]))
    for value in [-3, -sp.Rational(1, 2), 2]:
        equal((u*sigma*v.H).subs({sx.sympify(r): value}), sp.diag(value, 0))


@pytest.mark.parametrize('reference', [
    sp.Matrix([[sp.Symbol('a')]]),
    sp.diag(sp.Symbol('a', real=True), 1),
    sp.Matrix([[0, sp.Symbol('z', complex=True)], [0, 0]]),
])
def test_undecidable_rank_is_an_explicit_error(reference):
    with pytest.raises(NotImplementedError, match='rank|zero condition'):
        sx.ImmutableMatrix(reference).singular_value_decomposition()


@pytest.mark.parametrize('number', [sx.Float('0.125', 80), sx.Float(1.0)])
def test_float_inputs_are_not_presented_as_exact_svd(number):
    with pytest.raises(NotImplementedError, match='floating'):
        sx.Matrix([[number, 1], [0, 2]]).singular_value_decomposition()


@pytest.mark.parametrize('number', [sx.oo, -sx.oo, sx.zoo, sx.nan])
def test_nonfinite_inputs_are_rejected(number):
    with pytest.raises(ValueError, match='finite'):
        sx.ImmutableMatrix([[1, number]]).singular_value_decomposition()


def test_unresolved_quintic_spectrum_is_an_explicit_error():
    a = sx.eye(5)
    for i in range(4):
        a[i, i+1] = 1
    with pytest.raises(NotImplementedError, match='spectrum|eigenvalue'):
        a.singular_value_decomposition()


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL+1))
def test_results_pickle_and_sympy_interchange(protocol):
    a, *factors = check(sp.Matrix([[1, sp.I, 2], [sp.I, 1, 0]]))
    restored = pickle.loads(pickle.dumps(tuple(factors), protocol=protocol))
    for actual, expected in zip(restored, factors):
        assert actual == expected and isinstance(actual, sx.ImmutableMatrix)
        assert sx.ImmutableMatrix(actual.to_sympy()) == actual
    equal(restored[0]*restored[1]*restored[2].H, a)


def test_finite_scalar_codegen_and_independent_numpy_reconstruction():
    p = sx.Symbol('p', positive=True)
    matrix = p*sx.ImmutableMatrix([[1, sx.I, 2], [sx.I, 1, 0]])
    u, sigma, v = matrix.singular_value_decomposition()
    for generate in [sx.torchcode_function, sx.jaxcode_function]:
        source = generate('exact_svd', p, (u, sigma, v))
        ast.parse(source)
        assert 'symox' not in source and 'sympy' not in source
    sources = ['numpy.array(['+','.join(
        '['+','.join(sx.numpycode(factor[i, j]) for j in range(factor.cols))+']'
        for i in range(factor.rows))+'])' for factor in (u, sigma, v)]
    for point in [.1, 1., 2.5]:
        values = [np.asarray(eval(source, {'numpy': np, 'p': point}), dtype=np.complex128)
                  for source in sources]
        left, diagonal, right = values
        expected = point*np.asarray([[1, 1j, 2], [1j, 1, 0]], dtype=np.complex128)
        np.testing.assert_allclose(left.conj().T@left, np.eye(2), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(right.conj().T@right, np.eye(2), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(left@diagonal@right.conj().T, expected, rtol=1e-12, atol=1e-12)


def test_native_runtime_and_resource_rejections_precede_evaluation():
    subprocess.run([sys.executable, '-c', '''
import sys
sys.modules['sympy']=None
import symox as s
matrix=s.ImmutableMatrix([[1,s.I,2],[s.I,1,0]])
u,d,v=matrix.singular_value_decomposition()
assert (u*d*v.H-matrix).applyfunc(s.simplify)==s.zeros(2,3)
held=s.Pow(2,10**9,evaluate=False)
try:
    s.ImmutableMatrix([[held]]).singular_value_decomposition()
except NotImplementedError as error:
    assert 'resource limit' in str(error)
else:
    raise AssertionError('latent numeric power must be rejected before evaluation')
one=s.Integer(1)
y=s.Symbol('y')
raw=s._core._restore_arithmetic('Mul',s.Integer(2),[(y,one)],False)
for _ in range(80):
    raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
try:
    s.ImmutableMatrix([[raw]]).singular_value_decomposition()
except NotImplementedError as error:
    assert 'resource limit' in str(error)
else:
    raise AssertionError('latent mixed coefficients must be rejected before evaluation')
# Polynomial conversion stores a dense exponent vector per monomial; a width
# bound alone would miss its quadratic storage in the generator count.
wide=s.Add(*[s.Symbol('a'+str(i),positive=True) for i in range(2000)])
try:
    s.ImmutableMatrix([[wide]]).singular_value_decomposition()
except NotImplementedError as error:
    assert 'resource limit' in str(error)
else:
    raise AssertionError('polynomial scratch storage must be bounded before conversion')
'''], check=True, timeout=20)


@pytest.mark.parametrize('kind',['Symbol','Dummy'])
@pytest.mark.parametrize('flags',[{'finite':False},{'infinite':True},{'extended_real':True,'finite':False}])
def test_declared_nonfinite_symbols_rejected_before_spectral_work(kind,flags):
    z=getattr(sx,kind)('z',**flags)
    assert z.is_finite is False
    for matrix in [sx.ImmutableMatrix([[z,1]]),sx.ImmutableMatrix([[z],[1]])]:
        with pytest.raises(ValueError,match='finite'):
            matrix.singular_value_decomposition()


@pytest.mark.parametrize('flags',[{}, {'extended_real':True}])
def test_unknown_finiteness_does_not_reject_certifiable_symbolic_svd(flags):
    z=sx.Symbol('z',**flags)
    assert z.is_finite is None
    a=sx.ImmutableMatrix([[z,1]])
    u,d,v=a.singular_value_decomposition()
    assert (u.shape,d.shape,v.shape)==((1,1),(1,1),(2,1))
    points=[1.,-3.,0.] if flags else [1+2j,-3.,0.]
    for point in points:
        values=[np.asarray(f.subs({z:point}).to_sympy().evalf(),dtype=np.complex128) for f in [u,d,v]]
        left,diag,right=values
        np.testing.assert_allclose(left@diag@right.conj().T,np.asarray([[point,1]]),rtol=1e-13,atol=1e-13)
