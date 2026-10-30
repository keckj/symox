"""Native float64 array evaluation of lambdify (lambdarr): SymPy 1.14 values
and bit-for-bit differential checks against the established paths.

The established paths stay reachable with non-contiguous views of the same
data, which keep NumPy's ufunc evaluation of the printed text and the copying
evaluation of the compiled program; the native path must agree bit for bit."""
import numpy as np
import pytest

from symox import *
del test

x, y, z = symbols('x y z')
A = np.array([0.5, -2.0, 3.0, 0.25])


def bits(values):
    # bit patterns, with one NaN: the sign of a NaN depends on the operand
    # order that the compiled loops of NumPy and of the kernels use
    values = np.array(values, dtype=np.float64)
    values[np.isnan(values)] = np.nan
    return values.view(np.uint64).tolist()


def strided(a):
    # the same values in a non-contiguous view (established path)
    buffer = np.empty(a.shape + (2,), dtype=np.float64)
    view = buffer[..., 0]
    view[...] = a
    assert not view.flags.c_contiguous or a.size == 0
    return view


# SymPy 1.14 lambdify values at A, as float.hex.
SYMPY = [
    (x**2 + 3*x + 1, 'numpy',
     ['0x1.6000000000000p+1', '-0x1.0000000000000p+0', '0x1.3000000000000p+4', '0x1.d000000000000p+0']),
    (Piecewise((x**2, x < 1), (sin(x), True)), 'numpy',
     ['0x1.0000000000000p-2', '0x1.0000000000000p+2', '0x1.210386db6d55bp-3', '0x1.0000000000000p-4']),
    (exp(-x**2)*cos(x), 'numpy',
     ['0x1.5deebabf616cep-1', '-0x1.f383dc0a6fb92p-8', '-0x1.003815e93f37dp-13', '0x1.d206e984049a0p-1']),
    (sqrt(x**2 + 2*x + 5), 'numpy',
     ['0x1.4000000000000p+1', '0x1.1e3779b97f4a8p+1', '0x1.1e3779b97f4a8p+2', '0x1.2de32c6628741p+1']),
]


@pytest.mark.parametrize('expr,modules,expected', SYMPY)
def test_sympy_values_on_float64_arrays(expr, modules, expected):
    f = lambdify(x, expr, modules)
    result = f(A)
    assert type(result) is np.ndarray and result.dtype == np.float64 and result.shape == A.shape
    assert [float(v).hex() for v in result] == expected


def test_three_argument_polynomial_with_a_scalar():
    f = lambdify((x, y, z), x**2*y + 3*x*z - y*z**3 + 2, 'numpy')
    assert f(A, A[::-1].copy(), 0.5).tolist() == [2.78125, 10.625, -11.25, 2.34375]


KERNELS = [
    ((x,), sin(x)), ((x,), x**2 + 3*x + 1), ((x,), exp(-x**2)*cos(x)),
    ((x, y), sqrt(x**2 + y**2)), ((x, y, z), x**2*y + 3*x*z - y*z**3 + 2),
    ((x,), x**3 - x/3 + log(x + 3) + pi), ((x, y), x**Rational(3, 2) + y**-1 - E*x),
    ((x,), Integer(7)), ((x,), x), ((x, y), x*y**2 + 1),
]


def samples(shape, seed):
    rng = np.random.default_rng(seed)
    a = rng.uniform(-4, 4, size=shape)
    flat = a.reshape(-1)
    flat[:6] = [0.0, -0.0, np.inf, -np.inf, np.nan, 5e-324][:flat.size][:6]
    return a


@pytest.mark.parametrize('arguments,expr', KERNELS)
@pytest.mark.parametrize('shape', [(10,), (1000,), (3, 7), (0,)])
def test_numpy_kernels_match_the_ufunc_evaluation_bitwise(arguments, expr, shape):
    f = lambdify(arguments, expr, 'numpy')
    arrays = [samples(shape, i) for i in range(len(arguments))]
    result = f(*arrays)
    reference = f(*[strided(a) for a in arrays])
    assert result.shape == shape and result.dtype == np.float64
    assert result.flags.c_contiguous and result.flags.writeable and result.flags.owndata
    assert bits(result) == bits(reference)
    assert all(result is not a for a in arrays)
    # repeated calls see the current inputs
    arrays[0] += 0.5
    assert bits(f(*arrays)) == bits(f(*[strided(a) for a in arrays]))
    # scalar arguments (Python float and int) with the arrays
    if len(arguments) > 1:
        for scalar in (0.75, 2, -3):
            values = [scalar] + arrays[1:]
            assert bits(f(*values)) == bits(f(*([scalar] + [strided(a) for a in arrays[1:]])))


COMPILED = [
    ((x,), Piecewise((x**2, x < 1), (sin(x), True)), 'numpy'),
    ((x,), gamma(x), None),
    ((x, y), x**2*y + 3*x - y, None),
    ((x, y), Max(x, y) + besselj(0, x), None),
]


@pytest.mark.parametrize('arguments,expr,modules', COMPILED)
@pytest.mark.parametrize('shape', [(10,), (2, 5), (0,)])
def test_compiled_program_on_buffers_matches_the_copying_path(arguments, expr, modules, shape):
    f = lambdify(arguments, expr, modules)
    arrays = [samples(shape, 7 + i) for i in range(len(arguments))]
    result = f(*arrays)
    reference = f(*[strided(a) for a in arrays])
    assert result.shape == shape and result.dtype == np.float64 and result.flags.owndata
    assert bits(result) == bits(reference)
    if len(arguments) > 1:
        assert bits(f(arrays[0], 1.5)) == bits(f(strided(arrays[0]), 1.5))


def test_large_arrays_and_multiple_outputs():
    a = np.linspace(-3, 3, 10001)
    f = lambdify(x, exp(-x**2)*cos(x), 'numpy')
    assert bits(f(a)) == bits(f(strided(a)))
    g = lambdify(x, [x**2, Piecewise((x, x < 0), (sin(x), True))])
    result, reference = g(a), g(strided(a))
    assert isinstance(result, list) and len(result) == 2
    assert all(bits(r) == bits(s) for r, s in zip(result, reference))


def test_fallbacks_and_errors_are_unchanged():
    f = lambdify((x, y), sqrt(x**2 + y**2), 'numpy')
    # broadcasting shapes keep NumPy's evaluation
    assert f(np.ones((3, 1)), np.ones((1, 4))).shape == (3, 4)
    # float32 keeps the float64 VM precision
    r = f(np.ones(3, dtype=np.float32), np.ones(3, dtype=np.float32))
    assert r.dtype == np.float64
    # the compiled program keeps its equal-size requirement
    g = lambdify((x, y), Piecewise((x, x < y), (y, True)), 'numpy')
    with pytest.raises(ValueError, match='same size'):
        g(np.ones(3), np.ones(4))
    # an int beyond the float range keeps its error
    with pytest.raises(OverflowError):
        g(np.ones(3), 10**400)
    with pytest.raises(TypeError):
        f(np.ones(3))


def test_complex_arrays_keep_numpy_values():
    c = np.array([0.5 + 0.25j, -1.0 + 2.0j])
    f = lambdify(x, sin(x)*x + 1, 'numpy')
    assert f(c).tolist() == (np.sin(c)*c + 1).tolist()
    g = lambdify((x, y), x*y + y, 'numpy')
    assert g(c, 2.0).tolist() == (c*2.0 + 2.0).tolist()
    assert g(c, A[:2].copy()).tolist() == (c*A[:2] + A[:2]).tolist()


def test_complex_gamma_uses_the_generated_function():
    # SymPy 1.14 (SciPy's gamma); NumPy has no gamma for the object text
    f = lambdify(x, gamma(x))
    for _ in range(2):
        assert f(np.array([0.5 + 0.25j])).tolist() == [(1.3851135919886668-0.6731815357596996j)]
