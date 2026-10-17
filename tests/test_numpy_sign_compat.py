"""NumPy complex sign is consistent across the supported NumPy versions.

The unknown-complex printer intentionally differs from SymPy's direct
``numpy.sign`` spelling: NumPy 1.x has different complex semantics. Runtime
oracles and the simulated legacy ufunc verify the emitted expression instead.
"""
import warnings

import numpy as np
import pytest
import symox as sx


NUMPY2 = int(np.__version__.split('.')[0]) >= 2
PHASE_1_2 = complex(1 / np.sqrt(5.0), 2 / np.sqrt(5.0))
COMPLEX_NAN = complex(np.nan, np.nan)


class LegacySignNumpy:
    """Simulate NumPy 1.x complex sign while keeping the actual array APIs."""
    __version__ = "1.21.6"

    def __getattr__(self, name):
        return getattr(np, name)

    def sign(self, value):
        if np.iscomplexobj(value):
            real, imaginary = np.real(value), np.imag(value)
            return np.where(real != 0, np.sign(real), np.sign(imaginary)) + 0j
        return np.sign(value)


def run(source, values, numpy=np, argument='z'):
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        return eval(source, {'numpy': numpy}, {argument: values})


@pytest.mark.parametrize('numpy', [np, LegacySignNumpy()])
@pytest.mark.parametrize('dtype', [np.complex64, np.complex128])
def test_complex_sign_arrays_match_numpy2_including_zero_inf_and_nan(numpy, dtype):
    values = np.array([
        0, -0j, 1+2j, -3+4j, 2j, -2j,
        complex(np.inf, 0), complex(-np.inf, 2),
        complex(2, np.inf), complex(2, -np.inf),
        complex(np.inf, np.inf), complex(-np.inf, np.inf),
        complex(np.nan, 1), complex(1, np.nan),
        complex(np.inf, np.nan), complex(np.nan, np.inf),
    ], dtype=dtype).reshape(4, 4)
    source = sx.numpycode(sx.sign(sx.Symbol('z')))
    actual = run(source, values, numpy)
    expected = np.array([
        0j, 0j, PHASE_1_2, complex(-0.6, 0.8), 1j, -1j,
        1, -1, 1j, -1j,
        COMPLEX_NAN, COMPLEX_NAN, COMPLEX_NAN, COMPLEX_NAN,
        1, 1j,
    ], dtype=dtype).reshape(4, 4)
    np.testing.assert_allclose(actual, expected, equal_nan=True)
    if NUMPY2:
        np.testing.assert_allclose(actual, np.sign(values), equal_nan=True)
    assert actual.dtype == dtype
    assert actual.shape == values.shape


@pytest.mark.parametrize('value, golden', [
    (0, 0), (-2, -1), (2.0, 1),
    (float('inf'), 1), (float('-inf'), -1), (float('nan'), np.nan),
    (0j, 0j), (1+2j, PHASE_1_2), (-3+4j, complex(-0.6, 0.8)),
    (complex(float('inf'), 2), 1),
    (complex(float('nan'), 1), COMPLEX_NAN),
    (np.complex64(1+2j), PHASE_1_2), (np.complex128(1+2j), PHASE_1_2),
])
def test_scalar_and_runtime_real_dispatch_preserve_numpy_results(value, golden):
    source = sx.numpycode(sx.sign(sx.Symbol('z')))
    actual = run(source, value, LegacySignNumpy())
    expected = np.asarray(golden, dtype=np.asarray(value).dtype)[()]
    np.testing.assert_allclose(actual, expected, equal_nan=True)
    assert np.isscalar(actual)
    assert np.asarray(actual).dtype == np.asarray(expected).dtype


def test_subnormal_complex_division_uses_components_without_overflow():
    values = np.array([complex(1e-320, 1e-320), complex(-1e-320, 0), 0], dtype=np.complex128)
    source = sx.numpycode(sx.sign(sx.Symbol('z')))
    # NumPy 2 golden: the subnormal norm rounds before division, so this
    # component differs from sqrt(1/2). Preserve its exact binary value.
    component = float.fromhex('0x1.6a15d346a9f38p-1')
    expected = np.array([complex(component, component), -1, 0], dtype=np.complex128)
    actual = run(source, values)
    np.testing.assert_array_equal(actual, expected)
    if NUMPY2:
        np.testing.assert_array_equal(actual, np.sign(values))


def test_declared_real_arguments_keep_the_short_real_sign_path():
    real = sx.Symbol('r', real=True)
    source = sx.numpycode(sx.sign(real))
    assert source == 'numpy.sign(r)'
    values = np.array([-np.inf, -2, -0.0, 0.0, 2, np.inf, np.nan])
    np.testing.assert_equal(run(source, values, argument='r'), [-1, -1, 0, 0, 1, 1, np.nan])


def test_inline_binding_does_not_capture_user_symbols_or_repeat_argument():
    value = sx.Symbol('_sgn_value')
    other = sx.Symbol('_sgn_abs')
    expression = sx.sign(value)+other
    source = sx.numpycode(expression)
    actual = eval(source, {'numpy': LegacySignNumpy()}, {'_sgn_value': 1+2j, '_sgn_abs': 7})
    np.testing.assert_allclose(actual, PHASE_1_2+7)
    assert 'symox' not in source


def test_abs_derivative_of_complex_valued_function_has_correct_real_values():
    real = sx.Symbol('r', real=True)
    derivative = sx.Abs(real+sx.I).diff(real)
    source = sx.numpycode(derivative)
    values = np.array([-3., -1., 0., 1., 3.])
    actual = run(source, values, LegacySignNumpy(), argument='r')
    np.testing.assert_allclose(actual, values/np.sqrt(values**2+1), atol=1e-15)


@pytest.mark.skipif(not NUMPY2, reason='NumPy 2 fast path oracle')
def test_numpy2_calls_sign_once_without_running_complex_fallback():
    class ModernNumpy:
        __version__ = np.__version__
        calls = 0
        def sign(self, value):
            self.calls += 1
            return np.sign(value)
        def __getattr__(self, name):
            raise AssertionError(f'unnecessary fallback operation: {name}')
    module = ModernNumpy()
    values = np.array([0, 1+2j, complex(np.inf,2),complex(np.nan,1)])
    source = sx.numpycode(sx.sign(sx.Symbol('z')))
    np.testing.assert_equal(run(source, values, module), np.sign(values))
    assert module.calls == 1


@pytest.mark.parametrize('numpy', [np, LegacySignNumpy()])
def test_nontrivial_argument_is_evaluated_once(numpy):
    variable = sx.Symbol('z')
    argument = sx.Function('probe')(variable)
    source = sx.numpycode(sx.sign(argument), user_functions={'probe': 'probe'})
    calls = []
    def probe(value):
        calls.append(value)
        return value
    value = 1+2j
    actual = eval(source, {'numpy': numpy, 'probe': probe}, {'z': value})
    np.testing.assert_allclose(actual, PHASE_1_2)
    assert calls == [value]
