"""Native NVRTC execution, without an optional Python CUDA framework."""
import ctypes as ct
from concurrent.futures import ThreadPoolExecutor
import os
import subprocess
import sys

import numpy as np
import pytest
import symox as s
from symox.extra.codegen.cuda import _generate


def test_cuda_import_and_source_need_no_runtime():
    script = '''
import sys
import symox as s
from symox.extra.codegen.cuda import _generate
x = s.Symbol('x')
assert 'extern "C"' in _generate(x, x**2, 'float64', True, False)[0]
assert not {'torch','jax','cupy','cuda','sympy'}.intersection(sys.modules)
from symox.extra.codegen import cuda
assert cuda._DRIVER is None
'''
    subprocess.run([sys.executable, "-c", script], check=True, capture_output=True, text=True)


@pytest.mark.parametrize("settings", [dict(dtype="int64"), dict(cse="yes"), dict(fast_math=1)])
def test_source_option_validation(settings):
    x = s.Symbol("x")
    arguments = dict(dtype="float64", cse=True, fast_math=False)
    arguments.update(settings)
    with pytest.raises((TypeError, ValueError)):
        _generate(x, x, **arguments)


def test_cuda_argument_and_domain_validation():
    x, y = s.symbols("x y")
    with pytest.raises(ValueError):
        _generate((x, x), x, "float64", True, False)
    with pytest.raises(ValueError):
        _generate(x, x+y, "float64", True, False)
    with pytest.raises((ValueError, NotImplementedError)):
        _generate(x, s.I*x, "float64", True, False)
    with pytest.raises(TypeError):
        s.extra.cuda_lambdify(x, x, options="--use_fast_math")


@pytest.fixture(scope="module")
def runtime():
    try:
        kernel = s.extra.compile_cuda('extern "C" __global__ void smoke() {}', "smoke", cache_dir=False)
    except (ImportError, s.extra.CudaError) as error:
        if os.environ.get("SYMOX_CUDA_REQUIRE_RUNTIME") == "1":
            pytest.fail(str(error))
        pytest.skip(f"CUDA/NVRTC unavailable: {error}")
    kernel.close()


@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("cse", [False, True])
def test_broadcast_structured_piecewise(runtime, tmp_path, dtype, cse):
    x, y = s.symbols("x y")
    shared = s.sin(x+y)+s.cos(x*y)
    guarded = s.Piecewise((s.log(x)+1/x, x > 0), (x*x, True))
    f = s.extra.cuda_lambdify((x, (y,)), (s.Matrix([[shared**2, guarded], [1, x*y]]), [x > y, x+y]),
                       dtype=dtype, cse=cse, cache_dir=tmp_path)
    a = np.arange(-3, 3, dtype=dtype).reshape(2, 3).T[:, ::-1]
    b = np.array([1, 2], dtype=dtype)
    matrix, (boolean, total) = f(a, (b,))
    expected = np.empty(a.shape, dtype=dtype)
    np.log(a, where=a>0, out=expected)
    expected += np.divide(1, a, where=a>0, out=np.zeros_like(a))
    expected = np.where(a>0, expected, a*a)
    tolerance = 3e-6 if dtype == "float32" else 1e-13
    np.testing.assert_allclose(matrix[..., 0, 0], (np.sin(a+b)+np.cos(a*b))**2, rtol=tolerance, atol=tolerance)
    np.testing.assert_allclose(matrix[..., 0, 1], expected, rtol=tolerance, atol=tolerance)
    np.testing.assert_allclose(matrix[..., 1, 1], a*b)
    np.testing.assert_array_equal(boolean, a>b)
    assert boolean.dtype == np.bool_ and matrix.dtype == np.dtype(dtype)
    np.testing.assert_array_equal(total, a+b)
    assert f(np.empty((0, 2), dtype=dtype), (b,))[0].shape == (0, 2, 2, 2)
    with pytest.raises(TypeError):
        f(np.array([1j]), (b,))
    with pytest.raises(ValueError):
        f(a, b)
    f.close()
    with pytest.raises(RuntimeError, match="closed"):
        f(a, (b,))


def test_assumed_symbol_identity_and_constants(runtime, tmp_path):
    a, b = s.Symbol("same", positive=True), s.Symbol("same", negative=True)
    f = s.extra.cuda_lambdify((a, b), [a+2*b, s.true, s.Rational(1, 3)], cache_dir=tmp_path)
    actual = f(np.array([1., 2.]), -3.)
    np.testing.assert_array_equal(actual[0], [-5, -4])
    assert actual[1].dtype == np.bool_
    np.testing.assert_allclose(actual[2], 1/3)
    f.close()
    f = s.lambdify((), s.pi, "nvrtc", cache_dir=tmp_path)
    np.testing.assert_allclose(f(), np.pi, atol=1e-15)
    f.close()
    f = s.extra.cuda_lambdify((), s.Matrix(0, 0, []), cache_dir=tmp_path)
    assert f().shape == (0, 0)
    f.close()


def test_sparse_matrix_output(runtime, tmp_path):
    x = s.Symbol("x")
    expression = s.ImmutableSparseMatrix(2, 3, {(0, 1): s.sin(x), (1, 2): x*x})
    f = s.extra.cuda_lambdify(x, expression, cache_dir=tmp_path)
    try:
        values = np.array([-2., 0., 3.])
        result = f(values)
        expected = np.zeros((3, 2, 3))
        expected[:, 0, 1] = np.sin(values)
        expected[:, 1, 2] = values*values
        np.testing.assert_allclose(result, expected)
    finally:
        f.close()


def test_cache_context_threading_and_resource_lifetime(runtime, tmp_path):
    x = s.Symbol("x")
    first = s.extra.cuda_lambdify(x, s.sin(x), cache_dir=tmp_path)
    assert not first.kernel.cache_hit
    second = s.extra.cuda_lambdify(x, s.sin(x), cache_dir=tmp_path)
    assert second.kernel.cache_hit and first.kernel.cache_key == second.kernel.cache_key
    third = s.extra.cuda_lambdify(x, s.sin(x), dtype="float32", cache_dir=tmp_path)
    assert first.kernel.cache_key != third.kernel.cache_key
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(first, [np.arange(10.)+i for i in range(6)]))
    for i, result in enumerate(results):
        np.testing.assert_allclose(result, np.sin(np.arange(10.)+i))
    for function in (first, second, third):
        function.close()
        function.close()


def test_raw_kernel_launch_validation_and_compilation_log(runtime, tmp_path):
    with s.extra.compile_cuda('extern "C" __global__ void test(int n) {}', "test", cache_dir=tmp_path) as kernel:
        kernel.launch([ct.c_int(4)], grid=1, block=1)
        kernel.synchronize()
        with pytest.raises(TypeError):
            kernel.launch([4])
        with pytest.raises(ValueError):
            kernel.launch([ct.c_int(4)], block=0)
    with pytest.raises(s.extra.CudaError, match="nvrtcCompileProgram"):
        s.extra.compile_cuda("this is invalid cuda", "bad", cache_dir=tmp_path)
    with pytest.raises(s.extra.CudaError, match="cuModuleGetFunction"):
        s.extra.compile_cuda('extern "C" __global__ void test() {}', "missing", cache_dir=tmp_path)
    with pytest.raises(ValueError):
        s.extra.compile_cuda("", "bad", options=["--gpu-architecture=compute_50"])
