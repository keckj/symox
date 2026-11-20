"""cudacode logical operators, inline powers, reciprocal_sqrt and cudacode_function.

Compilation runs with NVRTC (``runtime`` fixture, CUDA driver required for
execution) and with hipcc when present (compilation only).
"""
import ctypes as ct
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import symox as s
from symox.extra import cudacode, cudacode_function

x, y, z = s.symbols("x y z", real=True)


# ---------------------------------------------------------------- printer
def test_logical_operators():
    assert cudacode(s.And(x < 1, y > 2)) == "((y > 2) && (x < 1))"
    assert cudacode(s.Or(x < 1, y > 2)) == "((y > 2) || (x < 1))"
    assert cudacode(s.Xor(x < 1, y > 2, evaluate=False)).count("!=") == 1
    assert "&&" in cudacode(s.Piecewise((x, s.And(x < 1, y > 2)), (y, True)))


def test_inline_integer_powers_of_small_bases():
    assert cudacode((x - y)**2) == "(x - y)*(x - y)"
    assert cudacode(1/(x + y)**3) == "1.0/((x + y)*(x + y)*(x + y))"
    assert cudacode((x - y)**5) == "pow(x - y, 5)"
    assert cudacode(s.sin(x)**2) == "pow(sin(x), 2)"
    assert cudacode((x - y)**2, inline_pow=0) == "pow(x - y, 2)"
    assert cudacode((x - y)**2, dtype="float") == "(x - y)*(x - y)"


def test_reciprocal_sqrt_is_opt_in():
    e = 1/s.sqrt(x + y) + 1
    assert cudacode(e) == "1 + rsqrt(x + y)"  # historical default
    assert cudacode(e, use_rsqrt=False) == "1 + pow(x + y, -1.0/2.0)"
    assert cudacode(e, use_rsqrt=False, reciprocal_sqrt=True) == "1 + 1.0/sqrt(x + y)"
    assert cudacode(e, dtype="float", use_rsqrt=False, reciprocal_sqrt=True) == "1 + 1.0F/sqrtf(x + y)"
    k = s.extra.cuda_kernel("k", [x, y], [e], use_rsqrt=False, reciprocal_sqrt=True)
    assert "1.0/sqrt(x + y)" in k
    assert "M_PI" in s.extra.cuda_kernel("k", [x], [s.pi*x])
    assert "M_PI" not in s.extra.cuda_kernel("k", [x], [s.pi*x], math_macros=False)


# ---------------------------------------------------------------- device function
def _pair_case():
    xs, ys, w = s.symbols("x0:3", real=True), s.symbols("y0:3", real=True), s.Symbol("w", real=True)
    d2 = sum((a - b)**2 for a, b in zip(xs, ys))
    outputs = [s.exp(-d2/(2*w**2)), s.Piecewise((1/s.sqrt(d2 - 1), d2 > 1), (0, True)),
               s.And(d2 < 4, s.Or(xs[0] > ys[0], s.Not(ys[1] < 0.5))), s.pi*xs[1] + s.Rational(1, 3)]
    return [list(xs), list(ys), w], outputs, (xs, ys, w)


@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("dialect", ["cuda", "hip"])
def test_device_function_structure(dtype, dialect):
    args, outputs, _ = _pair_case()
    src = cudacode_function("pair", args, outputs, dtype=dtype, dialect=dialect,
                            arg_names=["xi", "yj", "w"], result_type="pair_out",
                            result_fields=["o0", "o1", "o2", "o3"])
    T = "double" if dtype == "float64" else "float"
    assert f"__device__ __forceinline__ pair_out pair(const {T}* __restrict__ xi, const {T}* __restrict__ yj, const {T} w)" in src
    assert "    bool o2;" in src and f"    {T} o0;" in src
    for bad in ("M_PI", "rsqrt", "__expf", "#include"):
        assert bad not in src, bad
    assert not re.search(r"pow[f]?\((?:[^()]|\([^()]*\))*, -?[234]\)", src)  # integer powers inlined
    assert "1.0/sqrt" not in src  # reciprocal_sqrt is opt-in
    # Piecewise present: no shared temporary hoisted out of the branch
    assert "_sx_t" not in src


def test_device_function_scalar_cse_and_validation():
    src = cudacode_function("f", [x, y], s.exp(x + y)*s.sin(x + y), dtype="float32")
    assert "__device__ __forceinline__ float f(const float _sx_a0, const float _sx_a1)" in src
    assert "_sx_t0" in src
    src = cudacode_function("g", [[x, y, z]], x*z)
    assert "_sx_a0[1]" not in src and "_sx_a0[2]" in src
    assert "(void)_sx_a1" in cudacode_function("h", [x, y], 2*x)
    with pytest.raises(ValueError):
        cudacode_function("f", [x], y)
    with pytest.raises(ValueError):
        cudacode_function("f", [x, x], x)
    with pytest.raises(ValueError):
        cudacode_function("_sx_f", [x], x)
    with pytest.raises(ValueError):
        cudacode_function("f", [x], x, dialect="opencl")
    with pytest.raises(ValueError):
        cudacode_function("f", [x], x, result_type="R")
    with pytest.raises(TypeError):
        cudacode_function("f", [[x, [y]]], x)


# ---------------------------------------------------------------- compilation
def _wrapper(dtype, dialect):
    args, outputs, _ = _pair_case()
    T = "double" if dtype == "float64" else "float"
    src = cudacode_function("pair", args, outputs, dtype=dtype, dialect=dialect,
                            arg_names=["xi", "yj", "w"], result_type="pair_out",
                            result_fields=["o0", "o1", "o2", "o3"])
    head = "#include <hip/hip_runtime.h>\n" if dialect == "hip" else ""
    return head + src + f"""
extern "C" __global__ void eval_pairs(long long m, long long n, const {T}* x, const {T}* y, const {T} w, {T}* out) {{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= m * n) return;
    const pair_out r = pair(x + (k / n) * 3, y + (k % n) * 3, w);
    out[4 * k] = r.o0; out[4 * k + 1] = r.o1; out[4 * k + 2] = ({T})r.o2; out[4 * k + 3] = r.o3;
}}
"""


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
def test_device_function_runs_on_cuda(runtime, dtype):
    from symox.extra.codegen.cuda import _cleanup
    rng = np.random.default_rng(1)
    m, n, w = 5, 7, 0.8
    npd = np.float64 if dtype == "float64" else np.float32
    xa, ya = rng.normal(size=(m, 3)).astype(npd), rng.normal(size=(n, 3)).astype(npd)
    xa[0] = ya[0]  # d2 = 0: inactive branch would be 1/sqrt(-1)
    out = np.zeros((m, n, 4), dtype=npd)
    kernel = s.extra.compile_cuda(_wrapper(dtype, "cuda"), "eval_pairs", cache_dir=False)
    driver, allocations = kernel._driver, []
    try:
        with kernel._current():
            try:
                for array in (xa, ya, out):
                    pointer = ct.c_uint64()
                    driver.call("cuMemAlloc_v2", ct.byref(pointer), array.nbytes)
                    allocations.append(pointer)
                for array, pointer in zip((xa, ya), allocations):
                    driver.call("cuMemcpyHtoD_v2", pointer, ct.c_void_p(array.ctypes.data), array.nbytes)
                scalar = ct.c_double(w) if dtype == "float64" else ct.c_float(w)
                kernel.launch([ct.c_longlong(m), ct.c_longlong(n), allocations[0], allocations[1], scalar, allocations[2]],
                              grid=1, block=64)
                driver.call("cuCtxSynchronize")
                driver.call("cuMemcpyDtoH_v2", ct.c_void_p(out.ctypes.data), allocations[2], out.nbytes)
            finally:
                _cleanup([(driver.call, ("cuMemFree_v2", p)) for p in reversed(allocations)])
    finally:
        kernel.close()
    X, Y = xa.astype(np.float64)[:, None, :], ya.astype(np.float64)[None, :, :]
    d2 = ((X - Y)**2).sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        ref = np.stack([np.exp(-d2/(2*w*w)), np.where(d2 > 1, 1/np.sqrt(d2 - 1), 0.0),
                        ((d2 < 4) & ((X[..., 0] > Y[..., 0]) | ~(Y[..., 1] < 0.5))).astype(float),
                        np.pi*X[..., 1] + 1/3 + 0*d2], axis=-1)
    tol = 1e-5 if dtype == "float32" else 1e-13
    assert np.isfinite(out).all()
    np.testing.assert_allclose(out, ref, rtol=tol, atol=tol)


def _hipcc():
    for candidate in (shutil.which("hipcc"), "/opt/rocm/bin/hipcc"):
        if candidate and Path(candidate).exists():
            return candidate
    return None


@pytest.mark.skipif(_hipcc() is None, reason="hipcc unavailable")
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_device_function_compiles_with_hipcc(dtype, tmp_path):
    source = tmp_path / "pair.hip"
    source.write_text(_wrapper(dtype, "hip"))
    arch = os.environ.get("SYMOX_HIP_ARCH", "gfx1201")
    result = subprocess.run([_hipcc(), f"--offload-arch={arch}", "-O2", "-c", str(source), "-o", str(tmp_path / "pair.o")],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-4000:]
