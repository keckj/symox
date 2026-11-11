"""Operator precedence of the CUDA and SYCL printers, checked numerically.

A large family of products, quotients and negative powers (signs, compound
bases, n = 2..4, square roots, Piecewise) is printed with several option sets
and evaluated against NumPy ``lambdify`` in float64:

- CUDA source compiled with NVRTC and run on a GPU (``runtime`` fixture);
- the same CUDA source, and the SYCL source, compiled for the host with a C++
  shim (no GPU required) when a C++ compiler is available.

Regression: ``-x/s**2`` was printed ``-x/s*s`` (a negative coefficient lowered
the product precedence, leaving the inline denominator unparenthesized).
"""
import ctypes as ct
import itertools
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest
import symox as s
from symox import _core
from symox.extra import cudacode

x, y, a, b, c, d, w = SYMBOLS = s.symbols("x y a b c d w", real=True)


def corpus():
    bases = [w, x - y, c + d, 2 * w, w * d, c + d * w]
    out = []
    for base, n, sign in itertools.product(bases, (2, 3, 4), (1, -1)):
        p = base**n
        out += [
            sign / p, sign * x / p, sign * a * b / (c * p), sign * x / (p * d),
            sign * (x - y)**2 / (2 * p), x + sign / p, sign * x / p + y / c**2,
            s.exp(sign * x / p) / d, sign * x * p / c, sign * p, -p / (c * d),
            s.Piecewise((sign * x / p, x > 0), (a / p, True)),
            s.Piecewise((sign / p, a < b), (sign * x / (c * d**n), True)) * y,
            (sign * x / p) / (sign * a - 3), 1 / (1 + (sign * x / p)**2),
            sign * x / (p * s.sqrt(c)), sign * a / (p + c),
        ]
    out += [
        -x / s.sqrt(c), x / s.sqrt(c * d), -x * d**s.Rational(-3, 2), -1 / (x**2 * s.sqrt(c)),
        c**s.Rational(-1, 2) + 1, -c**s.Rational(-1, 2), (c + d)**s.Rational(-1, 2) * (-x),
        -x / (c * d), -x / c, -(x - y) / (c * d**2), a / b / c, -a / (b / c),
        x * w**s.Float(-1.0), -x * w**s.Float(-2.0), -x / (w**2 * c**3 * d**4),
        s.Max(x, y) / w**2 * -1, -s.Abs(x) / (x - y)**2, -x / s.exp(c)**2,
        s.Piecewise((1, -x / w**2 < -0.5), (0, True)), s.Piecewise((x, s.And(-a / c**2 < 0, x / d**3 > -1)), (y, True)),
        (x + a / w**2)**2, (x - y)**2 * (c + d)**-2, -(c + d)**-3 * x, x / (-(c + d))**3,
        2**(-x / w**2), s.sqrt(c + d / w**2) / -w**2, -x / w**2 / (c + d)**2, -(x / w**2) * (y / c**3),
    ]
    return out


EXPRS = corpus()
N = 64


def inputs():
    rng = np.random.default_rng(42)
    sign = np.where(rng.random(N) < 0.5, -1.0, 1.0)
    return {
        x: sign * rng.uniform(0.5, 2, N), y: rng.uniform(3, 4, N), a: -sign * rng.uniform(0.5, 2, N),
        b: rng.uniform(0.5, 2, N), c: rng.uniform(0.5, 2, N), d: rng.uniform(0.5, 2, N), w: rng.uniform(0.5, 2, N),
    }


def reference(dtype="double"):
    """float64 NumPy values on the inputs rounded to the tested dtype."""
    npd = np.float64 if dtype == "double" else np.float32
    data = inputs()
    values = [data[v].astype(npd).astype(np.float64) for v in SYMBOLS]
    out = []
    for e in EXPRS:
        r = s.lambdify(SYMBOLS, e, "numpy")(*values)
        out.append(np.broadcast_to(np.asarray(r, dtype=np.float64), (N,)))
    return np.array(out)


OPTIONS = [
    dict(dtype="double"),
    dict(dtype="double", use_rsqrt=False),
    dict(dtype="double", use_rsqrt=False, reciprocal_sqrt=True),
    dict(dtype="double", inline_pow=0),
    dict(dtype="double", cse=False),
    dict(dtype="float"),
    dict(dtype="float", fast_math=True),
    dict(dtype="float", use_rsqrt=False, reciprocal_sqrt=True, cse=False),
]


def tolerance(opts):
    if opts["dtype"] == "double":
        return 1e-12
    return 2e-3 if opts.get("fast_math") else 1e-5


def kernel_source(opts, exprs=None):
    exprs = EXPRS if exprs is None else exprs
    return 'extern "C" ' + s.extra.cuda_kernel("k", list(SYMBOLS), exprs, math_macros=False, **opts)


CHUNK = 256  # kernel parameter space is limited to 4 KiB on older GPUs


def check(values, opts, where):
    ref = reference(opts["dtype"])
    bad = []
    for e, got, want in zip(EXPRS, values, ref):
        if not np.allclose(got, want, rtol=tolerance(opts), atol=tolerance(opts)):
            bad.append(f"{where} {opts}: {e} -> {cudacode(e, **{k: v for k, v in opts.items() if k != 'cse'})}")
    assert not bad, "\n".join(bad[:20])


# ---------------------------------------------------------------- text
@pytest.mark.parametrize("expr,expected", [
    (-x / w**2, "-x/(w*w)"),
    (x / w**2, "x/(w*w)"),
    (-x / (x - y)**2, "-x/((x - y)*(x - y))"),
    (-x / y**3, "-x/(y*y*y)"),
    (-a * b / (c * d**2), "-a*b/(c*d*d)"),
    (-1 / (x * y**2), "-1/(x*y*y)"),
])
def test_negative_products_keep_inline_denominators_grouped(expr, expected):
    assert cudacode(expr) == expected


# ---------------------------------------------------------------- GPU
@pytest.fixture(scope="module")
def runtime():
    try:
        kernel = s.extra.compile_cuda('extern "C" __global__ void smoke() {}', "smoke", cache_dir=False)
    except (ImportError, s.extra.CudaError) as error:
        if os.environ.get("SYMOX_CUDA_REQUIRE_RUNTIME") == "1":
            pytest.fail(str(error))
        pytest.skip(f"CUDA/NVRTC unavailable: {error}")
    kernel.close()


@pytest.mark.parametrize("opts", OPTIONS, ids=lambda o: "-".join(f"{k}={v}" for k, v in o.items()))
def test_cuda_values_on_gpu(runtime, opts):
    npd = np.float64 if opts["dtype"] == "double" else np.float32
    data = inputs()
    arrays = [np.ascontiguousarray(data[v], dtype=npd) for v in SYMBOLS]
    rows = [_gpu_chunk(opts, EXPRS[i:i + CHUNK], arrays, npd) for i in range(0, len(EXPRS), CHUNK)]
    check(np.concatenate(rows).astype(np.float64), opts, "gpu")


def _gpu_chunk(opts, exprs, arrays, npd):
    from symox.extra.codegen.cuda import _cleanup
    outs = np.zeros((len(exprs), N), dtype=npd)
    kernel = s.extra.compile_cuda(kernel_source(opts, exprs), "k", cache_dir=False)
    driver, allocs = kernel._driver, []
    try:
        with kernel._current():
            try:
                for arr in arrays + list(outs):
                    p = ct.c_uint64()
                    driver.call("cuMemAlloc_v2", ct.byref(p), arr.nbytes)
                    allocs.append(p)
                for arr, p in zip(arrays, allocs):
                    driver.call("cuMemcpyHtoD_v2", p, ct.c_void_p(arr.ctypes.data), arr.nbytes)
                kernel.launch([ct.c_long(N), *allocs], grid=1, block=N)
                driver.call("cuCtxSynchronize")
                for row, p in zip(outs, allocs[len(arrays):]):
                    driver.call("cuMemcpyDtoH_v2", ct.c_void_p(row.ctypes.data), p, row.nbytes)
            finally:
                _cleanup([(driver.call, ("cuMemFree_v2", p)) for p in reversed(allocs)])
    finally:
        kernel.close()
    return outs


# ---------------------------------------------------------------- host C++
CXX = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")

CUDA_HOST_SHIM = r"""
#include <cmath>
#include <limits>
using namespace std;
struct dim3_ { long x; };
static dim3_ blockIdx = {0}, threadIdx = {0}, blockDim = {1}, gridDim = {1};
#define __global__
#define __device__
#define __forceinline__ inline
#define __restrict__
static inline double kp_rsqrt_(double v) { return 1.0 / sqrt(v); }
static inline float kp_rsqrtf_(float v) { return 1.0f / sqrtf(v); }
#define rsqrt(v) kp_rsqrt_(v)
#define rsqrtf(v) kp_rsqrtf_(v)
#define __expf(v) expf(v)
#define __logf(v) logf(v)
#define __sinf(v) sinf(v)
#define __cosf(v) cosf(v)
#define __tanf(v) tanf(v)
#define __powf(p, q) powf(p, q)
#define __fdividef(p, q) ((p) / (q))
"""


def _host_run(source, npd, nouts):
    tmp = tempfile.mkdtemp(prefix="symox_prec_")
    try:
        path, lib = os.path.join(tmp, "k.cpp"), os.path.join(tmp, "k.so")
        with open(path, "w") as fh:
            fh.write(source)
        r = subprocess.run([CXX, "-O1", "-std=c++17", "-shared", "-fPIC", "-o", lib, path], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
        so = ct.CDLL(lib)
        data = inputs()
        arrays = [np.ascontiguousarray(data[v], dtype=npd) for v in SYMBOLS]
        outs = np.zeros((nouts, N), dtype=npd)
        so.k.restype = None
        so.k(ct.c_long(N), *[arr.ctypes.data_as(ct.c_void_p) for arr in arrays],
             *[row.ctypes.data_as(ct.c_void_p) for row in outs])
        return outs.astype(np.float64)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@pytest.mark.skipif(CXX is None, reason="no host C++ compiler")
@pytest.mark.parametrize("opts", OPTIONS, ids=lambda o: "-".join(f"{k}={v}" for k, v in o.items()))
def test_cuda_values_on_host(opts):
    npd = np.float64 if opts["dtype"] == "double" else np.float32
    check(_host_run(CUDA_HOST_SHIM + kernel_source(opts), npd, len(EXPRS)), opts, "host")


SYCL_HOST_SHIM = r"""
#include <cmath>
#include <limits>
namespace sycl { using namespace std; }
"""


@pytest.mark.skipif(CXX is None, reason="no host C++ compiler")
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_sycl_printer_values_on_host(dtype):
    T = "double" if dtype == "float64" else "float"
    names = {v: f"v{i}" for i, v in enumerate(SYMBOLS)}
    params = ", ".join([f"const {T}* in{i}" for i in range(len(SYMBOLS))] + [f"{T}* out{k}" for k in range(len(EXPRS))])
    body = [f"    const {T} {names[v]} = in{i}[i];" for i, v in enumerate(SYMBOLS)]
    body += [f"    out{k}[i] = {_core._syclcode(e, dtype=dtype, symbol_names=names)};" for k, e in enumerate(EXPRS)]
    src = SYCL_HOST_SHIM + f'extern "C" void k(long n, {params}) {{\n  for (long i = 0; i < n; i++) {{\n' + "\n".join(body) + "\n  }\n}\n"
    npd = np.float64 if dtype == "float64" else np.float32
    check(_host_run(src, npd, len(EXPRS)), {"dtype": "double" if dtype == "float64" else "float"}, "sycl-host")
