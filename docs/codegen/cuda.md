# Native CUDA execution with NVRTC

[README](../../README.md) · [Code generation](overview.md)

`symox.extra.cuda_lambdify` compiles native CUDA source and executes it through the CUDA
Driver API. It requires the NVIDIA driver and NVRTC shared libraries, with no
PyTorch, JAX, CuPy or CUDA Python package dependency. Library loading is lazy.

```python
import numpy as np
import symox as sm

x, y = sm.symbols("x y")
f = sm.extra.cuda_lambdify((x, y), sm.sin(x) + y, dtype="float64")
result = f(np.arange(3.0)[:, None], np.arange(2.0))  # shape (3, 2)
f.close()
```

Host arrays broadcast and may be noncontiguous. Lists/tuples retain structure;
matrix and array dimensions follow broadcast batch dimensions. Results are
NumPy arrays, with synchronous transfers to and from the GPU. Boolean outputs,
float32/float64, safely renamed arguments, CSE and lazy Piecewise branches are
supported. Complex-valued expressions are not supported by this interface.
`lambdify(..., "cuda")` and `lambdify(..., "nvrtc")` use the same path.
`cuda_lambdify` also accepts `device`, `cse`, extra NVRTC `options` and
`fast_math=True` (adds `--use_fast_math`, off by default); the callable exposes
`kernel`, `__source__` and `close`.

`compile_cuda(source, name)` loads an `extern "C"` kernel for direct use with
device pointers and streams. Its `CudaKernel.launch` method accepts explicitly
typed ctypes arguments; callers own their buffers and asynchronous lifetimes.
Buffers and streams must belong to the selected device's primary CUDA context.
Raw launches are asynchronous; `synchronize()` and `close()` wait for that
primary context. Kernels support context managers. PTX is cached atomically under
`$XDG_CACHE_HOME/symox/nvrtc` (or `~/.cache/symox/nvrtc`), keyed by source,
compiler flags, NVRTC version and GPU architecture. Set `cache_dir=False` to
disable disk caching. Compilation errors include the NVRTC log.

`SYMOX_CUDA_LIBRARY` and `SYMOX_NVRTC_LIBRARY` can select library paths.
Run GPU integration tests with
`SYMOX_CUDA_REQUIRE_RUNTIME=1 python -m pytest tests/extra/test_cuda_runtime.py`.
See the [NVRTC documentation](https://docs.nvidia.com/cuda/nvrtc/).

## Device functions and printer options

`cudacode_function(name, args, exprs, ...)` returns a standalone
`__device__ __forceinline__` function for CUDA or HIP (`dialect="hip"`; the
source is the same CUDA-compatible C++ and compiles with NVRTC, nvcc, hipcc and
hiprtc). It needs no header and no math macro: constants such as `pi` are
printed as literals and `INFINITY`/`NAN` are defined only when absent.

```python
xs, ys = sm.symbols("x0:3", real=True), sm.symbols("y0:3", real=True)
d2 = sum((a - b)**2 for a, b in zip(xs, ys))
src = sm.extra.cudacode_function(
    "pair", [list(xs), list(ys)], [sm.exp(-d2), d2 < 1], dtype="float32",
    arg_names=["xi", "yj"], result_type="pair_out", result_fields=["value", "near"])
```

A Symbol argument is passed by value; a flat list of Symbols is passed as
`const T* __restrict__` and only used components are read. Several outputs
return a named struct (`bool` fields for boolean expressions). CSE is disabled
when a `Piecewise` is present, so an inactive branch is never evaluated.

`cudacode` prints `And`, `Or`, `Not`, `Xor`, `Implies`, `Equivalent` and `ITE`
as parenthesized C++ boolean expressions. Integer powers up to `inline_pow`
(default 4) of symbols and of small arithmetic bases (sums and products of at
most 16 symbols, numbers and inline powers) are written as products, for
example `(x - y)*(x - y)`; function calls are never repeated (`pow(sin(x), 2)`
is kept), and `inline_pow=0` restores `pow` everywhere. These products are
correctly rounded IEEE operations. `tests/extra/test_cuda_precedence.py`
compares printed products, quotients, negative powers and Piecewise conditions
with NumPy under several option sets, on a GPU through NVRTC and on the host
through a C++ shim.

`reciprocal_sqrt=True` (default `False`, effective with `use_rsqrt=False`)
prints `x**(-1/2)` as `1.0/sqrt(x)`: two correctly rounded operations instead of
`pow(x, -1.0/2.0)`, but the result at `x = -0.0` is `-inf` instead of `+inf`.
The approximate `rsqrt` intrinsic keeps its historical default in `cudacode`
and `cuda_kernel` (`use_rsqrt=True`); `cudacode_function` never enables fast
math and leaves `use_rsqrt` and `reciprocal_sqrt` off by default. `cuda_kernel` and `cuda_device_function`
also accept `math_macros=False` for NVRTC, which does not define `M_PI`.
