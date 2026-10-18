# Native CUDA execution with NVRTC

[README](../../README.md) · [Code generation](overview.md)

`cuda_lambdify` compiles native CUDA source and executes it through the CUDA
Driver API. It requires the NVIDIA driver and NVRTC shared libraries, with no
PyTorch, JAX, CuPy or CUDA Python package dependency. Library loading is lazy.

```python
import numpy as np
import symox as sm

x, y = sm.symbols("x y")
f = sm.cuda_lambdify((x, y), sm.sin(x) + y, dtype="float64")
result = f(np.arange(3.0)[:, None], np.arange(2.0))  # shape (3, 2)
f.close()
```

Host arrays broadcast and may be noncontiguous. Lists/tuples retain structure;
matrix and array dimensions follow broadcast batch dimensions. Results are
NumPy arrays, with synchronous transfers to and from the GPU. Boolean outputs,
float32/float64, safely renamed arguments, CSE and lazy Piecewise branches are
supported. Complex-valued expressions are not supported by this interface.
`lambdify(..., "cuda")` and `lambdify(..., "nvrtc")` use the same path.

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
`SYMOX_CUDA_REQUIRE_RUNTIME=1 python -m pytest tests/test_cuda_runtime.py`.
See the [NVRTC documentation](https://docs.nvidia.com/cuda/nvrtc/).
