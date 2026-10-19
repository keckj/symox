# Slang code generation

[README](../../README.md) · [Code generation](overview.md)

`symox.extra.slangcode` and the extensible `symox.extra.SlangCodePrinter` print real scalar expressions.
`slangcode_function` exports a complete Slang function with shared CSE. Scalar
results return directly; lists, tuples, matrices and arrays return a generated
struct whose `_sx_o_0`, `_sx_o_1`, ... fields follow depth-first row-major order.
Source export requires neither SlangPy nor a GPU:

```python
import symox as sm

x, y = sm.symbols("x y")
expr = sm.Matrix([[sm.sin(x*y)], [sm.sqrt(1+x*x+y*y)]])
source = sm.extra.slangcode_function("evaluate", (x, y), expr, dtype="float32")
```

For optional Python execution, install `symox[slang]`. SlangPy compiles and
vectorizes the scalar function for the chosen backend; no Torch or JAX runtime
is needed. CUDA and Vulkan are exercised by numerical tests, including software
Vulkan in CI. Other backend selections depend on the external runtime and device.

```python
import numpy as np
import slangpy as spy

device = spy.create_device(type=spy.DeviceType.cuda)  # or DeviceType.vulkan
f = sm.lambdify((x, y), expr, "slang", device=device)
inputs = spy.Tensor.from_numpy(device, np.linspace(0, 1, 1024, dtype=np.float32))
result = f(inputs, 0.3)       # GPU Tensor, shape (1024, 2, 1)
host_result = result.to_numpy()
f(inputs, 0.3, out=result)   # reuse the output allocation
```

Existing SlangPy tensors stay on device. NumPy inputs are uploaded; request
`return_type="numpy"` at construction for explicit host results. Lists and tuples
retain their structure, matrices and arrays have trailing expression dimensions,
and boolean outputs retain their logical type. `out=` reuses supplied tensors;
`_append_to=` records into a SlangPy command encoder for deferred submission.
Deferred execution cannot be combined with automatic NumPy output.

`slang_kernel` additionally exports an explicit `name_kernel` compute entry with
`StructuredBuffer` inputs, writable output buffers and a count. This lower-level
ABI uses equally sized contiguous scalar batches, one output buffer per expression
leaf; boolean outputs use uint32 values 0/1. The host controls allocation and
dispatch. SlangPy's higher-level adapter handles broadcasting and tensor strides.

The portable baseline is float32. Float64 source uses double literals and depends
on target capabilities; the Python adapter currently permits float64 execution
only with CUDA. CPU dispatch is not exposed by this adapter. Complex expressions
and unsupported special functions fail explicitly. CSE is disabled around
`Piecewise` graphs to keep inactive branches inside scalar short-circuit conditions.
Custom `user_functions` mappings name functions supplied by the surrounding Slang
module. Automatic differentiation is not part of this backend.

`python benches/bench_slang.py --help` describes the synchronized optical
benchmark, with source generation, module loading, first dispatch and warm
execution separated. It optionally compares against `torch.compile`. For one
local 65,536-ray float32 run, Slang CUDA took about 0.113 ms versus 0.095 ms for
Inductor; Slang is not assumed to outperform the framework compiler. Software
Vulkan inverse trigonometric builtins can also have lower numerical accuracy
than the tested NVIDIA implementation.

See the [Slang targets](https://docs.shader-slang.org/en/stable/external/slang/docs/user-guide/09-targets.html)
and [SlangPy vectorization](https://slangpy.shader-slang.org/en/latest/src/basics/firstfunctions.html)
documentation for external runtime details.
