# PyTorch code generation

[README](../../README.md) · [Code generation](overview.md)

Install `symox[torch]`, or install the PyTorch build appropriate for your GPU
separately. Importing symox does not import PyTorch. The PyTorch backend builds
tensor operations from the native expression tree without a SymPy dependency:

```python
import torch
import symox as sm

x, y = sm.symbols("x y")
expression = sm.Matrix([[sm.sin(x*y)], [sm.sqrt(1 + x*x + y*y)]])
f = sm.lambdify((x, y), expression, "torch")
xv = torch.linspace(0, 1, 1024, requires_grad=True)
result = f(xv, 0.3)                      # shape: (1024, 2, 1)
result.sum().backward()

compiled = sm.extra.torch_lambdify(
    (x, y), expression, compile=True,
    compile_options={"fullgraph": True},  # default backend: Inductor
)
result = compiled(xv, 0.3)
```

Tensor inputs determine the device and promoted dtype. Python scalar inputs
use PyTorch's default floating dtype when no tensor is supplied. Every call
also accepts `dtype=` and `device=`. Scalar outputs broadcast to the common
input shape. Lists and tuples retain their structure; symbolic matrices and
NumPy arrays become tensors with batch dimensions followed by expression
dimensions. Boolean outputs retain boolean dtype, and complex values remain
complex. Existing input tensors stay connected to autograd.
Constant complex scalar and matrix outputs also support CUDA Inductor with
`fullgraph=True`, including functions without inputs.

`torchcode_function("my_function", (x, y), expression)` returns a complete
Python module that needs only PyTorch at runtime. Generated callables expose
the full module through `__source__`; `inspect.getsource` retrieves the function
definition. `torchcode` and
`symox.extra.printing.pytorch.PyTorchPrinter` provide lower-level expression printing;
these expression snippets use the `_sx_tensor` helper included in exported
modules. Custom functions can be supplied through `user_functions`; source
export requires qualified import paths, whereas `torch_lambdify` also accepts
callables.

CSE shares work across outputs. For `Piecewise`, CSE is currently disabled to
keep expressions inside their branch masks. The printer neutralizes inactive
arguments of supported operations before evaluating branches, including
division, square roots and logarithms, so invalid inactive branches do not
poison gradients. Custom function mappings are responsible for their own
domains and autograd behavior. Unsupported symbolic functions fail during
generation instead of silently switching numerical backends.

Compilation is optional and its errors are propagated. The first call includes
tracing and compilation; shape, dtype or device changes may cause recompilation.
Use `compile_options` to pass PyTorch options such as `backend`, `dynamic`,
`fullgraph` or `mode`. Performance depends on workload and batch size: measure
warm execution separately from compilation, for example:

```bash
python benches/bench_pytorch.py --help
```

See PyTorch's [compilation documentation](https://docs.pytorch.org/docs/stable/generated/torch.compile.html)
and [autograd explanation of invalid operations](https://docs.pytorch.org/docs/stable/notes/autograd.html#division-by-zero-in-autograd).

For an indicative local measurement, the benchmark's default refraction kernel
with 65,536 float32 rays on an RTX 5090 took 0.549 ms in eager mode and 0.098 ms
with Inductor after compilation (5.6x). With gradients, it took 1.961 ms versus
0.552 ms (3.6x). The first compiled call took 2.4–4.2 seconds, excluding wrapper
setup. These are local measurements with PyTorch 2.11, not performance guarantees.
