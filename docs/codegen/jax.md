# JAX code generation

[README](../../README.md) · [Code generation](overview.md)

Install `symox[jax]` for the CPU backend, or install a suitable JAX accelerator
build separately. JAX also stays optional and is not imported by `import symox`.
The API uses the same input/output structure and native traversal as PyTorch:

```python
import jax
import jax.numpy as jnp
import symox as sm

x, y = sm.symbols("x y")
expr = sm.sin(x*y) + sm.sqrt(1 + x*x)
f = sm.lambdify((x, y), expr, "jax", jit=True)
values = f(jnp.linspace(0, 1, 1024), 0.3)
gradient = jax.grad(f, argnums=0)(0.2, 0.3)
batched = jax.vmap(f, in_axes=(0, None))(jnp.linspace(0, 1, 1024), 0.3)
```

`jax_lambdify` accepts `cse`, `user_functions`, `jit` and `jit_options`.
`jaxcode_function` exports a module importing only JAX; `jaxcode` and
`symox.printing.jax.JaxPrinter` expose expression printing. Generated functions
accept `dtype=` and `device=` (a JAX Device or a backend name such as `"cpu"`);
these keywords are static when JIT is enabled. Dtypes obey JAX's configuration,
including its default disabling of 64-bit arithmetic. No global configuration
is changed. `Piecewise` uses the same masking precautions as PyTorch, with CSE
disabled for conditional graphs. Use `block_until_ready` when timing JAX work.
`python benches/bench_jax.py --help` describes the reproducible JAX benchmark,
with compilation time and warm execution reported separately.

Both backends support elementary functions, relations, Boolean operations and
selected special functions. In particular, `loggamma` requires a provably
positive argument: the real `lgamma`/`gammaln` operations do not implement
SymPy's complex branches. Unsupported special functions raise
`NotImplementedError`. These backends extend the supported SymPy surface;
symox does not yet implement all of SymPy.

See JAX's [JIT documentation](https://docs.jax.dev/en/latest/jit-compilation.html)
and [FAQ on gradients through `where`](https://docs.jax.dev/en/latest/faq.html#gradients-contain-nan-where-using-where).
