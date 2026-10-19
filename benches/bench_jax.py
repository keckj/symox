"""Time generated JAX optical kernels, including first compilation separately.

Example: python benches/bench_jax.py --batch 65536 --backward --piecewise
Requires the optional jax dependency. No Torch, Poctix, or CuPy is imported.
"""

import argparse
import json
import statistics
import time

import jax
import jax.numpy as jnp
import numpy as np

import symox as sm


def optical_kernel(piecewise=False):
    x, y, z, angle, mu = sm.symbols("x y z angle mu")
    norm = sm.sqrt(x*x + y*y + z*z)
    dx = (x*sm.cos(angle) - z*sm.sin(angle))/norm
    dy = y/norm
    dz = (x*sm.sin(angle) + z*sm.cos(angle))/norm
    discriminant = 1 - mu*mu*(1-dz*dz)
    transmitted = (sm.Piecewise((sm.sqrt(discriminant), discriminant > 0), (0, True))
                   if piecewise else sm.sqrt(discriminant))
    return (x, y, z, angle, mu), sm.Matrix([[mu*dx], [mu*dy], [transmitted]])


def benchmark(device, batch, repeats, warmup, backward, dtype, piecewise=False):
    args, expression = optical_kernel(piecewise)
    target = jax.devices(device)[0]
    rng = np.random.default_rng(42)
    inputs = [jax.device_put(jnp.asarray(rng.standard_normal(batch), dtype=dtype), target)
              for _ in range(3)]
    inputs += [jax.device_put(jnp.asarray(.2, dtype=dtype), target),
               jax.device_put(jnp.asarray(.75, dtype=dtype), target)]
    inputs = tuple(inputs)

    def synchronize(value):
        return jax.tree_util.tree_map(lambda leaf: leaf.block_until_ready(), value)

    rows, reference = [], None
    for mode, cse, jit in [("eager", False, False), ("eager_cse", True, False),
                            ("jitted_cse", True, True)]:
        start = time.perf_counter()
        function = sm.extra.jax_lambdify(args, expression, cse=cse, jit=jit)
        if backward:
            def loss(*values):
                output = function(*values)
                return jnp.sum(output), output
            differentiated = jax.value_and_grad(loss, argnums=tuple(range(len(inputs))), has_aux=True)
            invoke_function = jax.jit(differentiated) if jit else differentiated

            def invoke():
                (_, value), gradients = invoke_function(*inputs)
                return value, gradients
        else:
            def invoke():
                return function(*inputs), ()
        codegen_ms = (time.perf_counter()-start)*1000
        synchronize(inputs)
        start = time.perf_counter()
        result = synchronize(invoke())
        first_call_ms = (time.perf_counter()-start)*1000
        if reference is None:
            reference = result
        else:
            tolerance = 2e-4 if dtype == jnp.float32 else 1e-10
            for actual, expected in zip(jax.tree_util.tree_leaves(result), jax.tree_util.tree_leaves(reference)):
                np.testing.assert_allclose(actual, expected, rtol=tolerance, atol=tolerance)
        for _ in range(warmup):
            synchronize(invoke())
        timings = []
        for _ in range(repeats):
            start = time.perf_counter()
            synchronize(invoke())
            timings.append((time.perf_counter()-start)*1000)
        rows.append({"mode": mode, "codegen_ms": codegen_ms, "first_call_ms": first_call_ms,
                     "warm_median_ms": statistics.median(timings), "warm_min_ms": min(timings)})
    eager_ms = rows[0]["warm_median_ms"]
    for row in rows:
        row["speedup_vs_eager"] = eager_ms/row["warm_median_ms"]
    return {"jax": jax.__version__, "device": str(target), "hardware": target.device_kind,
            "dtype": str(jnp.dtype(dtype)), "batch": batch, "backward": backward,
            "piecewise": piecewise, "repeats": repeats, "results": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cpu", help="JAX backend: cpu, gpu, or tpu")
    parser.add_argument("--batch", type=int, default=65536)
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--backward", action="store_true")
    parser.add_argument("--piecewise", action="store_true")
    parser.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    options = parser.parse_args()
    if options.batch < 1 or options.repeats < 1 or options.warmup < 0:
        parser.error("batch/repeats must be positive and warmup nonnegative")
    if options.dtype == "float64":
        jax.config.update("jax_enable_x64", True)
    print(json.dumps(benchmark(options.device, options.batch, options.repeats, options.warmup,
                               options.backward, getattr(jnp, options.dtype), options.piecewise), indent=2))
