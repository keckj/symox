"""Benchmark a generated Slang optical kernel on GPU without autograd.

Example: python benches/bench_slang.py --backend cuda --batch 65536
Source generation, loading/compilation, first dispatch and synchronized warm
execution are reported separately. Inputs are uploaded before timing; warm
calls include output allocation. Optional torch.compile comparison measures the
same math, shape, precision and GPU, with uploads outside its timed region.
"""

import argparse
import json
import os
import statistics
import time

os.environ.setdefault("TORCHINDUCTOR_COMPILE_THREADS", "1")

import numpy as np
import slangpy as spy
import symox as sm


def optical_kernel(piecewise=False):
    x, y, z, angle, mu = sm.symbols("x y z angle mu")
    norm = sm.sqrt(x*x+y*y+z*z)
    dx = (x*sm.cos(angle)-z*sm.sin(angle))/norm
    dy = y/norm
    dz = (x*sm.sin(angle)+z*sm.cos(angle))/norm
    discriminant = 1-mu*mu*(1-dz*dz)
    transmitted = (sm.Piecewise((sm.sqrt(discriminant), discriminant > 0), (0, True))
                   if piecewise else sm.sqrt(discriminant))
    return (x, y, z, angle, mu), sm.Matrix([[mu*dx], [mu*dy], [transmitted]])


def _timed_call(invoke, sync):
    sync()
    start = time.perf_counter()
    result = invoke()
    sync()
    return result, (time.perf_counter()-start)*1000


def _warm_times(invoke, sync, warmup, repeats):
    for _ in range(warmup):
        invoke()
    sync()
    return [_timed_call(invoke, sync)[1] for _ in range(repeats)]


def benchmark(backend, batch, repeats, warmup, dtype, torch_baseline=False, piecewise=False):
    device = spy.create_device(type=getattr(spy.DeviceType, backend), enable_debug_layers=False)
    arguments, expression = optical_kernel(piecewise)
    np_dtype = getattr(np, dtype)
    rng = np.random.default_rng(42)
    inputs = [rng.normal(size=batch).astype(np_dtype) for _ in range(3)]
    inputs += [np.full(batch, .2, dtype=np_dtype), np.full(batch, 1.3 if piecewise else .75, dtype=np_dtype)]
    uploaded = tuple(spy.Tensor.from_numpy(device, value) for value in inputs)
    device.wait()
    rows = []
    reference = None
    for cse in (False, True):
        start = time.perf_counter()
        source = sm.extra.slangcode_function("optical", arguments, expression, dtype=dtype, cse=cse)
        source_ms = (time.perf_counter()-start)*1000
        start = time.perf_counter()
        function = sm.extra.slang_lambdify(arguments, expression, device=device, dtype=dtype, cse=cse)
        load_ms = (time.perf_counter()-start)*1000
        result, first_ms = _timed_call(lambda: function(*uploaded), device.wait)
        actual = result.to_numpy()
        if reference is None:
            a, b, c, theta, eta = inputs
            norm = np.sqrt(a*a+b*b+c*c)
            dx = (a*np.cos(theta)-c*np.sin(theta))/norm
            dy = b/norm
            dz = (a*np.sin(theta)+c*np.cos(theta))/norm
            reference = np.stack([eta*dx, eta*dy, np.sqrt(np.maximum(1-eta*eta*(1-dz*dz), 0))], axis=-1)[..., None]
        # sqrt(discriminant) is ill-conditioned at total internal reflection;
        # float32 fusion changes a few rays close to that boundary by ~1e-5.
        tolerance = 3e-5 if dtype == "float32" else 1e-10
        np.testing.assert_allclose(actual, reference, rtol=tolerance, atol=tolerance)
        times = _warm_times(lambda: function(*uploaded), device.wait, warmup, repeats)
        rows.append({"mode": "slang_cse" if cse else "slang", "source_generation_ms": source_ms,
                     "source_bytes": len(source.encode()), "load_compile_ms": load_ms,
                     "first_dispatch_ms": first_ms, "warm_median_ms": statistics.median(times),
                     "warm_min_ms": min(times)})
    if torch_baseline:
        if backend != "cuda":
            raise ValueError("--torch-baseline requires --backend cuda to compare the same GPU API")
        import torch
        if not torch.cuda.is_available():
            raise RuntimeError("PyTorch CUDA unavailable")
        torch_inputs = tuple(torch.from_numpy(value).cuda() for value in inputs)
        sync = torch.cuda.synchronize
        start = time.perf_counter()
        function = sm.extra.torch_lambdify(arguments, expression, compile=True,
                                     compile_options={"fullgraph": True})
        generation_ms = (time.perf_counter()-start)*1000
        result, first_ms = _timed_call(lambda: function(*torch_inputs), sync)
        np.testing.assert_allclose(result.cpu().numpy(), reference, rtol=tolerance, atol=tolerance)
        times = _warm_times(lambda: function(*torch_inputs), sync, warmup, repeats)
        rows.append({"mode": "torch_compile", "source_generation_ms": generation_ms,
                     "load_compile_ms": 0, "first_dispatch_ms": first_ms,
                     "warm_median_ms": statistics.median(times), "warm_min_ms": min(times)})
    return {"slangpy": spy.__version__, "backend": backend, "device": str(device.info.adapter_name),
            "dtype": dtype, "batch": batch, "piecewise": piecewise, "repeats": repeats, "warmup": warmup,
            "warm_includes": "dispatch, synchronization and output allocation; inputs already on device",
            "outputs": "batch x 3 x 1", "results": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=["cuda", "vulkan"], default="cuda")
    parser.add_argument("--batch", type=int, default=65536)
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    parser.add_argument("--torch-baseline", action="store_true")
    parser.add_argument("--piecewise", action="store_true", help="Include total internal reflection; currently disables CSE")
    options = parser.parse_args()
    if options.batch < 1 or options.repeats < 1 or options.warmup < 0:
        parser.error("batch/repeats must be positive and warmup nonnegative")
    print(json.dumps(benchmark(options.backend, options.batch, options.repeats, options.warmup,
                               options.dtype, options.torch_baseline, options.piecewise), indent=2))
