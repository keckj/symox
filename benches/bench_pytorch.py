"""Benchmark exported optical kernels, reporting compile and warm times separately.

Example: python benches/bench_pytorch.py --device cuda --batch 65536 --backward
Requires the optional torch dependency. No Poctix or CuPy runtime is used.
"""

import argparse
import json
import os
import statistics
import time

# Avoid unnecessary parallel compiler processes on machines with many cores.
os.environ.setdefault("TORCHINDUCTOR_COMPILE_THREADS", "1")

import torch

import symox as sm


def optical_kernel(piecewise=False):
    """Rotate a ray direction, then refract it through a planar interface."""
    x, y, z, angle, mu = sm.symbols("x y z angle mu")
    norm = sm.sqrt(x*x + y*y + z*z)
    dx = (x*sm.cos(angle) - z*sm.sin(angle))/norm
    dy = y/norm
    dz = (x*sm.sin(angle) + z*sm.cos(angle))/norm
    discriminant = 1 - mu*mu*(1-dz*dz)
    transmitted = (sm.Piecewise((sm.sqrt(discriminant), discriminant > 0), (0, True))
                   if piecewise else sm.sqrt(discriminant))
    return (x, y, z, angle, mu), sm.Matrix([[mu*dx], [mu*dy], [transmitted]])


def benchmark(device, batch, repeats, warmup, backward, backend, dtype, piecewise=False):
    args, expression = optical_kernel(piecewise)
    rng = torch.Generator(device=device).manual_seed(42)
    inputs = [torch.randn(batch, device=device, dtype=dtype, generator=rng)
              for _ in range(3)]
    inputs += [torch.tensor(.2, device=device, dtype=dtype),
               torch.tensor(.75, device=device, dtype=dtype)]
    inputs = tuple(value.requires_grad_(backward) for value in inputs)

    def sync():
        if device.startswith("cuda"):
            torch.cuda.synchronize(device)

    rows = []
    reference = None
    for mode, cse, compile_ in [("eager", False, False), ("eager_cse", True, False),
                                ("compiled_cse", True, True)]:
        start = time.perf_counter()
        function = sm.torch_lambdify(args, expression, cse=cse, compile=compile_,
                                     compile_options={"backend": backend, "fullgraph": True}
                                     if compile_ else None)
        codegen_ms = (time.perf_counter()-start)*1000

        def invoke():
            value = function(*inputs)
            gradients = torch.autograd.grad(value.sum(), inputs) if backward else ()
            return value, gradients

        sync()
        start = time.perf_counter()
        result = invoke()
        sync()
        first_call_ms = (time.perf_counter()-start)*1000
        if reference is None:
            reference = result
        else:
            torch.testing.assert_close(result[0], reference[0])
            for actual, expected in zip(result[1], reference[1]):
                # Shared scalar parameters reduce many float32 gradients;
                # CSE and fusion change the addition order of these sums.
                tolerance = 1e-4 if dtype == torch.float32 else 1e-10
                torch.testing.assert_close(actual, expected, rtol=tolerance, atol=tolerance)
        for _ in range(warmup):
            invoke()
        sync()
        times = []
        for _ in range(repeats):
            start = time.perf_counter()
            invoke()
            sync()
            times.append((time.perf_counter()-start)*1000)
        rows.append({"mode": mode, "codegen_ms": codegen_ms,
                     "first_call_ms": first_call_ms, "warm_median_ms": statistics.median(times),
                     "warm_min_ms": min(times)})
    eager_ms = rows[0]["warm_median_ms"]
    for row in rows:
        row["speedup_vs_eager"] = eager_ms/row["warm_median_ms"]
    return {"torch": torch.__version__, "device": device,
            "hardware": torch.cuda.get_device_name(device) if device.startswith("cuda") else "cpu",
            "dtype": str(dtype), "batch": batch, "backward": backward, "piecewise": piecewise,
            "backend": backend, "repeats": repeats, "results": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--batch", type=int, default=65536)
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--backward", action="store_true")
    parser.add_argument("--piecewise", action="store_true", help="Include total internal reflection branch")
    parser.add_argument("--backend", default="inductor")
    parser.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    parser.add_argument("--threads", type=int, default=1)
    options = parser.parse_args()
    if options.batch < 1 or options.repeats < 1 or options.warmup < 0:
        parser.error("batch/repeats must be positive and warmup nonnegative")
    torch.set_num_threads(options.threads)
    print(json.dumps(benchmark(options.device, options.batch, options.repeats, options.warmup,
                               options.backward, options.backend, getattr(torch, options.dtype), options.piecewise),
                     indent=2))
