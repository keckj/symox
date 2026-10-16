"""Benchmarks of the printers, cse and code generation: symox vs SymPy.

    python benches/bench_printing.py
"""

from __future__ import annotations

import time

import sympy
from sympy.printing.numpy import NumPyPrinter

import symox


def timeit(fn, min_time=0.3, max_runs=20):
    fn()  # warm-up
    runs, total, best = 0, 0.0, float("inf")
    while total < min_time and runs < max_runs:
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        best = min(best, dt)
        total += dt
        runs += 1
    return best


def poly(m):
    x, y, z = m.symbols("x y z")
    return m.expand((x + y + z + 1) ** 15)


def nested(m):
    """A deep non-polynomial expression (derivatives of a rational/trig
    expression), typical of physics/engineering code generation."""
    x, y, z = m.symbols("x y z")
    e = m.sin(x * y + z) * m.exp(-x**2 / (1 + y**2)) / m.sqrt(1 + x**2 + z**2)
    return m.diff(e, x, 2) + m.diff(e, y, 2) + m.diff(e, x, y)


def jacobian(m):
    """Several outputs sharing many sub-expressions (cse use case)."""
    x, y, z, w = m.symbols("x y z w")
    f = m.exp(-(x - y) ** 2) * m.cos(z * w + x) / (1 + x**2 + y**2 + z**2 + w**2)
    out = []
    for a in (x, y, z, w):
        for b in (x, y, z, w):
            out.append(m.diff(f, a, b))
    return out


def same_expr(builder):
    """Build with symox and convert to SymPy so that both libraries print the
    very same expression."""
    xs = builder(symox)
    if isinstance(xs, list):
        return xs, [sympy.sympify(str(e)) for e in xs]
    return xs, sympy.sympify(str(xs))


BENCHES = []


def bench(name):
    def deco(f):
        BENCHES.append((name, f))
        return f

    return deco


for label, builder in (("poly15", poly), ("nested", nested)):
    ex, es = same_expr(builder)

    BENCHES.append((f"latex({label})", (lambda e=ex: symox.latex(e), lambda e=es: sympy.latex(e))))
    BENCHES.append((f"str({label})", (lambda e=ex: str(e), lambda e=es: str(e))))
    BENCHES.append((f"ccode({label})", (lambda e=ex: symox.ccode(e), lambda e=es: sympy.ccode(e))))
    BENCHES.append((f"pycode({label})", (lambda e=ex: symox.pycode(e), lambda e=es: sympy.pycode(e))))
    BENCHES.append(
        (f"numpy({label})", (lambda e=ex: symox.numpycode(e), lambda e=es: NumPyPrinter().doprint(e)))
    )
    BENCHES.append((f"pretty({label})", (lambda e=ex: symox.pretty(e), lambda e=es: sympy.pretty(e))))
    BENCHES.append((f"cse({label})", (lambda e=ex: symox.cse(e), lambda e=es: sympy.cse(e))))

jx, js = same_expr(jacobian)
BENCHES.append(("cse(jacobian 16 outputs)", (lambda: symox.cse(jx), lambda: sympy.cse(js))))
xs4 = symox.symbols("x y z w")
BENCHES.append(
    (
        "ccode_function(jacobian, cse)",
        (
            lambda: symox.ccode_function("f", xs4, jx, cse=True),
            lambda: [sympy.ccode(v) for v in sympy.cse(js)[1]],
        ),
    )
)
BENCHES.append(
    (
        "cuda_kernel(jacobian, cse)",
        (
            lambda: symox.cuda_kernel("k", xs4, jx, cse=True),
            lambda: [sympy.ccode(v) for v in sympy.cse(js)[1]],
        ),
    )
)


def main():
    print(f"{'benchmark':34s} {'symox':>12s} {'sympy':>12s} {'speedup':>9s}")
    for name, (fx, fs) in BENCHES:
        tx = timeit(fx)
        ts = timeit(fs, min_time=0.5, max_runs=5)
        print(f"{name:34s} {tx * 1e3:10.3f}ms {ts * 1e3:10.2f}ms {ts / tx:8.1f}x")


if __name__ == "__main__":
    main()
