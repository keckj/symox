"""Matrix benchmarks: symox vs SymPy.

    python benches/bench_matrices.py          # all benchmarks
    python benches/bench_matrices.py det8     # only those matching

Every benchmark also checks that both libraries print the same result.
"""

from __future__ import annotations

import random
import sys
import time

import symox
import sympy


def timeit(fn, min_time=0.3, max_runs=20, slow=1.0):
    t0 = time.perf_counter()
    r = fn()  # warm-up / first run
    first = time.perf_counter() - t0
    if first > slow:
        return first, r
    runs, total, best = 0, 0.0, first
    while total < min_time and runs < max_runs:
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        best = min(best, dt)
        total += dt
        runs += 1
    return best, r


def rand_int_matrix(m, n, seed, lo=-9, hi=9):
    rnd = random.Random(seed)
    return [[rnd.randint(lo, hi) for _ in range(n)] for _ in range(m)]


def bench_det8_symbolic(mod):
    x, y = mod.symbols("x y")
    M = mod.Matrix(8, 8, lambda i, j: (x if (i * 3 + j * 5) % 4 == 0 else y if (i + 2 * j) % 3 == 0 else 1)
                   * ((i * 7 + j * 3) % 5 - 2) + (x * y if i == j else 0))
    return lambda: M.det()


def bench_det5_generic(mod):
    syms = mod.symbols("a0:25")
    M = mod.Matrix(5, 5, list(syms))
    return lambda: M.det()


def bench_det20_int(mod):
    M = mod.Matrix(rand_int_matrix(20, 20, 1))
    return lambda: M.det()


def bench_inv30_rational(mod):
    rnd = random.Random(2)
    rows = [[mod.Rational(rnd.randint(-9, 9), rnd.randint(1, 9)) for _ in range(30)] for _ in range(30)]
    M = mod.Matrix(rows)
    return lambda: M.inv()


def bench_rref50_int(mod):
    M = mod.Matrix(rand_int_matrix(50, 50, 3))
    return lambda: M.rref()


def bench_rref_rankdef(mod):
    # 40 x 60 integer matrix of rank 30
    rnd = random.Random(4)
    A = [[rnd.randint(-5, 5) for _ in range(30)] for _ in range(40)]
    B = [[rnd.randint(-5, 5) for _ in range(60)] for _ in range(30)]
    C = [[sum(A[i][k] * B[k][j] for k in range(30)) for j in range(60)] for i in range(40)]
    M = mod.Matrix(C)
    return lambda: M.rref()


def bench_nullspace(mod):
    M = mod.Matrix(rand_int_matrix(20, 30, 5))
    return lambda: M.nullspace()


def bench_solve20(mod):
    M = mod.Matrix(rand_int_matrix(20, 20, 6))
    b = mod.Matrix([[i - 10] for i in range(20)])
    return lambda: M.solve(b)


def bench_charpoly12(mod):
    M = mod.Matrix(rand_int_matrix(12, 12, 7))
    x = mod.Symbol("x")
    return lambda: M.charpoly(x).as_expr()


def bench_inv4_symbolic(mod):
    x, y = mod.symbols("x y")
    M = mod.Matrix(4, 4, lambda i, j: x**((i + j) % 3) + y * (i - j) + (i + 1) * (j + 2) % 3)
    return lambda: M.inv()


def bench_inv3_generic(mod):
    syms = mod.symbols("a0:9")
    M = mod.Matrix(3, 3, list(syms))
    return lambda: M.inv()


def bench_rref_symbolic(mod):
    x, y = mod.symbols("x y")
    M = mod.Matrix(3, 5, lambda i, j: x**i + y * j + i * j)
    return lambda: M.rref()


def bench_eigen3(mod):
    M = mod.Matrix([[6, -2, 2], [-2, 3, -1], [2, -1, 3]])
    return lambda: M.eigenvects()


def bench_matmul_symbolic(mod):
    x, y = mod.symbols("x y")
    M = mod.Matrix(12, 12, lambda i, j: x**(i % 3) + y * j - i)
    return lambda: M * M


def bench_pow_symbolic(mod):
    x = mod.Symbol("x")
    M = mod.Matrix(4, 4, lambda i, j: x + i - j if (i + j) % 2 else 1)
    return lambda: M**6


def bench_jacobian(mod):
    xs = mod.symbols("x0:10")
    F = mod.Matrix([mod.sin(xs[i]) * xs[(i + 1) % 10] ** 2 + mod.exp(xs[i] * xs[(i + 3) % 10]) for i in range(10)])
    X = mod.Matrix(xs)
    return lambda: F.jacobian(X)


BENCHES = [
    ("det 8x8 symbolic (x, y)", bench_det8_symbolic),
    ("det 5x5 generic symbols", bench_det5_generic),
    ("det 20x20 integer", bench_det20_int),
    ("inv 30x30 rational", bench_inv30_rational),
    ("rref 50x50 integer", bench_rref50_int),
    ("rref 40x60 rank 30", bench_rref_rankdef),
    ("nullspace 20x30 integer", bench_nullspace),
    ("solve 20x20 integer", bench_solve20),
    ("charpoly 12x12 integer", bench_charpoly12),
    ("inv 4x4 polynomial", bench_inv4_symbolic),
    ("inv 3x3 generic symbols", bench_inv3_generic),
    ("rref 3x5 symbolic", bench_rref_symbolic),
    ("eigenvects 3x3", bench_eigen3),
    ("matmul 12x12 symbolic", bench_matmul_symbolic),
    ("power M**6 4x4 symbolic", bench_pow_symbolic),
    ("jacobian 10x10", bench_jacobian),
]


def main(filters):
    print(f"{'benchmark':28s} {'sympy':>11s} {'symox':>11s} {'speedup':>9s}  same")
    for name, b in BENCHES:
        if filters and not any(f in name or f in b.__name__ for f in filters):
            continue
        tx, rx = timeit(b(symox))
        ts, rs = timeit(b(sympy), min_time=0.5, max_runs=5)
        same = str(rx) == str(rs)
        print(f"{name:28s} {ts*1e3:9.2f}ms {tx*1e3:9.3f}ms {ts/tx:8.1f}x  {same}")


if __name__ == "__main__":
    main(sys.argv[1:])
