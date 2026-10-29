"""Before/after timings for symbolic derivative orders.

usage: PYTHONPATH=<tree>/python python bench.py <label> <mode>
mode "fresh": one call per case in this (fresh) process, cold caches.
mode "repeat": median of 7 rounds of 20 calls each, per case.
Each case prints label, mode, case, seconds, and a short value/exception tag.
"""
import statistics
import sys
import time

t0 = time.perf_counter()
from symox import *  # noqa: E402
IMPORT = time.perf_counter() - t0

x, y = symbols('x y')
n = Symbol('n', integer=True)
k = Symbol('k', integer=True)
f, g = Function('f'), Function('g')

CASES = {
    'diff_x_n': lambda: diff(x, (x, n)),
    'diff_xsin_n': lambda: diff(x*sin(x), (x, n)),
    'diff_fg_n': lambda: diff(f(x)*g(x), (x, n)),
    'held_diff_x': lambda: Derivative(f(x), (x, n)).diff(x),
    'held_subs_3': lambda: Derivative(f(x), (x, n)).subs(n, 3),
    'leibniz_subs_doit': lambda: diff(x*sin(x), (x, n)).subs(n, 3).doit(),
    'rewrite_exp': lambda: Derivative(sin(x), (x, n)).rewrite(exp),
    'sum_doit_derivs': lambda: Sum(Derivative(sin(x), (x, k)), (k, 0, 2)).doit(),
    'eq_parity': lambda: Eq(2*n, 1),
    'diff_literal_5': lambda: diff(sin(x)*exp(x), x, 5),
    'diff_literal_mixed': lambda: diff(x**3*y*sin(y), (x, 2), (y, 1)),
    'diff_sin_n': lambda: diff(sin(x), (x, n)),
}


def tag(call):
    try:
        value = call()
        return type(value).__name__
    except Exception as error:  # noqa: BLE001
        return 'EXC ' + type(error).__name__


def main():
    label, mode = sys.argv[1], sys.argv[2]
    print(f'{label}\t{mode}\timport\t{IMPORT:.6f}\t-')
    for name, call in CASES.items():
        if mode == 'fresh':
            start = time.perf_counter()
            result = tag(call)
            seconds = time.perf_counter() - start
        else:
            result = tag(call)
            rounds = []
            for _ in range(7):
                start = time.perf_counter()
                for _ in range(20):
                    tag(call)
                rounds.append((time.perf_counter() - start) / 20)
            seconds = statistics.median(rounds)
        print(f'{label}\t{mode}\t{name}\t{seconds:.6f}\t{result}')


main()
