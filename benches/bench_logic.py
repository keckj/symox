"""Boolean logic benchmarks: symox vs SymPy.

* ``satisfiable`` on random 3-SAT instances near the phase transition
  (clause/variable ratio 4.26), formulas built as expressions;
* ``simplify_logic`` of random functions of 6-8 variables (built as sums of
  minterms);
* ``to_cnf`` of a DNF and ``solveset``/inequality calls.

Run: python benches/bench_logic.py [--quick]
"""

from __future__ import annotations

import random
import sys
import time

import sympy
import symox
from sympy.logic.inference import satisfiable as sympy_satisfiable


def timeit(f, repeat=1):
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        r = f()
        best = min(best, time.perf_counter() - t0)
    return best, r


def random_3sat(rnd, n, m):
    return [[rnd.choice([1, -1]) * rnd.randint(1, n) for _ in range(3)] for _ in range(m)]


def build(mod, n, clauses):
    vs = mod.symbols(f"v1:{n + 1}")
    lit = lambda l: vs[l - 1] if l > 0 else mod.Not(vs[-l - 1])  # noqa: E731
    return mod.And(*[mod.Or(*[lit(l) for l in c]) for c in clauses])


def check_model(clauses, model, mod, n):
    vs = mod.symbols(f"v1:{n + 1}")
    val = {str(k): bool(v) for k, v in model.items()}
    for c in clauses:
        if not any(val.get(str(vs[abs(l) - 1]), False) == (l > 0) for l in c):
            return False
    return True


def bench_sat(quick):
    rnd = random.Random(42)
    sizes = [(20, 85), (40, 170), (60, 256)] if not quick else [(20, 85), (40, 170)]
    print("satisfiable() on random 3-SAT (ratio ~4.26):")
    for n, m in sizes:
        inst = [random_3sat(rnd, n, m) for _ in range(5 if quick else 10)]
        fs = [build(sympy, n, c) for c in inst]
        fx = [build(symox, n, c) for c in inst]
        ts, rs = timeit(lambda: [sympy_satisfiable(f) for f in fs])
        tx, rx = timeit(lambda: [symox.satisfiable(f) for f in fx], repeat=3)
        agree = all((a is False) == (b is False) for a, b in zip(rs, rx))
        valid = all(b is False or check_model(c, b, symox, n) for c, b in zip(inst, rx))
        nsat = sum(r is not False for r in rx)
        print(f"  n={n:3d} m={m:4d} x{len(inst)}: sympy {ts * 1e3:9.1f} ms  symox {tx * 1e3:8.2f} ms"
              f"  speedup {ts / tx:7.0f}x  (sat {nsat}/{len(inst)}, agree={agree}, models valid={valid})")
    # larger instances, symox only (SymPy takes too long), direct CNF
    for n, m in ([(200, 852), (300, 1278)] if not quick else [(150, 639)]):
        inst = [random_3sat(rnd, n, m) for _ in range(5)]
        fx = [build(symox, n, c) for c in inst]
        tx, rx = timeit(lambda: [symox.satisfiable(f) for f in fx])
        td, _ = timeit(lambda: [symox._core._sat_dimacs(n, c) for c in inst])
        print(f"  n={n:3d} m={m:4d} x5: symox {tx * 1e3:8.1f} ms (expression input), "
              f"{td * 1e3:7.1f} ms (DIMACS input), sat {sum(r is not False for r in rx)}/5")


def random_function(mod, rnd, n, density=0.5):
    vs = mod.symbols(f"a0:{n}")
    terms = []
    for r in range(2 ** n):
        if rnd.random() < density:
            lits = [vs[i] if (r >> (n - 1 - i)) & 1 else mod.Not(vs[i]) for i in range(n)]
            terms.append(mod.And(*lits))
    return mod.Or(*terms)


def bench_simplify(quick):
    print("simplify_logic() of random functions (sum of minterms):")
    for n in ([6, 7, 8] if not quick else [6, 7]):
        rnd1, rnd2 = random.Random(n), random.Random(n)
        k = 3
        fs = [random_function(sympy, rnd1, n) for _ in range(k)]
        fx = [random_function(symox, rnd2, n) for _ in range(k)]
        ts, rs = timeit(lambda: [sympy.simplify_logic(f) for f in fs])
        tx, rx = timeit(lambda: [symox.simplify_logic(f) for f in fx], repeat=3)
        same = sum(str(a) == str(b) for a, b in zip(rs, rx))
        sizes = [(sympy.count_ops(a), sympy.count_ops(sympy.sympify(str(b)))) for a, b in zip(rs, rx)]
        print(f"  n={n} x{k}: sympy {ts * 1e3:9.1f} ms  symox {tx * 1e3:8.2f} ms  speedup {ts / tx:7.0f}x"
              f"  (identical {same}/{k}, ops sympy/symox {sizes})")


def bench_misc():
    print("other:")
    a, b, c, d, e, f = sympy.symbols("a b c d e f")
    A, B, C, D, E, F = symox.symbols("a b c d e f")
    dnf_s = (a & b) | (c & d) | (e & f) | (a & ~c) | (b & ~e)
    dnf_x = (A & B) | (C & D) | (E & F) | (A & ~C) | (B & ~E)
    ts, _ = timeit(lambda: sympy.to_cnf(dnf_s), repeat=3)
    tx, _ = timeit(lambda: symox.to_cnf(dnf_x), repeat=3)
    print(f"  to_cnf(5-term DNF):           sympy {ts * 1e3:8.2f} ms  symox {tx * 1e3:7.3f} ms  speedup {ts / tx:6.0f}x")
    x = sympy.Symbol("x")
    X = symox.Symbol("x")
    from sympy.solvers.inequalities import solve_univariate_inequality as sui
    ts, _ = timeit(lambda: sui((x - 1) * (x - 2) * (x + 3) / (x - 5) >= 0, x), repeat=3)
    tx, _ = timeit(lambda: symox.solve_univariate_inequality((X - 1) * (X - 2) * (X + 3) / (X - 5) >= 0, X), repeat=3)
    print(f"  solve_univariate_inequality:  sympy {ts * 1e3:8.2f} ms  symox {tx * 1e3:7.3f} ms  speedup {ts / tx:6.0f}x")
    ts, _ = timeit(lambda: sympy.solveset(sympy.sin(2 * x) - sympy.Rational(1, 2), x, sympy.S.Reals), repeat=3)
    tx, _ = timeit(lambda: symox.solveset(symox.sin(2 * X) - symox.Rational(1, 2), X, symox.S.Reals), repeat=3)
    print(f"  solveset(sin(2x) - 1/2):      sympy {ts * 1e3:8.2f} ms  symox {tx * 1e3:7.3f} ms  speedup {ts / tx:6.0f}x")
    I1 = [sympy.Interval(i, i + sympy.Rational(3, 2)) for i in range(0, 60, 2)]
    I2 = [symox.Interval(i, i + symox.Rational(3, 2)) for i in range(0, 60, 2)]
    ts, _ = timeit(lambda: sympy.Complement(sympy.Interval(-5, 100), sympy.Union(*I1)), repeat=3)
    tx, _ = timeit(lambda: symox.Complement(symox.Interval(-5, 100), symox.Union(*I2)), repeat=3)
    print(f"  Interval - Union(30 intervals): sympy {ts * 1e3:6.2f} ms  symox {tx * 1e3:7.3f} ms  speedup {ts / tx:6.0f}x")


if __name__ == "__main__":
    quick = "--quick" in sys.argv
    bench_sat(quick)
    bench_simplify(quick)
    bench_misc()
