"""Native kernels against the reference bodies of symox itself.

For random forcing terms and characteristic polynomials, compares
``_undetermined_coefficients_match``, ``_solve_undetermined_coefficients``
and ``dsolve`` with the three kernels enabled and removed from ``_core``.
Prints one JSON line per case and a summary.

usage: PYTHONPATH=python python native_vs_reference.py N SEED
"""
import json
import random
import sys
import time
from contextlib import contextmanager

import symox
from symox import (Derivative, Function, Rational, S, Symbol, cos, dsolve, exp,
                   expand, sin, Add, Mul)
from symox import _core
from symox.solvers.ode.nonhomogeneous import (_undetermined_coefficients_match,
                                              _solve_undetermined_coefficients)

KERNELS = ('ode_uc_trialset', 'ode_uc_solve', 'ode_vop_solve')
x = Symbol('x')
f = Function('f')


@contextmanager
def reference():
    saved = {k: getattr(_core, k) for k in KERNELS}
    try:
        for k in KERNELS:
            delattr(_core, k)
        yield
    finally:
        for k, v in saved.items():
            setattr(_core, k, v)


def rat(rng, lo=-3, hi=3, den=(1, 1, 1, 2)):
    return Rational(rng.randint(lo, hi), rng.choice(den))


def monomial(rng, roots):
    k = rng.choice([0, 0, 1, 1, 2, 3])
    # often resonant with a characteristic root
    if roots and rng.random() < 0.5:
        a, b = rng.choice(roots)
    else:
        a = rng.choice([0, 0, 1, -1, 2, Rational(1, 2)])
        b = rng.choice([0, 0, 1, 2, 3])
    c = rat(rng)
    if c == 0:
        c = S.One
    t = c * x**k
    if a != 0:
        t *= exp(a*x)
    if b != 0:
        t *= rng.choice([sin, cos])(b*x)
    return t


def case(rng):
    # characteristic polynomial from roots (a +- i b) with multiplicities
    roots = []
    poly = S.One
    t = Symbol('t')
    order = 0
    while order < rng.choice([1, 2, 2, 3, 4]):
        if rng.random() < 0.3 and order <= 2:
            a, b = rng.choice([0, 1, -1]), rng.choice([1, 2])
            poly *= (t - a)**2 + b**2
            roots.append((a, b))
            order += 2
        else:
            a = rng.choice([0, 1, -1, 2, -2])
            poly *= (t - a)
            roots.append((a, 0))
            order += 1
    coeffs = symox.Poly(expand(poly), t).all_coeffs()[::-1]
    lead = rng.choice([1, 1, 2, Rational(1, 3)])
    lhs = Add(*[lead*c*Derivative(f(x), (x, i)) if i else lead*c*f(x)
                for i, c in enumerate(coeffs)])
    rhs = Add(*[monomial(rng, roots) for _ in range(rng.choice([1, 1, 2, 3]))])
    return lhs, rhs


def run(fn):
    t0 = time.perf_counter()
    try:
        out = fn()
    except Exception as e:  # noqa
        out = 'EXC %s: %s' % (type(e).__name__, e)
    return out, time.perf_counter() - t0


USED = {k: [0, 0] for k in KERNELS}


def count(name):
    kernel = getattr(_core, name)

    def run(*args):
        out = kernel(*args)
        USED[name][out is None] += 1
        return out
    setattr(_core, name, run)


def main():
    n, seed = int(sys.argv[1]), int(sys.argv[2])
    for k in KERNELS:
        count(k)
    rng = random.Random(seed)
    stats = {'n': 0, 'same': 0, 'diff': 0}
    for i in range(n):
        lhs, rhs = case(rng)
        eq = lhs - rhs
        hom = Add(eq, rhs)
        row = {'i': i, 'eq': str(eq)}
        m_nat, t1 = run(lambda: _undetermined_coefficients_match(-rhs, x, f(x), hom))
        with reference():
            m_ref, t2 = run(lambda: _undetermined_coefficients_match(-rhs, x, f(x), hom))
        row['match_same'] = m_nat == m_ref
        results = {}
        for hint in ('nth_linear_constant_coeff_undetermined_coefficients',
                     'nth_linear_constant_coeff_variation_of_parameters'):
            for simp in (False, True):
                if hint.endswith('parameters') and simp and i % 4:
                    continue
                key = '%s/%s' % (hint.split('_')[-1], simp)
                nat, tn = run(lambda: dsolve(eq, f(x), hint=hint, simplify=simp))
                with reference():
                    ref, tr = run(lambda: dsolve(eq, f(x), hint=hint, simplify=simp))
                same = str(nat) == str(ref) and (not hasattr(nat, 'rhs') or nat == ref)
                results[key] = {'same': same, 't_nat': round(tn, 4), 't_ref': round(tr, 4)}
                if not same:
                    results[key].update(nat=str(nat), ref=str(ref))
        row['hints'] = results
        ok = row['match_same'] and all(v['same'] for v in results.values())
        stats['n'] += 1
        stats['same' if ok else 'diff'] += 1
        print(json.dumps(row), flush=True)
    stats['kernels_native_declined'] = USED
    print(json.dumps({'summary': stats}))


if __name__ == '__main__':
    main()
