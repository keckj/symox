"""Random exp/primitive towers for is_log_deriv_k_t_radical_in_field and
parametric_log_deriv_heu (differential comparison of the tower kernels).

usage: random_towers.py OUT.json [COUNT] [SEED]

The backend is ``$PRDE_BACKEND_MODULE`` (``symox`` or ``sympy``); with
``PRDE_STOCK=1`` the symox native admission is disabled. The cases are
generated from the seed only (strings parsed by each backend), so every
backend sees the same calls. Each output row is ``[index, kind, case,
result, seconds, dummies, native]`` (``native``: number of native
tower answers during the call, 0 for SymPy).
"""
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import codec

S = codec.mod()
from importlib import import_module
prde = import_module(S.__name__ + '.integrals.prde')
risch = import_module(S.__name__ + '.integrals.risch')
if os.environ.get('PRDE_STOCK') == '1':
    import_module(S.__name__ + '.integrals._prde_native')._ROOT = None

OUT = sys.argv[1]
COUNT = int(sys.argv[2]) if len(sys.argv) > 2 else 200
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20261105

# Towers: (D expressions as strings in x, t0, t1, ..., generators).
TOWERS = [
    (['1', 't0'], 'exp(x)'),
    (['1', '1/x'], 'log(x)'),
    (['1', '2*t0'], 'exp(2x)'),
    (['1', '1/(x + 1)'], 'log(x+1)'),
    (['1', 't0', 't1/x'], 'exp(x), x**? exp(log x)'),
    (['1', '1/x', 't1'], 'log(x), exp(x)'),
    (['1', 't0', '1/x'], 'exp(x), log(x)'),
    (['1', '1/x', '1/(x*t0)'], 'log(x), log(log(x))'),
    (['1', 't0', '2*x*t1'], 'exp(x), exp(x**2)'),
    (['1', '(1 + t0**2)'], 'tan(x)'),
    (['1', 'x*t0'], 'exp(x**2/2)'),
]


def rpoly(rng, var, deg, coeff):
    return ' + '.join('(%s)*%s**%d' % (coeff(rng), var, i) for i in range(deg + 1))


def make(rng):
    D, label = rng.choice(TOWERS)
    k = len(D) - 1
    lower = ['x'] + ['t%d' % i for i in range(k - 1)]
    t = 't%d' % (k - 1)

    def coeff(r):
        c = r.choice([0, 0, 1, -1, 2, 3, -2])
        if r.random() < 0.3:
            c = '%d*%s' % (c, r.choice(lower))
        if r.random() < 0.15:
            c = '%s/(%s + %d)' % (c, r.choice(lower), r.choice([1, 2]))
        return str(c)

    kind = rng.choice(['log', 'log', 'log', 'heu'])
    mode = rng.random()
    spec = {'D': D, 'label': label, 'kind': kind}
    if kind == 'log':
        if mode < 0.55:
            # f = sum(r_i*Du_i/u_i) (+ n*Dt/t) (+ noise): a log-derivative
            # of a radical by construction, perturbed sometimes.
            us = [rpoly(rng, rng.choice([t] + lower), rng.choice([1, 1, 2]), coeff) for _ in range(rng.choice([1, 1, 2]))]
            rs = [rng.choice(['1', '2', '1/2', '-1', '3/2', '-1/3']) for _ in us]
            spec['logs'] = list(zip(us, rs))
            spec['noise'] = coeff(rng) if rng.random() < 0.3 else '0'
        else:
            spec['fa'] = rpoly(rng, t, rng.choice([0, 1, 2]), coeff)
            spec['fd'] = rpoly(rng, t, rng.choice([0, 1, 2]), coeff)
        spec['case'] = rng.choice(['auto'] * 8 + ['exp', 'primitive', 'base', 'tan', 'other_linear', 'bogus'])
    else:
        spec['fa'] = rpoly(rng, t, rng.choice([0, 1, 2]), coeff)
        spec['fd'] = rpoly(rng, t, rng.choice([0, 1, 2]), coeff)
        spec['wa'] = rpoly(rng, t, rng.choice([0, 1, 2]), coeff)
        spec['wd'] = rpoly(rng, t, rng.choice([0, 1]), coeff)
    return spec


def build(spec):
    D = spec['D']
    k = len(D) - 1
    syms = [S.Symbol('x')] + [S.Symbol('t%d' % i) for i in range(k)]
    ns = {str(s): s for s in syms}
    ex = lambda s: S.sympify(s, locals=ns)
    T = syms
    Dp = [S.Poly(ex(d), T[i]) for i, d in enumerate(D)]
    DE = risch.DifferentialExtension(extension={'D': Dp, 'T': T})
    t = T[-1]

    def deriv(e):
        return sum(S.diff(e, T[i])*ex(D[i]) for i in range(k + 1))

    if spec['kind'] == 'log':
        if 'logs' in spec:
            f = ex(spec['noise'])
            for u, r in spec['logs']:
                u = ex(u)
                if u == 0:
                    continue
                f += ex(r)*deriv(u)/u
            f = S.cancel(f)
            n, d = f.as_numer_denom()
            fa, fd = S.Poly(n, t), S.Poly(d, t)
        else:
            fa, fd = S.Poly(ex(spec['fa']), t), S.Poly(ex(spec['fd']), t)
        if fd.is_zero:
            return None
        return DE, ('is_log_deriv_k_t_radical_in_field', (fa, fd, DE, spec['case']))
    fa, fd, wa, wd = [S.Poly(ex(spec[n]), t) for n in ('fa', 'fd', 'wa', 'wd')]
    if fd.is_zero or wd.is_zero:
        return None
    return DE, ('parametric_log_deriv_heu', (fa, fd, wa, wd, DE))


def _native():
    # Native tower answers so far (symox only; 0 for SymPy).
    if S.__name__ != 'symox':
        return 0
    st = import_module('symox.integrals._prde_native').STATS
    return st['is_log_deriv_k_t_radical_in_field_tower'] + st['parametric_log_deriv_heu_tower']


rng = random.Random(SEED)
out = []
for idx in range(COUNT):
    spec = make(rng)
    try:
        built = build(spec)
    except Exception as e:
        out.append([idx, spec.get('kind'), spec, ['build-error', type(e).__name__], 0.0, 0])
        continue
    if built is None:
        continue
    DE, (name, args) = built
    d0 = S.Dummy().dummy_index
    n0 = _native()
    t0 = time.perf_counter()
    try:
        r = codec.value(S, getattr(prde, name)(*args))
    except Exception as e:
        r = ['raise', type(e).__name__, str(e)]
    dt = time.perf_counter() - t0
    out.append([idx, name, spec, r, dt, S.Dummy().dummy_index - d0 - 1, _native() - n0])
json.dump(out, open(OUT, 'w'))
print('ran', len(out))
if S.__name__ == 'symox':
    print('native', import_module('symox.integrals._prde_native').STATS)
