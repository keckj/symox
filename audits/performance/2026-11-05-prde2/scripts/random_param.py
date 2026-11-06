"""Random exp/primitive towers for the parametric Risch DE above the base
level: param_rischDE, limited_integrate, param_poly_rischDE and
prde_cancel_liouvillian (differential comparison of the native kernel).

usage: random_param.py OUT.json [COUNT] [SEED] [TIMEOUT]

The backend is ``$PRDE_BACKEND_MODULE`` (``symox`` or ``sympy``); with
``PRDE_STOCK=1`` the symox native admission is disabled. The cases are
generated from the seed only (strings parsed by each backend). Each row is
``[index, kind, spec, result, seconds, dummies, native]`` (``native``:
native tower answers of the four functions during the call, 0 for SymPy);
a call running longer than TIMEOUT seconds (default 20) is recorded as
``['timeout']``. ``$PRDE_ONLY`` (comma-separated indices) runs only those
draws (the generator still draws every case, so indices are stable).
"""
import json
import os
import random
import signal
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-11-05-prde', 'scripts'))
import codec

S = codec.mod()
from importlib import import_module
prde = import_module(S.__name__ + '.integrals.prde')
risch = import_module(S.__name__ + '.integrals.risch')
if os.environ.get('PRDE_STOCK') == '1':
    import_module(S.__name__ + '.integrals._prde_native')._ROOT = None

OUT = sys.argv[1]
COUNT = int(sys.argv[2]) if len(sys.argv) > 2 else 200
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20261106
TIMEOUT = int(sys.argv[4]) if len(sys.argv) > 4 else 20

TOWERS = [
    ['1', 't0'],               # exp(x)
    ['1', '1/x'],              # log(x)
    ['1', '2*t0'],             # exp(2x)
    ['1', '1/(x + 1)'],        # log(x + 1)
    ['1', 'x*t0'],             # exp(x**2/2)
    ['1', '1/x', 't1'],        # log(x), exp(x)
    ['1', 't0', '1/x'],        # exp(x), log(x)
    ['1', '1/x', '1/(x*t0)'],  # log(x), log(log(x))
    ['1', 't0', 't1/x'],
    ['1', '(1 + t0**2)'],      # tan(x): stock NotImplementedError paths
]
KINDS = ['param_rischDE', 'limited_integrate', 'param_poly_rischDE', 'prde_cancel_liouvillian']


def rpoly(rng, var, deg, coeff):
    return ' + '.join('(%s)*%s**%d' % (coeff(rng), var, i) for i in range(deg + 1))


def make(rng):
    D = rng.choice(TOWERS)
    k = len(D) - 1
    lower = ['x'] + ['t%d' % i for i in range(k - 1)]
    t = 't%d' % (k - 1)

    def coeff(r):
        c = r.choice([0, 1, 1, -1, 2, 3, -2])
        if r.random() < 0.3:
            c = '%d*%s' % (c, r.choice(lower))
        if r.random() < 0.15:
            c = '%s/(%s + %d)' % (c, r.choice(lower), r.choice([0, 1, 2]))
        return str(c)

    kind = rng.choice(KINDS)
    spec = {'D': D, 'kind': kind}
    small = lambda: rpoly(rng, t, rng.choice([0, 0, 1, 2]), coeff)
    if kind in ('param_rischDE', 'limited_integrate'):
        spec['fa'] = small()
        spec['fd'] = rng.choice(['1', 'x', t, '%s + %s' % (t, rng.choice(lower)), small()])
        spec['G'] = [(small(), rng.choice(['1', 'x', t, small()]))
                     for _ in range(rng.choice([1, 1, 2, 3]))]
    elif kind == 'param_poly_rischDE':
        spec['a'] = rng.choice(['1', '1', small()])
        spec['b'] = small()
        spec['q'] = [small() for _ in range(rng.choice([1, 2, 3]))]
        spec['n'] = rng.choice([-1, 0, 1, 2, 3])
    else:
        spec['b'] = rpoly(rng, t, 0, coeff)
        spec['q'] = [small() for _ in range(rng.choice([1, 2]))]
        spec['n'] = rng.choice([0, 1, 2])
    return spec


def build(spec):
    D = spec['D']
    k = len(D) - 1
    T = [S.Symbol('x')] + [S.Symbol('t%d' % i) for i in range(k)]
    ns = {str(s): s for s in T}
    ex = lambda s: S.sympify(s, locals=ns)
    DE = risch.DifferentialExtension(extension={'D': [S.Poly(ex(d), T[i]) for i, d in enumerate(D)], 'T': T})
    t = T[-1]
    P = lambda s: S.Poly(ex(s), t)
    kind = spec['kind']
    if kind in ('param_rischDE', 'limited_integrate'):
        fd = P(spec['fd'])
        G = [(P(a), P(d)) for a, d in spec['G']]
        if fd.is_zero or any(d.is_zero for _, d in G):
            return None
        return getattr(prde, kind), (P(spec['fa']), fd, G, DE)
    if kind == 'param_poly_rischDE':
        a = P(spec['a'])
        if a.is_zero:
            return None
        return prde.param_poly_rischDE, (a, P(spec['b']), [P(q) for q in spec['q']], spec['n'], DE)
    return prde.prde_cancel_liouvillian, (P(spec['b']), [P(q) for q in spec['q']], spec['n'], DE)


def _native():
    if S.__name__ != 'symox':
        return 0
    st = import_module('symox.integrals._prde_native').STATS
    return sum(st[k + '_tower'] for k in KINDS)


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


signal.signal(signal.SIGALRM, _alarm)
rng = random.Random(SEED)
ONLY = {int(i) for i in os.environ.get('PRDE_ONLY', '').split(',') if i}
out = []
for idx in range(COUNT):
    spec = make(rng)
    if ONLY and idx not in ONLY:
        continue
    try:
        built = build(spec)
    except Exception as e:
        out.append([idx, spec['kind'], spec, ['build-error', type(e).__name__], 0.0, 0, 0])
        continue
    if built is None:
        continue
    f, args = built
    d0 = S.Dummy().dummy_index
    n0 = _native()
    t0 = time.perf_counter()
    signal.alarm(TIMEOUT)
    try:
        r = codec.value(S, f(*args))
    except Timeout:
        r = ['timeout']
    except Exception as e:
        r = ['raise', type(e).__name__, str(e)]
    finally:
        signal.alarm(0)
    dt = time.perf_counter() - t0
    out.append([idx, spec['kind'], spec, r, dt, S.Dummy().dummy_index - d0 - 1, _native() - n0])
json.dump(out, open(OUT, 'w'))
print('ran', len(out))
if S.__name__ == 'symox':
    print('native', {k: v for k, v in import_module('symox.integrals._prde_native').STATS.items() if k.endswith('_tower')})
