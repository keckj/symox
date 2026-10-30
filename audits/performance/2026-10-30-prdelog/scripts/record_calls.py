"""Record the base-case calls of is_log_deriv_k_t_radical_in_field,
parametric_log_deriv_heu and limited_integrate made by risch_integrate on
the integrand corpus (symox, Python bodies), plus random base-case inputs.

usage: record_calls.py OUT.json
Each record is [kind, args] with every Poly given as [expr string, domain].
"""
import json
import os
import random
import signal
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from integrands import INTEGRANDS

import symox as S
from symox.integrals import prde, risch

x = S.Symbol('x')
records = []
seen = set()


def enc(p):
    return [str(p.as_expr()), str(p.get_domain())]


def add(kind, args):
    key = json.dumps([kind, args])
    if key not in seen:
        seen.add(key)
        records.append([kind, args])


def base(DE):
    return getattr(DE, 'case', None) == 'base' and DE.t == x


def wrap(name, encode):
    f = getattr(prde, name)

    def w(*a, **k):
        try:
            if encode(*a, **k) is not None:
                pass
        except Exception:
            pass
        return f(*a, **k)
    setattr(prde, name, w)
    if hasattr(risch, name):
        setattr(risch, name, w)


def e_log(fa, fd, DE, case='auto', z=None):
    if base(DE) and z is None and case in ('auto', 'base'):
        add('is_log_deriv', [enc(fa), enc(fd)])


def e_heu(fa, fd, wa, wd, DE, c1=None):
    if base(DE) and c1 is None:
        add('heu', [enc(fa), enc(fd), enc(wa), enc(wd)])


def e_lim(fa, fd, G, DE):
    if base(DE):
        add('limited', [enc(fa), enc(fd), [[enc(a), enc(d)] for a, d in G]])


wrap('is_log_deriv_k_t_radical_in_field', e_log)
wrap('parametric_log_deriv_heu', e_heu)
wrap('limited_integrate', e_lim)


def _alarm(*a):
    raise TimeoutError


signal.signal(signal.SIGALRM, _alarm)
for e in INTEGRANDS:
    signal.alarm(30)
    try:
        risch.risch_integrate(S.sympify(e, locals={'x': x}), x)
    except BaseException:
        pass
    signal.alarm(0)
n_corpus = len(records)

# Random base-case inputs: generic fractions and logarithmic derivatives.
rng = random.Random(20261030)


def rpoly(deg, lo=-4, hi=4):
    while True:
        c = [S.Rational(rng.randint(lo, hi), rng.choice([1, 1, 1, 2, 3])) for _ in range(deg + 1)]
        if c[0] != 0:
            return S.Poly(c, x)


def poly_s(p):
    return [str(p.as_expr()), str(p.get_domain())]


for _ in range(150):
    # f = sum(r_i * h_i'/h_i) (+ noise), the shape that succeeds.
    k = rng.randint(1, 3)
    num, den = S.Poly(0, x), S.Poly(1, x)
    for _ in range(k):
        h = rpoly(rng.randint(1, 2))
        r = S.Rational(rng.randint(-5, 5), rng.choice([1, 2, 3]))
        num, den = num * h + S.Poly(r, x) * h.diff(x) * den, den * h
    if rng.random() < 0.3:
        num = num + rpoly(rng.randint(0, 2)) * den * (1 if rng.random() < 0.5 else 0) + rpoly(0)
    if den.is_zero:
        continue
    add('is_log_deriv', [poly_s(num), poly_s(den)])
    add('is_log_deriv', [poly_s(rpoly(rng.randint(0, 3))), poly_s(rpoly(rng.randint(1, 3)))])
    wa, wd = (rpoly(rng.randint(0, 2)), rpoly(rng.randint(0, 2)))
    add('heu', [poly_s(num), poly_s(den), poly_s(wa), poly_s(wd)])
    m = rng.randint(1, 3)
    add('heu', [poly_s(num * wd + S.Poly(m, x) * wa * den), poly_s(den * wd), poly_s(wa), poly_s(wd)])
    G = [[poly_s(rpoly(rng.randint(0, 2))), poly_s(rpoly(rng.randint(1, 2)))] for _ in range(rng.randint(0, 2))]
    add('limited', [poly_s(num), poly_s(den), G])
    add('limited', [poly_s(rpoly(rng.randint(0, 3))), poly_s(rpoly(rng.randint(0, 3))), G])

with open(sys.argv[1], 'w') as fh:
    json.dump({'corpus': n_corpus, 'records': records}, fh, indent=0)
print('recorded', n_corpus, 'corpus calls and', len(records) - n_corpus, 'random calls')
