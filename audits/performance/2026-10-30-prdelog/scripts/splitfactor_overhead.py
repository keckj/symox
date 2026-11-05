"""Admission overhead of the native splitfactor (risch.splitfactor).

Times, per call, on splitfactor inputs recorded from risch_integrate:
  - risch.splitfactor (admission + kernel, or Python body),
  - _risch_splitfactor_native.reduce (admission + kernel),
  - the kernel alone (fields precomputed),
  - the admission stages: the authentication loops, the lookups loop and
    _prove_fields.
usage: splitfactor_overhead.py [OUT.json]
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from integrands import INTEGRANDS

import symox as S
from symox.integrals import risch
from symox.integrals import _risch_splitfactor_native as sf

x = S.Symbol('x')
root = sf._ROOT
assert root is not None
print('states', len(root['states']), 'lookups', len(root['lookups']), 'scopes', len(root['scopes']))

# Collect (p, DE state) pairs: snapshot the DE by rebuilding it per integrand
# and replaying the level when the call happens.
calls = []
orig = risch.splitfactor


def collect(p, DE, coefficientD=False, z=None):
    if coefficientD is False and z is None and len(calls) < 400:
        calls.append((p, DE, DE.level))
    return orig(p, DE, coefficientD, z)


import signal


def _alarm(*a):
    raise TimeoutError


signal.signal(signal.SIGALRM, _alarm)
risch.splitfactor = collect
for e in INTEGRANDS[:60]:
    signal.alarm(20)
    try:
        risch.risch_integrate(S.sympify(e, locals={'x': x}), x)
    except BaseException:
        pass
    signal.alarm(0)
risch.splitfactor = orig


def per_call(f, reps=200):
    t0 = time.perf_counter()
    for _ in range(reps):
        f()
    return (time.perf_counter() - t0) / reps


rows = []
for p, DE, level in calls[:120]:
    saved = DE.level
    # Restore the recorded level (DecrementLevel scopes have exited).
    while DE.level > level:
        DE.decrement_level()
    while DE.level < level:
        DE.increment_level()
    One = S.Poly(1, DE.t, domain=p.get_domain())
    full = per_call(lambda: risch.splitfactor(p, DE))
    native = sf.reduce(p, DE, One, False, None)
    red = per_call(lambda: sf.reduce(p, DE, One, False, None))
    reader = root['reader'][1]
    builtins = root['builtins']
    states = per_call(lambda: sf._states_match(root['states'], reader, builtins))

    def first_loop():
        for function, code, namespace, defaults, keywords, contents, attributes in root['states']:
            reader(function)
    loop1 = per_call(first_loop)

    def lookups():
        for scope, name, present, value in root['lookups']:
            (name in scope) is not present
            scope.get(name, builtins.get(name)) is not value
    look = per_call(lookups)
    fields = root['prove'](p, DE, sf._metadata, root['verify'], sf._lookup)
    prove = per_call(lambda: root['prove'](p, DE, sf._metadata, root['verify'], sf._lookup))
    if fields is not None:
        kern = per_call(lambda: root['kernel'](p, fields['d'], fields['t'], fields['T'], fields['D'], fields['level'], One))
    else:
        kern = None
    rows.append({'case': DE.case, 'deg': p.degree(DE.t), 'admitted': native is not None,
                 'splitfactor': full, 'reduce': red, 'states_match': states, 'auth_loop': loop1,
                 'lookups': look, 'prove_fields': prove, 'kernel': kern})
    while DE.level > saved:
        DE.decrement_level()
    while DE.level < saved:
        DE.increment_level()

adm = [r for r in rows if r['admitted']]
print('calls', len(rows), 'admitted', len(adm))


def mean(k, rs):
    v = [r[k] for r in rs if r[k] is not None]
    return sum(v) / len(v) * 1e6 if v else float('nan')


for k in ('splitfactor', 'reduce', 'auth_loop', 'states_match', 'lookups', 'prove_fields', 'kernel'):
    print('%-14s mean %8.1f us (admitted)  %8.1f us (all)' % (k, mean(k, adm), mean(k, rows)))
if len(sys.argv) > 1:
    json.dump({'states': len(root['states']), 'lookups': len(root['lookups']), 'rows': rows},
              open(sys.argv[1], 'w'), indent=0)
