"""Replay recorded prde calls in one backend and store the encoded results.

usage: replay.py IN.json OUT.json [KIND...]

The backend is ``$PRDE_BACKEND_MODULE`` (``symox`` or ``sympy``). With
``PRDE_STOCK=1`` the symox native admission is disabled (stock bodies).
Every result is encoded with ``codec.value``; an exception is stored as
``['raise', type name, message]``. Each output row is ``[index, kind,
result, seconds, dummies]`` (the number of Dummies the call created).
"""
import gzip
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import codec

S = codec.mod()
from importlib import import_module
prde = import_module(S.__name__ + '.integrals.prde')
if os.environ.get('PRDE_STOCK') == '1':
    nat = import_module(S.__name__ + '.integrals._prde_native')
    nat._ROOT = None
DEPOS = {'is_log_deriv_k_t_radical_in_field': 2, 'parametric_log_deriv_heu': 4,
         'parametric_log_deriv': 4, 'limited_integrate': 3, 'limited_integrate_reduce': 3,
         'param_rischDE': 3, 'param_poly_rischDE': 4, 'prde_cancel_liouvillian': 3}

data = json.load(gzip.open(sys.argv[1], 'rt') if sys.argv[1].endswith('.gz') else open(sys.argv[1]))
kinds = set(sys.argv[3:])
out = []
for idx, rec in enumerate(data['records']):
    name, args, kw, de = rec[0], rec[1], rec[2], rec[3]
    if args is None or (kinds and name not in kinds):
        continue
    dec = codec.Decoder(S)
    try:
        DE = dec.de(de)
        a = [dec.value(v) for v in args]
        a.insert(DEPOS[name], DE)
        k = {kk: dec.value(vv) for kk, vv in kw.items()}
    except Exception as e:
        out.append([idx, name, ['decode-error', repr(e)], 0.0])
        continue
    d0 = S.Dummy().dummy_index
    t0 = time.perf_counter()
    try:
        r = codec.value(S, getattr(prde, name)(*a, **k))
    except Exception as e:
        r = ['raise', type(e).__name__, str(e)]
    dt = time.perf_counter() - t0
    # Number of Dummies the call created (the stock bodies create them in a
    # fixed order; a native kernel must create as many).
    dummies = S.Dummy().dummy_index - d0 - 1
    out.append([idx, name, r, dt, dummies])
json.dump(out, open(sys.argv[2], 'w'))
print('replayed', len(out))
if S.__name__ == 'symox':
    print('native', import_module('symox.integrals._prde_native').STATS)
