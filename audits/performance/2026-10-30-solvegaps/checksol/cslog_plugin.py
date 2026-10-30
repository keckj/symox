"""pytest plugin: log every single-symbol checksol call made by SymPy 1.14
(expression, symbol, candidate, flags, result) to $CSLOG as JSON lines."""
import json
import os
import sys

import sympy
from sympy import srepr, Expr, Basic
import sympy.solvers.solvers as _ss

_orig = _ss.checksol
_out = open(os.environ['CSLOG'], 'a')
_depth = [0]


def _wrapped(f, symbol, sol=None, **flags):
    top = _depth[0] == 0
    _depth[0] += 1
    try:
        rec = None
        if top:
            try:
                if sol is not None:
                    d = {symbol: sol}
                elif isinstance(symbol, dict):
                    d = dict(symbol)
                else:
                    d = None
                if d is not None and len(d) == 1 and isinstance(f, Expr):
                    (k, v), = d.items()
                    if isinstance(k, sympy.Symbol) and isinstance(v, Basic):
                        rec = {'f': srepr(f), 'x': srepr(k), 'v': srepr(v),
                               'flags': {a: b for a, b in flags.items()
                                         if a in ('numerical', 'force', 'simplify', 'minimal')
                                         and isinstance(b, bool)}}
            except Exception:
                rec = None
        if sol is not None:
            res = _orig(f, symbol, sol, **flags)
        elif isinstance(symbol, dict):
            res = _orig(f, dict(symbol), **flags)
        else:
            res = _orig(f, symbol, sol, **flags)
        if rec is not None and len(rec['f']) + len(rec['v']) < 4000:
            rec['r'] = None if res is None else bool(res)
            rec['test'] = os.environ.get('PYTEST_CURRENT_TEST', '').split(' ')[0]
            _out.write(json.dumps(rec) + '\n')
            _out.flush()
        return res
    finally:
        _depth[0] -= 1


def _patch():
    for name, mod in list(sys.modules.items()):
        if getattr(mod, 'checksol', None) is _orig:
            setattr(mod, 'checksol', _wrapped)


_patch()


def pytest_collection_finish(session):
    _patch()
