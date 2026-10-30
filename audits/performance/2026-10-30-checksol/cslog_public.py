"""pytest plugin: log every top-level public checksol call made by SymPy
1.14's own tests (any equation form, any solution mapping, any flags) to
$CSLOG as JSON lines: srepr of the equation(s) and of the mapping, the
flags, the result or the exception, and the warnings."""
import json
import os
import sys
import warnings

import sympy
from sympy import srepr, Basic
import sympy.solvers.solvers as _ss

_orig = _ss.checksol
_out = open(os.environ['CSLOG'], 'a')
_depth = [0]


def _enc(v):
    if isinstance(v, Basic):
        return srepr(v)
    if isinstance(v, (list, tuple)):
        return [_enc(a) for a in v]
    if isinstance(v, (int, bool, float, str)) or v is None:
        return repr(v)
    raise TypeError(type(v))


def _wrapped(f, symbol, sol=None, **flags):
    top = _depth[0] == 0
    _depth[0] += 1
    try:
        rec = None
        if top:
            try:
                if sol is not None:
                    d = [(_enc(symbol), _enc(sol))]
                    form = 'pair'
                elif isinstance(symbol, dict):
                    d = [(_enc(k), _enc(v)) for k, v in symbol.items()]
                    form = 'dict'
                else:
                    d = None
                ff = list(f) if isinstance(f, (list, tuple)) else f
                rec = {'f': _enc(ff), 'flist': isinstance(f, (list, tuple)), 'sol': d, 'form': form,
                       'flags': {a: b for a, b in flags.items() if isinstance(b, bool)}}
            except Exception:
                rec = None
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            try:
                if sol is not None:
                    res = _orig(f, symbol, sol, **flags)
                elif isinstance(symbol, dict):
                    res = _orig(f, dict(symbol), **flags)
                else:
                    res = _orig(f, symbol, sol, **flags)
                err = None
            except Exception as e:
                res, err = None, e
        for x in w:
            warnings.warn_explicit(x.message, x.category, x.filename, x.lineno)
        if rec is not None and len(json.dumps(rec)) < 6000:
            rec['r'] = None if res is None else bool(res)
            rec['err'] = None if err is None else type(err).__name__
            rec['warn'] = len([x for x in w if 'could not verify' in str(x.message)])
            rec['test'] = os.environ.get('PYTEST_CURRENT_TEST', '').split(' ')[0]
            _out.write(json.dumps(rec) + '\n')
            _out.flush()
        if err is not None:
            raise err
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
