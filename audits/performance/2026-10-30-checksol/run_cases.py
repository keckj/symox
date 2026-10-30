"""Run the checksol cases with sympy or symox and write one JSON record per
case: result, exception type and message, warning messages and the
solution mapping after the call (SymPy simplifies it in place).
usage: run_cases.py sympy|symox cases.json out.json [--logged]"""
import json
import signal
import sys
import time
import warnings

lib = sys.argv[1]
if lib == 'symox':
    import symox as M
    from symox.solvers.solvers import checksol
else:
    import sympy as M
    from sympy.solvers.solvers import checksol
ns = dict(vars(M))
for name in 'x y z a b c p n m'.split():
    ns[name] = M.Symbol(name)
ns['f'] = M.Function('f')
if 'ExprCondPair' not in ns:
    # srepr of a Piecewise names its pairs
    if lib == 'symox':
        from symox.functions.elementary.piecewise import ExprCondPair
    else:
        from sympy.functions.elementary.piecewise import ExprCondPair
    ns['ExprCondPair'] = ExprCondPair
ns['m'] = M.Symbol('m')
try:
    from sympy.physics.units import m as _m  # noqa
except Exception:
    pass


class Timeout(Exception):
    pass


def _alarm(*a):
    raise Timeout()


signal.signal(signal.SIGALRM, _alarm)
cases = json.load(open(sys.argv[2]))
logged = '--logged' in sys.argv
out = []
for c in cases:
    rec = {'id': c['id']}
    try:
        if logged:
            fv = c['f']
            f = [eval(a, ns) for a in fv] if c['flist'] else eval(fv, ns)
            pairs = [(eval(k, ns), eval(v, ns)) for k, v in c['sol']]
            if c['form'] == 'pair':
                args = (f, pairs[0][0], pairs[0][1])
            else:
                args = (f, dict(pairs))
        else:
            f = eval(c['f'], ns)
            if c['sym'] == 'S:':
                d = eval(c['sol'], ns)
                args = (f, d)
            else:
                args = (f, eval(c['sym'], ns), eval(c['sol'], ns))
    except Exception as e:
        rec['setup'] = type(e).__name__
        out.append(rec)
        continue
    signal.alarm(60)
    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        try:
            r = checksol(*args, **c['flags'])
            rec['r'] = None if r is None else bool(r) if r in (True, False) else repr(r)
        except Timeout:
            rec['err'] = 'Timeout'
        except Exception as e:
            rec['err'] = type(e).__name__
            rec['msg'] = str(e)[:200]
        finally:
            signal.alarm(0)
    rec['t'] = time.perf_counter() - t0
    rec['warn'] = [str(x.message) for x in w if 'could not verify' in str(x.message)]
    if len(args) == 2 and isinstance(args[1], dict):
        rec['after'] = str(sorted((str(k), str(v)) for k, v in args[1].items()))
    out.append(rec)
json.dump(out, open(sys.argv[3], 'w'))
print(len(out), 'records')
