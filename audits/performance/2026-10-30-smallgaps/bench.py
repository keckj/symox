"""Timings of the smallgaps cases: median of repeated calls (seconds).

Usage: bench.py <module> <det_cases.py> <out.json>
For the baseline .so, run with a PYTHONPATH whose symox/_core.abi3.so is
the main build.
"""
import importlib
import json
import statistics
import sys
import time

modname, det_path, out = sys.argv[1:4]
m = importlib.import_module(modname)
cs = importlib.import_module(modname + '.solvers.solvers').checksol
units = importlib.import_module(modname + '.physics.units')
ns = {'__name__': 'det_cases_import'}
src = open(det_path).read().split('\nimport os\n')[0]
sys.argv = [det_path, modname, '/dev/null']
exec(compile(src, det_path, 'exec'), ns)
DET = dict(ns['cases']())


def med(f, n):
    ts = []
    for _ in range(n):
        t = time.perf_counter(); f(); ts.append(time.perf_counter() - t)
    return statistics.median(ts)


x, y = m.symbols('x y')
I, pi, log = m.I, m.pi, m.log
res = {}
for name in ('third4', 'n4_k0', 'n5_k0', 'n6_k0', 'n7_k0', 'n4_k3', 'n5_k3', 'n6_k3', 'prec50'):
    rows = DET[name]
    for meth in (None, 'bareiss'):
        def run():
            M = m.Matrix(rows)  # fresh matrix: no cached det
            return M.det() if meth is None else M.det(method=meth)
        res['det %s %s' % (name, meth)] = med(run, 5 if name.startswith(('n6', 'n7')) else 20)
r = m.erfinv(m.I**(1/m.pi))**3
res['checksol 17799'] = med(lambda: cs(-m.erf(x**(m.S(1)/3))**m.pi + m.I, x, r), 5)
res['solve 17799'] = med(lambda: m.solve(-m.erf(x**(m.S(1)/3))**m.pi + m.I, x), 3)
res['checksol units float'] = med(lambda: cs(units.cm*(x - 2), x, m.Float('2.0000000001')), 20)
res['checksol units 2cmx'] = med(lambda: cs(2*units.cm*x, x, 1), 20)
res['pow E branch 2I'] = med(lambda: (2*I)**(x/log(2*I)), 200)
res['pow E branch none'] = med(lambda: (x + 1)**(y/(x + 2)), 200)
res['pow plain'] = med(lambda: (x + 1)**(2*y + 1), 200)
json.dump(res, open(out, 'w'), indent=1, sort_keys=True)
for k, v in sorted(res.items()):
    print('%-28s %.6f' % (k, v))
