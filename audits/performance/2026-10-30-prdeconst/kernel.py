"""constant_system alone on the 20 generated cases of tests/test_native_prdeconst.py
(and SymPy's test_constant_system): native kernel vs the stock Python body
in the same process (the stock body is reached by making the admission
decline), and SymPy 1.14 on the same inputs.

usage: python kernel.py <symox|sympy> <tests dir>
prints: label case-index median-microseconds (of 7 runs of 20 calls)
"""
import sys, time, statistics, importlib
mod = sys.argv[1]
sys.path.insert(0, sys.argv[2])
S = importlib.import_module(mod)
prde = importlib.import_module(mod + '.integrals.prde')
risch = importlib.import_module(mod + '.integrals.risch')
PM = importlib.import_module(mod + '.polys.polymatrix').PolyMatrix
src = open(sys.argv[2] + '/test_native_prdeconst.py').read()
CASES = eval(src[src.index('CASES = [') + 8:src.index(']\n\n\ndef _expr') + 1])
x, t0, t1 = S.symbols('x t0 t1')
ns = {k: getattr(S, k) for k in ('Symbol', 'Integer', 'Rational', 'Add', 'Mul', 'Pow')}


def inputs(c):
    A = PM(c['rows'], c['cols'], [eval(a, ns) for a in c['A']], t1)
    u = PM(c['rows'], 1, [eval(b, ns) for b in c['u']], t1)
    if c['kind'] == 'base':
        DE = risch.DifferentialExtension(extension={'D': [S.Poly(1, x)]})
    else:
        d = {'t0': t0, '1/x': 1/x}[c['D'][1]]
        DE = risch.DifferentialExtension(extension={'D': [S.Poly(1, x), S.Poly(d, t0)]})
    return A, u, DE


def run(label):
    for i, c in enumerate(CASES):
        if c['inA'] == 'QQ' and c['inu'] == 'QQ':
            continue
        A, u, DE = inputs(c)
        prde.constant_system(A, u, DE)
        ts = []
        for _ in range(7):
            t = time.perf_counter()
            for _ in range(20):
                prde.constant_system(A, u, DE)
            ts.append((time.perf_counter() - t) / 20)
        print(label, i, '%.1f' % (statistics.median(ts) * 1e6))


if mod == 'sympy':
    run('sympy')
else:
    nat = importlib.import_module('symox.integrals._prde_native')
    run('native')
    nat.constant_system = lambda A, u, DE: None
    run('stock')
