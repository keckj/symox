"""is_deriv_k(fa, fd, DE) alone: median of 7 runs of 200 calls (us per call).

usage: python deriv_k.py <symox|sympy> [stock]
`stock` makes the symox admission decline (stock Python body, same process).
"""
import timeit, sys, importlib, statistics
S = importlib.import_module(sys.argv[1])
prde = importlib.import_module(sys.argv[1] + '.integrals.prde')
risch = importlib.import_module(sys.argv[1] + '.integrals.risch')
if sys.argv[2:] == ['stock']:
    prde._prde_native.is_deriv_k = lambda fa, fd, DE: None
Poly, symbols = S.Poly, S.symbols
x, t, t1, t2 = symbols('x t t1 t2')


def DE(D, exts, extargs):
    return risch.DifferentialExtension(extension={'D': D, 'exts': exts, 'extargs': extargs})


L2 = DE([Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)], [None, 'log', 'log'], [None, x, x + 1])
LE = DE([Poly(1, x), Poly(1/x, t1), Poly(t2, t2)], [None, 'log', 'exp'], [None, x, x])
LX = DE([Poly(1, x), Poly(1/x, t)], [None, 'log'], [None, x])
X = DE([Poly(1, x)], [None], [None])
cases = {
    'two logs (SymPy sym1)': (Poly(2*x**2 + 2*x, t2), Poly(1, t2), L2),
    'log+exp (SymPy sym2)': (Poly(x**2*t2**3, t2), Poly(1, t2), LE),
    'two logs, rational': (Poly(x**2, t2), Poly((x + 1)**5, t2), L2),
    'log, None': (Poly(x + 1, t), Poly(1, t), LX),
    'log, x**3': (Poly(x**3, t), Poly(1, t), LX),
    'no extension': (Poly(x, x), Poly(1, x), X),
}
for name, (fa, fd, de) in cases.items():
    n = 200
    best = statistics.median(timeit.repeat(lambda: prde.is_deriv_k(fa, fd, de), number=n, repeat=7)) / n * 1e6
    print('%-24s %9.1f us' % (name, best))
