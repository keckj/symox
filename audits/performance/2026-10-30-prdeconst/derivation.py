"""derivation(p, DE, coefficientD) alone: best of 5 runs of 2000 calls (us per call).

usage: python derivation.py <symox|sympy>
"""
import timeit, sys, importlib
S = importlib.import_module(sys.argv[1]); risch = importlib.import_module(sys.argv[1] + '.integrals.risch')
Poly, symbols = S.Poly, S.symbols
x, t, t0 = symbols('x t t0')
cases = {
 'tower ZZ(x,t0)': (risch.DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t0), Poly(t, t)]}), lambda: Poly(x*t**2 + t0*t + 1/x, t, domain='ZZ(x,t0)')),
 'exp ZZ[x]': (risch.DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]}), lambda: Poly(x*t**2 + t + 1, t)),
 'prim ZZ(x)': (risch.DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]}), lambda: Poly(x*t**2 + t/x, t, domain='ZZ(x)')),
}
for name, (DE, mk) in cases.items():
    p = mk()
    for cd in (False, True):
        n = 2000
        best = min(timeit.repeat(lambda: risch.derivation(p, DE, coefficientD=cd), number=n, repeat=5)) / n * 1e6
        print('%-16s coefficientD=%-5s %8.1f us' % (name, cd, best))
