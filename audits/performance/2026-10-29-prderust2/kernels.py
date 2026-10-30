"""Direct calls of the two new native entry points against their Python
bodies (``_ROOT = None``), same process, 7 rounds of 200 calls, median per call.
usage: PYTHONPATH=<worktree>/python python kernels.py"""
import statistics, time
import symox as sx
from symox.integrals import prde, _prde_native as native
from symox.integrals.risch import DifferentialExtension
from symox.polys.polymatrix import PolyMatrix as Matrix
x, t = sx.symbols('x t')
Poly, R = sx.Poly, sx.Rational
DE0 = DifferentialExtension(extension={'D': [Poly(1, x)]})
DE1 = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
G = [(Poly(x/3, x), Poly(2*x + 1, x)), (Poly(5, x), Poly(x**3, x))]
A = Matrix([[1, 2, R(1, 2)], [3, 4, 5], [R(2, 7), 0, 1]], t)
u = Matrix([[1], [R(2, 3)], [0]], t)
cases = [('limited_integrate_reduce base', lambda: prde.limited_integrate_reduce(Poly(R(1, 2)*x + 3, x), Poly(3*x**2 - x, x), G, DE0)),
         ('constant_system primitive 3x3', lambda: prde.constant_system(A, u, DE1))]
root = native._ROOT
for name, f in cases:
    out = {}
    for label, r in (('native', root), ('python', None)):
        native._ROOT = r
        f()
        rounds = []
        for _ in range(7):
            t0 = time.perf_counter()
            for _ in range(200):
                f()
            rounds.append((time.perf_counter() - t0) / 200)
        out[label] = statistics.median(rounds) * 1e6
    native._ROOT = root
    print('%-32s native %8.1f us  python %8.1f us  ratio %.3f' % (name, out['native'], out['python'], out['native'] / out['python']))
