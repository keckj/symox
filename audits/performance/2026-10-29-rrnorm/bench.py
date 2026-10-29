"""Before/after timings of the public calls touched by the rrnorm task."""
import sys, time
from symox import *
x, y, z = symbols('x y z')
p = Symbol('p', positive=True)
cases = [
    ('Poly.gcd RR[y]', lambda: Poly(1.5*x*z + 1.5*y, x, z).gcd(Poly(3.0*x*z + 3.0*y, x, z))),
    ('Poly.lcm RR[y]', lambda: Poly(0.5*x**2*z - 0.5*y**2, x, z).lcm(Poly(1.5*x*z + 1.5*y, x, z))),
    ('Poly.cofactors RR[y]', lambda: Poly(1.5*x*z + 1.5*y, x, z).cofactors(Poly(3.0*x*z + 3.0*y, x, z))),
    ('gcd expr RR[y]', lambda: gcd(1.5*x*z + 1.5*y, 3.0*x*z + 3.0*y, x, z)),
    ('Poly.sqf_list RR[y]', lambda: Poly(1.5*(x*z + y)**2*(x + 1), x, z).sqf_list()),
    ('Poly.integrate RR[y]', lambda: Poly(1.5*x*z + 1.5*y, x, z).integrate(x)),
    ('polylog sign', lambda: polylog(3, -Rational(1, 2)).is_extended_negative),
    ('hyperexpand lerchphi', lambda: hyperexpand(hyper([S.Half]*3 + [1], [Rational(3, 2)]*3, Rational(1, 4)))),
    ('solve exp tower +1', lambda: solve(exp(exp(exp(x))) + 1, x)),
    ('solve exp(exp(x)+1)-4', lambda: solve(exp(exp(x) + 1) - 4, x)),
    ('simplify log', lambda: simplify(log(-1 + 2*log(2)))),
    ('integrate atan(p*x)/x**2', lambda: integrate(atan(p*x)/x**2, (x, 1, oo))),
]
for name, f in cases:
    best = None
    for _ in range(3):
        t = time.perf_counter()
        try:
            r = f()
        except Exception as e:
            r = '%s: %s' % (type(e).__name__, e)
        dt = time.perf_counter() - t
        best = dt if best is None else min(best, dt)
    print('%-28s %9.4fs  %s' % (name, best, str(r)[:90]), flush=True)
