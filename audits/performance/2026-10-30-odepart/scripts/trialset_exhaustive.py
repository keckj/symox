"""Native trial sets against the reference body on every monomial
c*x**k*exp(a*x)*T(b*x), k <= 4, for eq_homogeneous S.Zero and for
several homogeneous operators (resonant and not), plus sums of two
monomials. Prints the number of cases, mismatches and native declines.
usage: PYTHONPATH=python python trialset_exhaustive.py
"""
import itertools, json
from symox import Rational, S, Symbol, Function, exp, sin, cos, Derivative, Add, expand, Poly
from symox import _core
from symox.solvers.ode.nonhomogeneous import _undetermined_coefficients_match as ucm

x, t = Symbol('x'), Symbol('t')
f = Function('f')
K = ('ode_uc_trialset',)


def ref(*a):
    saved = _core.ode_uc_trialset
    del _core.ode_uc_trialset
    try:
        return ucm(*a)
    finally:
        _core.ode_uc_trialset = saved


def op(poly):
    c = Poly(expand(poly), t).all_coeffs()[::-1]
    return Add(*[ci*Derivative(f(x), (x, i)) if i else ci*f(x) for i, ci in enumerate(c)])


ops = [op(p) for p in [(t - 1)**2, t**2 + 1, (t**2 + 1)**2, t**3, (t + 2)*(t - Rational(1, 2)),
                       (t - 1)*((t - 1)**2 + 4), t - 3, (t**2 - 2)]]
monos = []
for k, a, b, tr in itertools.product(range(5), [0, 1, -2, Rational(1, 2)], [0, 1, 2, Rational(1, 2)], [sin, cos]):
    if b == 0 and tr is cos:
        continue
    m = Rational(3, 2)*x**k
    if a:
        m *= exp(a*x)
    if b:
        m *= tr(b*x)
    monos.append(m)
n = bad = declined = 0
cases = [(m, None) for m in monos] + [(m, o) for m in monos for o in ops]
cases += [(m1 + m2, o) for m1, m2 in zip(monos[::7], monos[3::7]) for o in (None,) + tuple(ops[:3])]
for expr, o in cases:
    args = (expr, x) if o is None else (expr, x, f(x), o)
    nat = _core.ode_uc_trialset(expr, x, f(x) if o is not None else x, o)
    if nat is None:
        declined += 1
    a, b = ucm(*args), ref(*args)
    n += 1
    if a != b:
        bad += 1
        print('MISMATCH', expr, o, a, b)
print(json.dumps({'cases': n, 'mismatches': bad, 'native_declined': declined}))
