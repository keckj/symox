"""Differential corpus for the ratint2 lot: RootSum constructions, doit and
evalf, ratint/integrate and ratint_logpart with algebraic and parametric
coefficients (QQ<a>[p], QQ<a>(t)), apart(full=True), EX elements.
Usage: python diffcorpus.py sympy|symox OUT.jsonl

Each case is run once for warm-up, then timed as the best of three runs
(capped at about 10 s of total work per case); a case without result
within 60 s is recorded as CaseTimeout. Extra arguments select cases."""
import json, sys, time

mode, out = sys.argv[1], sys.argv[2]
if mode == 'sympy':
    import sympy as M
    from sympy.integrals.rationaltools import ratint_logpart, ratint
else:
    import symox as M
    from symox.integrals.rationaltools import ratint_logpart, ratint
names = ('symbols', 'sqrt', 'Poly', 'QQ', 'S', 'integrate', 'I', 'Rational', 'RootSum',
         'Lambda', 'log', 'exp', 'sin', 'apart', 'EX', 'Function', 'PurePoly', 'cbrt', 'diff')
(symbols, sqrt, Poly, QQ, S, integrate, I, Rational, RootSum, Lambda, log, exp, sin, apart,
 EX, Function, PurePoly, cbrt, diff) = (getattr(M, n) for n in names)

x, y, z, a, b, t = symbols('x y z a b t')
p = symbols('p', positive=True)
q = symbols('q')
f = Function('f')
s2, s3 = sqrt(2), sqrt(3)
K = QQ.algebraic_field(s2)


def rs(*args, **kw):
    return lambda: RootSum(*args, **kw)


def lp(num, den, dom):
    return lambda: ratint_logpart(Poly(num, x, domain=dom), Poly(den, x, domain=dom), x)


def kp(num, den):
    def run():
        dom = QQ.algebraic_field(s2)[p]
        return ratint_logpart(Poly(num, x, domain=dom), Poly(den, x, domain=dom), x)
    return run


def kt(num, den):
    def run():
        dom = QQ.algebraic_field(s2).frac_field(t)
        return ratint_logpart(Poly(num, x, domain=dom), Poly(den, x, domain=dom), x)
    return run


cases = [
    # RootSum constructor
    ('RootSum cubic id', rs(x**3 + x + 1)),
    ('RootSum cubic Lambda 1/(x+1)', rs(x**3 + x + 1, Lambda(x, 1/(x + 1)))),
    ('RootSum cubic Lambda x**2', rs(x**3 + x + 1, Lambda(x, x**2))),
    ('RootSum cubic Lambda x/(x**2+2)', rs(x**3 + x + 1, Lambda(x, x/(x**2 + 2)))),
    ('RootSum cubic log', rs(x**3 + x + 1, Lambda(x, log(x)))),
    ('RootSum cubic y*log(z-y)', rs(x**3 + x + 1, Lambda(y, y*log(z - y)))),
    ('RootSum cubic 2+3*log', rs(x**3 + x + 1, Lambda(y, 2 + 3*y*log(z - y)))),
    ('RootSum cubic exp', rs(x**3 + x + 1, Lambda(x, exp(x)))),
    ('RootSum cubic Function f', rs(x**3 + x + 1, f)),
    ('RootSum constant Lambda', rs(x**3 + x + 1, Lambda(x, 7))),
    ('RootSum quartic a param 1/(x-a)', rs(x**4 + a*x + 1, Lambda(x, 1/(x - a)))),
    ('RootSum quintic sum x**3', rs(x**5 - x + 1, Lambda(x, x**3))),
    ('RootSum quintic 1/x', rs(x**5 - x + 1, Lambda(x, 1/x))),
    ('RootSum product quad*cubic log', rs((x**2 + 1)*(x**3 + x + 1), Lambda(x, log(x)))),
    ('RootSum quadratic=True log', rs(x**2 + 3*x + 1, Lambda(x, log(x)), quadratic=True)),
    ('RootSum quadratic=False log', rs(x**2 + 3*x + 1, Lambda(x, log(x)))),
    ('RootSum linear factor log', rs((x - 2)*(x**3 + x + 1), Lambda(x, log(x)))),
    ('RootSum repeated factor log', rs((x**3 + x + 1)**2, Lambda(x, log(x)))),
    ('RootSum auto=False', rs(x**3 + x + 1, Lambda(x, 1/(x + 1)), auto=False)),
    ('RootSum scaled 4*x**3+2', rs(4*x**3 + 2*x + 3, Lambda(x, log(x)))),
    ('RootSum integer basis 8', rs(x**3 + 8*x + 64, Lambda(x, log(x)))),
    ('RootSum integer basis rational', rs(x**3 + 8*x + 64, Lambda(x, 1/(x + 1)))),
    ('RootSum symbolic coeff a', rs(x**3 + a*x + a**3, Lambda(x, log(x)))),
    ('RootSum symbolic coeff rational', rs(x**3 + a*x + a**3, Lambda(x, x**2))),
    ('RootSum param p', rs(x**3 + p*x + 1, Lambda(x, log(x - p)))),
    ('RootSum param rational p', rs(x**3 + p*x + 1, Lambda(x, 1/(x - p)))),
    ('RootSum sqrt2 coeff log', rs(x**3 + s2*x + 1, Lambda(x, log(x)))),
    ('RootSum sqrt2 coeff rational', rs(x**3 + s2*x + 1, Lambda(x, x**2))),
    ('RootSum x gen keyword', lambda: RootSum(y**3 + y + 1, Lambda(x, log(x)), y)),
    ('RootSum fraction coeff', rs(x**3 + x/2 + Rational(1, 3), Lambda(x, log(x)))),
    ('RootSum quadratic=True param', rs(x**2 + a*x + 1, Lambda(x, log(x)), quadratic=True)),
    ('RootSum multivariate error', rs(x**3 + y*x + 1 + z*y)),
    ('RootSum bad func error', rs(x**3 + x + 1, Lambda((x, y), x + y))),
    ('RootSum args', lambda: RootSum(x**3 + x + 1, Lambda(x, log(x))).args),
    ('RootSum hash eq', lambda: RootSum(x**3 + x + 1, Lambda(x, log(x))) == RootSum(y**3 + y + 1, Lambda(y, log(y)))),
    ('RootSum doit cubic', lambda: RootSum(x**3 + x + 1, Lambda(x, log(x))).doit()),
    ('RootSum doit quartic', lambda: RootSum(x**4 + 1, Lambda(x, log(x))).doit()),
    ('RootSum doit quintic stays', lambda: RootSum(x**5 - x + 1, Lambda(x, log(x))).doit()),
    ('RootSum evalf cubic', lambda: RootSum(x**3 + x + 1, Lambda(x, log(x - 1))).evalf()),
    ('RootSum evalf quintic', lambda: RootSum(x**5 - x + 1, Lambda(x, log(x + 2))).evalf(30)),
    ('RootSum evalf symbolic', lambda: RootSum(x**3 + a*x + 1, Lambda(x, log(x))).evalf()),
    ('RootSum diff', lambda: diff(RootSum(x**3 + x + 1, Lambda(y, y*log(z - y))), z)),
    ('RootSum subs', lambda: RootSum(x**3 + x + 1, Lambda(y, log(z - y))).subs(z, 2)),
    # ratint / integrate
    ('integrate 1/(x**3+x+1) real=False', lambda: ratint(1/(x**3 + x + 1), x, real=False)),
    ('integrate 1/(x**5+x+1)', lambda: integrate(1/(x**5 - x + 1), x)),
    ('integrate x**2/(x**5-x+1)', lambda: integrate(x**2/(x**5 - x + 1), x)),
    ('integrate 1/(x**3+p*x+1)', lambda: integrate(1/(x**3 + p*x + 1), x)),
    ('integrate 1/(x**2+sqrt(2)*x+p)', lambda: integrate(1/(x**2 + s2*x + p), x)),
    ('integrate 1/(x**3+sqrt(2)*x+1)', lambda: integrate(1/(x**3 + s2*x + 1), x)),
    ('integrate 1/(x**2+sqrt(2)*p)', lambda: integrate(1/(x**2 + s2*p), x)),
    ('integrate x/(x**4+a)', lambda: integrate(x/(x**4 + a), x)),
    ('integrate 1/(x**4+a*x+1)', lambda: integrate(1/(x**4 + a*x + 1), x)),
    ('ratint 1/(x**3+a) real=False', lambda: ratint(1/(x**3 + a), x, real=False)),
    ('logpart extension=True param', lambda: ratint_logpart(Poly(x, x, extension=True), Poly(x**2 + s2*x + p, x, extension=True), x)),
    ('logpart QQ<s2> 1/(x**3+s2*x+1)', lp(1, x**3 + s2*x + 1, K)),
    ('logpart QQ<s2>[p] 1/(x**2+s2*p)', kp(1, x**2 + s2*p)),
    ('logpart QQ<s2>[p] x/(x**2+s2*x+p)', kp(x, x**2 + s2*x + p)),
    ('logpart QQ<s2>[p] 1/(x**2+p)', kp(1, x**2 + p)),
    ('logpart QQ<s2>[p] 1/(x**3+s2*p)', kp(1, x**3 + s2*p)),
    ('logpart QQ<s2>[p] 1/((x-p)(x+s2))', kp(1, (x - p)*(x + s2))),
    ('logpart QQ<s2>(t) 1/(x**2+s2*t)', kt(1, x**2 + s2*t)),
    ('logpart QQ<s2>(t) x/(x**2+s2*x+t)', kt(x, x**2 + s2*x + t)),
    ('logpart QQ<s2>(t) 1/(x**3+t*x+s2)', kt(1, x**3 + t*x + s2)),
    ('logpart QQ<s2>(t) 1/((x-t)(x+s2))', kt(1, (x - t)*(x + s2))),
    ('logpart QQ<s2>(t) (x+1)/(x**2+t*x+s2)', kt(x + 1, x**2 + t*x + s2)),
    # apart(full=True)
    ('apart full 1/(x**3+x+1)', lambda: apart(1/(x**3 + x + 1), x, full=True)),
    ('apart full doit 1/(x**3+x+1)', lambda: apart(1/(x**3 + x + 1), x, full=True).doit()),
    ('apart full x/(x**4+1)', lambda: apart(x/(x**4 + 1), x, full=True)),
    ('apart full 1/(x**2+1)**2', lambda: apart(1/(x**2 + 1)**2, x, full=True)),
    ('apart full 1/((x-1)(x**3+2))', lambda: apart(1/((x - 1)*(x**3 + 2)), x, full=True)),
    ('apart full 1/(x**5-x+1)', lambda: apart(1/(x**5 - x + 1), x, full=True)),
    ('apart full a/(x**3+a)', lambda: apart(a/(x**3 + a), x, full=True)),
    ('apart full 1/(x**2+sqrt(2))', lambda: apart(1/(x**2 + s2), x, full=True)),
    # EX elements
    ('EX sqrt2*(1+sqrt2)', lambda: EX(s2*(1 + s2))),
    ('EX (x+1)**2', lambda: EX((x + 1)**2)),
]

import signal


class CaseTimeout(Exception):
    pass


def _alarm(*_):
    raise CaseTimeout('no result within %d s' % LIMIT)


LIMIT = 60
signal.signal(signal.SIGALRM, _alarm)
only = sys.argv[3:]

with open(out, 'w') as fh:
    for name, fn in cases:
        if only and name not in only:
            continue
        signal.alarm(LIMIT)
        try:
            r = fn(); text = str(r); err = None
        except Exception as e:
            text = None; err = f'{type(e).__name__}: {str(e)[:200]}'
        finally:
            signal.alarm(0)
        best = None
        if err is None:
            total = 0.0
            for _ in range(3):
                t0 = time.perf_counter(); fn(); dt = time.perf_counter() - t0
                best = dt if best is None else min(best, dt); total += dt
                if total > 10: break
        fh.write(json.dumps({'case': name, 'result': text, 'error': err, 'seconds': best}) + '\n')
        fh.flush()
        print(name, '|', err or f'{len(text)} chars', '|', best, flush=True)
