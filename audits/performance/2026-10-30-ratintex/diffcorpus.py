"""Differential corpus for the ratintex lot: ratint over QQ<sqrt(2)>,
log_to_atan over EX, real=False RootSum, roots of decomposable and cubic
polynomials.  Usage: python diffcorpus.py sympy|symox OUT.jsonl

Each case is run once for warm-up, then timed as the best of three runs
(capped at about 10 s of total work per case)."""
import json, sys, time

mode, out = sys.argv[1], sys.argv[2]
if mode == 'sympy':
    import sympy as M
    from sympy.integrals.rationaltools import ratint_logpart, ratint, log_to_atan
else:
    import symox as M
    from symox.integrals.rationaltools import ratint_logpart, ratint, log_to_atan
from_names = ('symbols', 'sqrt', 'Poly', 'QQ', 'S', 'roots', 'integrate', 'I', 'Rational')
symbols, sqrt, Poly, QQ, S, roots, integrate, I, Rational = (getattr(M, n) for n in from_names)

x = symbols('x')
a = symbols('a')
p = symbols('p', positive=True)
s2 = sqrt(2)
K = QQ.algebraic_field(s2)
K3 = QQ.algebraic_field(sqrt(3))

def lp(num, den, dom):
    return lambda: ratint_logpart(Poly(num, x, domain=dom), Poly(den, x, domain=dom), x)

cases = [
    # ratint_logpart over QQ<a> (TypeError before)
    ('logpart QQ<s2> 1/(x**2+s2)', lp(1, x**2 + s2, K)),
    ('logpart QQ<s2> 1/((x-s2)(x+1))', lp(1, (x - s2)*(x + 1), K)),
    ('logpart QQ<s2> x/(x**3+s2)', lp(x, x**3 + s2, K)),
    ('logpart QQ<s2> 1/(x**3-s2*x+1)', lp(1, x**3 - s2*x + 1, K)),
    ('logpart QQ<s2> (x+s2)/(x**4+s2*x**2+1)', lp(x + s2, x**4 + s2*x**2 + 1, K)),
    ('logpart QQ<s3> 1/(x**2-s3*x+1)', lp(1, x**2 - sqrt(3)*x + 1, K3)),
    ('logpart QQ<s2> extension=True', lambda: ratint_logpart(Poly(1, x, extension=True), Poly(x**2 + s2, x, extension=True), x)),
    ('logpart QQ<s2>/QQ<s3> 1/(x**2+s3)', lambda: ratint_logpart(Poly(1, x, domain=K), Poly(x**2 + sqrt(3), x, domain=K3), x)),
    ('logpart QQ<s2>/QQ<s3> (x+s2)/(x**2+s3*x+1)', lambda: ratint_logpart(Poly(x + s2, x, domain=K), Poly(x**2 + sqrt(3)*x + 1, x, domain=K3), x)),
    ('logpart QQ<s2>/QQ<2**(1/3)> x/(x**3+2**(1/3))', lambda: ratint_logpart(Poly(x, x, domain=K), Poly(x**3 + 2**Rational(1, 3), x, domain=QQ.algebraic_field(2**Rational(1, 3))), x)),
    ('ratint tuple QQ<s2>', lambda: ratint((Poly(1, x, domain=K), Poly(x**2 + s2, x, domain=K)), x)),
    # log_to_atan over EX and public integrals that reach it
    ('log_to_atan QQ/EX doctest', lambda: log_to_atan(Poly(x + S(1)/2, x, domain='QQ'), Poly(sqrt(3)/2, x, domain='EX'))),
    ('log_to_atan EX non-constant', lambda: log_to_atan(Poly(x**2 + s2*x, x, domain='EX'), Poly(x + s2, x, domain='EX'))),
    ('log_to_atan EX deg2/deg1', lambda: log_to_atan(Poly(x**2 + s2, x, domain='EX'), Poly(s2*x + 1, x, domain='EX'))),
    ('integrate 1/(x**4+2)', lambda: integrate(1/(x**4 + 2), x)),
    ('integrate 1/(x**2+sqrt(2))', lambda: integrate(1/(x**2 + s2), x)),
    ('integrate 1/(x**3+sqrt(2))', lambda: integrate(1/(x**3 + s2), x)),
    ('integrate 1/(x**2+p)', lambda: integrate(1/(x**2 + p), x)),
    ('integrate 1/(x**5+1)', lambda: integrate(1/(x**5 + 1), x)),
    ('integrate x/(x**5+1)', lambda: integrate(x/(x**5 + 1), x)),
    ('integrate 1/(x**3-2)', lambda: integrate(1/(x**3 - 2), x)),
    ('integrate 1/(x**3+x+1)', lambda: integrate(1/(x**3 + x + 1), x)),
    ('integrate 1/(x**3+2*x-1)', lambda: integrate(1/(x**3 + 2*x - 1), x)),
    ('integrate 1/((x**2+1)*(x**3+2))', lambda: integrate(1/((x**2 + 1)*(x**3 + 2)), x)),
    # real=False RootSum
    ('ratint real=False x/(x**3-a*x+1)', lambda: ratint(x/(x**3 - a*x + 1), x, real=False)),
    ('ratint real=False 1/(x**4+x+a)', lambda: ratint(1/(x**4 + x + a), x, real=False)),
    ('ratint real=False 1/(x**5+1)', lambda: ratint(1/(x**5 + 1), x, real=False)),
    ('ratint real=False 1/(x**2+2)', lambda: ratint(1/(x**2 + 2), x, real=False)),
    # roots
    ('roots x**8-2*x**4+x**2-3', lambda: roots(x**8 - 2*x**4 + x**2 - 3, x)),
    ('roots (x**4-2*x**2-3)**2', lambda: roots((x**4 - 2*x**2 - 3)**2, x)),
    ('roots octic in x**2+1', lambda: roots((x**2 + 1)**4 - 2*(x**2 + 1)**2 + (x**2 + 1) - 3, x)),
    ('roots case 20 x**4+2*x**2+3', lambda: roots(x**4 + 2*x**2 + 3, x)),
    ('roots x**6+x**2+1', lambda: roots(x**6 + x**2 + 1, x)),
    ('roots x**4+x**2+1', lambda: roots(x**4 + x**2 + 1, x)),
    ('roots x**4-x**2+3', lambda: roots(x**4 - x**2 + 3, x)),
    ('roots x**8+x**4+3', lambda: roots(x**8 + x**4 + 3, x)),
    ('roots x**4+a*x**2+1', lambda: roots(x**4 + a*x**2 + 1, x)),
    ('roots x**6-3*x**3+1', lambda: roots(x**6 - 3*x**3 + 1, x)),
]

with open(out, 'w') as fh:
    for name, f in cases:
        try:
            r = f(); text = str(r); err = None
        except Exception as e:
            text = None; err = f'{type(e).__name__}: {str(e)[:200]}'
        best = None
        if err is None:
            total = 0.0
            for _ in range(3):
                t = time.perf_counter(); f(); dt = time.perf_counter() - t
                best = dt if best is None else min(best, dt); total += dt
                if total > 10: break
        fh.write(json.dumps({'case': name, 'result': text, 'error': err, 'seconds': best}) + '\n')
        fh.flush()
        print(name, '|', err or f'{len(text)} chars', '|', best, flush=True)
