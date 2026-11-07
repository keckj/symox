"""Differential corpus for the algpoly lot: composite domains over an
algebraic ground, QQ<a>[t] and QQ<a>(t), their elements and polynomials,
construct_domain(extension=True), unification, EX elements and the
ratint_logpart cases left by the ratint2 lot.

Usage: python diffcorpus.py sympy|symox OUT.jsonl [case names...]

Each record keeps the printed result (or the exception class and message)
and a `value`, the result as an expression (as_expr applied through
tuples and lists), used by compare.py to tell a representation difference
from a wrong value."""
import json, sys, time

mode, out = sys.argv[1], sys.argv[2]
if mode == 'sympy':
    import sympy as M
    from sympy.integrals.rationaltools import ratint_logpart
    from sympy.polys.constructor import construct_domain
else:
    import symox as M
    from symox.integrals.rationaltools import ratint_logpart
    from symox.polys.constructor import construct_domain
names = ('symbols', 'sqrt', 'Poly', 'QQ', 'ZZ', 'EX', 'I', 'Rational', 'cancel', 'factor', 'gcd',
         'resultant', 'sqf_list', 'factor_list', 'discriminant')
(symbols, sqrt, Poly, QQ, ZZ, EX, I, Rational, cancel, factor, gcd, resultant, sqf_list, factor_list,
 discriminant) = (getattr(M, n) for n in names)

x, y, t, u = symbols('x y t u')
p = symbols('p', positive=True)
s2, s3 = sqrt(2), sqrt(3)
A = QQ.algebraic_field(s2)
R = A[t]           # QQ<sqrt(2)>[t]
F = A.frac_field(t)  # QQ<sqrt(2)>(t)


def PR(e, *g):
    return Poly(e, *(g or (x,)), domain=R)


def PF(e, *g):
    return Poly(e, *(g or (x,)), domain=F)


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


f1 = lambda: PR(s2*t*x**2 + (t + s2)*x + 2)
f2 = lambda: PR((x + s2*t)*(x - 1))
f3 = lambda: PR((x + s2*t)**2*(x + t))
g1 = lambda: PF(x**2/(s2*t) + x/t + 1)
g2 = lambda: PF((x + s2/t)*(t*x - 1))

cases = [
    # domains and elements
    ('R str', lambda: R),
    ('F str', lambda: F),
    ('R element', lambda: R(s2*t + 1)),
    ('R element repr', lambda: repr(R(s2*t**2 + 2))),
    ('R element sqrt2/2', lambda: R(s2/2*t)),
    ('R mul', lambda: R(s2*t + 1)*R(t - s2)),
    ('R add', lambda: R(s2*t) + R(t + s2)),
    ('R pow', lambda: R(s2*t + 1)**3),
    ('R exquo', lambda: R.exquo(R(2*t**2 - 4), R(s2*t + 2))),
    ('R exquo fails', lambda: R.exquo(R(t), R(t + 1))),
    ('R gcd', lambda: R.gcd(R(2*t**2 - 4), R(s2*t + 2))),
    ('R gcd monomial', lambda: R.gcd(R(2*s2*t**2), R(4*t))),
    ('R lcm', lambda: R.lcm(R(2*t**2 - 4), R(s2*t + 2))),
    ('R from 1/t fails', lambda: R(1/t)),
    ('R is_negative', lambda: R.is_negative(R(-s2*t + 1))),
    ('F element', lambda: F(t/(s2*t + 1))),
    ('F element inverse', lambda: F(1/(s2*t))),
    ('F element sqrt2/t', lambda: F(s2/t)),
    ('F element scaled', lambda: F(t/(2*t + 2))),
    ('F element sum', lambda: F(1/t + s2/(t + 1))),
    ('F element power', lambda: F((s2*t + 1)**2/t**3)),
    ('F add', lambda: F(1/t) + F(s2/(t + 1))),
    ('F sub', lambda: F(t) - F(1/(s2*t))),
    ('F mul', lambda: F(t/(s2*t + 2))*F((t + s2)/t)),
    ('F div', lambda: F(t)/F(s2*t + 2)),
    ('F div 2', lambda: F(s2*t + 1)/F(2*t)),
    ('F div by zero', lambda: F(t)/F(0)),
    ('F pow', lambda: F(t/(s2*t + 1))**-2),
    ('F gcd', lambda: F.gcd(F((s2*t + 1)/t), F(1/(s2*t)))),
    ('F lcm', lambda: F.lcm(F((s2*t + 1)/t), F(1/(s2*t)))),
    ('F numer denom', lambda: (F.numer(F(t/(2*t + 2))), F.denom(F(t/(2*t + 2))))),
    ('F from x fails', lambda: F(x)),
    ('F of R element', lambda: F.convert(R(s2*t), R)),
    # polynomials over QQ<a>[t]
    ('PR str', f1),
    ('PR domain', lambda: f1().domain),
    ('PR as_expr', lambda: f1().as_expr()),
    ('PR add', lambda: f1() + f2()),
    ('PR mul', lambda: f1()*f2()),
    ('PR pow', lambda: f2()**2),
    ('PR LC', lambda: f1().LC()),
    ('PR diff', lambda: f1().diff(x)),
    ('PR eval', lambda: f1().eval(x, s2)),
    ('PR exquo', lambda: (f2()*f1()).exquo(f2())),
    ('PR div', lambda: f1().div(f2())),
    ('PR rem', lambda: f3().rem(f2())),
    ('PR pdiv', lambda: f1().pdiv(f2())),
    ('PR prem', lambda: f3().prem(f2())),
    ('PR gcd', lambda: f3().gcd(f2())),
    ('PR gcd 2', lambda: PR(2*x**2 - 4*t**2).gcd(PR(s2*x + 2*t))),
    ('PR lcm', lambda: f3().lcm(f2())),
    ('PR cofactors', lambda: f3().cofactors(f2())),
    ('PR cancel', lambda: f3().cancel(f2())),
    ('PR cancel signs', lambda: PR(-s2*x - 2*t).cancel(PR(-x**2 + 2*t**2))),
    ('PR content', lambda: PR(2*s2*t*x + 4*t).content()),
    ('PR primitive', lambda: PR(2*s2*t*x + 4*t**2).primitive()),
    ('PR resultant', lambda: f1().resultant(f2())),
    ('PR resultant includePRS', lambda: f1().resultant(f3(), includePRS=True)),
    ('PR discriminant', lambda: f1().discriminant()),
    ('PR sqf_list', lambda: f3().sqf_list()),
    ('PR sqf_part', lambda: f3().sqf_part()),
    ('PR factor_list', lambda: f3().factor_list()),
    ('PR factor_list 2', lambda: PR(x**2 - 2*t**2).factor_list()),
    ('PR to_field', lambda: f1().to_field()),
    ('PR monic fails', lambda: f1().monic()),
    ('PR clear_denoms', lambda: f1().clear_denoms()),
    ('PR multivariate gcd', lambda: PR((x + s2*y*t)*(x - y), x, y).gcd(PR((x + s2*y*t)*(y + t), x, y))),
    ('PR set_domain EX', lambda: f1().set_domain(EX)),
    # polynomials over QQ<a>(t)
    ('PF str', g1),
    ('PF domain', lambda: g1().domain),
    ('PF as_expr', lambda: g1().as_expr()),
    ('PF add', lambda: g1() + g2()),
    ('PF mul', lambda: g1()*g2()),
    ('PF LC', lambda: g2().LC()),
    ('PF monic', lambda: g2().monic()),
    ('PF div', lambda: g1().div(g2())),
    ('PF rem', lambda: (g1()*g2() + PF(x/t)).rem(g2())),
    ('PF exquo', lambda: (g1()*g2()).exquo(g2())),
    ('PF gcd', lambda: (g1()*g2()).gcd(g2()*PF(x + t))),
    ('PF lcm', lambda: g1().lcm(g2())),
    ('PF cofactors', lambda: (g1()*g2()).cofactors(g2())),
    ('PF cancel', lambda: (g1()*g2()).cancel(g2()*PF(x + t))),
    ('PF content', lambda: PF(2*s2*x/t + 4/(t + 1)).content()),
    ('PF primitive', lambda: PF(2*s2*x/t + 4/(t + 1)).primitive()),
    ('PF clear_denoms', lambda: g1().clear_denoms()),
    ('PF resultant', lambda: g1().resultant(g2())),
    ('PF discriminant', lambda: g2().discriminant()),
    ('PF sqf_list', lambda: (g2()**2*PF(x + t)).sqf_list()),
    ('PF factor_list', lambda: (g2()*PF(x**2 - 2*t**2)).factor_list()),
    ('PF to_ring fails', lambda: g1().to_ring()),
    ('PF invert', lambda: PF(x + s2).invert(PF(x**2 + t))),
    ('PF gcdex', lambda: PF(x + s2).gcdex(PF(x**2 + t))),
    ('PF set_domain EX', lambda: g1().set_domain(EX)),
    # construction with extension=True and unification
    ('construct_domain ext poly', lambda: construct_domain([s2*t, 1], extension=True)),
    ('construct_domain ext frac', lambda: construct_domain([s2/t, 1], extension=True)),
    ('construct_domain ext field', lambda: construct_domain([s2*t, 1], extension=True, field=True)),
    ('Poly extension ring', lambda: Poly(s2*x*t + 1, x, extension=True)),
    ('Poly extension frac', lambda: Poly(x/(s2*t) + 1, x, extension=True)),
    ('Poly extension domain', lambda: Poly(s2*x*t + s3, x, extension=True).domain),
    ('unify R with ZZ[y]', lambda: PR(x*t) + Poly(x*y, x, domain='ZZ[y]')),
    ('unify R with QQ', lambda: PR(x*t) + Poly(x/2, x)),
    ('unify R with A', lambda: PR(x*t) * Poly(s2*x, x, domain=A)),
    ('unify R with F', lambda: PR(x*t) + PF(x/t)),
    ('unify R with ZZ(y)', lambda: PR(x*t) * Poly(x/y, x, domain='ZZ(y)')),
    ('unify R generator clash', lambda: PR(x*t) + Poly(x*t, x, t)),
    ('unify R with QQ<sqrt3>', lambda: PR(x*t) + Poly(s3*x, x, extension=True)),
    ('unify R with QQ_I', lambda: PR(x*t) + Poly(I*x/2, x)),
    ('unify F with ZZ_I[y]', lambda: PF(x/t) * Poly(I*x*y, x, domain='ZZ_I[y]')),
    ('unify R with RR', lambda: PR(x*t) + Poly(0.5*x, x)),
    ('unify F with EX', lambda: PF(x/t) + Poly(x*sqrt(y), x, domain=EX)),
    ('unify R with QQ<sqrt3>(y)', lambda: PR(x*t) + Poly(s3*x/y, x, domain=QQ.algebraic_field(s3).frac_field(y))),
    ('domain unify R F', lambda: R.unify(F)),
    ('domain unify R ZZ[y]', lambda: R.unify(ZZ[y])),
    # public functions
    ('cancel extension', lambda: cancel((x**2 - 2*t**2)/(x - s2*t), x, extension=True)),
    ('factor extension param', lambda: factor(x**2 - 2*t**2, x, extension=s2)),
    ('gcd extension param', lambda: gcd(x**2 - 2*t**2, x + s2*t, x, extension=True)),
    ('resultant extension param', lambda: resultant(x**2 + s2*t, x - t, x, extension=True)),
    ('sqf_list extension param', lambda: sqf_list((x + s2*t)**2*(x + 1), x, extension=True)),
    ('factor_list extension param', lambda: factor_list(x**2 - 2*t**2, x, extension=True)),
    ('discriminant extension param', lambda: discriminant(x**2 + s2*t*x + 1, x, extension=True)),
    # EX elements
    ('EX sqrt2*(1+sqrt2)', lambda: EX(s2*(1 + s2))),
    ('EX (x+1)**2', lambda: EX((x + 1)**2)),
    ('EX add', lambda: EX(x) + EX(1/x)),
    ('EX mul numbers', lambda: EX(2)*EX(Rational(1, 3))),
    ('EX mul', lambda: EX(x + 1)*EX(x - 1)),
    ('EX div', lambda: EX(x**2 - 1)/EX(x + 1)),
    ('EX pow', lambda: EX(s2 + 1)**2),
    ('EX sub zero', lambda: EX(0) - EX(x + 1)),
    ('EX add zero', lambda: EX(s2*(1 + s2)) + EX(0)),
    ('EX convert', lambda: EX.convert(s2*(1 + s2))),
    # ratint_logpart cases left by ratint2
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
]


def value(r):
    """The result as a SymPy-parsable expression string, or None."""
    try:
        if isinstance(r, (tuple, list)):
            parts = [value(v) for v in r]
            return None if None in parts else '(' + ', '.join(parts) + ',)'
        if hasattr(r, 'as_expr') and callable(r.as_expr):
            return str(r.as_expr())
        if hasattr(r, 'ex'):
            return str(r.ex)
        if isinstance(r, (bool, int)):
            return str(r)
        if hasattr(r, 'is_Number') or hasattr(r, 'free_symbols'):
            return str(r)
    except Exception:
        return None
    return None


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
        val = None
        try:
            t0 = time.perf_counter()
            r = fn(); text = str(r); err = None; val = value(r)
            dt = time.perf_counter() - t0
        except (Exception, BaseException) as e:  # a native panic is a BaseException
            if isinstance(e, KeyboardInterrupt):
                raise
            text = None; err = f'{type(e).__name__}: {str(e)[:200]}'; dt = None
        finally:
            signal.alarm(0)
        fh.write(json.dumps({'case': name, 'result': text, 'value': val, 'error': err, 'seconds': dt}) + '\n')
        fh.flush()
        print(name, '|', err or text[:150], flush=True)
