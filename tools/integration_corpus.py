"""Differential corpus of definite integrals and transforms against SymPy.

Each case is evaluated by SymPy 1.14 and by symox in separate persistent
worker processes (one backend at a time, never in parallel), with a hard
per-case timeout enforced by the parent: a stuck worker is killed and
restarted. Results are compared in the parent with SymPy:

* SAME        structurally or symbolically identical results;
* EQUIV       different form, numerically equal at sample parameters;
* MISSING     symox keeps the integral, SymPy evaluates it (a gap);
* EXTRA       symox evaluates, SymPy does not; checked against mpmath
              quadrature at sample parameters (EXTRA-OK, EXTRA-WRONG or
              EXTRA-UNCHECKED);
* DIFF        both evaluate but differ numerically (a wrong result unless
              SymPy is the one that is wrong; quadrature arbitrates);
* COND        same values, different convergence conditions;
* ERROR       an exception in symox where SymPy returned a result or a
              different exception class;
* TIMEOUT     the symox worker exceeded the per-case limit.

Usage: python tools/integration_corpus.py [--timeout S] [--only SUBSTR]
       [--reference-cache FILE] [--json OUT]
Run it under a memory cap (systemd-run --user --scope -p MemoryMax=4G ...).
"""
from __future__ import annotations

import argparse
import json
import os
import select
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# (setup, call). The setup defines symbols; the call is evaluated to give the
# result. Names come from `from <backend> import *`.
PARAMS = {
    'none': "s,a,b,n,k,w=symbols('s a b n k w')",
    'pos': "s,a,b,k,w=symbols('s a b k w',positive=True); n=Symbol('n',positive=True,integer=True)",
    'real': "s,a,b,k,w=symbols('s a b k w',real=True); n=Symbol('n',integer=True)",
    'cplx': "s,a,b,k,w=symbols('s a b k w',complex=True); n=Symbol('n',integer=True)",
    'nz': "s,a,b,k,w=symbols('s a b k w',nonzero=True); n=Symbol('n',integer=True,nonnegative=True)",
}
BASE = "x,t,y=symbols('x t y'); xp=Symbol('xp',positive=True); tr=Symbol('tr',real=True)"

LAPLACE_KERNELS = [
    '1', 't', 't**2', 't**3', 'sqrt(t)', '1/sqrt(t)', 't**a', 'exp(-t)', 'exp(-2*t)',
    't*exp(-t)', 'exp(a*t)', 'exp(-a*t)', 'sin(t)', 'cos(t)', 'sin(a*t)', 'cos(a*t)',
    'sinh(t)', 'cosh(t)', 't*sin(t)', 'exp(-t)*sin(t)', 'Heaviside(t-1)', 'Heaviside(t-a)',
    '1-Heaviside(t-1)', 'Heaviside(t)-Heaviside(t-2)', 't*Heaviside(t-1)', 'DiracDelta(t-1)',
    'log(t)', 'cos(sqrt(t))', 'erf(t)', 'besselj(0,t)', 'sin(t)/t', 'exp(-t**2)', '1/(t+1)',
]
TRANSFORM_CASES = [
    ('laplace_transform({f}, t, s)', f) for f in [
        'cos(sqrt(t))', 'sin(sqrt(t))', 'sqrt(t)*exp(-t)', 'besselj(0,t)', 'besselj(1,t)',
        'erf(t)', 'erfc(sqrt(t))', 'log(t)', 't**a', 'exp(-t**2)', '1/(t+1)', 'Si(t)',
        'atan(t)', 'sin(t)/t', 'exp(-a*t)*cos(b*t)', 'Heaviside(t-a)*sin(t)', 'exp(-t)*t**a',
        'cosh(sqrt(t))', 'exp(-1/t)', 'besseli(0,t)', 'sinh(t)/t',
    ]
] + [
    ('inverse_laplace_transform({f}, s, t)', f) for f in [
        '1/sqrt(s)', 'exp(-sqrt(s))', 'log(s)/s', '1/(s**2+1)**2', 'exp(-a*s)/s', 'atan(1/s)',
        '1/sqrt(s**2+1)', 'exp(-1/s)/s', 's/(s**2+a**2)', '1/(s*(s+1))',
    ]
] + [
    ('fourier_transform({f}, x, k)', f) for f in [
        'exp(-x**2)', 'exp(-Abs(x))', '1/(1+x**2)', 'Heaviside(x)*exp(-x)', 'exp(-a*x**2)',
        'x*exp(-x**2)', 'sin(x)/x', 'Heaviside(1-Abs(x))',
    ]
] + [
    ('mellin_transform({f}, x, s)', f) for f in [
        'exp(-x)', '1/(1+x)', 'sin(x)', 'cos(x)', 'exp(-x**2)', 'log(1+x)', 'besselj(0,x)',
        'Heaviside(1-x)', 'atan(x)', '1/(1+x)**2',
    ]
] + [
    ('hankel_transform({f}, x, k, 0)', f) for f in ['exp(-x)', '1/x', 'exp(-x**2)', 'exp(-a*x)/x']
]

DEFINITE = [
    # improper, no parameters
    ('exp(-x)', '(x,0,oo)'), ('x*exp(-x)', '(x,0,oo)'), ('1/(1+x**2)', '(x,-oo,oo)'),
    ('1/x**2', '(x,1,oo)'), ('1/sqrt(x)', '(x,0,1)'), ('log(x)', '(x,0,1)'),
    ('sin(x)/x', '(x,0,oo)'), ('sin(x)**2/x**2', '(x,0,oo)'), ('exp(-x**2)', '(x,-oo,oo)'),
    ('x**2*exp(-x**2)', '(x,-oo,oo)'), ('1/x', '(x,0,1)'), ('1/x', '(x,-1,1)'),
    ('1/(x-1)', '(x,0,2)'), ('1/x**2', '(x,-1,1)'), ('tan(x)', '(x,0,pi)'),
    ('1/sin(x)', '(x,0,pi)'), ('sec(x)**2', '(x,0,pi)'), ('1/(x**2-1)', '(x,0,2)'),
    ('log(x)/(1+x**2)', '(x,0,oo)'), ('x/(exp(x)-1)', '(x,0,oo)'), ('exp(-x)*sin(x)', '(x,0,oo)'),
    ('exp(-x)*cos(x)', '(x,0,oo)'), ('sin(x)', '(x,0,oo)'), ('cos(x)', '(x,0,oo)'),
    ('sin(x**2)', '(x,0,oo)'), ('cos(x)/x', '(x,1,oo)'), ('exp(I*x)', '(x,0,oo)'),
    ('1/(1+x**4)', '(x,0,oo)'), ('x**2/(1+x**4)', '(x,0,oo)'), ('1/sqrt(1-x**2)', '(x,-1,1)'),
    ('log(1+x)/x', '(x,0,1)'), ('exp(-x)/sqrt(x)', '(x,0,oo)'), ('1/(x*(1+x))', '(x,1,oo)'),
    ('exp(-x)*log(x)', '(x,0,oo)'), ('1/cosh(x)', '(x,-oo,oo)'), ('x/sinh(x)', '(x,0,oo)'),
    ('exp(x)', '(x,-oo,0)'), ('exp(x)', '(x,0,oo)'), ('x**3', '(x,-oo,oo)'), ('x', '(x,-oo,oo)'),
    ('1/(x**2+2*x+2)', '(x,-oo,oo)'), ('exp(-x)*Heaviside(x-1)', '(x,0,oo)'),
    ('exp(-x)*Heaviside(x-1)', '(x,-oo,oo)'), ('Heaviside(x)*exp(-x**2)', '(x,-oo,oo)'),
    ('DiracDelta(x-1)*exp(-x)', '(x,0,oo)'), ('exp(-Abs(x))', '(x,-oo,oo)'),
    ('besselj(0,x)', '(x,0,oo)'), ('exp(-x)*besselj(0,x)', '(x,0,oo)'), ('erf(x)*exp(-x)', '(x,0,oo)'),
    ('exp(-x)*exp(-2*x)', '(x,0,oo)'), ('x**2*exp(-3*x)', '(x,0,oo)'),
    ('1/(1+x)**2', '(x,0,oo)'), ('x/(1+x)**3', '(x,0,oo)'), ('atan(x)/x**2', '(x,1,oo)'),
    ('exp(-x)*x**Rational(1,3)', '(x,0,oo)'), ('cos(x)*exp(-x**2)', '(x,-oo,oo)'),
    # with parameters
    ('exp(-s*x)', '(x,0,oo)'), ('x*exp(-s*x)', '(x,0,oo)'), ('x**2*exp(-s*x)', '(x,0,oo)'),
    ('exp(-a*x)', '(x,0,oo)'), ('exp(a*x)', '(x,0,oo)'), ('exp(a*x)', '(x,-oo,0)'),
    ('x**a', '(x,0,1)'), ('x**a', '(x,1,oo)'), ('x**n', '(x,0,1)'), ('x**a*exp(-x)', '(x,0,oo)'),
    ('1/(x**2+a**2)', '(x,0,oo)'), ('1/(x**2+a**2)', '(x,-oo,oo)'), ('1/(x**2+a)', '(x,0,oo)'),
    ('exp(-a*x**2)', '(x,-oo,oo)'), ('exp(-a*x**2)', '(x,0,oo)'), ('cos(a*x)*exp(-x)', '(x,0,oo)'),
    ('sin(a*x)/x', '(x,0,oo)'), ('cos(s*x)', '(x,0,oo)'), ('sin(s*x)', '(x,0,oo)'),
    ('exp(-s*x)*sin(x)', '(x,0,oo)'), ('exp(-s*x)*cos(x)', '(x,0,oo)'),
    ('exp(-s*x)*sin(a*x)', '(x,0,oo)'), ('exp(-a*x)*exp(-b*x)', '(x,0,oo)'),
    ('exp(-x)*exp(-s*x)', '(x,0,oo)'), ('x*exp(-x)*exp(-s*x)', '(x,0,oo)'),
    ('exp(-s*x)*Heaviside(x-1)', '(x,0,oo)'), ('exp(-s*x)*Heaviside(x-a)', '(x,0,oo)'),
    ('(1-Heaviside(x-1))*exp(-s*x)', '(x,0,oo)'), ('exp(-s*x)*DiracDelta(x-1)', '(x,0,oo)'),
    ('exp(-s*x)*DiracDelta(x)', '(x,-1,oo)'), ('exp(-s*x)*Heaviside(x-5)', '(x,-1,oo)'),
    ('exp(-s*x)*sinh(x)', '(x,0,oo)'), ('exp(-s*x)*cosh(a*x)', '(x,0,oo)'),
    ('exp(-s*x)*sqrt(x)', '(x,0,oo)'), ('exp(-s*x)/sqrt(x)', '(x,0,oo)'),
    ('exp(-s*x)*x**a', '(x,0,oo)'), ('exp(-s*x)*log(x)', '(x,0,oo)'),
    ('exp(-s*x**2)', '(x,-oo,oo)'), ('exp(-(x-a)**2)', '(x,-oo,oo)'),
    ('x**(s-1)*exp(-x)', '(x,0,oo)'), ('x**(s-1)/(1+x)', '(x,0,oo)'),
    ('1/(1+a*x)', '(x,0,1)'), ('1/(x-a)', '(x,0,1)'), ('1/(x+a)', '(x,0,1)'),
    ('exp(-a*Abs(x))', '(x,-oo,oo)'), ('a/(a**2+x**2)', '(x,-oo,oo)'),
    ('log(a*x)', '(x,0,1)'), ('x*exp(-a*x**2)', '(x,0,oo)'), ('exp(-x**2/a)', '(x,-oo,oo)'),
    ('besselj(0,a*x)*exp(-x)', '(x,0,oo)'), ('exp(I*a*x)*exp(-x)', '(x,0,oo)'),
    ('cos(a*x)/(1+x**2)', '(x,0,oo)'), ('exp(-k*x)*cos(w*x)', '(x,0,oo)'),
    ('exp(-s*x)*x**n', '(x,0,oo)'), ('exp(-a*x)', '(x,b,oo)'), ('exp(-s*x)', '(x,a,oo)'),
    ('Heaviside(x-a)', '(x,0,1)'), ('DiracDelta(x-a)', '(x,0,1)'),
    ('sin(x)**2', '(x,0,a)'), ('1/x', '(x,a,1)'), ('x**a', '(x,b,1)'),
]


def build_cases():
    cases = []
    for assumption in PARAMS:
        for f, lim in DEFINITE:
            cases.append(dict(id=f'int[{assumption}] {f} {lim}', setup=PARAMS[assumption],
                              call=f'integrate({f}, {lim})', integrand=f, limits=lim))
        for f in LAPLACE_KERNELS:
            g = f'({f})*exp(-s*t)'
            cases.append(dict(id=f'lap-int[{assumption}] {f}', setup=PARAMS[assumption],
                              call=f'integrate({g}, (t,0,oo))', integrand=g, limits='(t,0,oo)'))
        for template, f in TRANSFORM_CASES:
            call = template.format(f=f)
            cases.append(dict(id=f'tr[{assumption}] {call}', setup=PARAMS[assumption], call=call))
            if template.startswith('laplace') or template.startswith('fourier') or template.startswith('mellin'):
                cases.append(dict(id=f'tr-nc[{assumption}] {call}', setup=PARAMS[assumption],
                                  call=call[:-1] + ', noconds=True)'))
    return cases


WORKER = r'''
import sys, json, traceback, warnings
warnings.simplefilter('ignore')
backend = sys.argv[1]
if backend == 'sympy':
    import sympy as M
else:
    import symox as M
import sympy
ns0 = {}
exec('from %s import *' % backend, ns0)
def conv(v):
    if isinstance(v, tuple):
        return sympy.Tuple(*[conv(i) for i in v])
    if backend == 'symox' and hasattr(v, 'to_sympy'):
        return v.to_sympy()
    return sympy.sympify(v)
for line in sys.stdin:
    case = json.loads(line)
    ns = dict(ns0)
    try:
        exec(BASE, ns)
        exec(case['setup'], ns)
        value = eval(case['call'], ns)
        out = dict(ok=True, srepr=sympy.srepr(conv(value)), text=str(value))
    except BaseException as e:
        out = dict(ok=False, error=type(e).__name__, text=str(e)[:300])
    sys.stdout.write(json.dumps(out) + '\n')
    sys.stdout.flush()
'''


class Worker:
    def __init__(self, backend):
        self.backend = backend
        self.proc = None

    def start(self):
        env = dict(os.environ)
        code = 'BASE=%r\n' % BASE + WORKER
        self.proc = subprocess.Popen([sys.executable, '-c', code, self.backend], stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, env=env)

    def run(self, case, timeout):
        if self.proc is None or self.proc.poll() is not None:
            self.start()
        self.proc.stdin.write(json.dumps(case) + '\n')
        self.proc.stdin.flush()
        ready, _, _ = select.select([self.proc.stdout], [], [], timeout)
        if not ready:
            self.proc.kill()
            self.proc.wait()
            self.proc = None
            return dict(ok=False, error='TIMEOUT', text='')
        line = self.proc.stdout.readline()
        if not line:
            self.proc = None
            return dict(ok=False, error='CRASH', text='')
        return json.loads(line)

    def close(self):
        if self.proc and self.proc.poll() is None:
            self.proc.stdin.close()
            self.proc.wait(timeout=10)


def evaluated(e):
    import sympy as sp
    return not e.has(sp.Integral, sp.LaplaceTransform, sp.InverseLaplaceTransform,
                     sp.FourierTransform, sp.MellinTransform, sp.InverseMellinTransform,
                     sp.HankelTransform) if isinstance(e, sp.Basic) else True


def samples(expr, setup):
    """Parameter samples consistent with the case's assumptions."""
    import sympy as sp
    syms = sorted(expr.free_symbols - {sp.Symbol('x'), sp.Symbol('t')}, key=str)
    tries = [sp.Rational(7, 5), sp.Rational(5, 2), sp.Rational(13, 10) + sp.I / 3, sp.Rational(3, 1)]
    result = []
    for value in tries:
        mapping = {}
        for i, sym in enumerate(syms):
            v = value + sp.Rational(i, 7)
            if sym.is_integer:
                v = sp.Integer(2 + i)
            elif sym.is_real and not v.is_real:
                v = sp.re(v)
            mapping[sym] = v
        result.append(mapping)
    return result


def branch_value(e, mapping):
    """Value of e at the sample: the Piecewise branch whose condition holds."""
    import sympy as sp
    v = e.subs(mapping)
    if isinstance(v, sp.Piecewise):
        for expr, cond in v.args:
            try:
                c = cond if cond in (sp.true, sp.false) else cond.subs(mapping)
            except TypeError:
                return None
            if c == sp.true:
                return expr
            if c != sp.false:
                try:
                    if bool(sp.N(c)):
                        return expr
                    continue
                except TypeError:
                    return None
        return None
    return v


def numeric(e):
    import sympy as sp
    if e is None or not evaluated(e) or e.has(sp.AccumBounds):
        return None
    try:
        v = complex(sp.N(e, 30))
    except (TypeError, ValueError):
        if e in (sp.oo, -sp.oo, sp.zoo):
            return e
        return None
    return v


def quadrature(case, mapping):
    import sympy as sp
    import mpmath
    if 'integrand' not in case:
        return None
    ns = {}
    exec('from sympy import *', ns)
    exec(BASE, ns)
    exec(case['setup'], ns)
    f = eval(case['integrand'], ns)
    var, lo, hi = eval(case['limits'], ns)
    f = f.subs(mapping)
    lo = sp.sympify(lo).subs(mapping)
    hi = sp.sympify(hi).subs(mapping)
    if f.has(sp.DiracDelta) or f.free_symbols - {var}:
        return None
    g = sp.lambdify(var, f, 'mpmath')
    try:
        mpmath.mp.dps = 20
        a = mpmath.mpf(-mpmath.inf) if lo == -sp.oo else mpmath.mpmathify(complex(lo)).real
        b = mpmath.inf if hi == sp.oo else mpmath.mpmathify(complex(hi)).real
        v, err = mpmath.quad(g, [a, b], error=True, maxdegree=8)
        if abs(err) > 1e-6 * max(1, abs(v)) or not mpmath.isfinite(v):
            return None
        return complex(v)
    except Exception:
        return None


def close(u, v):
    if isinstance(u, complex) and isinstance(v, complex):
        return abs(u - v) <= 1e-7 * max(1.0, abs(u), abs(v))
    return u == v


_NAMESPACE = None


def srepr_namespace():
    global _NAMESPACE
    if _NAMESPACE is None:
        import sympy as sp
        import sympy.integrals.risch  # noqa: F401  (NonElementaryIntegral)
        ns = dict(sp.__dict__)
        stack = [sp.Basic]
        while stack:
            cls = stack.pop()
            ns.setdefault(cls.__name__, cls)
            stack.extend(cls.__subclasses__())
        _NAMESPACE = ns
    return _NAMESPACE


def compare(case, ref, nat):
    import sympy as sp
    if not nat['ok']:
        if nat['error'] == 'TIMEOUT':
            return 'TIMEOUT', ''
        if not ref['ok'] and ref['error'] == nat['error']:
            return 'SAME', 'same exception'
        return 'ERROR', f"{nat['error']}: {nat['text'][:120]}"
    if not ref['ok']:
        return ('EXTRA-UNCHECKED' if ref['error'] != 'TIMEOUT' else 'REFTIMEOUT'), f"reference {ref['error']}"
    r = sp.sympify(eval(ref['srepr'], srepr_namespace()))
    n = sp.sympify(eval(nat['srepr'], srepr_namespace()))
    if r == n:
        return 'SAME', ''
    if isinstance(r, sp.Tuple) and isinstance(n, sp.Tuple) and len(r) == len(n):
        head = compare_values(case, r[0], n[0])
        if head[0] in ('SAME', 'EQUIV') and tuple(r[1:]) != tuple(n[1:]):
            return 'COND', head[1] + ' aux: %s vs %s' % (tuple(r[1:]), tuple(n[1:]))
        return head
    return compare_values(case, r, n)


def compare_values(case, r, n):
    import sympy as sp
    if r == n:
        return 'SAME', ''
    re_, ne_ = evaluated(r) or isinstance(r, sp.Piecewise), evaluated(n) or isinstance(n, sp.Piecewise)
    if not ne_ and not re_:
        return 'SAME', 'both unevaluated'
    if not ne_:
        return 'MISSING', ''
    # numeric comparison at parameter samples
    verdicts = []
    for mapping in samples(sp.Tuple(r, n), case['setup']):
        nv = numeric(branch_value(n, mapping))
        rv = numeric(branch_value(r, mapping)) if re_ else None
        qv = quadrature(case, mapping) if (rv is None or nv is None or not close(nv, rv)) else None
        if nv is None:
            if rv is not None:
                verdicts.append('native-held')
            continue
        if rv is not None and close(nv, rv):
            verdicts.append('eq')
        elif qv is not None:
            verdicts.append('quad-ok' if close(nv, qv) else 'quad-wrong')
        elif rv is not None:
            verdicts.append('ne')
        else:
            verdicts.append('unchecked')
    if 'quad-wrong' in verdicts or ('ne' in verdicts and 'quad-ok' not in verdicts):
        return ('DIFF' if re_ else 'EXTRA-WRONG'), ','.join(verdicts)
    if not re_:
        return ('EXTRA-OK' if 'quad-ok' in verdicts else 'EXTRA-UNCHECKED'), ','.join(verdicts)
    if 'native-held' in verdicts and 'eq' not in verdicts:
        return 'COND', ','.join(verdicts)
    if isinstance(r, sp.Piecewise) != isinstance(n, sp.Piecewise) or (
            isinstance(r, sp.Piecewise) and [c for _, c in r.args] != [c for _, c in n.args]):
        return 'COND', ','.join(verdicts)
    return 'EQUIV', ','.join(verdicts)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--timeout', type=float, default=60)
    p.add_argument('--only', default='')
    p.add_argument('--reference-cache', default=str(Path(os.environ.get('TMPDIR', '/tmp')) / 'integration_corpus_reference.json'))
    p.add_argument('--json', default='')
    p.add_argument('--show', default='MISSING,EXTRA-WRONG,EXTRA-UNCHECKED,DIFF,COND,ERROR,TIMEOUT')
    args = p.parse_args()
    cases = [c for c in build_cases() if args.only in c['id']]
    cache_path = Path(args.reference_cache)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    todo = [c for c in cases if c['id'] not in cache]
    if todo:
        w = Worker('sympy')
        for i, c in enumerate(todo):
            cache[c['id']] = w.run(c, args.timeout * 3)
            if i % 20 == 19:
                cache_path.write_text(json.dumps(cache))
        w.close()
        cache_path.write_text(json.dumps(cache))
    w = Worker('symox')
    results = []
    counts = {}
    for c in cases:
        start = time.time()
        nat = w.run(c, args.timeout)
        elapsed = time.time() - start
        try:
            verdict, note = compare(c, cache[c['id']], nat)
        except Exception as e:  # comparison failure is reported, not hidden
            verdict, note = 'COMPARE-FAIL', f'{type(e).__name__}: {e}'[:200]
        counts[verdict] = counts.get(verdict, 0) + 1
        results.append(dict(id=c['id'], verdict=verdict, note=note, time=round(elapsed, 2),
                            reference=cache[c['id']].get('text', cache[c['id']].get('error')),
                            native=nat.get('text', nat.get('error'))))
        if verdict in args.show.split(','):
            print(f"{verdict:16} {c['id']}\n    ref: {results[-1]['reference'][:200]}\n    nat: {results[-1]['native'][:200]}  {note}", flush=True)
    w.close()
    print(json.dumps(counts, sort_keys=True))
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=1))


if __name__ == '__main__':
    main()
