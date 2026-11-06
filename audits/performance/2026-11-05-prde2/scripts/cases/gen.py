"""Run the regression cases in one backend (argv[1]: sympy|symox) and print
encoded results (srepr, domains, Dummy count)."""
import sys, importlib
M = importlib.import_module(sys.argv[1])
prde = importlib.import_module(sys.argv[1] + '.integrals.prde')
risch = importlib.import_module(sys.argv[1] + '.integrals.risch')
Poly, S, Dummy, srepr = M.Poly, M.S, M.Dummy, M.srepr
x, t, t0, t1 = M.symbols('x t t0 t1')
DEx = risch.DifferentialExtension

def enc(v):
    if isinstance(v, Poly):
        return ('P', srepr(v.as_expr()), tuple(srepr(g) for g in v.gens), str(v.domain))
    if isinstance(v, (list, tuple)):
        return tuple(enc(i) for i in v)
    if hasattr(v, 'rows') and hasattr(v, 'domain'):
        return ('M', v.rows, v.cols, tuple(srepr(v[i, j].as_expr()) for i in range(v.rows) for j in range(v.cols)), str(v.domain))
    if v is None:
        return None
    return srepr(v)

def P(e, g, **kw):
    return Poly(M.sympify(e), g, **kw)

def cases():
    D1 = lambda: DEx(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    Dexp = lambda: DEx(extension={'D': [Poly(1, x), Poly(-t, t)]})
    Dt = lambda: DEx(extension={'D': [Poly(1, x), Poly(t, t)]})
    D3 = lambda: DEx(extension={'D': [Poly(1, x), Poly(1/x, t0), Poly(1/(x*t0), t1)]})
    yield 'cancel_prim', lambda: prde.prde_cancel_liouvillian(P('-1/(x - 1)', t), [P('-x + 1', t), P('1', t)], 1, D1())
    yield 'cancel_exp', lambda: prde.prde_cancel_liouvillian(Poly(0, t, domain='QQ[x]'), [Poly(1, t, domain='QQ(x)')], 0, Dexp())
    yield 'risch_exp', lambda: prde.param_rischDE(P('-t - 1', t, field=True), P('t + x', t, field=True),
        [(Poly(t + x, t, domain='ZZ(x)'), Poly(1, t, domain='QQ')), (Poly(0, t, domain='QQ'), Poly(1, t, domain='QQ'))], Dt())
    yield 'risch_prim', lambda: prde.param_rischDE(P('1', t), P('x', t), [(P('1', t), P('x*t', t)), (P('t', t), P('1', t))], D1())
    yield 'poly_prim', lambda: prde.param_poly_rischDE(P('1', t), P('0', t), [P('t', t), P('1/x', t)], 2, D1())
    yield 'poly_exp_large', lambda: prde.param_poly_rischDE(P('1', t), P('t', t), [P('t**2', t), P('t', t)], 2, Dt())
    yield 'poly_spde', lambda: prde.param_poly_rischDE(P('t', t), P('1', t), [P('t**2 + x', t), P('t', t)], 2, Dt())
    yield 'poly_neg', lambda: prde.param_poly_rischDE(P('t', t), P('1', t), [P('t**2 + x', t), P('t', t)], -1, Dt())
    yield 'li_prim', lambda: prde.limited_integrate(P('1', t), P('x', t), [(P('1', t), P('x', t))], D1())
    yield 'li_prim2', lambda: prde.limited_integrate(P('t', t), P('x', t), [(P('1', t), P('x', t))], D1())
    yield 'li_exp', lambda: prde.limited_integrate(P('t', t), P('1', t), [(P('1', t), P('x', t))], Dt())
    yield 'li_none', lambda: prde.limited_integrate(P('1', t), P('t', t), [(P('t', t), P('1', t))], Dt())
    yield 'li_tower3', lambda: prde.limited_integrate(P('1', t1), P('x*t0', t1), [(P('1', t1), P('x', t1))], D3())
    yield 'risch_tower3', lambda: prde.param_rischDE(P('0', t1), P('1', t1), [(P('1', t1), P('x*t0', t1)), (P('t1', t1), P('x', t1))], D3())

def more():
    Dt = lambda: DEx(extension={'D': [Poly(1, x), Poly(t, t)]})
    D1 = lambda: DEx(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    G2 = lambda: [(P('1', t), P('1', t)), (P('t', t), P('x', t))]
    yield 'wn_root1', lambda: prde.param_rischDE(P('t + 1', t), P('t + x', t), G2(), Dt())
    yield 'wn_symroot', lambda: prde.param_rischDE(P('x*t', t), P('t + x', t), G2(), Dt())
    yield 'wn_strip', lambda: prde.param_rischDE(P('x*t + 1', t), P('x*t', t), G2(), Dt())
    yield 'wn_prim_root', lambda: prde.param_rischDE(P('2/x', t), P('t + x', t), G2(), D1())
    yield 'li_exp_den', lambda: prde.limited_integrate(P('1', t), P('t + 1', t), [(P('t', t), P('t + 1', t))], Dt())
    yield 'li_prim_den', lambda: prde.limited_integrate(P('1', t), P('x*t + x', t), [(P('1', t), P('x*t + x', t))], D1())
    yield 'li_prim_base_den', lambda: prde.limited_integrate(P('t', t), P('x - 1', t), [(P('1', t), P('x - 1', t))], D1())
    yield 'cancel_exp_den', lambda: prde.prde_cancel_liouvillian(Poly(0, t, domain='QQ(x)'), [Poly(1/(x + 1), t, domain='QQ(x)'), Poly(t/(x-1), t, domain='QQ(x)')], 1, DEx(extension={'D': [Poly(1, x), Poly(-t, t)]}))
    yield 'poly_prim_bsmall', lambda: prde.param_poly_rischDE(P('1', t), P('0', t), [P('t/(x+1)', t), P('1/(x**2-1)', t), P('t**2', t)], 3, D1())

sys.path.insert(0, '/tmp/claude-1000/-home-keckj-Documents-symox/70ac2d6f-3a62-411a-a6c7-c410b2a23124/scratchpad/prde2')
from more2 import more2
only = set(sys.argv[2:])
import itertools
for name, f in itertools.chain(cases(), more(), more2(prde, P, DEx, Poly, x, t)):
    if only and name not in only:
        continue
    d0 = Dummy().dummy_index
    try:
        r = enc(f())
    except Exception as e:
        r = ('raise', type(e).__name__, str(e))
    n = Dummy().dummy_index - d0 - 1
    print(repr((name, r, n)))
if sys.argv[1] == 'symox':
    print(importlib.import_module('symox.integrals._prde_native').STATS, file=sys.stderr)
