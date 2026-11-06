import sys, json
from importlib import import_module
M = sys.argv[1]
S = import_module(M)
prde = import_module(M + '.integrals.prde')
risch = import_module(M + '.integrals.risch')
CASES = [
 # [kind, D, args]
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '2', ['t0 + x', 'x*t0 - 1'], 2]],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '0', ['t0**2 + 1', '2*t0'], 2]],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '2*t0', ['t0 + x'], 1]],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', '-3*t0', ['t0 + x'], 1]],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['1', 't0**2 + 1', ['t0**3', 'x*t0'], 1]],
 ['param_poly_rischDE', ['1', '1 + t0**2'], ['2', '-x', ['x*t0', '1'], 1]],
 ['param_poly_rischDE', ['1', 't0**3 + 1'], ['1', 't0', ['t0**2 + x', 't0'], 2]],
 ['param_poly_rischDE', ['1', 't0**2'], ['1', '1', ['t0', 'x'], 1]],
 ['param_poly_rischDE', ['1', 't0 + 1'], ['1', '1', ['t0', 'x'], 1]],
 ['param_poly_rischDE', ['1', '1/x', '1 + t1**2'], ['1', '1', ['t1 + t0', 'x'], 1]],
 ['param_rischDE', ['1', '1 + t0**2'], ['1', 't0', [['1', 't0'], ['t0', '1']]]],
 ['param_rischDE', ['1', '1 + t0**2'], ['2*t0', 't0**2 + 1', [['1', 't0**2 + 1']]]],
 ['param_rischDE', ['1', '1 + t0**2'], ['x', '1', [['t0', '1'], ['1', '1']]]],
 ['param_rischDE', ['1', 't0**2'], ['1', 't0', [['1', 't0']]]],
 ['param_rischDE', ['1', 't0 + 1'], ['1', '1', [['t0', '1']]]],
 ['limited_integrate', ['1', '1 + t0**2'], ['2*t0', 't0**2 + 1', [['1', '1']]]],
 ['prde_no_cancel_b_small', ['1', '1 + t0**2'], ['3', ['t0**2 + x', 'x*t0'], 2]],
 ['prde_no_cancel_b_small', ['1', 't0**3 + 1'], ['t0', ['t0**2 + x', '2'], 1]],
 ['prde_no_cancel_b_large', ['1', '1 + t0**2'], ['t0**2', ['t0**3 + x', 't0'], 1]],
]
def build(D, args, kind):
    k = len(D) - 1
    T = [S.Symbol('x')] + [S.Symbol('t%d' % i) for i in range(k)]
    ns = {str(s): s for s in T}
    ex = lambda s: S.sympify(s, locals=ns)
    DE = risch.DifferentialExtension(extension={'D': [S.Poly(ex(d), T[i]) for i, d in enumerate(D)], 'T': T})
    P = lambda s: S.Poly(ex(s), T[-1])
    if kind == 'param_poly_rischDE':
        a, b, q, n = args
        return (P(a), P(b), [P(x) for x in q], n, DE), 2
    if kind in ('param_rischDE', 'limited_integrate'):
        fa, fd, G = args
        return (P(fa), P(fd), [(P(x), P(y)) for x, y in G], DE), 2
    b, q, n = args
    return (P(b), [P(x) for x in q], n, DE), 1
def ser(o):
    if isinstance(o, S.Poly):
        return ['Poly', S.srepr(o.as_expr()), [S.srepr(g) for g in o.gens], str(o.domain)]
    if isinstance(o, (list, tuple)):
        return [ser(x) for x in o]
    if hasattr(o, 'to_Matrix'):
        o = o.to_Matrix()
    return S.srepr(o)
def run(kind, D, args):
    a, qi = build(D, args, kind)
    try:
        r = getattr(prde, kind)(*a)
        out = ser(r)
        if kind.startswith('prde_no_cancel'):
            out = [out, ser(a[qi])]
        return out
    except Exception as e:
        return ['raise', type(e).__name__, str(e)]
res = []
for kind, D, args in CASES:
    run(kind, D, args)
    d0 = S.Dummy().dummy_index
    r = run(kind, D, args)
    res.append([r, S.Dummy().dummy_index - d0 - 1])
if M == 'symox':
    print(json.dumps(import_module('symox.integrals._prde_native').STATS), file=sys.stderr)
json.dump(res, sys.stdout)
