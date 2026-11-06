"""mform.py MODULE INDICES OUT [SOURCE]: public forms of the M of
param_poly_rischDE calls (to_Matrix srepr, domain, getitem) of SOURCE
(default raw/nl_sympy.json.gz), run from the worktree root."""
import sys, json, gzip, importlib
mod = sys.argv[1]
S = importlib.import_module(mod)
prde = importlib.import_module(mod + '.integrals.prde')
risch = importlib.import_module(mod + '.integrals.risch')
src = sys.argv[4] if len(sys.argv) > 4 else 'audits/performance/2026-11-06-prde3/raw/nl_sympy.json.gz'
A = {r[0]: r for r in json.load(gzip.open(src, 'rt') if src.endswith('.gz') else open(src))}
out = {}
for k in map(int, sys.argv[2].split(',')):
    spec = A[k][2]; D = spec['D']; n = len(D) - 1
    T = [S.Symbol('x')] + [S.Symbol('t%d' % i) for i in range(n)]
    ns = {str(s): s for s in T}
    ex = lambda s: S.sympify(s, locals=ns)
    DE = risch.DifferentialExtension(extension={'D': [S.Poly(ex(d), T[i]) for i, d in enumerate(D)], 'T': T})
    P = lambda s: S.Poly(ex(s), T[-1])
    h, M = prde.param_poly_rischDE(P(spec['a']), P(spec['b']), [P(q) for q in spec['q']], spec['n'], DE)
    out[k] = [S.srepr(e) for e in M.to_Matrix()] + [str(M.domain), repr(M)[:0]]
    out[k].append([S.srepr(M[i, j].as_expr()) for i in range(M.rows) for j in range(M.cols)])
json.dump(out, open(sys.argv[3], 'w'))
