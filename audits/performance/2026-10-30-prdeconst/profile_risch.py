"""cProfile of 5 repeated risch_integrate calls (after one warm-up call).

usage: python profile_risch.py <symox|sympy> <integrand> <function prefix>...
prints the cumulative and own time of each matching function and its main callees.
"""
import sys, cProfile, pstats, importlib, time
mod=sys.argv[1]; expr=sys.argv[2]
S=importlib.import_module(mod); risch=importlib.import_module(mod+'.integrals.risch')
x=S.Symbol('x'); f=S.sympify(expr, locals={'x':x})
risch.risch_integrate(f,x)
pr=cProfile.Profile(); pr.enable()
for _ in range(5): risch.risch_integrate(f,x)
pr.disable()
st=pstats.Stats(pr)
def short(k): return '%s:%d(%s)'%(k[0].split('/')[-1],k[1],k[2])
callers={}
for k,(cc,nc,tt,ct,cl) in st.stats.items():
    for c,v in cl.items(): callers.setdefault(c,[]).append((k,v))
for target in sys.argv[3:]:
    for k in st.stats:
        if short(k).startswith(target):
            print('==',short(k),'cum %.3f tot %.3f'%(st.stats[k][3],st.stats[k][2]))
            for callee,v in sorted(callers.get(k,[]),key=lambda a:-a[1][3])[:8]:
                print('   %-60s n=%d cum=%.3f'%(short(callee),v[1],v[3]))
