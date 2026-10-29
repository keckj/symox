import os, sys, time
exec("from %s import *" % os.environ.get("LIB","sympy"))
x=symbols('x')
kind=sys.argv[1]; n=int(sys.argv[2]); e=x
if kind=='mixed':
    for _ in range(n): e=exp(e)+1
    eq=e-5
elif kind=='exp':
    for _ in range(n): e=exp(e)
    eq=e-2
elif kind=='expp1':
    for _ in range(n): e=exp(e)
    eq=e+1
t=time.perf_counter(); r=solve(eq,x); print(kind, n, '%.3f'%(time.perf_counter()-t), str(r)[:200], flush=True)
