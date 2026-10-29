import os, sys, time
exec("from %s import *" % os.environ.get("LIB","sympy"))
x=symbols('x')
n=int(sys.argv[1])
t=time.perf_counter(); r=solve(x**n*exp(x)-2,x); print(n, '%.2f'%(time.perf_counter()-t), len(r), str(r)[:120], flush=True)
