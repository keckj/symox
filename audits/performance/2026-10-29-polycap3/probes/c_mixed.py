import os, sys, time
exec("from %s import *" % os.environ.get("LIB","sympy"))
x=symbols('x')
n=int(sys.argv[1]); e=x
for _ in range(n): e=exp(e)+1
t=time.perf_counter(); r=solve(e-5,x); print(n, '%.2f'%(time.perf_counter()-t), str(r)[:150], flush=True)
