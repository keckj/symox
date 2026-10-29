import time, sys, hashlib
import os
exec("from %s import *" % os.environ.get("LIB","sympy"))
x=symbols('x')
n=int(sys.argv[1])
t=time.perf_counter()
r=factor_list(x**n+x+1)
t1=time.perf_counter()-t
print(n, "time %.2f"%t1, [(Poly(p,x).degree(),e) for p,e in r[1]], hashlib.sha256(str(r).encode()).hexdigest()[:16], flush=True)
