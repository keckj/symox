import time, sys; sys.set_int_max_str_digits(0)
import os
exec("from %s import *" % os.environ.get("LIB","sympy"))
x=symbols('x')
f=expand((x+1)**200); g=expand((x-2)**150+1)
t=time.perf_counter()
r,prs=resultant(f,g,x,includePRS=True)
t1=time.perf_counter()-t; print("t",t1,flush=True)
import hashlib
h=hashlib.sha256(str(prs).encode()).hexdigest()[:16]
print("time %.2f len %d bits %d hash %s rhash %s"%(t1,len(prs),int(r).bit_length(),h,hashlib.sha256(str(r).encode()).hexdigest()[:16]))
print([Poly(p,x).degree() for p in prs][:10], sum(len(str(p)) for p in prs))
