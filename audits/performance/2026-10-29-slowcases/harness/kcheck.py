import sys, time, resource
from symox import *
x=symbols('x'); k=int(sys.argv[1]); f=1/(x**5+x+1)**k
t=time.perf_counter(); r=integrate(f,x); dt=time.perf_counter()-t
print('k',k,'time %.2f'%dt,'maxrss MB',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024,'len',len(str(r)),'RootSum' in str(r), flush=True)
d=diff(r,x)-f
for pt in [Rational(1,3),Rational(-5,7),Integer(2),Rational(7,2)]:
    print(pt, N(d.subs(x,pt),30))
open(sys.argv[2],'w').write(srepr(r))
