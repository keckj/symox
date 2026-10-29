import os, sys, time, resource
exec("from %s import *" % os.environ.get("LIB","sympy"))
x=symbols('x')
k=int(sys.argv[1]) if len(sys.argv)>1 else 3
t=time.perf_counter(); r=integrate(1/(x**5+x+1)**k, x); dt=time.perf_counter()-t
print('k',k,'time %.2f'%dt, 'maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024, 'len', len(str(r)), flush=True)
open(os.environ.get('OUT','/dev/null'),'w').write(srepr(r))
