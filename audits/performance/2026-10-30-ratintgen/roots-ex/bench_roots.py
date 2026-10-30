import sys,time,json
from symox import *
x=symbols('x'); s2,s3,s5=sqrt(2),sqrt(3),sqrt(5)
cases=[('integrate','1/(x**3+x+1)'),('integrate','1/(x**5+1)'),('integrate','x/(x**5+1)'),('integrate','1/(x**3-2)'),
 ('roots','sqrt(3)*x**2 - 2*x + sqrt(2)'),('roots','x**4 - sqrt(3)*x**2 - 1'),('roots','x**4 + 2*x**2 + 3'),('roots','x**4 - 5*x**2 + 6'),('roots','sqrt(2)*x**4 + 3*x**2 + sqrt(5)')]
out={}
for kind,e in cases:
    f=sympify(e); ts=[]
    for i in range(5):
        clear_cache() if 'clear_cache' in globals() else None
        t0=time.perf_counter(); integrate(f,x) if kind=='integrate' else roots(f,x,multiple=True); ts.append(time.perf_counter()-t0)
    ts.sort(); out[kind+'('+e+')']=round(ts[2]*1e3,2)
print(json.dumps(out))
