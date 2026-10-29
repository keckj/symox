import os, sys, time, json
from symox import *
x,y=symbols('x y'); a=Symbol('a',positive=True); n=symbols('n')
def nest(fn,k,e):
    for _ in range(k): e=fn(e)
    return e
cases=[
 'exp(x)-2','exp(x)+1','exp(2*x)-3*exp(x)+2','exp(exp(x))-2','exp(exp(x))+1','exp(exp(x)+1)-5',
 'exp(exp(exp(x)))+1','exp(exp(exp(x)))-2','exp(x**2)-3','exp(sin(x))-2','exp(-x)+1','exp(-exp(x))-2',
 'exp(exp(-x))-3','exp(2*exp(x))-4','exp(exp(2*x))-4','exp(exp(x)/2)-3','2**exp(x)-8','2**(2**x)-16',
 'exp(log(x)+1)-2','log(exp(x)+1)-2','log(log(x))-1','exp(x)*x-1','exp(exp(x))*exp(x)-1','exp(exp(x)+x)-1',
 'exp(a*exp(x))-2','exp(y*exp(x))-2','exp(exp(x)+y)-y','exp(x+exp(1))-3','exp(I*exp(x))-1','exp(exp(I*x))-2',
 'exp(exp(x))-exp(2)','exp(exp(x)**2)-4','exp(sqrt(exp(x)))-3','sqrt(exp(exp(x)))-2','exp(exp(x))**2-4',
 'exp(exp(x)-1)-exp(1)','exp(3*exp(x))-8','exp(exp(x)/3)+1','cosh(exp(x))-2','exp(exp(x))+exp(-exp(x))-3',
 'exp(log(2)*exp(x))-4','exp(exp(x))-x','exp(exp(exp(x)+1)+1)-5','exp(exp(exp(exp(x))))+1','exp(x)-y',
 'exp(exp(x))-y','exp(exp(y*x))-2','exp(exp(x+y))-3','3**exp(x)-9','exp(2**x)-3','exp(exp(x))-1',
]
out={}
for c in cases:
    t=time.perf_counter()
    try: r=srepr(solve(eval(c),x))
    except Exception as e: r='ERR %s: %s'%(type(e).__name__,str(e)[:100])
    out[c]=(r,time.perf_counter()-t)
json.dump(out,open(sys.argv[1],'w'),indent=0)
