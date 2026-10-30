import sys,time
lib=sys.argv[1]
if lib=='sympy':
    from sympy import *
else:
    sys.path.insert(0,'python'); from symox import *
x=Symbol('x'); p=Symbol('p',positive=True); r=Symbol('r',real=True)
for s in ["frac(x).series(x,oo)","frac(1/x).series(x)","frac(1/x).series(x,0,3)","frac(x).series(x,-oo)","frac(-1/x).series(x)","frac(1/x).nseries(x,0,3)","frac(1/p).series(p)","frac(x**2+x).series(x,oo)","frac(1/x).series(x,0,0)","(2*frac(1/x)).series(x)","(frac(1/x)+1).series(x)","(frac(1/x)+x).series(x)","frac(1/x).series(x,dir='-')","frac(exp(1/x)).series(x)", "frac(1/x+1).series(x)",
"floor(1/p).series(p)","ceiling(1/p).series(p)","floor(p).series(p,oo)","floor(1/p+p).series(p)","floor(1/p).nseries(p,0,3)","floor(1/x).nseries(x,0,3)","ceiling(1/x).nseries(x,0,3)","floor(1/sin(p)).series(p)","floor(1/r).series(r)","floor(1/x).series(x)","floor(1/p).series(p,0,0)", "floor(-1/p).series(p)"]:
    t=time.perf_counter()
    try: v=eval(s)
    except Exception as e: v="EXC "+type(e).__name__+" "+str(e)[:60]
    print(s,"->",v)
