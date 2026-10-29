import sys, time, signal, json
from symox import *
from symox.integrals.rationaltools import ratint
x=symbols('x')
cases=["1/(x**3+x+1)","1/(x**3-2)","1/(x**4+1)","1/(x**5+1)","1/(x**6-2)","1/(x**8+1)","x**2/(x**4+x+1)"]
for n in range(3,9):
    for a in (1,-1,2,-2,3,4,-4,9,16,1/S(4)):
        cases.append(f"1/(x**{n}+({a}))")
cases+=["x/(x**4+1)","x**2/(x**4+1)","x**2/(x**6+1)","(x**2+1)/(x**4+1)","(x**2-1)/(x**4+1)","1/(x**4+x**2+1)",
 "x/(x**4+x**2+1)","1/(x**6+x**3+1)","1/(x**5-x+1)","1/((x**2+1)*(x**2+2))","1/((x**2+1)*(x**2+4)*(x-1))",
 "(3*x**3+2*x+5)/(x**4+3*x**2+1)","1/(x**4+2*x**2+5)","1/(x**4-x**2+1)","x**3/(x**4+4)","1/(x**3+3*x+1)",
 "1/(x**4+6*x**2+1)","(x+1)/(x**3-x**2+x-1)","1/(x**4+x**3+x**2+x+1)","1/(x**6-1)","1/(x**6+8)","x**4/(x**6+1)",
 "1/(x**4+x+1)","1/(x**3-3*x+1)","1/(x**4-4*x**2+2)","1/(x**4+8*x**2+4)","(2*x+1)/(x**4+2*x**3+3*x**2+2*x+2)"]
def handler(*a): raise TimeoutError
signal.signal(signal.SIGALRM, handler)
out={}
for w in cases:
    f=sympify(w)
    signal.alarm(60)
    try:
        t0=time.perf_counter(); r=ratint(f,x); dt=time.perf_counter()-t0
        out[w]=(srepr(r),dt)
    except TimeoutError:
        out[w]=("TIMEOUT",60)
    except Exception as e:
        out[w]=("ERR "+type(e).__name__+" "+str(e)[:100],0)
    signal.alarm(0)
json.dump(out,open(sys.argv[1],'w'))
print(len(out))
