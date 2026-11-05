import sys,time
lib=sys.argv[1]
if lib=='sympy':
    from sympy import *
else:
    sys.path.insert(0,'python'); from symox import *
x,y=symbols('x y')
for s in ["O(3**x,(x,oo))","O(x**x,(x,oo))","O(factorial(x),(x,oo))","O(2**x,(x,oo))","O(exp(1/x))","O(log(log(x)))","O(x*log(x),(x,oo))","O(x*sin(x),(x,oo))","O(sin(x),(x,oo))","O(exp(-x),(x,-oo))","O(exp(x),(x,oo))", "O(x**2*exp(x),(x,oo))","O(log(x)+1)","O(exp(1/x)*x)", "O(x+log(x),(x,oo))"]:
    t=time.perf_counter()
    try: r=eval(s)
    except Exception as e: r="EXC "+type(e).__name__+" "+str(e)[:80]
    print(s,"->",r, "%.3f"%(time.perf_counter()-t))
