import sys,time
lib=sys.argv[1]
if lib=='sympy':
    from sympy import *
else:
    sys.path.insert(0,'python'); from symox import *
x,y=symbols('x y')
for s in ["O(exp(1/x))+1","O(sin(1/x))","O(sin(x)**2,(x,oo))","O(x*sin(x)**2,(x,oo))","O(cos(x)/x,(x,oo))","O(exp(1/x)*(x+1))",
 "O(1/sin(x),(x,oo))","O(sin(x)+cos(x),(x,oo))","O(exp(1/x))+O(exp(2/x))","O(x*exp(1/x))+O(exp(1/x))","O(exp(1/x))+O(1/x)",
 "O(sin(x),(x,oo))+O(1/x,(x,oo))","O(exp(-1/x))+O(x)","O(x)+O(exp(-1/x))","O(x+exp(-1/x))","O(sin(2*x+1)*x**2,(x,oo))",
 "O(x**2*exp(1/x)+x)","O(exp(1/x)/x+exp(1/x))","O(sin(x)*exp(x),(x,oo))","O(log(x)*sin(1/x))","O(exp(1/x)+x**-20, (x, 0))+O(x**-30)",
 "O(x**3+x**x,(x,oo))","O(x*factorial(x)+x**2,(x,oo))","O(y*exp(1/x),x)+O(x,x)","O(x*sin(x)+x**2,(x,oo))","O(tan(x),(x,oo))",
 "O(sin(x)**-2,(x,oo))","O(sin(1/x)+x)","O(sin(1/x)/x+1)"]:
    t=time.perf_counter()
    try: r=eval(s)
    except Exception as e: r="EXC "+type(e).__name__+" "+str(e)[:80]
    print(s,"->",r)
