import sys,time
lib=sys.argv[1]
if lib=='sympy':
    from sympy import *
else:
    sys.path.insert(0,'python'); from symox import *
x,y=symbols('x y')
for s in ["O(x + exp(1/x))","O(exp(1/x)+x**-20)","O(exp(1/x))+O(x)","O(exp(2/x))+O(x**3)","O(exp(2/x))+O(x**-3)","O(exp(1/x)+1)","O(log(x)+log(log(x)))",
 "O(exp(-x)+1/x**4,(x,-oo))","O(2**x + factorial(x), (x, oo))","O(2**x + x**x + factorial(x), (x, oo))","O(x + factorial(x), (x, oo))","O(2**x + x**x, (x, oo))","O(x**2 + x**x, (x, oo))","O(x**x + x**(-2), (x, oo))","O(2**x + 3**x, (x, oo))","O(x*log(x) + sin(x), (x, oo))","O(x*sin(x) + 1, (x, oo))",
 "O(x**2,(x,oo)).subs(x,cos(x))","O(x*y).subs(x,cos(x))"]:
    t=time.perf_counter()
    try: r=eval(s)
    except Exception as e: r="EXC "+type(e).__name__+" "+str(e)[:80]
    print(s,"->",r, "%.3f"%(time.perf_counter()-t))
