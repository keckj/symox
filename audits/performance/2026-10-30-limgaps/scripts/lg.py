import sys
lib=sys.argv[1]
if lib=='sympy':
    from sympy import *
else:
    sys.path.insert(0,'python'); from symox import *
x=Symbol('x'); y=Symbol('y',positive=True)
for s in ["loggamma(x)","x*loggamma(x)","loggamma(x)/x","loggamma(x)+x","2*loggamma(x)","loggamma(x)**2","1/loggamma(x)","loggamma(2*x)","loggamma(x)+loggamma(x+1)","loggamma(-x)","exp(loggamma(x))","loggamma(x)-x*log(-x)","loggamma(x**2)","loggamma(-x**2)", "y*loggamma(x)", "loggamma(x)+3"]:
    try:
        r=limit(eval(s),x,-oo)
    except Exception as e:
        r="EXC "+type(e).__name__
    print(s, "->", r)
