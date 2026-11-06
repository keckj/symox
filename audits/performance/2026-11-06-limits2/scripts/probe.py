import sys, json, signal
if sys.argv[1]=='sympy': from sympy import *
else: from symox import *
x=Symbol('x'); y=Symbol('y'); w=Symbol('w'); a=Symbol('a'); n=Symbol('n', integer=True)
p=Symbol('p', positive=True); k=Symbol('k', nonnegative=True); r=Symbol('r', real=True)
cases = [
 "(I/2)**x, x, oo", "(2*I)**x, x, oo", "(-2)**x, x, oo", "I**x, x, oo", "(1+I)**x, x, oo", "(I/2)**(x**2), x, oo",
 "(I/2)**(-x), x, oo", "2*(I/2)**x, x, oo", "x*(I/2)**x, x, oo", "(I/2)**x + 1, x, oo", "(I/2)**(1/x), x, 0, '+'",
 "a**x, x, oo", "p**x, x, oo", "(-p)**x, x, oo", "exp(I*x)/x, x, oo", "exp(a*x), x, oo", "exp(r*x), x, oo",
 "exp(x*y)*x, x, oo", "x**w, x, oo", "x**w, x, 0, '+'", "x**w, x, 0, '-'", "2*x**w, x, oo", "x**w*log(x), x, oo",
 "(x+1)**w - x**w, x, oo", "x**(w+1)/x**w, x, oo", "x**p, x, oo", "x**k, x, 0, '+'", "x**r, x, 0, '+'",
 "sin(x)**w, x, 0, '+'", "(1+1/x)**(w*x), x, oo", "x**y*exp(-x), x, oo", "exp(x)**w, x, oo", "exp(w*x), x, oo",
 "exp(x*(y+I)), x, oo", "(y*x)**n, x, oo", "x**n, x, oo", "x**n, x, -oo", "exp(I*y*x), x, oo", "log(x)**w, x, oo",
 "x**(1/w), x, oo", "(x**2+1)**w/x**(2*w), x, oo", "(n+1)**k/((n+1)**(k+1)-n**(k+1)), n, oo",
 "x**a/(x**a+1), x, oo", "x**I, x, oo", "x**(1+I)/x, x, oo", "(-x)**w, x, oo", "(3/x)**w, x, 0, '+'",
 "exp(-x)*x**w, x, oo", "exp(x)*x**w, x, oo", "Mod(x, a), x, a", "x**y, x, 1", "(x-1)**w, x, 1, '+'",
 "(x-1)**w*sin(x-1), x, 1, '+'", "(x-1)**w/(x-1), x, 1, '-'",
]
out={}
def h(*a): raise TimeoutError()
signal.signal(signal.SIGALRM, h)
for c in cases:
    signal.alarm(60)
    try:
        v=str(eval("limit(%s)" % c))
    except Exception as ex:
        v="EXC %s %s" % (type(ex).__name__, str(ex)[:80])
    signal.alarm(0)
    out[c]=v
json.dump(out, open(sys.argv[2],'w'), indent=0)
