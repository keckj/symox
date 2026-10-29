import os, sys, time
exec("from %s import *" % os.environ.get("LIB","sympy"))
x=symbols('x'); f=Function('f')
def nest(fn, n, e):
    for _ in range(n): e = fn(e)
    return e
C = {
 'exp10': lambda: solve(nest(exp, 10, x) - 2, x),
 'log10': lambda: solve(nest(log, 10, x) - 1, x),
 'sqrt10': lambda: solve(nest(sqrt, 10, x) - 2, x),
 'mixed10': lambda: solve(nest(lambda e: exp(e)+1, 10, x) - 5, x),
 'ratroot3': lambda: solve(expand((x - 2**41)*(x**2 + x + 1)), x),
 'ratroot5': lambda: solve(expand((x - 2**41)*(x**4 - 2)), x),
 'ratroot3q': lambda: solve(expand((3*x - 2**41)*(x**2 + 1)), x),
 'roots3': lambda: roots(expand((x - 2**41)*(x**2 + x + 1)), x),
 'ode3': lambda: dsolve(f(x).diff(x,3) - 2**41*f(x).diff(x,2) + f(x).diff(x) - 2**41*f(x)),
 'ode5': lambda: dsolve(f(x).diff(x,5) - 2**41*f(x).diff(x,4) - 2*f(x).diff(x) + 2**42*f(x)),
 'lam4096': lambda: solve(x**4096*exp(x) - 2, x),
 'lam4097': lambda: solve(x**4097*exp(x) - 2, x),
 'lam5001': lambda: solve(x**5001*exp(x) - 3, x),
 'lamlog4097': lambda: solve(x**4097 - exp(x), x),
}
k=sys.argv[1]
t=time.perf_counter()
try:
    r=C[k]()
    s=str(r); s = s if len(s)<300 else s[:300]+'...(%d)'%len(s)
except Exception as e:
    s='%s: %s'%(type(e).__name__, str(e)[:200])
print(k, '%.2fs'%(time.perf_counter()-t), s, flush=True)
