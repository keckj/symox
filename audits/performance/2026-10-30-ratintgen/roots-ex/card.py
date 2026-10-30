import sys,time
if sys.argv[1]=='sympy':
    from sympy import *
else:
    from symox import *
x=symbols('x')
for e in ['1/(x**3-2)','1/(x**3+x+1)','1/(x**3+2*x-1)','1/(x**5+1)','x/(x**5+1)','1/(x**3+2)','x/(x**3-2)']:
    f=sympify(e); t0=time.perf_counter(); r=integrate(f,x); dt=time.perf_counter()-t0
    print('%8.1fms'%(dt*1e3), e, '=>', r)
