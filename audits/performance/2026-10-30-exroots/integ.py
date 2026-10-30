import sys,time
import re as _re
if sys.argv[1]=='sympy':
    from sympy import *
else:
    from symox import *
x=symbols('x')
for e in ['1/(x**3-2)','1/(x**3+x+1)','1/(x**3+2*x-1)','1/(x**5+1)']:
    f=sympify(e); t0=time.perf_counter(); r=integrate(f,x); dt=time.perf_counter()-t0
    s=srepr(r); s=_re.sub(r"Dummy\('(\w+)', dummy_index=\d+\)", r"Dummy('\1')", s)
    print('%8.1fms'%(dt*1e3), e, '=>', r)
    print('SREPR', e, s)
