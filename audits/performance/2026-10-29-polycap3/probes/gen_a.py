import os, sys, hashlib
sys.set_int_max_str_digits(0)
exec("from %s import *" % os.environ.get("LIB","sympy"))
x,y=symbols('x y')
cases=[
 ("x**4 + 1", "x**2 + 1"),
 ("x - 1", "x**3 + 2"),
 ("x**3 + 2", "x - 1"),
 ("x**3 - 2*x", "x**5 + 3"),
 ("2*x**3 + 3*x - 7", "4*x**2 - 6"),
 ("x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5", "3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21"),
 ("x**6 - 1", "x**3 - 1"),
 ("(x - 1)**3*(x + 2)", "(x - 1)**2*(x**2 + 5)"),
 ("7*x**5 - 3", "2*x**2"),
 ("x**10 - 2*x**3 + 1", "x**7 + 5*x - 11"),
]
for a,b in cases:
    f=expand(sympify(a)); g=expand(sympify(b))
    r,prs=resultant(f,g,x,includePRS=True)
    print(repr((a,b,str(r),[str(p) for p in prs])))
for n,m in [(40,30),(60,45)]:
    f=expand((x+1)**n); g=expand((x-2)**m+1)
    r,prs=resultant(f,g,x,includePRS=True)
    print(n,m,len(prs),hashlib.sha256(str(prs).encode()).hexdigest(),hashlib.sha256(str(r).encode()).hexdigest())
    s=subresultants(f,g,x)
    print("sub",hashlib.sha256(str(s).encode()).hexdigest())
