import sys
if sys.argv[1]=='sympy':
    from sympy import *
else:
    from symox import *
x,a,b=symbols('x a b')
s5,s2,s3=sqrt(5),sqrt(2),sqrt(3)
P=[ (Rational(3125,4) + 3125*s5/4)*x**4 + (-125*s5/4 - Rational(125,2))*x**2 + 5*s5/64 + Rational(75,64),
    (Rational(3125,4) - 3125*s5/4)*x**4 + (-Rational(125,2) + 125*s5/4)*x**2 + Rational(75,64) - 5*s5/64,
    x**2 + s2*x + 1, s3*x**2 - 2*x + s2, x**2/s2 - s3, x**4 + s2*x**2 + 1, x**4 - s3*x**2 - 1, x**3 + s2*x + 1,
    x**2 + (1+s2)*x + s3, x**2 - s5, s2*x - s3, x**4+x**2*s2/2+Rational(1,7), x**2+a*x+b, x**4+a*x**2+1,
    x**3-2*x, x**4-5*x**2+6, x**4+x**3-x, x**3+x+1, x**4+1, x**4+2*x**2+3, x**2+x*I+s2, x**4 - 2*s2*x**2 + 1,
    (x**2 + s2)*(x**2 - s3), x**3 + s2*x**2, x**4*s2 + 3*x**2 + s5]
for p in P:
    try: r=roots(p,x); m=roots(p,x,multiple=True)
    except Exception as e: r='ERR '+type(e).__name__; m=''
    print(str(p),'=>',r,'|',m)
