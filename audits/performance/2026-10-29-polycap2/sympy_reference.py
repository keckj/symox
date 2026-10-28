import time
from sympy import *
x,y=symbols('x y')
R=Rational
print("r1", resultant(x**3/2 + x - R(1,3), R(3,5)*x**2 - 7, x))
print("r2", resultant(3*x**4 - 2*x + 5, 7*x**3 + x**2 - 4, x))
print("r3", resultant(x**2-1, x-1, x), resultant(2*x+1, 4*x**2-1, x))
print("r4", resultant(x**2+y, x*y-1, x))
print("prs", resultant(x**3 - 2*x + 1, 2*x**2 + 3, x, includePRS=True))
print("poly", Poly(x**5+3*x+1, x).resultant(Poly(x**2-2, x)), Poly(x**5+3*x+1, x).resultant(Poly(x**2-2, x), includePRS=True))
print("poly_qq", Poly(x**5/3+3*x+1, x).resultant(Poly(x**2-R(1,2), x)))
print("disc", discriminant(x**2+3*x+1), discriminant(x**3/2-x+R(1,7)), discriminant(5*x**5 - 3*x**2 + 7), discriminant(4*x**4 - 1, x), discriminant(x**2-2*x+1))
print("discpoly", Poly(x**3-2,x).discriminant(), Poly(x**3/3-2,x).discriminant())
f=expand((x+1)**60); g=expand((x-2)**45+1)
print("big", resultant(f,g,x) == (3**45-1)**60)
print("type", type(resultant(x**2+1, x-1, x)), type(Poly(x**2+1).resultant(Poly(x-1))))
print("polys", resultant(x**2+1, x-1, x, polys=True), type(resultant(x**2+1, x-1, x, polys=True)))
c=Poly(cyclotomic_poly(105,x),x).all_coeffs(); print("c105", c)
for n in [20001, 9699690//19//17]:
    t=time.time(); p=Poly(cyclotomic_poly(n,x),x); cs=p.all_coeffs(); print("cyc",n,time.time()-t, p.degree(), max(abs(v) for v in cs), sum(int(v)*(i+1) for i,v in enumerate(cs)) % (2**61-1))
for e in [x**4100 + x**2, (x**3+x)**1400 if False else x**6000 + 2*x**3 + 1, expand((x**2+x+1)**3)**1]:
    t=time.time(); print("dec", time.time()-t, decompogen(e, x))
t=time.time(); print("dec8", decompogen(x**8194 + x**2, x), time.time()-t)
