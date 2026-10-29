import os
exec("from %s import *" % os.environ["LIB"])
x=symbols('x')
def mixed(n):
    e=x
    for _ in range(n): e=exp(e)+1
    return e-5
def nest(n):
    e=x
    for _ in range(n): e=exp(e)
    return e
cases=[exp(exp(exp(x)))+1, mixed(3), mixed(4), nest(4)-2, nest(5)+1, nest(3)-x*0-3, exp(exp(x))-exp(2)]
for c in cases: print(solve(c,x))
for e in [exp(sqrt(2)), exp(I*pi/3), exp(log(3)), exp(exp(1)), exp(2*I), sin(sqrt(3)), sin(pi/5), cos(log(2)), log(sqrt(2)+1), log(1+exp(1)), sqrt(exp(1)), Pow(pi,Rational(1,3)), Pow(sqrt(2)+1,Rational(1,2)), exp(1+I*pi)]:
    print(e, e.is_algebraic)
