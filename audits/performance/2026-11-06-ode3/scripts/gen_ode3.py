import sys
from sympy import *
x, y = symbols('x y')
f = Function('f')
fx = f(x); d1 = fx.diff(x); d2 = fx.diff(x, 2); d3 = fx.diff(x, 3)
ENV = dict(globals())
CASES = [
 ('fx*(d1**2 - 1)', 'nth_algebraic'),
 ('d2 - d1/x', 'nth_algebraic'),
 ('d1*d2', 'nth_algebraic'),
 ('d3 + d2 - d1', None),
 ('(d1 - x)*(d2 + 1)', 'nth_algebraic'),
 ('d2*(fx - 1)', 'nth_algebraic'),
 ('(fx**2 - 4)*(d1 + fx)', 'factorable'),
 ('fx*(d1 - x)', 'factorable'),
 ('d1**2 - 1', 'factorable'),
 ('(d1 - y)**2', 'factorable'),
 ('fx**2*d1', 'factorable'),
 ('x**2*d2 + x*d1 + (x**2 - 4)*fx', '2nd_linear_bessel'),
 ('4*x**2*d2 + 8*x*d1 + (x - 3)*fx', '2nd_linear_bessel'),
 ('(x - 1)**2*d2 + (x - 1)*d1 + ((x - 1)**2 - 1)*fx', '2nd_linear_bessel'),
 ('x**2*d2 + x*d1 + (y**2*x**2 - 1)*fx', '2nd_linear_bessel'),
 ('x*d2 + d1 + x*fx', '2nd_linear_bessel'),
 ('d2 + fx/x', '2nd_linear_bessel'),
]
for src, hint in CASES:
    eq = eval(src, ENV)
    c = classify_ode(eq, fx)
    out = [src, hint, str(c)]
    if hint:
        try: out.append(str(dsolve(eq, fx, hint=hint)))
        except Exception as e: out.append('ERR %s' % type(e).__name__)
    print(repr(tuple(out)) + ',')
