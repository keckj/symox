from symox import *
x, a, b = symbols('x a b')
f = Function('f')
g = Function('g')
df = f(x).diff(x)
EQS = [
    (df - (x**2 + 1) - f(x)**2, f(x)),
    (df - 2*f(x)**2/x - 3/x**2, f(x)),
    (2*df + 3*f(x)**2 - 1/x**2, f(x)),
    (x*df - f(x)**2 - 3*f(x) + 4, f(x)),
    (df + f(x)**2 - f(x)/x + 1/x**2, f(x)),
    (df - f(x)**2 - a*x, f(x)),
    (df - a*f(x)**2 - b, f(x)),
    (df - x*f(x)**2 + f(x) - 2/x, f(x)),
    ((x**2 + 1)*df + f(x)**2 - x, f(x)),
    (df - exp(x)*f(x)**2 - 1, f(x)),
    (df - sqrt(2)*f(x)**2 - x, f(x)),
    (df - f(x)**2 - f(x)*sin(x), f(x)),
    (df - f(x)**3 - 1, f(x)),
    (df - f(x)**2 - f(2*x), f(x)),
    (df - Rational(1, 3)*f(x)**2 + x/(x**2 - 1), f(x)),
    (df*x**2 - x**2*f(x)**2 - 2, f(x)),
    (df + f(x)**2 - 2/(x**2 - 2*x), f(x)),
    (df - f(x)**2 + I*x, f(x)),
    (df - f(x)**2 - 1.5*x, f(x)),
    (df - f(x)**2 - oo, f(x)),
    (df - f(x)**2*g(x) - x, f(x)),
    (df**2 - f(x)**2, f(x)),
    (df - (f(x) + x)**2, f(x)),
    (df - f(x)**2/(x + f(x)) - 1, f(x)),
    (df - 3*f(x) - 4*f(x)**2 + 7/x**3, f(x)),
    (df + (x**2 + 2)*f(x)**2/(x - 1) - x, f(x)),
    (df - Integral(x, x)*f(x)**2, f(x)),
    (df - log(2)*f(x)**2 - x, f(x)),
    (df - pi*f(x)**2 - x, f(x)),
    (Eq(df, f(x)**2 + x**3), f(x)),
]
