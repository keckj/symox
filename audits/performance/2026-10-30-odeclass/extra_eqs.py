from symox import *
x, y, a, b = symbols('x y a b')
f = Function('f')
fx = f(x); df = fx.diff(x)
SRC = [
 "df - x*fx", "-df + 2*fx", "3*df - 2*x*fx**2", "-Rational(3,2)*df + x*fx - 1",
 "exp(x)*df - exp(fx)", "3*exp(x)*exp(fx)*df - x", "exp(x + fx)*df - 1", "E*df - fx",
 "df - exp(-x**2) + 2*x*fx", "df + fx - exp(x)*fx**2", "-df - exp(x)*fx + x",
 "sqrt(2)*df - fx", "-sqrt(2)*x*df + fx**2", "2**x*df - fx", "-x*2**x*df + 3*fx",
 "df*(fx + x) - fx + x", "(x*fx + 1)*df + fx**2", "df**2 - fx", "df*fx**2 + x*fx",
 "df - fx**Rational(1,3)*x", "df + a*fx - b*fx**a", "df - fx**2 - fx/x - 1/x**2",
 "2*df + 3*fx**2 + fx/x - 5/x**2", "df - log(fx)*x", "x*df - fx + x*exp(fx/x)",
 "df + fx*cos(x) - sin(2*x)/2", "df - (fx + x)/(fx - x)", "(3*x + 2*fx + 1)*df + x - fx",
 "df - x/fx + fx/x", "df + 2*fx/x - x**2*fx**3", "df - exp(x)*exp(fx)", "df*exp(-fx) - x",
 "df + x*fx - x*fx**2*exp(x)", "fx*df - x*exp(x**2 - fx**2)", "df - a", "df - fx", "df",
 "x*df - 2*x - 2*fx", "df - (1 - fx**2)", "(1 + x**2)*df - 1 - fx**2", "-2*df + 4*x*fx",
 "df - 2*x*exp(x)*fx/(exp(x) + 1)", "cos(x)*df + sin(x)*fx - 1", "df + fx/(x*log(x)) - 1/x",
]
EQS = [(eval(s), fx) for s in SRC]
