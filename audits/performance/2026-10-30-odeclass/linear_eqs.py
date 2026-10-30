from symox import *
x, y, a, b, k = symbols('x y a b k')
f = Function('f'); g = Function('g')
fx = f(x); d1 = fx.diff(x); d2 = fx.diff(x, 2); d3 = fx.diff(x, 3); d4 = fx.diff(x, 4)
RHS = ["0", "exp(x)", "exp(-2*x)", "x**2 + 1", "sin(x)", "cos(3*x + 1)", "x*exp(x)*sin(x)",
       "sin(x)*cos(x)", "sin(x)**2", "log(x)", "1/x", "tan(x)", "exp(x**2)", "2**x", "3**(x/2 + 1)",
       "exp(x)*exp(2*x)", "sinh(x)", "cosh(2*x)*x**3", "exp(y*x)", "sin(a*x + b)", "x*y + 1",
       "sqrt(x)", "exp(x)/(x + 1)", "x**y", "sin(x**2)", "exp(sqrt(2)*x)", "exp(I*x)", "y",
       "sin(y)*cos(y)*exp(x)", "sin(y)*exp(x)", "exp(x*y + x*a)", "g(x)", "Derivative(g(x), x)",
       "sin(x)*exp(x) + x**2*cos(x)", "1.5*exp(x)", "exp(x)*log(2)", "(x + 1)**2", "exp(2)*exp(x)",
       "abs(x)", "Heaviside(x)", "x**Rational(3, 2)", "exp(-x)*sin(2*x)", "pi**x", "exp(x)**2",
       "sqrt(exp(x))", "exp(E*x)", "exp(x)*exp(y)", "atan(x)", "exp(sin(y))*exp(cos(y))*x"]
LHS = ["d2 + f(x)", "d2 - 3*d1 + 2*f(x)", "d3 - d1", "d4 + 2*d2 + f(x)", "d1 + 2*f(x)",
       "a*d2 + b*f(x)", "2*d2 + k*d1", "d2 + x*f(x)", "x**2*d2 + x*d1 + f(x)",
       "x**2*d2 - 2*x*d1 + 2*f(x)", "x**3*d3 + x*d1 - f(x)", "2*x**2*d2 + 3*x*d1 - f(x)",
       "x*d1 - f(x)", "x**2*d2 + f(x)**2", "d2*f(x) + d1", "d2 + 1.5*f(x)", "d2 + I*f(x)",
       "x**2*d2 + 3*d1", "(x**2 + 1)*d2 + x*d1", "x**2*(d2 + d1/x) - f(x)", "y*x**2*d2 + y*x*d1 + f(x)",
       "exp(x)*d2 + f(x)", "d2 + sin(f(x))", "d2 + d1**2", "x**2*d2 + x*d1 + (x**2 - 1)*f(x)",
       "d2 + 2*d1 + f(x)", "d3 + 3*d2 + 3*d1 + f(x)", "-d2 - f(x)", "x**2*d2/2 + x*d1/2 - f(x)",
       "x*d2 + d1", "x**2*d2 + x*d1 + a*f(x)"]
SRC = [l + " - (" + r + ")" for l in LHS for r in RHS[:25]] + ["d2 + f(x) - (" + r + ")" for r in RHS[25:]] + \
      ["x**2*d2 + x*d1 + f(x) - (" + r + ")" for r in RHS + ["log(x)**2*x", "x**3*y", "1/log(x)", "sin(log(x))", "x*log(x)", "log(x)/x**2", "x**Rational(1, 2)", "a*x**2 + b/x"]]
EQS = [(eval(s), fx) for s in SRC]
