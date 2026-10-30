"""Gap cases of the limgaps lot: python cases.py sympy|symox (from the worktree root)."""
import sys, time, signal
lib = sys.argv[1]
if lib == 'sympy':
    from sympy import *
    from sympy.calculus.accumulationbounds import AccumBounds
else:
    sys.path.insert(0, 'python')
    from symox import *
x, y = symbols('x y')
p = Symbol('p', positive=True)
CASES = [
 "limit(((x**(x+1) + (x+1)**x)/x**(x+1))**x, x, oo)",
 "frac(x).series(x, oo)",
 "frac(x).series(x, oo, 3)",
 "limit(factorial(x)/x, x, -oo)",
 "limit(factorial(x), x, -oo)",
 "limit(loggamma(x), x, -oo)",
 "limit(1 + O(1/x), x, 1)",
 "O(exp(1/x) + 1)",
 "O(log(x) + log(log(x)))",
 "O(x**2, (x, oo)).subs(x, cos(x))",
 "O(x**2, (x, 1)).subs(x, cos(x))",
 "O(x*y).subs(x, cos(x))",
 "O(1, x).diff(x)",
 "limit(Abs(x + y)**2, x, 0)",
 "floor(1/p).series(p, 0)",
 "exp(floor(x)).series(x, oo)",
 "limit(cos(x)*floor(x)/x, x, oo)",
 "floor(1/x).nseries(x, 0, 3)",
 "sin(Piecewise((x, x>0), (-x, True))).series(x, 0, 4)",
 "limit(gamma(x)/x, x, -oo)",
 "limit(1/factorial(x), x, -oo)",
 "limit(x*loggamma(x), x, -oo)",
 "limit(sin(x)**2, x, oo)",
 "limit(cos(x)**3 + 1, x, oo)",
 "limit(sin(x)**2*cos(x), x, oo)",
]
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TimeoutError()))
for c in CASES:
    signal.alarm(60)
    t = time.perf_counter()
    try:
        r = repr(eval(c))
    except Exception as e:
        r = 'EXC %s: %s' % (type(e).__name__, str(e)[:80])
    signal.alarm(0)
    print('%-60s -> %s   [%.1f ms]' % (c, r, 1000*(time.perf_counter()-t)))
