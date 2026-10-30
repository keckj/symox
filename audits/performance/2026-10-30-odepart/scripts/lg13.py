import sys, time
L = __import__(sys.argv[1])
x = L.Symbol('x'); f = L.Function('f'); fx = f(x)
t = time.perf_counter(); r = L.dsolve(fx.diff(x) + L.cos(x)*fx - L.exp(2*x), fx); print('%.1fs' % (time.perf_counter() - t), r, file=sys.stderr)
