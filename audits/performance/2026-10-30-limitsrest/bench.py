"""Limits and series public-call benchmark.

usage: python bench.py sympy|symox CALL [repeats]
Prints one JSON line: first-call time and the median of the repeated calls
(seconds), plus the str of the result (for the cross-check). SymPy caches
results, so the repeated calls of both libraries are warm calls."""
import importlib, json, statistics, sys, time, warnings

lib = importlib.import_module(sys.argv[1])
x, w, a, b = lib.symbols('x w a b')
L, S, oo, I = lib.limit, lib.series, lib.oo, lib.I
CALLS = {
    'lim_pow4': lambda: L(2**(-x)*lib.sqrt(4**(x + 1) + 1), x, oo),
    'lim_pow4b': lambda: L(4**x/2**(2*x) + 1/x, x, oo),
    'lim_ei': lambda: L((lib.Ei(x - lib.exp(-lib.exp(x))) - lib.Ei(x))*lib.exp(-x)*lib.exp(lib.exp(x))*x, x, oo),
    'lim_ugamma': lambda: L(lib.uppergamma(x, 1)/lib.gamma(x), x, oo),
    'lim_stirling': lambda: L(lib.log(lib.factorial(x)) - (x*lib.log(x) - x + lib.log(2*lib.pi*x)/2), x, oo),
    'lim_sinx': lambda: L(lib.sin(x)/x, x, 0),
    'ser_abs_cdir': lambda: lib.Abs(x).series(x, 0, cdir=-1),
    'ser_log_cdir': lambda: lib.log(x).series(x, 0, cdir=-1),
    'ser_abs_cdirI': lambda: lib.Abs(x).series(x, 0, cdir=I),
    'ser_nonsym': lambda: S(lib.exp(a*b), a*b, n=8),
    'ser_logx': lambda: lib.log(-lib.log(x)).nseries(x, logx=-x),
    'ser_sin': lambda: lib.sin(x).series(x, 0, 10),
    # round 2
    'lim_xlog3': lambda: L((1 + x**lib.log(3))**(1/x), x, 0),
    'lim_expsqrt3': lambda: L((lib.exp(x**lib.sqrt(3)) - 1)/x, x, 0),
    'ns_floor_cdir': lambda: (x*lib.floor(x)).nseries(x, 0, 3, cdir=-1),
    'ns_log_cdir': lambda: lib.log(x**2 - x**3).nseries(x, 0, 3, cdir=-1),
    'lim_absxy': lambda: L((lib.Abs(x + a) - lib.Abs(x - a))/(2*x), x, 0),
    'lim_fact_neg': lambda: L(lib.factorial(x), x, -oo),
    'lim_order_1': lambda: L(1 + lib.O(x**2), x, 1),
    'ord_subs_cos': lambda: lib.O(x**2).subs(x, lib.cos(x)),
    'ord_diff': lambda: lib.O(x**2).diff(x),
    'ord_loglog': lambda: lib.Order(lib.log(x) + 1/lib.log(x)),
    'ser_floor_oo': lambda: lib.floor(x).series(x, oo),
}
name = sys.argv[2]
repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 5
f = CALLS[name]
warnings.simplefilter('ignore')
t0 = time.perf_counter()
try:
    out = str(f())
except Exception as e:  # noqa
    out = 'ERR:%s:%s' % (type(e).__name__, e)
first = time.perf_counter() - t0
times = []
for _ in range(repeats):
    t0 = time.perf_counter()
    try:
        f()
    except Exception:  # noqa
        pass
    times.append(time.perf_counter() - t0)
print(json.dumps({'call': name, 'lib': sys.argv[1], 'first': first, 'median': statistics.median(times),
                  'result': out[:2000]}))
