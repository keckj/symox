"""Extra fps cases beyond test_formal.py, replayed with SymPy or symox.

extra.py LIB OUT

Each case is an expression string evaluated with ``from LIB import *`` and
``x, y, a, n`` symbols (``n`` integer). Recorded per case: the str of the
value (``EXC <class>`` on exception) and the elapsed time. For a series
value, ``truncate(8)`` is recorded instead of the series node; a top-level
``fps(...)`` call that returns its input is recorded as ``EXC HELD <str>``
(not a truth: SymPy failing the same way gives SAME). Each case runs
under a 30 s alarm. The output keys and format match ``tl.py`` so that
``compare.py`` classifies them (truth = value recorded without exception).
"""
import ast, json, signal, sys, time, warnings
warnings.simplefilter("ignore")
sys.setrecursionlimit(10000)

LIB, OUT = sys.argv[1:3]

CASES = [
    # elementary functions
    "fps(exp(x))", "fps(sin(x))", "fps(cos(x))", "fps(tan(x))", "fps(log(1 + x))",
    "fps(atan(x))", "fps(asin(x))", "fps(acos(x))", "fps(sinh(x))", "fps(cosh(x))",
    "fps(atanh(x))", "fps(asinh(x))", "fps(sqrt(1 + x))", "fps((1 + x)**Rational(1, 3))",
    "fps(exp(2*x))", "fps(exp(-x**2))", "fps(sin(x)**2)", "fps(sin(x)*cos(x))",
    "fps(erf(x))", "fps(Si(x))", "fps(cos(sqrt(x)))", "fps(exp(x) - 1)",
    "fps(log(1 - x)/x)", "fps(1/sqrt(1 - x**2))", "fps(x*exp(x)**2)",
    # rational functions
    "fps(1/(1 - x))", "fps(1/(1 + x**2))", "fps((x + 1)/(x**2 - 5*x + 6))",
    "fps(x/(1 - x - x**2))", "fps(1/(x**2 + x + 1))", "fps(1/(x**2 + x + 1), full=True)",
    "fps(1/(1 - x)**3)", "fps((2*x + 3)/(x**3 - 1))", "fps(1/(x**2 + 1), full=True)",
    "fps(1/(x**4 + 1), full=True)", "fps(x**3/(1 - x))", "fps(1/(1 - x) + x**2)",
    # products
    "fps(exp(x)*sin(x))", "fps(x*exp(x))", "fps(x**2*cos(x))", "fps(exp(x)*log(1 + x))",
    "fps(exp(x)/(1 - x))", "fps(sin(x)/x)", "fps((1 + x)*log(1 + x))",
    "fps(exp(x)).product(fps(sin(x))).truncate(8)",
    "fps(exp(x)).product(fps(cos(x)), x, 10).truncate(8)",
    "fps(sin(x)).product(fps(1/(1 - x))).truncate(8)",
    "fps(log(1 + x)).product(fps(exp(-x))).truncate(6)",
    # compositions
    "fps(exp(sin(x)))", "fps(sin(sin(x)))", "fps(log(cos(x)))", "fps(exp(exp(x)))",
    "fps(exp(x)).compose(fps(sin(x))).truncate(8)",
    "fps(cos(x)).compose(fps(log(1 + x))).truncate(6)",
    "fps(1/(1 - x)).compose(fps(x*exp(x))).truncate(6)",
    "fps(exp(x)).inverse().truncate(8)", "fps(cos(x)).inverse().truncate(8)",
    "fps(1 + x + x**2).inverse().truncate(6)",
    # operations
    "(fps(exp(x)) + fps(sin(x))).truncate(8)", "(fps(exp(x)) - fps(log(1 + x))).truncate(8)",
    "(3*fps(sin(x))).truncate(8)", "(fps(sin(x)) + 1).truncate(8)", "(-fps(cos(x))).truncate(8)",
    "fps(sin(x)).diff(x).truncate(8)", "fps(exp(x)).integrate().truncate(8)",
    "fps(log(1 + x)).integrate(x).truncate(8)", "fps(atan(x)).diff(x).truncate(8)",
    "fps(sin(x)).polynomial(8)", "fps(exp(x)).infinite", "fps(cos(x)).as_leading_term(x)",
    "fps(sin(x))[0:6]", "fps(exp(x)).zero_coeff()",
    # nonzero points
    "fps(sin(x), x0=1)", "fps(exp(x), x0=2)", "fps(log(x), x0=1)", "fps(1/x, x0=1)",
    "fps(x**2 + 1, x0=1)", "fps(cos(x), x0=pi)", "fps(1/(1 - x), x0=2)", "fps(atan(x), x0=1)",
    "fps(exp(x), x0=-1, dir=-1)",
    # x0 = oo, -oo
    "fps(exp(1/x), x0=oo)", "fps(1/(1 + x), x0=oo)", "fps(log(1 + 1/x), x0=oo)",
    "fps(x/(x + 1), x0=oo)", "fps(sin(1/x), x0=-oo)", "fps(atan(x), x0=oo)",
    "fps(log(1 + x), x0=oo)", "fps(1/(x**2 + 1), x0=-oo)",
    # dir
    "fps(sin(x), dir=-1)", "fps(log(1 + x), dir=-1)", "fps(exp(x), x0=1, dir=-1)",
    "fps(1/(1 - x), dir=-1)", "fps(sin(x), dir='-')", "fps(sin(x), dir='+')",
    "fps(sin(x), dir=2)", "fps(x, dir='-')", "compute_fps(sin(x), x, dir='-')[2]",
    # hyper / rational / order / full
    "fps(sin(x), hyper=False)", "fps(log(1 + x), hyper=False)", "fps(atan(x), hyper=False)",
    "fps(exp(x), rational=False)", "fps(log(1 + x), rational=False)",
    "fps(sin(x), hyper=False, rational=False)", "fps(exp(x), order=1)", "fps(sin(x), order=1)",
    "fps(sin(x), order=2)", "fps(cos(x)**2, order=2)", "fps(log(1 + x**2), full=True)",
    "fps(atan(x), full=True)", "fps(1/(x**3 - 2), full=True)",
    # symbolic
    "fps(x**n*sin(x))", "fps(exp(a*x), x)", "fps(1/(1 - a*x), x)", "fps(sin(a*x), x)",
    "fps(x**a*exp(x), x)", "fps(y*x)", "fps(x + y, x)", "fps(3)", "fps(x)", "fps(x**3 + 2*x)",
    # failures (returned unchanged)
    "fps(1/x)", "fps(sqrt(x))", "fps(log(x))", "fps(exp(1/x))", "fps(abs(x))", "fps(gamma(x))",
]


class Alarm(Exception):
    pass


def _alarm(*a):
    raise Alarm()


signal.signal(signal.SIGALRM, _alarm)


def main():
    ns = {}
    exec("from %s import *\nx, y, a = symbols('x y a')\nn = symbols('n', integer=True)" % LIB, ns)
    try:
        exec("from %s.series.formal import compute_fps, FormalPowerSeries" % LIB, ns)
    except ImportError:
        pass
    fps_class = ns['FormalPowerSeries']
    rec = {}
    for i, case in enumerate(CASES):
        t = time.perf_counter()
        signal.alarm(30)
        ok = True
        try:
            v = eval(case, ns)
            top = ast.parse(case, mode='eval').body
            held = (isinstance(top, ast.Call) and getattr(top.func, 'id', None) == 'fps'
                    and not isinstance(v, fps_class))
            if isinstance(v, fps_class):
                v = v.truncate(8)
            if isinstance(v, (list, tuple)):
                r = type(v).__name__ + "(" + ", ".join(map(str, v)) + ")"
            else:
                r = str(v)
            if held:
                # fps returned its input: recorded like a decline
                r, ok = "EXC HELD " + r, False
        except Alarm:
            r, ok = "TIMEOUT", False
        except BaseException as e:
            r, ok = "EXC " + type(e).__name__, False
        finally:
            signal.alarm(0)
        rec["extra:%03d %s" % (i, case)] = [r, ok, round(time.perf_counter() - t, 4), i]
    json.dump(rec, open(OUT, "w"), indent=0, sort_keys=True)


main()
