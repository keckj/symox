"""Native formal power series (``sympy/series/formal.py``) against live SymPy.

Each case is evaluated with SymPy 1.14 and with symox; the str of the value
(``truncate(8)`` for a series, the exception class name on failure) must
agree. A top-level ``fps`` call that returns its input must do so in both.
"""
import pytest
import sympy as sp
import symox as sx

CASES = [
    # elementary functions
    "fps(exp(x))", "fps(sin(x))", "fps(cos(x))", "fps(tan(x))", "fps(log(1 + x))",
    "fps(atan(x))", "fps(asin(x))", "fps(acos(x))", "fps(sinh(x))", "fps(cosh(x))",
    "fps(atanh(x))", "fps(asinh(x))", "fps(sqrt(1 + x))", "fps((1 + x)**Rational(1, 3))",
    "fps(exp(-x**2))", "fps(sin(x)**2)", "fps(erf(x))", "fps(cos(sqrt(x)))",
    "fps(log(1 - x)/x)", "fps(1/sqrt(1 - x**2))",
    # rational functions
    "fps(1/(1 - x))", "fps(1/(1 + x**2))", "fps((x + 1)/(x**2 - 5*x + 6))",
    "fps(x/(1 - x - x**2))", "fps(1/(x**2 + x + 1), full=True)", "fps(1/(1 - x)**3)",
    "fps(1/(x**4 + 1), full=True)", "fps(x**3/(1 - x))",
    # products and compositions
    "fps(exp(x)*sin(x))", "fps(x**2*cos(x))", "fps(exp(x)/(1 - x))", "fps(sin(x)/x)",
    "fps(exp(sin(x)))", "fps(log(cos(x)))",
    "fps(exp(x)).product(fps(sin(x)), x, 4).truncate(4)",
    "fps(exp(x)).product(fps(sin(x)), x, 2).truncate(2)",
    "fps(exp(x)).compose(fps(sin(x))).truncate(8)",
    "fps(1/(1 - x)).compose(fps(x*exp(x))).truncate(6)",
    "fps(exp(x)).inverse().truncate(8)", "fps(cos(x)).inverse().truncate(8)",
    # operations
    "(fps(exp(x)) + fps(sin(x))).truncate(8)", "(3*fps(sin(x))).truncate(8)",
    "(fps(sin(x)) + 1).truncate(8)", "(-fps(cos(x))).truncate(8)",
    "fps(sin(x)).diff(x).truncate(8)", "fps(exp(x)).integrate().truncate(8)",
    "fps(sin(x)).polynomial(8)", "fps(exp(x)).infinite", "fps(cos(x)).as_leading_term(x)",
    "fps(sin(x))[0:6]", "fps(exp(x)).zero_coeff()", "integrate(fps(sin(x)), (x, 0, 1))",
    "integrate(fps(sin(x)), x).function", "str(fps(sin(x)))",
    # nonzero points, infinity, direction
    "fps(sin(x), x0=1)", "fps(exp(x), x0=2)", "fps(log(x), x0=1)", "fps(1/x, x0=1)",
    "fps(cos(x), x0=pi)", "fps(atan(x), x0=1)",
    "fps(exp(1/x), x0=oo)", "fps(1/(1 + x), x0=oo)", "fps(sin(1/x), x0=-oo)",
    "fps(atan(x), x0=oo)", "fps(sin(x), dir=-1)", "fps(log(1 + x), dir=-1)",
    "fps(sin(x), dir='-')", "fps(sin(x), dir=2)",
    # options
    "fps(sin(x), hyper=False)", "fps(atan(x), hyper=False)", "fps(log(1 + x), rational=False)",
    "fps(sin(x), hyper=False, rational=False)", "fps(sin(x), order=1)",
    "fps(log(1 + x**2), full=True)", "fps(atan(x), full=True)",
    # symbolic, held and failing inputs
    "fps(x**n*sin(x))", "fps(x**n*atan(x), x, oo)", "fps(x**n*exp(x), x)",
    "fps(x**(n - 2)*sin(x) + x**n*exp(x), x)", "fps(exp(a*x), x)", "fps(1/(1 - a*x), x)", "fps(y*x)",
    "fps(x + y, x)", "fps(3)", "fps(x)", "fps(1/x)", "fps(sqrt(x))", "fps(log(x))",
    "fps(abs(x))", "fps(x + y)",
    # module functions
    "rational_algorithm(1/(1 - x), x, k)", "rational_algorithm(log(1 + x), x, k)",
    "rational_independent([sin(x), cos(x), x*sin(x)], x)",
    "next(simpleDE(asin(x), x, f))", "exp_re(-f(x) + Derivative(f(x), x), r, k)",
    "hyper_re(Derivative(f(x), x, 2) + f(x), r, k)",
    "rsolve_hypergeometric(exp(x), x, S(1), k + 1, k, 1)",
    "rsolve_hypergeometric(sin(x), x, S(-1), (k + 1)*(k + 2), k, 2)",
    "rsolve_hypergeometric(log(1 + x), x, -k, k + 1, k, 1)",
    "solve_de(exp(x), x, -f(x) + Derivative(f(x), x), 1, f, k)",
    "hyper_algorithm(sin(x), x, k)",
]


def _run(lib, case):
    names = {}
    exec("from %s import *\nfrom %s.series.formal import *\n"
         "x, y, a = symbols('x y a')\nn = symbols('n', integer=True)\n"
         "k = Dummy('k')\nf, r = Function('f'), Function('r')" % (lib, lib), names)
    names['k'] = names['symbols']('k')
    try:
        value = eval(case, names)
    except Exception as error:
        return 'EXC ' + type(error).__name__
    if isinstance(value, names['FormalPowerSeries']):
        return 'FPS ' + str(value.truncate(8))
    return str(value)


@pytest.mark.parametrize('case', CASES)
def test_same_as_reference(case):
    assert _run('symox', case) == _run('sympy', case)


@pytest.mark.parametrize('case', [
    "fps(exp(x)).product(fps(sin(x))).truncate(8)",
    "fps(exp(x)).product(fps(cos(x)), x, 10).truncate(8)",
    "fps(sin(x)).product(fps(1/(1 - x))).truncate(8)",
])
def test_product_value_as_reference(case):
    # Deliberate deviation: the reference convolves the coefficients by a
    # symbolic FFT whose radicals it leaves unsimplified; the native product
    # convolves them exactly. The coefficients are the same numbers.
    names = dict(vars(sp), x=sp.Symbol('x'))
    expected = eval(case, names).removeO()
    actual = eval(case, dict(vars(sx), x=sx.Symbol('x'))).removeO().to_sympy()
    difference = sp.expand(expected - actual)
    assert all(abs(sp.N(difference.coeff(names['x'], i), 50)) < 1e-40 for i in range(10))
    assert all(c.is_Rational for c in sp.Poly(actual, names['x']).coeffs())


def test_series_node_round_trip():
    x = sx.Symbol('x')
    series = sx.fps(sx.sin(x))
    assert series.func(*series.args) == series
    # Each fps has its own Dummy index: compare the printed forms.
    assert str(series.to_sympy()) == str(sp.fps(sp.sin(sp.Symbol('x'))))
    assert series.subs(x, sx.Symbol('y')) == series


@pytest.mark.parametrize('source,point', [
    ("x**2 + 1", 1), ("(1 + x)**2", 0), ("x*(x + 1)", 0), ("(x - 2)**3", 0),
    ("(1 + x)**2*(x - 3)", 2), ("x**3 + 2*x", 0),
])
def test_polynomial_keeps_taylor_coefficients(source, point):
    # Deliberate deviation: the reference takes the coefficients (and the
    # constant term) of the unexpanded polynomial, so fps((1 + x)**2) gives
    # O(x**4). The native series has the Taylor coefficients of the input.
    X = sp.Symbol('x')
    actual = sx.fps(eval(source, dict(vars(sx), x=sx.Symbol('x'))), sx.Symbol('x'), point)
    truncated = actual.truncate(6).to_sympy()
    assert sp.expand(truncated.removeO() - sp.sympify(source, locals={'x': X})) == 0
    assert truncated.getO() == sp.Order((X - point)**6, (X, point))


def test_symbolic_order_keeps_the_terms_it_claims():
    # Deliberate deviation: the reference truncate drops x**(n + 3)/2 but
    # still claims O(x**(n + 4)); the native one keeps every term below the
    # order term at a finite point.
    x, n = sx.Symbol('x'), sx.Symbol('n', integer=True)
    assert sx.fps(x**(n + 1)*sx.exp(x), x).truncate(3) == \
        x**(n + 3)/2 + x**(n + 2) + x**(n + 1) + sx.O(x**(n + 4), (x, 0))
