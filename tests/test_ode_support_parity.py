"""Core behaviors the ported SymPy ODE solvers rely on, checked against SymPy."""
import pytest
import sympy as sp
import symox as sx


def run(module, source):
    x, y, a, b = module.symbols('x y a b')
    xp = module.Symbol('xp', positive=True)
    names = vars(module) | {
        'x': x, 'y': y, 'a': a, 'b': b, 'xp': xp,
        'f': module.Function('f'), 'g': module.Function('g'),
        'u': module.Function('u'),
    }
    return eval(source, names)


@pytest.mark.parametrize('source', [
    # Integrals over an applied function and function-aware has()
    'Integral(x*f(x), f(x)).doit()',
    '(sin(f(x)) + x).has(f(x))',
    # Substituting a lower derivative inside a higher one
    'Derivative(f(x), x, 2).subs(Derivative(f(x), x), y)',
    'Derivative(f(x, y), x, 2, y).subs(Derivative(f(x, y), x), g(x))',
    # Powers with symbolic exponents
    'diff(x**a, x)',
    'diff(x**(a*x), x)',
    # cse keeps limit tuples intact
    'cse([Integral(exp(x)*x, x) + exp(x), Integral(exp(x)*x, x)])',
    'cse(Integral(f(x)*x, (x,)) + x)',
    # powsimp when -b and b coincide
    'powsimp(zoo**x*(-zoo)**y)',
    # log expansion of numeric and assumption-dependent arguments
    'expand(log(2*x))',
    'expand_log(log(8*x/3))',
    'expand_log(log(x**2*y), force=True)',
    'expand_log(log(-2*xp))',
    'expand_log(log(xp**a))',
    'expand_log(log(exp(xp)*xp**2))',
    # logcombine
    'logcombine(log(x) + 2*log(y), force=True)',
    'logcombine(a*log(x) - log(y))',
    'logcombine(log(xp) - log(2))',
    'logcombine(log(x) + log(y) - log(x*y) + exp(log(x) + log(2)), force=True)',
    'logcombine(Integral(log(x) + log(y), x), force=True)',
    # roots of an equation and solve fallbacks
    'roots(Eq(x**2, 4), x)',
    'solve(exp(f(x)) - y, f(x))',
    'solve([], [])',
    # Polynomials over RR and float coefficients
    'sorted(map(str, Poly(x*y + a, x).free_symbols))',
    'factor(2.5*x**2 - 10.0)',
    'factor_list(0.1*x**2 + 0.3*x + 0.2)',
    'Poly(x**2 + 2.0*x + 1.0).sqf_list()',
    'factor(1.5*x*y + 3.0*y)',
    'roots(x + 51.05*y, x)',
    'roots(2.0*x**2 - y, x, multiple=True)',
    # Powers of float multiples
    'sqrt(0.5*y)', '(0.5*y)**Rational(1, 3)', 'sqrt(-2.0*y)', '(2.5*y**2)**Rational(1, 2)',
    # exp facts, substitution and powsimp used by the solvers
    '(exp(-x).is_zero, exp(x).is_nonzero, exp(a).is_finite)',
    'exp(a + x).subs(exp(a), y)',
    'exp(-2*x + a).subs(exp(-x), b)',
    '(4**x - 3*y + 2).subs(2**x, y)',
    'powsimp(y*exp(-a*x))',
    'solve([a - x, b - y*exp(-x)], [x, y], check=False)',
    # collect rebuilds exp(c*t) as exp(t)**c
    'collect(exp(x**3/3)*y + 1, y)',
    'collect(a*exp(3*x/2) + a, a)',
    'collect(exp(x**3/3)*f(x).diff(x) - 1, f(x))',
    # Symbolic quadratics, radicals and Lambert forms met by implicit solutions
    'solve(x**2 + a*x + 1, x)',
    'solve(x**2 - 2*a*x + a**2 - 4*b**2, x)',
    'solve(sqrt(y + 1) - 2*sqrt(y - 1), y)',
    'solve(sqrt(y + 1)/sqrt(y - 1) - a, y)',
    'solve(a + x*exp(-y) + y*exp(-y), y)',
    'solve(LambertW(x) - a, x)',
    'solve(LambertW(a*exp(-x)) - y, a)',
    'solve(log(y - a)/6 + log(y + a)/6 - Rational(1, 2), y)',
    # Non-symbol targets are recast by substitution; held derivatives and
    # integrals of the target are masked
    'solve(Derivative(f(x), x) - f(x) - x, f(x))',
    'solve(exp(x) + exp(-x) - 2, exp(x))',
    'solve(a + x - Integral(g(y), y), y)',
    'solve(a + x - Derivative(g(y), y), y)',
    # Sign of a cancelled denominator and combined exponentials
    'simplify((a + x)/(a - x))',
    'simplify((a**2 + x**2)/(a**2 - x**2))',
    'simplify(x/(1 - x))',
    'simplify((a + b)/(b - a))',
    'simplify(b/(b*x - a*y))',
    'simplify(b/(a*y - b*x))',
    'powsimp(3*exp(6*a)*exp(3*x**2))',
    'simplify(3*exp(6*a)*exp(3*x**2) + 1)',
    # Equal exponents merged by simplify's powsimp step
    'simplify(sqrt(2)*sqrt(a + b))',
    'simplify(-sqrt(3)*x*sqrt(a + b)/3)',
    'solve(y**2 - 2*a - 2*b, y)',
    'solve(f(x)**2/2 - a - x, f(x))',
    # terms_gcd options used by constantsimp
    'terms_gcd(x - x**4/12, clear=False)',
    'terms_gcd(a*(x - x**4/12) + b*(1 - x**3/6) + O(x**6), clear=False, deep=True, expand=False)',
    'terms_gcd(x/2 + y/4, clear=False)',
    'terms_gcd(Eq(2*x + 2, 4*y), deep=True)',
    'terms_gcd(x*sin(2*y + 2*x*y)**2 + x, deep=True, clear=False)',
    'Poly((x + 1)**2 + x, expand=False).gens',
    # simplify keeps the cancelled fraction on a count_ops tie, so the
    # factorable method reaches its fixed point as in the reference
    'simplify(-2*x/(x + cos(y)) - y/(x + cos(y)))',
    'solve(x*f(x).diff(x) + 2*x + f(x) + cos(f(x))*f(x).diff(x), f(x).diff(x))',
    'classify_ode(f(x).diff(x) - (-2*x - f(x))/(x + cos(f(x))), f(x))',
    # simplify applies factor_terms(sign=False) and distributes a hollow
    # number*(sum) product, as the reference does
    'simplify(x - x*f(x))',
    'simplify((x - x*y)/(1 + y))',
    'simplify(2*x + 2*x*y)',
    'simplify(sin(2*x + 2) + 2*x*y)',
    'solve(f(x).diff(x) + x*f(x) - x, f(x).diff(x))',
    'classify_ode(f(x).diff(x) + x*f(x) - x, f(x))',
    # a product nested in an unevaluated product is parenthesized
    'Mul._from_args((Rational(1, 8), I*exp(I*x)))',
    'terms_gcd(Eq(2*x + 2, 4*y))',
    'terms_gcd(0.5*x + 1.5*x**2)',
    'terms_gcd(a*exp(x/3) + x*exp(-x), clear=False, deep=True, expand=False)',
    'Poly(a*exp(x/3) + x, expand=False).gens',
    'Poly((a + b*x)*exp(-x), expand=False).clear_denoms(convert=True)',
    'terms_gcd(x**2 + I*x**2*exp(I*x)/8, clear=False, deep=True, expand=False)',
    'terms_gcd(x/2 + (a + x)*exp(x)/3, clear=False, deep=True, expand=False)',
    # Solvers and classification end to end
    'pdsolve((x + y)*u(x, y).diff(x) + u(x, y).diff(y), u(x, y))',
    'dsolve(f(x).diff(x) - f(x), f(x), simplify=False)',
    'dsolve(f(x).diff(x, 2) - x*f(x), f(x))',
    'dsolve(x**2*f(x).diff(x, 2) + x*f(x).diff(x) + (x**2 - 4)*f(x), f(x))',
    'dsolve(f(x).diff(x) - 2*f(x), f(x), ics={f(0): 3})',
    'dsolve(f(x).diff(x, 2) + f(x), f(x), ics={f(0): 1, f(x).diff(x).subs(x, 0): 2})',
    'classify_ode(f(x).diff(x, 2) - x*f(x), f(x))',
    'dsolve(f(x).diff(x, 2) + x*f(x), f(x), hint="2nd_power_series_ordinary", n=6)',
    'dsolve(f(x).diff(x, 2) + x*f(x), f(x), hint="2nd_power_series_ordinary", n=10)',
    'dsolve(f(x)*f(x).diff(x) - 1, f(x))',
    'classify_ode(x**2*f(x).diff(x, 2) + x*f(x).diff(x) + (x**2 - 4)*f(x), f(x))',
])
def test_matches_reference(source):
    assert str(run(sx, source)) == str(run(sp, source))


def test_frobenius_series_is_a_true_solution():
    # SymPy 1.14 calls series(ser, x=x0, n=n+1), which leaves ser unchanged and
    # yields wrong coefficients for this equation.  Symox expands in x, so the
    # residual only starts beyond the requested order.
    x = sx.Symbol('x')
    f = sx.Function('f')
    eq = x*(x - 1)*f(x).diff(x, 2) + (sx.Rational(3, 2) - 2*x)*f(x).diff(x) + 2*f(x)
    sol = sx.dsolve(eq, f(x), hint='2nd_power_series_regular')
    C1, C2 = sx.symbols('C1 C2')
    body = sol.rhs.removeO()
    for k, v in ((C1, 1), (C2, 0)), ((C1, 0), (C2, 1)):
        y = body.subs({k[0]: k[1], v[0]: v[1]})
        res = sx.expand(x*(x - 1)*y.diff(x, 2) + (sx.Rational(3, 2) - 2*x)*y.diff(x) + 2*y)
        lowest = min(term.as_coeff_exponent(x)[1] for term in sx.Add.make_args(res)) if res else sx.oo
        assert lowest >= 4
    assert sol.rhs.removeO() == C2*(1 - 4*x/3) + C1*x**sx.Rational(5, 2)*(
        5*x**2/56 + 3*x/14 + 1)


def test_initial_values_may_name_constants():
    # SymPy 1.14 confuses a user symbol C1 in the initial data with the
    # generated constant ("Couldn't solve for initial conditions").  Symox
    # keeps the user symbol distinct.
    x, C1 = sx.symbols('x C1')
    f = sx.Function('f')
    eq = f(x).diff(x) - f(x)
    sol = sx.dsolve(eq, f(x), ics={f(0): C1})
    assert sol == sx.Eq(f(x), C1*sx.exp(x))
    assert sx.checkodesol(eq, sol, f(x)) == (True, 0)


@pytest.mark.parametrize('source', [
    # SymPy divides by b**2 + c**2 = 0 here and returns False
    'u(x, y).diff(x) + I*u(x, y).diff(y)',
    'u(x, y).diff(x) + I*u(x, y).diff(y) + 2*u(x, y) - x',
    # SymPy substitutes x*y inside an integral over x, capturing the variable
    'x*u(x, y).diff(x) - y*u(x, y).diff(y) + y**2*u(x, y) - g(x, y)',
])
def test_pdsolve_returns_checked_solutions(source):
    x, y = sx.symbols('x y')
    u = sx.Function('u')
    eq = run(sx, source)
    sol = sx.pdsolve(eq, u(x, y))
    assert isinstance(sol, sx.Eq)
    assert sx.checkpdesol(eq, sol, u(x, y)) == (True, 0)


def _series_residual_order(eq, f, x, y):
    """Lowest power of x in the residual of y, with log(x) held as a symbol."""
    L = sx.Symbol('L')
    res = sx.expand(eq.subs(f(x), y).doit().subs(sx.log(x), L))
    if res == 0:
        return sx.oo
    return min(t.as_coeff_exponent(x)[1]
               for c in sx.Poly(res, L).all_coeffs() if c != 0
               for t in sx.Add.make_args(sx.expand(c)))


@pytest.mark.parametrize('source, lowest', [
    # Bessel of order 0: double indicial root 0
    ('x**2*f(x).diff(x, 2) + x*f(x).diff(x) + x**2*f(x)', 6),
    # Bessel of order 1: roots 1 and -1, inconsistent recurrence at resonance
    ('x**2*f(x).diff(x, 2) + x*f(x).diff(x) + (x**2 - 1)*f(x)', 6),
    ('x**2*f(x).diff(x, 2) + 3*x*f(x).diff(x) + (1 - x)*f(x)', 6),
    ('x*f(x).diff(x, 2) - f(x)', 5),
])
def test_frobenius_logarithmic_second_solution(source, lowest):
    # SymPy 1.14 returns only the first Frobenius series in these cases; the
    # second element of the basis carries log(x).
    x = sx.Symbol('x')
    f = sx.Function('f')
    C1, C2 = sx.symbols('C1 C2')
    eq = run(sx, source)
    sol = sx.dsolve(eq, f(x), hint='2nd_power_series_regular', n=6)
    reference = sp.dsolve(run(sp, source), sp.Function('f')(sp.Symbol('x')),
                          hint='2nd_power_series_regular', n=6)
    body = sol.rhs.removeO()
    assert sol.rhs.getO() == sx.Order(x**6)
    # constantsimp may take a common factor out of both terms ((C1*A +
    # C2*B)/x), so read the basis from the derivatives in the constants.
    parts = [body.diff(C) for C in (C1, C2)]
    assert sum(part.has(sx.log) for part in parts) == 1
    first = [part for part in parts if not part.has(sx.log)][0]
    assert str(first) == str(reference.rhs.removeO().diff(sp.Symbol('C1')))
    for part in parts:
        assert _series_residual_order(eq, f, x, part) >= lowest
