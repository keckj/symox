"""Forms of SymPy 1.14's simplify sequence and of its discrete probability
sums (expected values computed with SymPy 1.14)."""
from symox import Rational, cos, exp, lambdify, nonlinsolve, simplify, sin, symbols
from symox.stats import Geometric, P, Poisson


def test_simplify_returns_reference_sequence_forms():
    x = symbols('x')
    assert str(simplify((x + 5)/(16*x**2 + 28*x + 6))) == '(x + 5)/(2*(8*x**2 + 14*x + 3))'
    assert str(simplify((sin(x) + cos(x))**2)) == '2*sin(x + pi/4)**2'


def test_bounded_discrete_probabilities_sum_the_pmf_directly():
    p, t = symbols('p t', positive=True)
    assert str(P(Geometric('X', p) <= 4)) == 'p*(-p + (1 - p)**3 + (1 - p)**2 + 2)'
    assert P(Poisson('X', 3*t) < 4) == (Rational(9, 2)*t**3 + Rational(9, 2)*t**2 + 3*t + 1)*exp(-3*t)
    assert P(Poisson('X', 2) < 3) == 5*exp(-2)


def test_system_solvers_keep_native_normal_forms():
    # The quartic numerator system is still solved completely: its
    # restriction factors keep the normal form whose exact zero test decides.
    x, y = symbols('x y')
    solutions = nonlinsolve([x**3 - 1/y, x + y - 1], [x, y])
    assert len(solutions) == 4
    import mpmath
    with mpmath.workdps(80):
        residuals = lambdify([x, y], [x**3 - 1/y, x + y - 1], 'mpmath')
        for a, b in solutions:
            values = [lambdify([], c, 'mpmath')() for c in (a, b)]
            assert all(abs(r) < mpmath.mpf(10)**-60 for r in residuals(*values))
