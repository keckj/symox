"""Wave 9 solint-rest2: precision, log-singular leading terms, Frullani-type
integrals, SingularityFunction integrals, ratint and sec/csc/cot solveset.

Expected values are hard-coded from SymPy 1.14 unless a comment says the
SymPy result is mathematically wrong (symox then keeps the correct value).
"""
from symox import *
del test

from symox.integrals.rationaltools import ratint

x, y, a, b, c = symbols('x y a b c')


def _cancelled():
    r = CRootOf(7*x**5 - 7*x**3 + 1, 1)
    return root(r**15, 3) - root(r**15, 5) + Rational(1, 7)


def _is_scaled_zero(value):
    # SymPy prints -0.e-110: a zero Float with no significant digit. The sign
    # and the exact exponent are rounding noise and are not reproduced.
    # SymPy's scaled zero compares unequal to 0 (and is_zero is False).
    return str(value).lstrip('-').startswith('0.e-')


def test_cancelled_crootof_sum_is_a_scaled_zero():
    f = _cancelled()
    assert _is_scaled_zero(f.n(2))
    assert _is_scaled_zero(f.n(15))
    assert _is_scaled_zero((2*f).n(2))
    assert _is_scaled_zero(re(f).n(2))
    assert _is_scaled_zero((x + f).subs(x, 0).n(2))
    assert Abs(f).n(2) == 0
    assert im(f).n(2) == 0


def test_cancelled_crootof_sum_nested_in_functions():
    f = _cancelled()
    # SymPy prints 6.1e-238 (rounding noise); a scaled zero is the exact
    # information available.
    assert _is_scaled_zero((f**2).n(2))
    assert _is_scaled_zero((pi*f).n(2))
    assert _is_scaled_zero((f*(f + 1)).n(2))
    assert str(exp(f).n(2)) == '1.0'
    assert str((f + 1).n(2)) == '1.0'
    assert str(cos(1 + f).n(2)) == '0.54'


def test_log_singular_leading_terms():
    assert Ei(x).as_leading_term(x) == log(x) + EulerGamma
    assert Ci(x).as_leading_term(x) == log(x) + EulerGamma
    assert Chi(x).as_leading_term(x) == log(x) + EulerGamma
    # SymPy raises NotImplementedError for expint; E1(x) ~ -log(x) - EulerGamma.
    assert expint(1, x).as_leading_term(x) == -log(x) - EulerGamma


def test_log_singular_limits():
    assert limit(Ei(-x) - Ei(-2*x), x, 0, '+') == -log(2)
    assert limit(Ci(x) - Ci(2*x), x, 0, '+') == -log(2)
    assert limit(Chi(x) - Chi(2*x), x, 0, '+') == -log(2)


def test_frullani_integrals():
    assert integrate((exp(-x) - exp(-2*x))/x, (x, 0, oo)) == log(2)
    assert integrate((cos(x) - cos(2*x))/x, (x, 0, oo)) == log(2)


def test_singularity_function_over_the_real_line():
    assert integrate(SingularityFunction(x, 0, -1), (x, -oo, oo)) == 1
    assert integrate(SingularityFunction(x, 2, -1)*x**2, (x, -oo, oo)) == 4
    assert integrate(SingularityFunction(x, 0, -2), (x, -oo, oo)) == 0
    assert integrate(SingularityFunction(x, 0, -1), (x, 0, oo)) == 0
    assert integrate(SingularityFunction(x, 0, -1), (x, -oo, 5)) == 1
    assert integrate(SingularityFunction(x, a, 0), (x, -oo, 5)) == SingularityFunction(5, a, 1)


def test_singularity_function_divergent_integrals():
    # SymPy 1.14 returns 0 for both, which is wrong: <x - a>^0 is a unit
    # step and <x - 1>^1 a ramp, so the integrals over (-oo, oo) diverge.
    assert integrate(SingularityFunction(x, a, 0), (x, -oo, oo)) == oo
    assert integrate(SingularityFunction(x, 1, 1), (x, -oo, oo)) == oo


def test_ratint_cubic_with_one_real_root():
    # SymPy returns a real closed form with a log of the real root, a log of
    # the quadratic factor and an atan; symox's radical form differs.
    result = ratint(1/(x**3 + x + 1), x)
    assert not result.has(I)
    check = diff(result, x) - 1/(x**3 + x + 1)
    for point in (Rational(3, 7), -3, 5):
        assert abs(complex(N(check.subs(x, point), 30))) < 1e-20
    assert abs(N(result.subs(x, 2) - result.subs(x, 1), 30)
               - Rational(184135444388838876317217752915, 10**30)) < 1e-25


def test_reciprocal_trig_solveset_forms():
    def lattice(*bodies):
        sets = ['ImageSet(Lambda(_n, %s), Integers)' % body for body in bodies]
        return sets[0] if len(sets) == 1 else 'Union(%s)' % ', '.join(sets)
    cases = [
        (sec(2*x) - 3, S.Complexes, lattice('_n*pi - asec(3)/2', '_n*pi + asec(3)/2')),
        (sec(2*x) - 3, S.Reals, lattice('_n*pi - asec(3)/2', '_n*pi + asec(3)/2')),
        (sec(x) - 3, S.Complexes, lattice('2*_n*pi - asec(3)', '2*_n*pi + asec(3)')),
        (sec(x) - 3, S.Reals, lattice('2*_n*pi - asec(3)', '2*_n*pi + asec(3)')),
        (csc(x) - 3, S.Complexes, lattice('2*_n*pi + acsc(3)', '2*_n*pi - acsc(3) + pi')),
        (csc(x) - 3, S.Reals, lattice('2*_n*pi + acsc(3)', '2*_n*pi - acsc(3) + pi')),
        (sec(x) - Rational(1, 2), S.Complexes, lattice('2*_n*pi - asec(1/2)', '2*_n*pi + asec(1/2)')),
        (sec(x) - Rational(1, 2), S.Reals, 'EmptySet'),
        (sec(x) - I, S.Complexes, lattice('2*_n*pi - asec(I)', '2*_n*pi + asec(I)')),
        (cot(x) - 3, S.Complexes, lattice('_n*pi + acot(3)')),
        (csc(2*x + 1) + 2, S.Complexes, lattice('_n*pi - 1/2 + 7*pi/12', '_n*pi - 1/2 + 11*pi/12')),
        (sec(x) + 1, S.Complexes, lattice('2*_n*pi + pi')),
        (sec(x) - 1, S.Reals, lattice('2*_n*pi')),
        (sec(x) - sqrt(2), S.Complexes, lattice('2*_n*pi + pi/4', '2*_n*pi + 7*pi/4')),
        (sec(3*x) - 2, S.Reals, lattice('2*_n*pi/3 + pi/9', '2*_n*pi/3 + 5*pi/9')),
    ]
    for equation, domain, expected in cases:
        assert str(solveset(equation, x, domain)) == expected, (equation, domain)


def test_sec_solve_and_symbolic_systems():
    assert set(solve(sec(2*x) - 3, x)) == {
        -I*log(-sqrt(Rational(1, 3) - 2*sqrt(2)*I/3)),
        -I*log(Rational(1, 3) - 2*sqrt(2)*I/3)/2,
        -I*log(-sqrt(Rational(1, 3) + 2*sqrt(2)*I/3)),
        -I*log(Rational(1, 3) + 2*sqrt(2)*I/3)/2}
    assert solve([a*x + b*y - 1, c*x - y], [x, y]) == {x: 1/(a + b*c), y: c/(a + b*c)}
    assert nonlinsolve([a*x**2 - y, x - b], [x, y]) == FiniteSet((b, a*b**2))
