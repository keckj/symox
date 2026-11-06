"""Wave 6 ode-solve leftovers; expected values come from SymPy 1.14."""
from symox import (I, Integer, Mod, Rational, S, acos, cos, exp, log, sin,
                   pi, roots, solve, solveset, sqrt, symbols, tanh)

x, y, C1, a, n = symbols('x y C1 a n')


# -- solveset: Mod equations over the integers (_solve_modular) -----------

def test_solveset_modular_integers():
    assert str(solveset(Mod(x, 3) - 1, x, S.Integers)) == \
        'ImageSet(Lambda(_n, 3*_n + 1), Integers)'
    assert str(solveset(3 - Mod(5*x - 8, 7), x, S.Integers)) == \
        'ImageSet(Lambda(_n, 7*_n + 5), Integers)'
    assert str(solveset(Mod(2*x, 7) - 3, x, S.Integers)) == \
        'ImageSet(Lambda(_n, 7*_n + 5), Integers)'
    assert str(solveset(Mod(x**2, 160) - 9, x, S.Integers)) == (
        'Union(ImageSet(Lambda(_n, 160*_n + 3), Integers), '
        'ImageSet(Lambda(_n, 160*_n + 13), Integers), '
        'ImageSet(Lambda(_n, 160*_n + 67), Integers), '
        'ImageSet(Lambda(_n, 160*_n + 77), Integers), '
        'ImageSet(Lambda(_n, 160*_n + 83), Integers), '
        'ImageSet(Lambda(_n, 160*_n + 93), Integers), '
        'ImageSet(Lambda(_n, 160*_n + 147), Integers), '
        'ImageSet(Lambda(_n, 160*_n + 157), Integers))')
    assert str(solveset(Mod(7**x, 41) - 15, x, S.Integers)) == \
        'ImageSet(Lambda(_n, 40*_n + 3), Naturals0)'


def test_solveset_modular_empty_and_unsupported():
    # proven absence of solutions
    assert solveset(7 - Mod(x, 5), x, S.Integers) is S.EmptySet
    assert solveset(Mod(x, 3) - I, x, S.Integers) is S.EmptySet
    # other domains
    assert str(solveset(Mod(x, 3) - 1, x, S.Naturals0)) == \
        'Intersection(ImageSet(Lambda(_n, 3*_n + 1), Integers), Naturals0)'
    assert str(solveset(Mod(x, 3) - 1, x, S.Reals)) == \
        'ConditionSet(x, Eq(Mod(x, 3) - 1, 0), Reals)'


def test_modular_pow():
    assert pow(Integer(3), 4, 5) == 1
    assert pow(Integer(3), -1, 5) == 2
    assert pow(x, 2, 3) == Mod(x**2, 3)


# -- solveset: radical equations (_solve_radical) -------------------------

def test_solveset_radicals():
    assert solveset(sqrt(x + 1) - sqrt(x) - 1, x, S.Reals) == {0}
    assert solveset(sqrt(x + 1) - sqrt(x) - 1, x) == {0}
    assert solveset(sqrt(x) + sqrt(x - 1) - 1, x, S.Reals) == {1}
    assert solveset(sqrt(17*x - sqrt(x**2 - 5)) - 7, x, S.Reals) == {3}
    assert solveset(sqrt(x) - 2, x, S.Reals) == {4}
    assert solveset(sqrt(x) + 2, x, S.Reals) is S.EmptySet


def test_solveset_real_cube_root_intentional():
    # SymPy returns {-8} through the real cube root, but with the principal
    # power (-8)**(1/3) = 1 + sqrt(3)*I, so -8 is not a solution (SymPy's
    # own Complexes result is EmptySet).
    assert solveset(x**Rational(1, 3) + 2, x, S.Reals) is S.EmptySet
    assert ((-8)**Rational(1, 3) + 2).expand(complex=True) != 0


def test_set_dummy_eq():
    r = solveset(Mod(x, 3) - 1, x, S.Integers)
    assert r.dummy_eq(solveset(Mod(x, 3) - 1, x, S.Integers))
    assert not r.dummy_eq(solveset(Mod(x, 5) - 1, x, S.Integers))


# -- solve: _tsolve steps -------------------------------------------------

def test_solve_symbolic_exponent():
    assert solve(x**y - 2, x) == [2**(1/y)]
    assert solve(x**n - 1, x) == [1]
    assert solve(-2*y/sqrt(y**3) - (C1 - x), y) == \
        [4/(C1**2 - 2*C1*x + x**2)]


def _same(u, v, subs):
    return abs(complex((u - v).subs(subs).evalf(30))) < 1e-20


def test_solve_logcombine():
    # rational_riccati_01; SymPy: [-sqrt(2)/tanh(sqrt(2)*(C1 - x))]
    eq = x + sqrt(2)*log(y - sqrt(2))/4 - sqrt(2)*log(y + sqrt(2))/4 - C1
    sol = solve(eq, y)
    assert len(sol) == 1
    ref = -sqrt(2)/tanh(sqrt(2)*(C1 - x))
    for pt in ({x: Rational(1, 3), C1: Rational(-2, 7)},
               {x: Rational(5, 2), C1: Rational(1, 9)}):
        assert _same(sol[0], ref, pt)
    # almost_lin_01: three cube roots, SymPy's form: the logcombine step
    # combines the logs of y only, then `_invert` gives
    # exp(3*C1*x + 3)**(1/x) (equal to exp(3*C1 + 3/x) for real C1 and x)
    eq = log(y + 1)/3 + log(y**2 - y + 1)/3 - 1/x - C1
    sol = solve(eq, y)
    w = (exp(3*C1*x + 3)**(1/x) - 1)**Rational(1, 3)
    assert set(sol) == {w, (-1 - sqrt(3)*I)*w/2, (-1 + sqrt(3)*I)*w/2}
    for s in sol:
        assert _same(((y + 1)*(y**2 - y + 1)).subs(y, s), exp(3*C1 + 3/x),
                     {x: Rational(3, 2), C1: Rational(1, 5)})


def test_exptrigsimp_fractional_power_intentional():
    # SymPy's exptrigsimp turns (exp(a) - 1)**(1/3) into
    # (-1)**(1/3)*(1 - exp(a))**(1/3), which is wrong for a > 0.
    from symox import exptrigsimp, simplify
    e = (-1 - sqrt(3)*I)*(exp(a) - 1)**Rational(1, 3)/2
    for r in (exptrigsimp(e), simplify(e)):
        assert _same(r, e, {a: 2})


# -- roots(..., trig=True) ------------------------------------------------

def test_roots_cubic_trig():
    eq = -x**3 + 2*x**2 + 3*x - 2
    t = acos(8*sqrt(13)/169)/3
    assert roots(eq, trig=True, multiple=True) == [
        Rational(2, 3) + 2*sqrt(13)*cos(t)/3,
        -2*sqrt(13)*sin(-t + pi/6)/3 + Rational(2, 3),
        -2*sqrt(13)*cos(-t + pi/3)/3 + Rational(2, 3)]
    assert roots(x**3 - 3*x + 1, trig=True) == {
        2*cos(2*pi/9): 1, 2*cos(4*pi/9): 1, -2*cos(pi/9): 1}
    # factorable or one real root: the radical forms
    assert roots(x**3 - x, trig=True) == {-1: 1, 1: 1, 0: 1}
    assert roots(x**3 + x + 1, trig=True) == roots(x**3 + x + 1)
