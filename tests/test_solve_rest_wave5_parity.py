"""Wave 5 solve/solveset/inequality parity with SymPy 1.14.

Expected values are the str() of the SymPy 1.14 results (errors as
"E:<type>: <message prefix>"). Two values deliberately differ (documented in
docs/algorithmic-parity.md): SymPy returns only (pi/4, 3*pi/4) for
cos(2*x) < 0 on [0, 2*pi], which misses (5*pi/4, 7*pi/4), and its
minsolve_linear_system misses solutions with a single nonzero value.
"""
import pytest


def _namespace():
    ns = {}
    exec('from symox import *', ns)
    exec('from symox.solvers.inequalities import *', ns)
    exec('from symox.solvers.solveset import *', ns)
    exec("x, y, z, a, b, c, t = symbols('x y z a b c t')\n"
         "xr, yr = symbols('xr yr', real=True)\n", ns)
    return ns


_NS = {}


CASES = [
    ('periodicity(sin(x)**2, x)',
     '2*pi'),
    ('periodicity(sin(x)*cos(x), x)',
     'pi'),
    ('periodicity(tan(x)**2 + sin(x), x)',
     '2*pi'),
    ('periodicity(sin(2*x) + cos(3*x), x)',
     '2*pi'),
    ('periodicity(Abs(sin(x)), x)',
     'pi'),
    ('periodicity(exp(sin(x)), x)',
     '2*pi'),
    ('periodicity(sin(x) + x, x)',
     'None'),
    ('periodicity(Mod(x, 3), x)',
     '3'),
    ('periodicity(sin(x)**2 + cos(x)**2, x)',
     '0'),
    ('periodicity(sin(x/3) + cos(x/2), x)',
     '12*pi'),
    ('solveset(sin(x)**2 > S(1)/4, x, Interval(0, 2*pi))',
     'Union(Interval.open(pi/6, 5*pi/6), Interval.open(7*pi/6, 11*pi/6))'),
    ('solveset(cos(2*x) < 0, x, Interval(0, 2*pi))',
     'Union(Interval.open(pi/4, 3*pi/4), Interval.open(5*pi/4, 7*pi/4))'),
    ('solveset(tan(x) > 1, x, Interval(0, pi))',
     'Interval.open(pi/4, pi/2)'),
    ('solve_univariate_inequality(sin(x)**2 < S(1)/2, x, relational=False)',
     'Union(Interval.Ropen(0, pi/4), Interval.open(3*pi/4, pi))'),
    ('solve([x + y - 2, x - y + z - 1], [x, y, z], particular=True)',
     '{x: 2, y: 0, z: -1}'),
    ('solve([x + y - 2, x - y + z - 1], [x, y, z], particular=True, quick=True)',
     '{x: 2, y: 0, z: -1}'),
    # SymPy's minsolve never tries a single nonzero value and returns
    # [{x: 3/2, y: 0, z: 3/2}]; y = 3 alone solves the system
    ('solve([x + 2*y + 3*z - 6, 2*x + y - 3], particular=True, dict=True)',
     '[{x: 0, y: 3, z: 0}]'),
    ('solve([x + y - 2], [x, y], quick=True)',
     'E:ValueError: when using `quick`, `particular` should be True'),
    ('solve(x + y - 2, [x, y], particular=True)',
     "E:ValueError: \nThe 'particular/quick' flag is usually used with systems of equations.\nEither pass your equation in a list or consider using a solver like\n`diophanti"),
    ('solve(x**2 - y, exclude=[y])',
     '[-sqrt(y), sqrt(y)]'),
    ('solve([x + y - 3, x - y - a], exclude=[a])',
     '{x: a/2 + 3/2, y: 3/2 - a/2}'),
    ('solve([x + y - 3, x - y - 1], [x, y], manual=True)',
     '[(2, 1)]'),
    ('solve([x**2 - y, x + y - 2], [x, y], manual=True)',
     '[(-2, 4), (1, 1)]'),
    ('solve(x**2 - 4, x, numerical=True)',
     '[-2, 2]'),
    ('solve(x**2 - 4, x, minimal=True)',
     '[-2, 2]'),
    ('solve(x**2 - 4, x, positive=True)',
     '[-2, 2]'),
    ('solve(Abs(xr) - 2, xr)',
     '[-2, 2]'),
    ('solve(Abs(xr - 1) - xr, xr)',
     '[1/2]'),
    ('solve(Abs(x) - 2, x)',
     'E:NotImplementedError: solving Abs(x) when the argument is not real or imaginary.'),
    ('solveset(x + Abs(x - 1) - 1, x, S.Reals)',
     'Interval(-oo, 1)'),
    ('solveset(x + Abs(x - 1) - 1, x, Interval(0, oo))',
     'Interval(0, 1)'),
    ('solveset(Abs(x) + Abs(x - 1) - 1, x, S.Reals)',
     'Interval(0, 1)'),
    ('solveset(Abs(x**2 - 4) - x, x, S.Reals)',
     '{-1/2 + sqrt(17)/2, 1/2 + sqrt(17)/2}'),
    ('solveset(Abs(x - 1) - Abs(x + 1), x, S.Reals)',
     '{0}'),
    ('solveset(Abs(x) - y, x, S.Reals)',
     'ConditionSet(x, Contains(y, Interval(0, oo)), {-y, y})'),
    ('solveset(Abs(sin(x)) - S.Half, x, S.Reals)',
     'Union(ImageSet(Lambda(_n, 2*_n*pi + pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 5*pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 7*pi/6), Integers), ImageSet(Lambda(_n, 2*_n*pi + 11*pi/6), Integers))'),
    ('solveset(Abs(cos(x)) - 1, x, Interval(0, 2*pi))',
     '{0, pi, 2*pi}'),
    ('solveset(Abs(x) + 1, x, S.Reals)',
     'EmptySet'),
    ('solveset(Abs(x) - 2, x)',
     'E:ValueError: \nAbsolute values cannot be inverted in the complex domain.'),
    ('solveset(Max(x, 1) - 3, x, S.Reals)',
     '{3}'),
    ('solveset(Min(x, 1) - Max(x, 0), x, S.Reals)',
     'Interval(0, 1)'),
    ('solveset(Piecewise((x - 1, x > 0), (x + 1, True)), x, S.Reals)',
     '{-1, 1}'),
    ('solveset(Abs(x - 1) < Abs(x), x, S.Reals)',
     'Interval.open(1/2, oo)'),
    ('solve_univariate_inequality(Max(x, 1) > 2, x)',
     '(2 < x) & (x < oo)'),
    ('solve_univariate_inequality(Min(x**2, 4) <= 1, x, relational=False)',
     'Interval(-1, 1)'),
    ('solve_univariate_inequality(Max(x, 1) < 0, x)',
     'False'),
    ('solve_univariate_inequality(Max(x, 1) >= 1, x)',
     'E:NotImplementedError: \nThe inequality, Max(1, x) >= 1, cannot be solved using\nsolve_univariate_inequality.'),
    ('solve_univariate_inequality(floor(x) > 1, x)',
     'E:NotImplementedError: \nThe inequality, floor(x) > 1, cannot be solved using\nsolve_univariate_inequality.'),
    ('reduce_inequalities(sqrt(x) > 1, x)',
     '1 < x'),
    ('reduce_inequalities(exp(x) < 1, x)',
     'x < 0'),
    ('reduce_inequalities(x**2 > 1, x)',
     '((-oo < x) & (x < -1)) | ((1 < x) & (x < oo))'),
    ('reduce_inequalities(log(x) > 0, x)',
     '1 < x'),
]


def _run(src):
    if not _NS:
        _NS.update(_namespace())
    try:
        out = str(eval(src, _NS))
    except Exception as e:  # noqa: BLE001
        out = 'E:' + type(e).__name__ + ': ' + str(e)[:150]
    return out.replace('symox', 'sympy')


@pytest.mark.parametrize('src,expected', CASES)
def test_solve_rest_wave5(src, expected):
    assert _run(src) == expected


def test_complex_abs_raises():
    from symox import Abs, S, solveset, symbols
    x = symbols('x')
    with pytest.raises(ValueError, match='Absolute values cannot be inverted'):
        solveset(Abs(x) - 2, x, S.Complexes)
