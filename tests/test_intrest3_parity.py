"""Integration leftovers of wave 11: guards on a symbolic frequency, held
values where a pole may lie on the interval, the degenerate imaginary poles
of x/(x**2 + a**2), and the held simple delta of deltaintegrate. The
expected values are those of SymPy 1.14 unless a comment says otherwise."""
import pytest
from symox import *
del test

x = Symbol('x')
a = Symbol('a')
ar = Symbol('a', real=True)
ap = Symbol('a', positive=True)
anz = Symbol('a', nonzero=True)
br = Symbol('b', real=True)
b = Symbol('b')
n = Symbol('n', integer=True)


@pytest.mark.parametrize('f, expected', [
    (sin(a*x)**2, Piecewise(((a*x/2 - sin(a*x)*cos(a*x)/2)/a, Ne(a, 0)), (0, True))),
    (sin(a*x)*cos(a*x), Piecewise((sin(a*x)**2/(2*a), Ne(a, 0)), (0, True))),
    (cos(a*x)**2, Piecewise(((a*x/2 + sin(a*x)*cos(a*x)/2)/a, Ne(a, 0)), (x, True))),
    (sin(a*x)**3, Piecewise(((cos(a*x)**3/3 - cos(a*x))/a, Ne(a, 0)), (0, True))),
    (sin(a*x)**2*cos(a*x), Piecewise((sin(a*x)**3/(3*a), Ne(a, 0)), (0, True))),
    (sin(a*x)**4, Piecewise(((3*a*x/8 - sin(a*x)**3*cos(a*x)/4 - 3*sin(a*x)*cos(a*x)/8)/a, Ne(a, 0)), (0, True))),
    (tan(a*x)**2, Piecewise(((-a*x + sin(a*x)/cos(a*x))/a, Ne(a, 0)), (0, True))),
    (sec(a*x)**2, Piecewise((sin(a*x)/(a*cos(a*x)), Ne(a, 0)), (x, True))),
])
def test_symbolic_frequency_keeps_the_guard(f, expected):
    assert integrate(f, x) == expected


@pytest.mark.parametrize('f', [sin(a*x)**2*cos(a*x)**2, sin(2*a*x + 1)**2])
def test_symbolic_frequency_guard_other_forms(f):
    # SymPy's own forms differ (other reduction formulas); the value is
    # checked: a Piecewise whose a = 0 branch is f(0)*x and whose generic
    # branch differentiates back to f.
    r = integrate(f, x)
    assert isinstance(r, Piecewise)
    (generic, cond), (at_zero, _) = r.args
    assert cond == Ne(a, 0)
    assert at_zero == f.subs(a, 0)*x
    assert simplify(generic.diff(x) - f) == 0


def test_definite_symbolic_frequency_at_zero():
    assert integrate(sin(n*x)**2, (x, 0, pi)) == Piecewise((pi/2, Ne(n, 0)), (0, True))
    assert integrate(sin(ap*x)**2, x) == (ap*x/2 - sin(ap*x)*cos(ap*x)/2)/ap


def test_degenerate_imaginary_poles():
    # the pole at x = 0 only reaches the interval at a = 0, where the
    # value is not finite
    assert integrate(x/(x**2 + ar**2), (x, 0, 1)) == -log(ar**2)/2 + log(ar**2 + 1)/2
    assert integrate(1/(x**2 + ar**2), (x, 0, 1)) == atan(1/ar)/ar
    assert integrate(x/(x**2 + ar**2), (x, 1, 2)) == -log(ar**2 + 1)/2 + log(ar**2 + 4)/2
    # a complex a may put a pole on the interval (a = I/2): held
    assert isinstance(integrate(x/(x**2 + a**2), (x, 0, 1)), Integral)


def test_improper_rational_with_undecided_pole_is_held():
    # SymPy returns 1/a, 1/a and 1/(2*a**2): wrong for a < 0, where a pole
    # lies on the interval and the integral diverges
    assert isinstance(integrate(1/x**2, (x, ar, oo)), Integral)
    assert isinstance(integrate(1/(x + ar)**2, (x, 0, oo)), Integral)
    assert isinstance(integrate(1/x**3, (x, anz, oo)), Integral)
    # decided cases are unchanged
    assert integrate(1/x**2, (x, ap, oo)) == 1/ap
    assert integrate(1/(x**2 + ar**2), (x, 0, oo)) == Piecewise(
        (pi/(2*ar), (Abs(arg(ar)) < pi/2) | ((Abs(arg(ar)) < pi/2) & Ne(Abs(arg(ar)), pi/2))),
        (Integral(1/(ar**2 + x**2), (x, 0, oo)), True))


def test_deltaintegrate_keeps_the_slope():
    from symox.integrals.deltafunctions import deltaintegrate
    # SymPy expands to Heaviside(x - b/a)/Abs(a); the expansion is held
    # while a may vanish, and the reference fallback Heaviside(a*x - b)
    # dropped the 1/a
    r = deltaintegrate(DiracDelta(ar*x - br), x)
    assert r == Piecewise((Heaviside(ar*x - br)/ar, Ne(ar, 0)), (x*DiracDelta(br), True))
    assert integrate(DiracDelta(anz*x - br), x) == Heaviside(anz*x - br)/anz
    assert deltaintegrate(DiracDelta(x - 1), x) == Heaviside(x - 1)
