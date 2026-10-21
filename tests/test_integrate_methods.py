"""The risch, manual and heurisch options of integrate and Integral.doit.

The expected results were produced by the reference SymPy 1.14 on the same
inputs. ``risch=True`` runs only the Risch algorithm and holds a proven
non-elementary part as a NonElementaryIntegral, ``manual=True`` and
``heurisch=True`` run the reference integrator with that method alone, and
a False option removes the method from the default chain. The default
chain reports a NonElementaryIntegral only from a Risch proof.
"""
import pytest
from symox import *
from symox.integrals.risch import NonElementaryIntegral

# the star import brings the library test runner, which is not a test here
del test

x, y = symbols('x y')
a = Symbol('a', positive=True)

VALUES = [
    ('integrate(1/(x**2+a), x, heurisch=True)',
     'Mul', 'atan(x/sqrt(a))/sqrt(a)'),
    ('integrate(1/(x**2+y), x, heurisch=True)',
     'Add', '-sqrt(-1/y)*log(x - y*sqrt(-1/y))/2 + sqrt(-1/y)*log(x + y*sqrt(-1/y))/2'),
    ('integrate(sin(x)**3*cos(x)**2, x, heurisch=True)',
     'Add', 'cos(x)**5/5 - cos(x)**3/3'),
    ('integrate(x*exp(x), x, heurisch=True)',
     'Add', 'x*exp(x) - exp(x)'),
    ('integrate(exp(-x**2), x, heurisch=True)',
     'Mul', 'sqrt(pi)*erf(x)/2'),
    ('integrate(x*sin(x*y), x, heurisch=True)',
     'Piecewise', 'Piecewise((-x*cos(x*y)/y + sin(x*y)/y**2, Ne(y, 0)), (0, True))'),
    ('integrate(sin(x)*cos(x), x, manual=True)',
     'Mul', '-cos(x)**2/2'),
    ('integrate(x*exp(x), x, manual=True)',
     'Add', 'x*exp(x) - exp(x)'),
    ('integrate(exp(x)/x, x, manual=True)',
     'Ei', 'Ei(x)'),
    ('integrate(x**2*sqrt(1-x**2), x, manual=True)',
     'Add', 'sqrt(1 - x**2)*(x**3/4 - x/8) + asin(x)/8'),
    ('integrate(1/(x**5+x+1), x, manual=True)',
     'Add', '-log(x**2 + x + 1)/14 + 5*sqrt(3)*atan(2*sqrt(3)*(x + 1/2)/3)/21 - 4*Integral(x/(x**3 - x**2 + 1), x)/7 + Integral(x**2/(x**3 - x**2 + 1), x)/7 + 5*Integral(1/(x**3 - x**2 + 1), x)/7'),
    ('integrate(x*exp(x), (x, 0, 1), manual=True)',
     'One', '1'),
    ('integrate(sin(x), (x, 0, pi), heurisch=True)',
     'Integer', '2'),
    ('integrate(exp(-x), (x, 0, oo), manual=True)',
     'One', '1'),
    ('integrate(exp(-x**2), x, risch=True)',
     'NonElementaryIntegral', 'Integral(exp(-x**2), x)'),
    ('integrate(x**x, x, risch=True)',
     'NonElementaryIntegral', 'Integral(x**x, x)'),
    ('integrate(exp(x)/x, x, risch=True)',
     'NonElementaryIntegral', 'Integral(exp(x)/x, x)'),
    ('integrate((x**2-1)*sin(x), x, risch=True)',
     'Integral', 'Integral((x - 1)*(x + 1)*sin(x), x)'),
    ('integrate(5, x, risch=True)',
     'Mul', '5*x'),
    ('integrate(x + x**x, x)',
     'NonElementaryIntegral', 'Integral(x + x**x, x)'),
    ('integrate(x**x, (x, 0, 1))',
     'NonElementaryIntegral', 'Integral(x**x, (x, 0, 1))'),
    ('integrate(exp(-x**2)*log(x), x)',
     'Add', '-x*hyper((1/2, 1/2), (3/2, 3/2), -x**2) + sqrt(pi)*log(x)*erf(x)/2'),
    ('integrate(sin(x)/x, x, risch=False)',
     'Si', 'Si(x)'),
    ('integrate(exp(x**2), x, risch=False, meijerg=False)',
     'Mul', 'sqrt(pi)*erfi(x)/2'),
    ('integrate(x*sqrt(1+x), x, manual=False)',
     'Add', '2*x**2*sqrt(x + 1)/5 + 2*x*sqrt(x + 1)/15 - 4*sqrt(x + 1)/15'),
    ('integrate(1/(x*log(x)), x, heurisch=False)',
     'log', 'log(log(x))'),
    ('Integral(exp(x)*exp(exp(x)), x).doit(risch=True)',
     'exp', 'exp(exp(x))'),
    ('Integral(sin(x)*exp(x), x).doit(manual=True)',
     'Add', 'exp(x)*sin(x)/2 - exp(x)*cos(x)/2'),
    ('integrate(exp(x**2)*log(x), x)',
     'Add', '-x*hyper((1/2, 1/2), (3/2, 3/2), x**2) + sqrt(pi)*log(x)*erfi(x)/2'),
    ('integrate(atan(x)*x, x, manual=True)',
     'Add', 'x**2*atan(x)/2 - x/2 + atan(x)/2'),
    ('integrate(sqrt(x**2+1), x, heurisch=True)',
     'Add', 'x*sqrt(x**2 + 1)/2 + asinh(x)/2'),
    ('integrate((x**2-1)*sin(x)/log(x), x)',
     'Integral', 'Integral((x - 1)*(x + 1)*sin(x)/log(x), x)'),
    ('integrate(y*exp(x)/log(x), x)',
     'Mul', 'y*Integral(exp(x)/log(x), x)'),
    ('integrate(exp(x**2), x, meijerg=False)',
     'Mul', 'sqrt(pi)*erfi(x)/2'),
    ('integrate(1/(x**2+y), x)',
     'Add', '-sqrt(-1/y)*log(x - y*sqrt(-1/y))/2 + sqrt(-1/y)*log(x + y*sqrt(-1/y))/2'),
    ('integrate(sin(x)*exp(x**2), x)',
     'Integral', 'Integral(exp(x**2)*sin(x), x)'),
    ('integrate(x**x*log(x) + x**x, x)',
     'Pow', 'x**x'),
    ('integrate(1/log(x) + x, x)',
     'Add', 'x**2/2 + li(x)'),
    ('integrate(exp(-x**2)/(x+1), x)',
     'NonElementaryIntegral', 'Integral(exp(-x**2)/(x + 1), x)'),
    # an Integral integrand is denested into the limits, a
    # NonElementaryIntegral one is not
    ('integrate(Integral(y, x), x)', 'Mul', 'x**2*y/2'),
    ('integrate(Integral(x*y, x), y, risch=True)', 'Mul', 'x**2*y**2/4'),
    ('integrate(Integral(y, x), (x, 0, 1), risch=False)', 'Mul', 'y/2'),
    ('integrate(Integral(y, x), x, manual=True)', 'Mul', 'x**2*y/2'),
    ('integrate(Integral(y, x), x, heurisch=True)', 'Mul', 'x**2*y/2'),
    ('integrate(integrate(exp(-x**2)/(x+1), x), y)', 'Mul', 'y*Integral(exp(-x**2)/(x + 1), x)'),
    ('integrate(Integral(sin(x), (x, 0, y)), y)', 'Add', 'y - sin(y)'),
]

ERRORS = [
    ('integrate(x, (x, 0, 1), risch=True)', ValueError, 'risch=True is only allowed for indefinite integrals.'),
    ('integrate(Integral(y, (x, 1, 2)), x, risch=True)', ValueError, 'risch=True is only allowed for indefinite integrals.'),
    ('Integral(Integral(y, (x, 1, 2)), x).doit(risch=True)', ValueError, 'risch=True is only allowed for indefinite integrals.'),
    ('integrate(x, x, risch=True, manual=True)', ValueError, 'At most one of manual, meijerg, risch, heurisch can be True'),
    ("integrate(x, x, conds='bad')", ValueError, 'conds must be one of "separate", "piecewise", "none", got: bad'),
]


@pytest.mark.parametrize('call, kind, expected', VALUES)
def test_method_options_match_the_reference(call, kind, expected):
    result = eval(call)
    assert (type(result).__name__, str(result)) == (kind, expected)


@pytest.mark.parametrize('call, error, message', ERRORS)
def test_method_option_errors_match_the_reference(call, error, message):
    with pytest.raises(error) as info:
        eval(call)
    assert str(info.value).startswith(message)


def test_only_a_risch_proof_gives_a_nonelementary_integral():
    proved = integrate(exp(-x**2), x, risch=True)
    assert isinstance(proved, NonElementaryIntegral)
    # an equal integral built elsewhere is not marked
    assert not isinstance(Integral(exp(-x**2), x), NonElementaryIntegral)
    # undecided by every method: a plain unevaluated integral
    undecided = integrate(sin(x)*exp(x**2), x)
    assert undecided.is_Integral and not isinstance(undecided, NonElementaryIntegral)
    # the Risch algorithm does not handle the trigonometric case
    held = integrate((x**2 - 1)*sin(x), x, risch=True)
    assert held.is_Integral and not isinstance(held, NonElementaryIntegral)


def test_nonelementary_integral_doit_without_risch_stays_marked():
    i = NonElementaryIntegral(exp(-x**2)/(x + 1), x)
    assert isinstance(i.doit(risch=False), NonElementaryIntegral)
