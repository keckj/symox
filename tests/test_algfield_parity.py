"""Algebraic fields, root objects and inexact domains compared with SymPy 1.14.

Every expected value below was produced by SymPy 1.14 and is hard-coded, so
the test does not need the reference package.
"""
import pytest
import symox as sm
from symox import (CRootOf, I, Poly, Rational as R, apart, cancel, div,
                   factor, factor_terms, gcd, gcd_terms, lcm, simplify, sqrt,
                   together)
from symox.polys.polyerrors import PolynomialError

x, y = sm.symbols('x y')
r = CRootOf(x**3 - x - 1, 0)


# Polynomial operations with a root object coefficient (extension=True).
# SymPy infers the generator x and reads CRootOf(x**3 - x - 1, 0) as an
# algebraic coefficient although its defining polynomial uses x.
@pytest.mark.parametrize('expr, expected', [
    (lambda: gcd(x**3 - x - 1, x - r, extension=True),
     'x - CRootOf(x**3 - x - 1, 0)'),
    (lambda: cancel((x**3 - x - 1)/(x - r), extension=True),
     'x**2 + x*CRootOf(x**3 - x - 1, 0) - 1 + CRootOf(x**3 - x - 1, 0)**2'),
    (lambda: div(x**3 - x - 1, x - r, extension=True),
     '(x**2 + x*CRootOf(x**3 - x - 1, 0) - 1 + CRootOf(x**3 - x - 1, 0)**2, 0)'),
    (lambda: lcm(x**2 - r**2, x - r, extension=True),
     'x**2 - CRootOf(x**3 - x - 1, 0)**2'),
    (lambda: factor(x**3 - x - 1, extension=r),
     '(x - CRootOf(x**3 - x - 1, 0))*(x**2 + x*CRootOf(x**3 - x - 1, 0) - 1 + CRootOf(x**3 - x - 1, 0)**2)'),
    (lambda: gcd(y**3 - y - 1, y - r, extension=True),
     'y - CRootOf(x**3 - x - 1, 0)'),
    (lambda: Poly(x - r, extension=True),
     "Poly(x - CRootOf(x**3 - x - 1, 0), x, domain='QQ<CRootOf(x**3 - x - 1, 0)>')"),
])
def test_root_object_coefficients(expr, expected):
    assert str(expr()) == expected


def test_explicit_generator_inside_root_object_raises():
    # with x given explicitly SymPy refuses the root object as a coefficient
    with pytest.raises(PolynomialError, match='contains an element of the set of generators'):
        gcd(x**3 - x - 1, x - r, x, extension=True)


# Exact sign and zero decisions of numeric expressions with root objects:
# evalf gives the sign, minimal_polynomial decides a zero.
@pytest.mark.parametrize('expr, expected', [
    (lambda: (r - 2).is_negative, True),
    (lambda: (r - 1).is_positive, True),
    (lambda: (r**2 - 2).is_negative, True),
    (lambda: (r**3 - r - 2).is_negative, True),
    (lambda: (r - sqrt(2)).is_positive, False),
    (lambda: (r + r).is_positive, True),
    (lambda: (r**3 - r - 1).is_zero, True),
    (lambda: (r**3 - r - 1).is_extended_positive, False),
    (lambda: (r**3 - r - 1).is_nonnegative, True),
    (lambda: (CRootOf(x**5 - x - 1, 0)**5 - CRootOf(x**5 - x - 1, 0) - 1).is_zero, True),
    (lambda: (CRootOf(x**3 - 2, 0) - 2**R(1, 3)).is_zero, True),
    (lambda: CRootOf(x**3 - x - 1, 1).is_positive, False),
    (lambda: (CRootOf(x**3 - x - 1, 1) - CRootOf(x**3 - x - 1, 2)).is_zero, None),
    (lambda: (r**3 - r - 1).equals(0), True),
    (lambda: (CRootOf(x**3 - 2, 0) - 2**R(1, 3)).equals(0), True),
])
def test_root_object_signs(expr, expected):
    assert expr() is expected


@pytest.mark.parametrize('expr, expected', [
    (lambda: simplify(r**3 - r - 1), '0'),
    (lambda: simplify(CRootOf(x**5 - x - 1, 0)**5 - CRootOf(x**5 - x - 1, 0) - 1), '0'),
    (lambda: simplify(CRootOf(x**3 - 2, 0) - 2**R(1, 3)), '0'),
    (lambda: simplify(r**3 - r - 2), '-2 - CRootOf(x**3 - x - 1, 0) + CRootOf(x**3 - x - 1, 0)**3'),
    (lambda: simplify((sqrt(2) + sqrt(3))**2 - 5 - 2*sqrt(6)), '0'),
    (lambda: sm.Max(r, sqrt(2)), 'sqrt(2)'),
    (lambda: sm.Min(r, 1), '1'),
    (lambda: r < 2, 'True'),
    (lambda: r > sqrt(2), 'False'),
    (lambda: sm.sign(r - 2), '-1'),
    (lambda: sm.Abs(r - 2), '2 - CRootOf(x**3 - x - 1, 0)'),
    (lambda: sm.floor(r*10), '13'),
    (lambda: sm.ceiling(CRootOf(x**5 - x - 1, 0)*100), '117'),
    (lambda: sm.Piecewise((1, r > 1), (0, True)), '1'),
])
def test_root_object_decisions(expr, expected):
    assert str(expr()) == expected


@pytest.mark.parametrize('expr, expected', [
    (lambda: (r - 2).evalf(30), '-0.675282042755253974039091145522'),
    (lambda: (r**2 + CRootOf(x**5 - x - 1, 0)).evalf(25), '2.922181644508111444305555'),
    (lambda: (CRootOf(x**3 - x - 1, 1) + 1).evalf(20), '0.33764102137762698702 - 0.5622795120623012439*I'),
    (lambda: (r*sm.pi).evalf(30), '4.16172420255857193206273483720'),
    (lambda: sm.N(sqrt(r), 40), '1.150963925257758035680601218461375863667'),
    (lambda: (CRootOf(x**3 - x - 1, 1)*CRootOf(x**3 - x - 1, 2)).evalf(20), '0.75487766624669276005'),
])
def test_root_object_evalf(expr, expected):
    assert str(expr()) == expected


# Extended sign queries of numbers: Expr._eval_is_extended_positive_negative
# evaluates the value and proves a radical zero with minimal_polynomial.
@pytest.mark.parametrize('expr, attr, expected', [
    ((sqrt(2) + sqrt(3))**2 - 5 - 2*sqrt(6), 'is_extended_positive', False),
    ((sqrt(2) + sqrt(3))**2 - 5 - 2*sqrt(6), 'is_extended_nonnegative', True),
    ((1 + 2**R(1, 3))**3 - 3 - 3*2**R(1, 3) - 3*2**R(2, 3), 'is_extended_negative', False),
    (sqrt(2) - 1, 'is_extended_positive', True),
    (sm.pi - 3, 'is_extended_positive', True),
    (sm.E - 3, 'is_extended_negative', True),
    (I + 1, 'is_extended_positive', False),
    (sqrt(2)*I, 'is_extended_positive', False),
    (x - 1, 'is_extended_positive', None),
    (sm.zoo + 1, 'is_extended_positive', False),
    (sm.oo - 1, 'is_extended_positive', True),
    (sm.log(2) - 1, 'is_extended_nonpositive', True),
    (sm.sin(1) - sm.cos(1), 'is_extended_positive', True),
])
def test_extended_sign_of_numbers(expr, attr, expected):
    assert getattr(expr, attr) is expected


# Integer._eval_power collects perfect roots prime by prime.
@pytest.mark.parametrize('expr, expected', [
    (lambda: 12**R(1, 6), '2**(1/3)*3**(1/6)'),
    (lambda: (-12)**R(1, 6), '(-3)**(1/6)*2**(1/3)'),
    (lambda: (-48)**R(1, 6), '(-3)**(1/6)*2**(2/3)'),
    (lambda: (-72)**R(1, 4), '(-1)**(1/4)*2**(3/4)*sqrt(3)'),
    (lambda: (-12)**R(5, 6), '2*(-3)**(5/6)*2**(2/3)'),
    (lambda: (-12)**R(-1, 6), '-(-3)**(5/6)*2**(2/3)/6'),
    (lambda: (-R(12, 5))**R(1, 6), '(-3)**(1/6)*2**(1/3)*5**(5/6)/5'),
    (lambda: 300**R(1, 6), '10**(1/3)*3**(1/6)'),
    (lambda: 12**R(5, 6), '2*2**(2/3)*3**(5/6)'),
    (lambda: 12**R(-1, 6), '2**(2/3)*3**(5/6)/6'),
    (lambda: 72**R(1, 4), '2**(3/4)*sqrt(3)'),
    (lambda: R(12, 5)**R(1, 6), '2**(1/3)*3**(1/6)*5**(5/6)/5'),
    (lambda: 18**R(1, 4), '2**(1/4)*sqrt(3)'),
    (lambda: 720**R(1, 6), '2**(2/3)*3**(1/3)*5**(1/6)'),
    (lambda: 48**R(1, 6), '2**(2/3)*3**(1/6)'),
    (lambda: 12**R(7, 6), '12*2**(1/3)*3**(1/6)'),
    (lambda: 12**R(1, 4), 'sqrt(2)*3**(1/4)'),
    (lambda: 6**R(1, 6), '6**(1/6)'),
    (lambda: 2**R(7, 3), '4*2**(1/3)'),
    (lambda: (-8)**R(1, 6), '(-1)**(1/6)*sqrt(2)'),
    (lambda: (-27)**R(1, 3), '3*(-1)**(1/3)'),
    (lambda: (-108)**R(1, 3), '3*(-1)**(1/3)*2**(2/3)'),
    (lambda: (-12)**R(1, 2), '2*sqrt(3)*I'),
    (lambda: (2**60*3**7)**R(1, 6), '3072*3**(1/6)'),
])
def test_integer_root_extraction(expr, expected):
    assert str(expr()) == expected


# A float content is a common factor of gcd_terms (Term reads
# as_coeff_mul with rational=True).
@pytest.mark.parametrize('expr, expected', [
    (lambda: factor_terms(2.0*x + 2.0*y), '2.0*(x + y)'),
    (lambda: gcd_terms(2.0*x + 2.0), '2.0*(x + 1)'),
    (lambda: together(1.0*x + 1.0), '1.0*(x + 1)'),
    (lambda: together(-0.5/(1.0*x + 1.0)), '-0.5/(x + 1)'),
    (lambda: factor_terms(-1.0*x - 1.0), '-1.0*(x + 1)'),
    (lambda: factor_terms(2.0*x + 4.0), '2.0*x + 4.0'),
    (lambda: factor_terms(2.5*x*y + 2.5*x), '2.5*x*(y + 1)'),
    (lambda: factor_terms(x*1.5 + 3*y), '1.5*x + 3*y'),
    (lambda: factor(1.0*x + 1.0), '1.0*(x + 1)'),
])
def test_float_content(expr, expected):
    assert str(expr()) == expected


# apart with Float coefficients decomposes over RR by undetermined
# coefficients instead of converting the floats to binary rationals.
@pytest.mark.parametrize('expr, expected', [
    (lambda: apart(1/(x**2 - 1.0), x), '-0.5/(x + 1) + 0.5/(x - 1)'),
    (lambda: apart(1.0/(x**2 - 3*x + 2), x), '-1.0/(x - 1) + 0.5/(0.5*x - 1.0)'),
    (lambda: apart(1/(x**2 - 1.5), x), '0.666666666666667/(0.666666666666667*x**2 - 1.0)'),
    (lambda: apart(2.5*x/(x**2 - 1), x), '1.25/(x + 1) + 1.25/(x - 1)'),
    (lambda: apart((x + 1.0)/(x*(x + 2)), x), '0.25/(0.5*x + 1.0) + 0.5/x'),
    (lambda: apart(x**3/(x**2 - 1.0), x), '1.0*x + 0.5/(x + 1) + 0.5/(x - 1)'),
    (lambda: apart(1/(x**2 + x + 0.3), x), '1.0/(1.0*x**2 + 1.0*x + 0.3)'),
    (lambda: apart(1/((x - 1.0)**2*(x + 2)), x),
     '-0.111111111111111/(x - 1) + 0.333333333333333/(x - 1)**2 + 0.0555555555555556/(0.5*x + 1.0)'),
    (lambda: apart(1.0/(x**2 + 1), x), '1.0/(x**2 + 1)'),
    (lambda: apart(0.5/(x - 1), x), '0.5/(x - 1)'),
])
def test_apart_inexact(expr, expected):
    assert str(expr()) == expected
