"""Real conjugate pairs must survive compact principal-root presentation.

A complex-log primitive can have the correct derivative locally and still
jump by a period integral. Verify exact derivatives, period increments and
both limits at every type of apparent pole of the half-angle substitution.
"""
import pytest
import sympy as sp
import symox as sx


def test_closed_unknown_function_and_parameters_keep_their_original_coefficients():
    from symox.integrals.heurisch import _cartesian_numeric_coefficients
    z, a = sx.symbols('z a')
    unknown = sx.Function('f')(1)
    for polynomial in (z+unknown+sx.I, z+sx.sqrt(a+sx.I)):
        assert _cartesian_numeric_coefficients(polynomial, (z,)) == polynomial


@pytest.fixture(params=[(1, 3), (1, 8), (2, 3)])
def primitive(request):
    b, c = request.param
    x = sx.Symbol('x')
    result = sx.integrate(1/(b+c*sx.cos(x)**2), x)
    return b, c, sx.to_sympy(result)


def test_exact_derivative_and_period_increments(primitive):
    b, c, result = primitive
    x = sp.Symbol('x')
    t = sp.Symbol('t', real=True)
    floors = result.atoms(sp.floor)
    assert floors and not result.has(sp.I, sp.log)
    # On each interval between substitution poles, the floors are constant.
    # tan(x/2)=t and dt/dx=(1+t**2)/2 turn the derivative identity rational.
    smooth = result.xreplace(dict.fromkeys(floors, sp.S.Zero))
    smooth = smooth.subs(sp.tan(x/2), t)
    actual = sp.diff(smooth, t)*(1+t**2)/2
    expected = 1/(b+c*((1-t**2)/(1+t**2))**2)
    assert sp.cancel(actual-expected) == 0
    growth = 2*sp.pi/sp.sqrt(b*(b+c))
    assert sp.simplify(result.subs(x, 2*sp.pi)-result.subs(x, 0)-growth) == 0
    assert sp.simplify(result.subs(x, 4*sp.pi)-result.subs(x, -2*sp.pi)-3*growth) == 0


@pytest.mark.parametrize('period', [-1, 0, 1])
def test_continuous_extension_at_half_angle_poles(primitive, period):
    b, c, result = primitive
    x = sp.Symbol('x')
    point = (2*period+1)*sp.pi
    left = sp.limit(result, x, point, dir='-')
    right = sp.limit(result, x, point, dir='+')
    assert left.is_finite is right.is_finite is True
    assert sp.simplify(left-right) == 0
    # The neighboring two removable poles are one full period apart.
    next_value = sp.limit(result, x, point+2*sp.pi, dir='+')
    assert sp.simplify(next_value-right-2*sp.pi/sp.sqrt(b*(b+c))) == 0


def test_one_global_constant_on_positive_and_negative_periods(primitive):
    b, c, result = primitive
    x = sp.Symbol('x')
    scale = sp.sqrt(b*(b+c))
    # An independent tangent substitution at period pi gives this globally
    # unwrapped reference, with different apparent poles from tan(x/2).
    reference = (sp.atan(sp.sqrt(sp.Rational(b, b+c))*sp.tan(x))
                 + sp.pi*sp.floor((x-sp.pi/2)/sp.pi))/scale
    constant = sp.simplify(result.subs(x, 0)-reference.subs(x, 0))
    for period in (-3, -1, 0, 1, 3):
        for phase in (sp.Rational(1, 3), sp.Rational(7, 3)):
            point = 2*period*sp.pi+phase
            difference = (result-reference-constant).subs(x, point).evalf(55)
            assert abs(difference) < sp.Float('1e-45', 55)
