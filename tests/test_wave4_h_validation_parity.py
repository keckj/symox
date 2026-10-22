"""Wave-4 records where SymPy 1.14 is wrong and symox keeps the correct value.

Each test checks the symox value and an independent invariant (definition,
mpmath at high precision, equilibrium or normalization), never the reference.
"""
import mpmath
import pytest

import symox as S
from symox import Min, MatrixSymbol, Rational, Symbol, sqrt, symbols


def _num(value, dps=30):
    return complex(S.N(value, dps))


def test_diagonalof_length_is_min_of_dimensions():
    n = Symbol('n', integer=True)
    length = S.DiagonalOf(MatrixSymbol('x', n, 3)).diagonal_length
    assert length == Min(3, n)
    # A 2 x 3 matrix has a diagonal of length 2, not 3.
    assert length.subs(n, 2) == 2
    assert S.DiagonalOf(MatrixSymbol('x', 2, 3)).diagonal_length == 2


def test_hn1_rewrite_hankel1_uses_half_integer_order():
    nu = Symbol('nu', integer=True)
    z = Symbol('z')
    value = S.hn1(nu, z).rewrite(S.hankel1).subs({nu: 1, z: 2})
    with mpmath.workdps(40):
        zz = mpmath.mpf(2)
        expected = mpmath.sqrt(mpmath.pi / (2 * zz)) * mpmath.hankel1(1.5, zz)
        closed = -mpmath.exp(1j * zz) * (zz + 1j) / zz**2
        assert abs(expected - closed) < mpmath.mpf(10)**-35
        expected = complex(expected)
    assert abs(_num(value) - expected) < 1e-25
    assert abs(_num(S.hn1(1, 2)) - expected) < 1e-25


def test_hn2_rewrite_hankel2_uses_half_integer_order():
    nu = Symbol('nu', integer=True)
    z = Symbol('z')
    form = S.hn2(nu, z).rewrite(S.hankel2)
    value = form.subs({nu: 2, z: Rational(7, 5)})
    with mpmath.workdps(40):
        zz = mpmath.mpf(7) / 5
        scale = mpmath.sqrt(mpmath.pi / (2 * zz))
        expected = scale * (mpmath.besselj(2.5, zz) - 1j * mpmath.bessely(2.5, zz))
        assert abs(expected - scale * mpmath.hankel2(2.5, zz)) < mpmath.mpf(10)**-35
        expected = complex(expected)
    assert abs(_num(value) - expected) < 1e-25
    assert abs(_num(S.hn2(2, Rational(7, 5))) - expected) < 1e-25


def test_cable_with_inconsistent_fixed_load_points_has_no_solution():
    from symox.physics.continuum_mechanics.cable import Cable
    # With vertical loads only the horizontal tension H is constant. The
    # node equilibria give H = 3/(1.37 - 0.63) at Z and H = 8/(0.63 + 2/3)
    # at X, which differ, so no equilibrium exists for these fixed points.
    slopes = [(7.26 - 10) / 2, (6 - 7.26) / 2, (10 - 6) / 6]
    h_z = 3 / (slopes[1] - slopes[0])
    h_x = 8 / (slopes[2] - slopes[1])
    assert abs(h_z - h_x) > 2
    c = Cable(('A', 0, 10), ('B', 10, 10))
    c.apply_load(-1, ('Z', 2, 7.26, 3, 270))
    c.apply_load(-1, ('X', 4, 6, 8, 270))
    with pytest.raises(ValueError):
        c.solve()


def test_cable_distributed_load_reactions_are_in_equilibrium():
    from symox.physics.continuum_mechanics.cable import Cable
    c = Cable(('A', 0, 40), ('B', 100, 20))
    c.apply_load(0, ('X', 850))
    c.solve(58.58)
    loads = {str(k): float(v) for k, v in c.reaction_loads.items()}
    # Parabolic cable with lowest point at x0: 40 - 20 = w/(2H)*(x0**2 - (100 - x0)**2).
    w, x0 = 850.0, 58.58
    H = w * (x0**2 - (100 - x0)**2) / (2 * 20)
    assert abs(loads['R_A_x'] + loads['R_B_x']) < 1e-6
    assert abs(abs(loads['R_A_x']) - H) < 1e-3
    assert abs(loads['R_A_y'] + loads['R_B_y'] - w * 100) < 1e-3
    assert abs(loads['R_A_y'] - w * x0) < 1e-3
    assert abs(loads['R_B_y'] - w * (100 - x0)) < 1e-3


def test_gaussian_conj_matches_abcd_propagation():
    from symox.physics.optics import gaussian_conj
    s_out, z_out, m = gaussian_conj(3, 2, 5)
    # q at the lens for a waist at distance s before it, then a thin lens f.
    q = 3 + 2j
    q_out = 1 / (1 / q - 1 / 5)
    assert abs(float(s_out) - (-q_out.real)) < 1e-12
    assert z_out == Rational(25, 4)
    assert abs(float(z_out) - q_out.imag) < 1e-12
    assert abs(float(m) - (q_out.imag / 2) ** 0.5) < 1e-12


def test_beta_noncentral_cdf_is_one_beyond_support():
    from symox.stats import BetaNoncentral, cdf
    value = cdf(BetaNoncentral('x', 1, 1, 1), evaluate=False)(2).doit()
    assert value == 1


def _check_cdf(cdf_expr, var, pdf, lower, upper, outside_low, outside_high, inside):
    for point in outside_low:
        assert abs(_num(cdf_expr.subs(var, point))) < 1e-25
    for point in outside_high:
        assert abs(_num(cdf_expr.subs(var, point)) - 1) < 1e-25
    with mpmath.workdps(30):
        for point in inside:
            expected = mpmath.quad(pdf, [lower, point])
            assert abs(_num(cdf_expr.subs(var, point)) - complex(expected)) < 1e-20
        assert abs(mpmath.quad(pdf, [lower, upper]) - 1) < mpmath.mpf(10)**-25


def test_bounded_pareto_cdf_is_bounded():
    from symox.stats import BoundedPareto, cdf
    L, H, x = symbols('L H x', positive=True)
    form = cdf(BoundedPareto('X', 2, L, H))(x).subs({L: 1, H: 3})
    pdf = lambda t: 2 * t**-3 / (1 - mpmath.mpf(1) / 9)
    _check_cdf(form, x, pdf, 1, 3, [Rational(1, 2)], [3, 4, 10],
               [Rational(3, 2), 2, Rational(5, 2)])


def test_power_function_cdf_is_bounded():
    from symox.stats import PowerFunction, cdf
    a, b = symbols('a b', real=True)
    z = Symbol('z')
    form = cdf(PowerFunction('X', 2, a, b))(z).subs({a: 1, b: 4})
    pdf = lambda t: 2 * (t - 1) / 9
    _check_cdf(form, z, pdf, 1, 4, [0, Rational(1, 2)], [4, 5],
               [2, 3, Rational(7, 2)])


def test_reciprocal_cdf_is_bounded():
    from symox.stats import Reciprocal, cdf
    a, b, x = symbols('a b x', positive=True)
    form = cdf(Reciprocal('R', a, b))(x).subs({a: 1, b: 5})
    pdf = lambda t: 1 / (t * mpmath.log(5))
    _check_cdf(form, x, pdf, 1, 5, [Rational(1, 2)], [5, 7],
               [2, 3, Rational(9, 2)])


def test_levy_cdf_is_zero_below_location():
    from symox.stats import Levy, cdf
    mu = Symbol('mu', real=True)
    c = Symbol('c', positive=True)
    z = Symbol('z')
    form = cdf(Levy('x', mu, c))(z).subs({mu: 1, c: 2})
    for point in (Rational(1, 2), 0, -3):
        value = _num(form.subs(z, point))
        assert value == 0
    with mpmath.workdps(30):
        pdf = lambda t: mpmath.sqrt(2 / (2 * mpmath.pi)) * mpmath.exp(-2 / (2 * (t - 1))) / (t - 1) ** 1.5
        for point in (2, 5):
            expected = mpmath.quad(pdf, [1, 1.5, point])
            assert abs(_num(form.subs(z, point)) - complex(expected)) < 1e-15


def test_logit_normal_cdf_is_real_outside_unit_interval():
    from symox.stats import LogitNormal, cdf
    mu = Symbol('mu', real=True)
    s = Symbol('s', positive=True)
    z = Symbol('z')
    form = cdf(LogitNormal('x', mu, s))(z).subs({mu: 0, s: 1})
    pdf = lambda t: mpmath.exp(-mpmath.log(t / (1 - t))**2 / 2) / (mpmath.sqrt(2 * mpmath.pi) * t * (1 - t))
    _check_cdf(form, z, pdf, 0, 1, [-1, Rational(-1, 2)], [1, 2],
               [Rational(1, 10), Rational(1, 2), Rational(9, 10)])


def test_dirichlet_two_component_marginal_is_beta():
    from symox.stats import MultivariateBeta, marginal_distribution
    a1, a2 = symbols('a1 a2', positive=True)
    x = Symbol('x')
    C = MultivariateBeta('C', a1, a2)
    marginal = marginal_distribution(C, C[0])(x)
    beta_pdf = x**(a1 - 1) * (1 - x)**(a2 - 1) / S.beta(a1, a2)
    for values in ({a1: 2, a2: 3}, {a1: Rational(1, 2), a2: Rational(5, 2)}):
        for point in (Rational(1, 5), Rational(1, 2), Rational(4, 5)):
            got = _num(marginal.subs(values).subs(x, point))
            want = _num(beta_pdf.subs(values).subs(x, point))
            assert abs(got - want) < 1e-25
    with mpmath.workdps(30):
        f = S.lambdify(x, marginal.subs({a1: 2, a2: 3}), 'mpmath')
        assert abs(mpmath.quad(f, [0, 1]) - 1) < mpmath.mpf(10)**-20


def test_gue_density_is_normalized():
    from symox.stats import GaussianUnitaryEnsemble, density
    G = GaussianUnitaryEnsemble('U', 2)
    X = MatrixSymbol('X', 2, 2)
    value = density(G)(X)
    assert S.simplify(value - 2 * S.exp(-S.Trace(X**2)) / S.pi**2) == 0
    # Tr H**2 = h11**2 + h22**2 + 2*(a**2 + b**2) for H = [[h11, a+ib], [a-ib, h22]].
    with mpmath.workdps(30):
        weight = (mpmath.sqrt(mpmath.pi) ** 2) * (mpmath.sqrt(mpmath.pi / 2) ** 2)
        assert abs(weight - mpmath.pi**2 / 2) < mpmath.mpf(10)**-25
        assert abs(weight * 2 / mpmath.pi**2 - 1) < mpmath.mpf(10)**-25


def test_implicit_region_without_regular_point_raises():
    from symox.abc import x, y, z
    from symox.vector import ImplicitRegion
    equation = (x + y + z)**4
    # The gradient 4*(x + y + z)**3 vanishes wherever the equation holds.
    for derivative in (S.diff(equation, v) for v in (x, y, z)):
        assert S.factor(derivative) == 4 * (x + y + z)**3
    with pytest.raises(NotImplementedError):
        ImplicitRegion((x, y, z), equation).regular_point()
    point = ImplicitRegion((x, y), x**2 + (y - 3)**2 - 4)  # circle: regular everywhere
    p = point.regular_point()
    assert S.simplify(p[0]**2 + (p[1] - 3)**2 - 4) == 0
