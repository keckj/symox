"""Multiprecision special functions follow SymPy's mpmath evaluation."""

import mpmath
import pytest
import sympy as sp

import symox as sm


def build(module, source):
    namespace = {name: getattr(module, name) for name in dir(module) if not name.startswith("_")}
    namespace["R"] = module.Rational
    return eval(source, namespace)


SAME_AS_REFERENCE = [
    ("meijerg([[], []], [[R(1,3)], []], R(2,3))", 50),
    ("meijerg([[1],[]],[[R(1,2)],[0]], 3)", 50),
    ("meijerg([[1,1],[]],[[1],[0]], R(1,2))", 40),
    ("hyper([1,2],[3],R(1,2))", 50),
    ("hyper([1,1],[2],3)", 40),
    ("hyper([R(1,2)],[R(3,2),R(3,2)],-4)", 40),
    ("besselj(2, 3)", 50),
    ("besselj(R(1,3), -2-I)", 40),
    ("bessely(R(1,3), 2)", 50),
    ("besseli(1, 5)", 50),
    ("besselk(2, R(1,2))", 50),
    ("hankel1(R(1,3), 2)", 30),
    ("jn(2, 3)", 30),
    ("yn(3, R(1,2))", 30),
    ("airyai(2)", 50),
    ("airybi(-3)", 50),
    ("airyaiprime(1)", 50),
    ("airybiprime(1)", 50),
    ("loggamma(2+3*I)", 30),
    ("gamma(1+I)", 30),
    ("digamma(R(1,3))", 30),
    ("polygamma(2, R(1,2))", 40),
    ("polygamma(R(1,2), 3)", 30),
    ("erf(2+I)", 30),
    ("erfc(30)", 50),
    ("erfc(100)", 20),
    ("erfi(2)", 50),
    ("expint(2, 3)", 50),
    ("Ei(-1)", 40),
    ("li(3)", 50),
    ("Li(5)", 30),
    ("Si(2)", 50),
    ("Ci(2)", 50),
    ("Shi(2)", 50),
    ("Chi(2)", 50),
    ("fresnels(2)", 50),
    ("fresnelc(2)", 50),
    ("zeta(3)", 50),
    ("zeta(2, R(1,3))", 50),
    ("dirichlet_eta(3)", 50),
    ("polylog(3, -2)", 40),
    ("lerchphi(R(1,2), 2, R(1,3))", 30),
    ("elliptic_k(3)", 40),
    ("elliptic_e(R(1,3))", 40),
    ("elliptic_f(1, R(1,2))", 50),
    ("uppergamma(R(1,3), 100)", 40),
    ("harmonic(R(1,2), 3)", 30),
    ("bernoulli(R(1,2), R(1,3))", 30),
    ("euler(R(1,2))", 30),
    ("Ynm(2, 1, R(1,2), R(1,3))", 30),
    ("appellf1(1, 1, 1, 2, R(1,2), R(1,3))", 30),
    ("LambertW(3+I)", 30),
]


@pytest.mark.parametrize("source, digits", SAME_AS_REFERENCE)
def test_special_values_match_the_reference(source, digits):
    expected = build(sp, source).evalf(digits)
    actual = build(sm, source).evalf(digits)
    assert str(actual) == str(expected)


def test_meijer_g_gives_fifty_correct_digits():
    value = sm.meijerg([[], []], [[sm.Rational(1, 3)], []], sm.Rational(2, 3)).evalf(50)
    assert str(value) == "0.44851116544806341805071002975596652025901404786260"


# Conditioned arguments: SymPy rounds the argument once and loses digits.
# The native evaluation refines until two precisions agree.
CONDITIONED = [
    ("elliptic_k(1 - R(1,10**30))", lambda: mpmath.ellipk(1 - mpmath.mpf(10) ** -30)),
    ("zeta(1+R(1,10**20))", lambda: mpmath.zeta(1 + mpmath.mpf(10) ** -20)),
    ("erfinv(1-R(1,10**30))", lambda: mpmath.erfinv(1 - mpmath.mpf(10) ** -30)),
    ("elliptic_pi(2, R(1,2))", lambda: mpmath.ellippi(2, mpmath.mpf(1) / 2)),
    ("LambertW(-exp(-1)+R(1,10**20))", lambda: mpmath.lambertw(-mpmath.exp(-1) + mpmath.mpf(10) ** -20)),
]


@pytest.mark.parametrize("source, reference", CONDITIONED)
def test_conditioned_arguments_keep_the_requested_digits(source, reference):
    digits = 30
    value = sm.to_sympy(build(sm, source).evalf(digits))
    with mpmath.workdps(digits + 60):
        exact = mpmath.mpmathify(reference())
        approximation = value._to_mpmath(mpmath.mp.prec)
        assert abs(approximation - exact) <= abs(exact) * mpmath.mpf(10) ** (1.5 - digits)


@pytest.mark.parametrize("source, digits", [
    ("meijerg([[], []], [[R(1,3)], []], exp_polar(3*I*pi)*2)", 30),
    ("meijerg([[], []], [[R(1,3)], []], exp_polar(-5*I*pi/2)*R(1,2))", 30),
])
def test_meijer_g_on_a_polar_branch(source, digits):
    expected = build(sp, source).evalf(digits)
    actual = build(sm, source).evalf(digits)
    assert abs(sm.to_sympy(actual) - expected) < sp.Float(10) ** (1 - digits) * abs(expected)
    assert str(sp.re(sm.to_sympy(actual))) == str(sp.re(expected))


def test_exp_polar_evaluates_on_the_principal_sheet_only():
    assert str(sm.exp_polar(sm.I * sm.pi).evalf()) == str(sp.exp_polar(sp.I * sp.pi).evalf())
    assert str(sm.exp_polar(3 * sm.I).evalf()) == str(sp.exp_polar(3 * sp.I).evalf())
    assert str(sm.exp_polar(4 * sm.I).evalf()) == str(sp.exp_polar(4 * sp.I).evalf())


@pytest.mark.parametrize("source", [
    "erfinv(R(27,10)).evalf()",
    "erfinv(Float(2.7))",
    "erfinv(2+I).evalf()",
    "zeta(Float(1.0))",
    "lerchphi(2,1,0).evalf()",
    "loggamma(Float(-2.0))",
])
def test_mpmath_errors_propagate_as_in_the_reference(source):
    with pytest.raises(ValueError):
        build(sp, source)
    with pytest.raises(ValueError):
        build(sm, source)


def test_assumption_queries_do_not_raise_evaluation_errors():
    assert build(sm, "erfinv(R(27,10)).is_positive") is build(sp, "erfinv(R(27,10)).is_positive")


FLOAT_ARGUMENTS = [
    "airyai(Float('13.1'))", "airybi(Float('-1.3'))", "li(Float('0.3'))", "li(Float('-1.3'))",
    "elliptic_e(Float('2.7'))", "besselj(2, Float('2.7'))", "zeta(Float('0.3'))", "Ei(Float('-1.3'))",
    "polylog(2, Float('1.0'))", "besselj(2, Float('3.0', 50))", "zeta(Float('3', 50))",
    "erf(Float('0.5', 40))", "gamma(Float('0.3', 50))",
    "erf(Float('1.14')*I)", "gamma(2.0*I)", "erf(1.0*pi)", "besselj(1, 1 + 2.0*I)",
    "lowergamma(R(1,2), -Float('1.3'))",
]


@pytest.mark.parametrize("source", FLOAT_ARGUMENTS)
def test_float_arguments_evaluate_like_the_reference(source):
    assert str(build(sm, source)) == str(build(sp, source))


@pytest.mark.parametrize("source", [
    "sqrt(-Float('1.3'))", "(-Float('1.3'))**R(3,2)", "(-Float('1.3'))**R(1,3)", "(-Float('1.3'))**R(2,3)",
    "(-Float('1.3'))**Float('0.5')", "Integer(-2)**Float('0.5')", "(-Float('1.3'))**Float('-0.5')",
    "Float('-2.0')**R(1,4)", "Float('-1.3', 30)**R(1,3)", "(-Float('1.3'))**Float('0.3')",
])
def test_negative_float_powers_follow_the_reference(source):
    assert str(build(sm, source)) == str(build(sp, source))
