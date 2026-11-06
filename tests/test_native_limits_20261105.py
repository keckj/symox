"""Regression tests of the limits lot of 2026-11-05 (native Gruntz API,
floats of Limit.doit, Bessel functions at an imaginary infinity, the
NotImplementedError of the reference, directed coefficients)."""
import pytest

from symox import (
    I, Limit, Mod, PoleError, Rational, S, Symbol, besseli, besselj, besselk, sign as _sign,
    bessely, exp, log, oo, pi, sin, sinh, sqrt, symbols, tanh, zoo,
)
from symox import limit
from symox.series.gruntz import (
    SubsSet, compare, gruntz, limitinf, mrv, mrv_leadterm, rewrite, sign,
)

x, y, n, c, a = symbols('x y n c a')


def test_gruntz_api_compare_sign_limitinf():
    assert compare(exp(x), x**5, x) == ">"
    assert compare(exp(x**2), exp(x)**2, x) == ">"
    assert compare(exp(x), exp(x + exp(-x)), x) == "="
    assert compare(2, 1, x) == ">"
    assert sign(-x, x) == -1
    assert sign(exp(-x), x) == 1
    assert sign(1/x, x) == 1
    assert limitinf(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x) == -1
    assert limitinf(x**2/(x + 1), x) is oo


def test_gruntz_api_mrv_rewrite_leadterm():
    s, e = mrv(exp(x + exp(-x)), x)
    assert isinstance(s, SubsSet)
    assert set(s.keys()) == {exp(-x), exp(x + exp(-x))}
    w = Symbol('w')
    s, e = mrv(exp(x), x)
    f, logw = rewrite(e, s, x, w)
    assert f == 1/w and logw == -x
    assert mrv_leadterm(exp(x)/(exp(x) + 1) - 1, x) == (-1, 1)
    with pytest.raises(TypeError):
        rewrite(exp(x), {}, x, w)


def test_gruntz_entry():
    assert gruntz(exp(x)*(exp(1/x + exp(-x)) - exp(1/x)), x, oo) == 1
    assert gruntz(sin(x)/x, x, 0) == 1
    assert gruntz(1/x, x, 0, '-') is -oo
    with pytest.raises(NotImplementedError):
        gruntz(exp(x*y), x, oo)


def test_limit_floats_are_rationalized():
    # Limit.doit: nsimplify(e) when e has floats
    assert limit(x**0.5, x, 16) == 4
    assert Limit(x**0.5, x, 16).doit() == 4
    assert limit(x**0.25, x, 16) == 2
    assert limit((1 + 0.5/x)**x, x, oo) == exp(S.Half)


def test_bessel_at_imaginary_infinity():
    assert limit(bessely(1, x), x, I*oo) is -oo
    assert limit(bessely(1, x), x, -I*oo) is -oo
    assert limit(bessely(0, x), x, I*oo) == oo*I
    assert limit(besselk(1, x), x, I*oo) == 0
    assert limit(besselk(0, x), x, -I*oo) == 0
    assert limit(besselj(1, x), x, I*oo) == oo*I
    assert limit(besseli(1, x), x, I*oo) == 0


def test_reference_not_implemented_errors():
    for f in (lambda: limit(exp(x*y), x, oo),
              lambda: limit(exp(-x*y), x, oo),
              lambda: limit(log(n)**(n/log(n))/c**n, n, oo),
              lambda: limit(exp(I*x)*sin(pi*x), x, oo),
              lambda: limit(exp((2 - n)*x), x, oo),
              lambda: Limit(exp(x), x, zoo).doit(),
              lambda: Limit(x**2/(x + 1), x, zoo).doit(),
              lambda: limit(x**2/(x + 1), x, zoo)):
        with pytest.raises(NotImplementedError):
            f()


def test_limit_direction_admission():
    with pytest.raises(TypeError):
        Limit(x, x, 0, dir=0)
    with pytest.raises(ValueError):
        Limit(x, x, 0, dir='0')


def test_directed_infinity_of_unit_coefficient():
    # limitinf: c0 matching I*a (c0/I free of I) gives c0*oo
    assert limit(x**Rational(-2, 3), x, 0, '-') == -oo*(-1)**Rational(1, 3)
    assert limit((-x)**Rational(-2, 3), x, 0, '+') == -oo*(-1)**Rational(1, 3)


def test_hyperbolic_form_is_kept():
    assert limit(x/((x - 6)*sinh(tanh(0.03*x)) + tanh(x) - 0.5), x, oo) == 1/sinh(1)
    assert limit(sinh(tanh(x)), x, oo) == sinh(1)


def test_unit_power_sign_coefficient():
    # issue 8634: the coefficient exp(I*pi*n) of (-t)**n is written (-1)**n
    m = Symbol('m', integer=True, positive=True)
    assert limit(x**m, x, -oo) == oo*_sign((-1)**m)


def test_symbolic_exponent_leading_term():
    # Limit.doit leadterm step: the sign of a symbolic exponent decides
    w = Symbol('w')
    p = Symbol('p', positive=True)
    m = Symbol('m', integer=True)
    assert limit(x**(-pi), x, 0, '-') == -oo*(-1)**(1 - pi)
    assert limit(a/x**p, x, 0, '-') == -oo*(-1)**(1 - p)*_sign(a)
    assert limit(a/x**p, x, 0, '+') == oo*_sign(a)
    assert limit(2/(x - 1)**p, x, 1, '-') == -oo*(-1)**(1 - p)
    assert limit(exp(x)/x**p, x, 0, '-') == -oo*(-1)**(1 - p)
    assert limit(x**p, x, 0) == 0
    for e, z0, d, msg in ((x**w, 0, '+', 'w'), (x**w, 0, '-', 'w'),
                          (x**(3 - w), 0, '+', '3 - w'), (y*x**w, 0, '+', 'w'),
                          (sin(x)*x**w, 0, '+', 'w + 1'), (x**I, 0, '+', 'I'),
                          ((x - 2)**w, 2, '-', 'w'), (x**m, oo, '+', '-m'),
                          (2*x**m, oo, '+', '-m'), ((x + 1)**m, oo, '+', '-m')):
        with pytest.raises(NotImplementedError, match='Not sure of sign of %s$' % msg.replace('+', r'\+')):
            limit(e, x, z0, d)


def test_symbolic_exponent_at_infinity():
    # gruntz cannot decide the sign of the exponent; the heuristics rebuild
    # a lone power and hold a product
    w = Symbol('w')
    assert limit(x**w, x, oo) == oo**w
    assert limit(x**(3 - w), x, oo) == oo**(3 - w)
    assert limit((x + 1)**w, x, oo) == oo**w
    assert limit(log(x)**w, x, oo) == oo**w
    assert limit(x**w, x, -oo) == (-oo)**w
    assert limit(x**(I*w), x, oo) == oo**(I*w)
    for e in (2*x**w, y*x**w, x**w*log(x), x**w/(x + 1), 2*x**I):
        assert isinstance(limit(e, x, oo), Limit)
    with pytest.raises(NotImplementedError, match=r'Result depends on the sign of sign\(w\)$'):
        limit(exp(x)**w, x, oo)
    with pytest.raises(NotImplementedError, match=r'Result depends on the sign of sign\(w\)$'):
        limit(x**w + 1, x, oo)


def test_engine_value_is_kept_over_a_held_symbolic_power():
    # the held/error outcome of the symbolic-power steps stands only when
    # the native engine does not decide (issue 22334: 1/(k + 1))
    m = Symbol('m', positive=True)
    for k in (Symbol('k', positive=True), Symbol('k', nonnegative=True)):
        e = (m + 1)**k/((m + 1)**(k + 1) - m**(k + 1))
        assert (limit(e, m, oo) - 1/(k + 1)).simplify() == 0
    # a held Limit.doit takes the oscillating route, not the whole limit()
    j = Symbol('j', integer=True)
    with pytest.raises(NotImplementedError):
        Limit(j/(j + 1) + (I/2)**j, j, oo).doit()


def test_constant_base_pow_heuristics():
    # Limit.pow_heuristics: b**p with b free of x is exp(limit(p*log(b)))
    assert limit((I/2)**x, x, oo) == 0
    assert limit((I/2)**(x**2), x, oo) == 0
    assert limit((2*I)**x, x, oo) is zoo
    assert limit((-2)**x, x, oo) is zoo
    assert limit((1 + I)**x, x, oo) is zoo
    assert limit((I/2)**(-x), x, oo) is zoo
    assert limit((I/2)**(1/x), x, 0, '+') == 0
    assert limit(Rational(-1, 2)**x, x, oo) == 0
    assert str(limit(a**x, x, oo)) == 'exp(oo*sign(log(a)))'
    # products and sums go to gruntz, which raises
    with pytest.raises(NotImplementedError, match="Result depends on the sign"):
        limit(2*(I/2)**x, x, oo)


def test_powsimp_before_gruntz():
    # Limit.doit merges the product with powsimp before pow_heuristics and
    # gruntz
    b = Symbol('b')
    assert limit(2**x*(I/4)**x, x, oo) == 0
    assert limit(3**x*Rational(-1, 4)**x, x, oo) == 0
    for e in (exp(b*x)*exp(-x), x*exp(b*x)*exp(-x), exp(b*x - x)):
        with pytest.raises(NotImplementedError, match=r"sign of sign\(b - 1\)"):
            limit(e, x, oo)


def test_expanded_parametric_coefficient():
    # limitinf reads the sign of the leading coefficient factor by factor
    # (the reference coefficient (b - 1)/x, expanded to b/x - 1/x here)
    b = Symbol('b')
    xp = Symbol('xp', positive=True)
    assert limitinf((b*xp - xp)/log(xp), xp) == oo*_sign(b - 1)
    assert limitinf(xp*(b - 1)/log(xp), xp) == oo*_sign(b - 1)
