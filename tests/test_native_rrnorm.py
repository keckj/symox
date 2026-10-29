"""SymPy 1.14 parity for the RR[y] gcd/sqf_list scale, the numeric sign of
functions evaluated by evalf, and log forms in solve/simplify.

Every expected value below is hard-coded from SymPy 1.14.
"""
from symox import (Abs, I, Poly, Rational, S, Symbol, exp, exp_polar, gcd,
                   hyper, hyperexpand, lcm, log, pi, polylog, simplify, solve,
                   sqf_list, symbols, besselj, Ei, li, lerchphi, zeta, asin,
                   hyper as _hyper)

x, y, z = symbols('x y z')


def _s(v):
    return str(v)


def test_rr_ring_gcd_lcm_cofactors_scale():
    p1 = Poly(1.5*x*z + 1.5*y, x, z)
    p2 = Poly(3.0*x*z + 3.0*y, x, z)
    assert _s(p1.gcd(p2)) == "Poly(1.5*x*z + 1.5*y, x, z, domain='RR[y]')"
    assert _s(p1.lcm(p2)) == "Poly(3.0*x*z + 3.0*y, x, z, domain='RR[y]')"
    assert _s(p1.cofactors(p2)) == ("(Poly(1.5*x*z + 1.5*y, x, z, domain='RR[y]'), "
                                    "Poly(1.0, x, z, domain='RR[y]'), Poly(2.0, x, z, domain='RR[y]'))")
    assert _s(gcd(1.5*x*z + 1.5*y, 3.0*x*z + 3.0*y, x, z)) == "1.5*x*z + 1.5*y"
    assert _s(lcm(1.5*x*z + 1.5*y, 3.0*x*z + 3.0*y, x, z)) == "3.0*x*z + 3.0*y"
    f = Poly(0.5*x**2*z - 0.5*y**2, x, z)
    g = Poly(1.5*x*z + 1.5*y, x, z)
    assert _s(f.gcd(g)) == "Poly(0.5, x, z, domain='RR[y]')"
    assert _s(f.lcm(g)) == ("Poly(1.5*x**3*z**2 + 1.5*y*x**2*z - 1.5*y**2*x*z - 1.5*y**3, "
                            "x, z, domain='RR[y]')")
    assert _s(Poly(2.0*x*z*y + 4.0*y, x, z).gcd(Poly(3.0*x*z*y + 6.0*y**2, x, z))) == \
        "Poly(1.0*y, x, z, domain='RR[y]')"


def test_rr_ring_sqf_list():
    assert _s(Poly(1.5*x*z + 1.5*y*x, x, z).sqf_list()) == \
        "(1.00000000000000, [(Poly(1.5*x*z + 1.5*y*x, x, z, domain='RR[y]'), 1)])"
    assert _s(Poly(1.5*(x*z + y)**2*(x + 1), x, z).sqf_list()) == (
        "(1.00000000000000, [(Poly(1.0*x + 1.0, x, z, domain='RR[y]'), 1), "
        "(Poly(1.0*x*z + 1.0*y, x, z, domain='RR[y]'), 2)])")
    assert _s(Poly(3.0*y*(x*z + y)**2, x, z).sqf_list()) == \
        "(1.0*y, [(Poly(1.0*x*z + 1.0*y, x, z, domain='RR[y]'), 2)])"


def test_rr_ring_integrate_domain():
    p1 = Poly(1.5*x*z + 1.5*y, x, z)
    assert _s(p1.integrate(x)) == "Poly(0.75*x**2*z + 1.5*y*x, x, z, domain='RR(y)')"
    assert _s(p1.integrate(z)) == "Poly(0.75*x*z**2 + 1.5*y*z, x, z, domain='RR(y)')"
    assert _s(Poly(1.5*x + y, x).integrate(x)) == "Poly(0.75*x**2 + 1.0*y*x, x, domain='RR(y)')"


def test_numeric_extended_sign():
    e = polylog(3, -Rational(1, 2))
    assert (e.is_extended_negative, e.is_extended_positive, e.is_extended_real,
            e.is_real, e.is_negative, e.is_finite, e.is_zero) == \
        (True, False, True, None, None, None, False)
    e = polylog(3, Rational(1, 2))
    assert (e.is_extended_positive, e.is_extended_negative, e.is_negative) == (True, False, False)
    e = besselj(1, 4)
    assert (e.is_extended_negative, e.is_extended_nonpositive, e.is_real) == (True, True, None)
    assert Ei(-1).is_extended_negative is True
    assert li(Rational(1, 2)).is_extended_negative is True
    assert lerchphi(Rational(1, 2), 2, 3).is_extended_positive is True
    e = asin(2)
    assert (e.is_positive, e.is_negative) == (False, False)
    e = polylog(2, 3)
    assert (e.is_extended_negative, e.is_extended_positive, e.is_extended_real) == (None, None, None)
    e = _hyper([1, 2], [3], Rational(1, 2))
    assert (e.is_extended_negative, e.is_extended_positive) == (None, None)
    assert zeta(Rational(1, 2)).is_extended_negative is True


def test_hyperexpand_lerchphi_abs():
    # sympy/simplify/tests/test_hyperexpand.py::test_lerchphi, line 623
    lhs = hyperexpand(hyper([S.Half, S.Half, S.Half, 1], [Rational(3, 2)]*3, Rational(1, 4)))
    rhs = Abs(-polylog(3, exp_polar(I*pi)/2) + polylog(3, S.Half))
    assert _s(lhs) == "-polylog(3, -1/2) + polylog(3, 1/2)"
    assert _s(rhs) == "-polylog(3, -1/2) + polylog(3, 1/2)"
    assert lhs == rhs


def test_nested_exp_solve_log_forms():
    def tower(n, shift=0):
        e = x
        for _ in range(n):
            e = exp(e) + shift
        return e
    assert _s(solve(tower(2) + 1, x)) == "[log(I*pi)]"
    assert _s(solve(tower(3) + 1, x)) == "[log(log(I*pi))]"
    assert _s(solve(tower(1) + 1, x)) == "[I*pi]"
    assert _s(solve(tower(2, 1) - 5, x)) == "[log(-1 + 2*log(2))]"
    assert _s(solve(tower(3, 1) - 5, x)) == "[log(1 - log(-1 + log(4))) + I*pi]"
    assert _s(solve(exp(x) - log(4), x)) == "[log(2*log(2))]"
    assert _s(solve(exp(x) - 2*I, x)) == "[log(2) + I*pi/2]"


def test_simplify_log_reference_forms():
    assert _s(simplify(log(-1 + 2*log(2)))) == "log(-1 + 2*log(2))"
    assert _s(simplify(log(2*log(2)))) == "log(2*log(2))"
    assert _s(simplify(2*log(2))) == "log(4)"
    assert _s(simplify(log(I*pi))) == "log(I*pi)"
