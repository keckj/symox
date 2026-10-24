"""Analytic continuation on and near cuts, branch points and poles.

Expected values are SymPy 1.14 results.  Where symox intentionally differs
(SymPy loses the value or raises next to a pole), the value is checked
against mpmath at higher precision instead.
"""
import mpmath
import pytest

from symox import *  # noqa: F401,F403
from symox import (I, Float, N, Rational, S, Symbol, elliptic_k, erfinv, exp_polar, gamma,
                   loggamma, meijerg, pi, zeta, zoo, oo, nan)

try:
    del test  # noqa: F821
except NameError:
    pass

x = Symbol('x')


def _parts(value):
    re, im = value.as_real_imag()
    return [mpmath.mpf(str(re)), mpmath.mpf(str(im))]


def _close(value, expected, dps):
    with mpmath.workdps(dps + 20):
        got = _parts(value)
        want = _parts(S(expected))
        scale = max(abs(w) for w in want) or 1
        tol = mpmath.mpf(10) ** (2 - dps)
        return all(abs(g - w) <= tol * scale for g, w in zip(got, want))


# (expression, SymPy 1.14 value at 30 digits)
NEAR_CUT = [
    ("root(-2 + I/10**30, 3)", "0.629960524947436582383605303639 + 1.09112363597172140356007261419*I"),
    ("root(-2 - I/10**30, 5, -1)", "-0.354967313104630125990361345903 - 1.09247705577745372665759105996*I"),
    ("(-2 - I/10**30)**(Rational(1, 3) + I)", "27.3471127046682687025489243607 - 10.1081869797889872729116946463*I"),
    ("asinh(2*I + Rational(1, 10**30))", "1.31695789692481670862504634731 + 1.57079632679489661923132169164*I"),
    ("asinh(2*I - Rational(1, 10**30))", "-1.31695789692481670862504634731 + 1.57079632679489661923132169164*I"),
    ("acosh(-2 - I/10**30)", "1.31695789692481670862504634731 - 3.14159265358979323846264338328*I"),
    ("atanh(2 + I/10**30)", "0.549306144334054845697622618461 + 1.57079632679489661923132169164*I"),
    ("atanh(2 - I/10**30)", "0.549306144334054845697622618461 - 1.57079632679489661923132169164*I"),
    ("acoth(Rational(1, 2) - I/10**30)", "0.549306144334054845697622618461 + 1.57079632679489661923132169164*I"),
    ("asech(-Rational(1, 2) + I/10**30)", "1.31695789692481670862504634731 - 3.14159265358979323846264338328*I"),
    ("acot(I/2 - Rational(1, 10**30))", "-1.57079632679489661923132169164 - 0.549306144334054845697622618461*I"),
    ("LambertW(-exp(-1) - I/10**30, -1)", "-3.08884301561304385595708671677 - 7.46148928565425455690611661218*I"),
    ("polylog(3, 5 - I/10**30)", "4.80534410296559047335545852397 - 4.06881863619625864550808108922*I"),
    ("lerchphi(3 + I/10**30, 2, Rational(1, 2))", "4.30830591009261429471108281005 + 1.99266227072611637422301777311*I"),
    ("hyper([Rational(1, 2), Rational(1, 2)], [1], 2 - I/10**30)",
     "0.834626841674073186281429732799 - 0.834626841674073186281429732799*I"),
    ("besselk(Rational(1, 3), -2 - I/10**30)", "0.0582724806480826243794713144524 + 6.88294639550414978305930387339*I"),
    ("bessely(1, -2 + I/10**30)", "0.107032431540937546888370772278 - 1.15344961551374677440489648454*I"),
    ("expint(Rational(1, 3), -2 + I/10**30)", "-4.44762891495433434096582625718 - 0.738755044800734611294076736024*I"),
    ("Ei(-2 - I/10**30)", "-0.048900510708061119567239835228 - 3.14159265358979323846264338328*I"),
    ("li(-2 + I/10**30)", "0.0355322759135604677869013496219 + 3.73515048255125463977638213891*I"),
    ("elliptic_k(2 - I/10**30)", "1.31102877714605990523241979495 - 1.31102877714605990523241979495*I"),
    ("elliptic_f(1, 3 - I/10**30)", "1.00107738045610623607965958638 - 0.728453936819176796102799766312*I"),
    ("legendre(Rational(1, 2), -2 - I/10**30)", "-0.142611928112603630781604769358 - 1.32913816218535775563582817623*I"),
    ("loggamma(-Rational(201, 2) - I/10**30)", "-364.900968309427351822756570463 + 317.300858012569117084726981711*I"),
    ("zeta(Rational(1, 2), -Rational(3, 2))", "-0.604898643421630370247265914236 - 2.23071014330082108153411674911*I"),
    ("gamma(-3 + I/10**20)", "-0.209352944738633412121136873875 + 16666666666666666666.6666666667*I"),
    ("LambertW(Rational(1, 10**40)) - Rational(1, 10**40)", "-1.00000000000000000000000000000e-80"),
]


@pytest.mark.parametrize("src, expected", NEAR_CUT)
def test_near_cut_values_match_sympy(src, expected):
    value = N(eval(src), 30)
    assert _close(value, expected, 30), (src, value)


# Float arguments on a singularity: SymPy maps mpmath inf/nan to oo, -oo,
# zoo or nan (Expr._from_mpmath, gamma intlike, loggamma is_integer).
SINGULAR = [
    ("erfinv(1.0)", oo), ("erfinv(-1.0)", -oo), ("elliptic_k(1.0)", oo), ("li(1.0)", -oo),
    ("gamma(-3.0)", zoo), ("gamma(0.0)", zoo), ("gamma(-100.0)", zoo), ("gamma(Float(-3, 30))", zoo),
    ("N(digamma(-3.0))", nan), ("N(polygamma(1, -3.0))", nan),
    ("atanh(1.0)", oo), ("atanh(-1.0)", -oo), ("acoth(1.0)", oo), ("acoth(-1.0)", -oo),
    ("atan(1.0*I)", oo*I), ("atan(-1.0*I)", -oo*I), ("acot(1.0*I)", -oo*I),
    ("loggamma(0.0)", oo), ("N(expint(1, 0.0))", oo), ("besselk(0, 0.0)", oo), ("bessely(1, 0.0)", zoo),
    ("Ei(0.0)", -oo), ("elliptic_pi(1.0, 0.5)", oo), ("polylog(1, 1.0)", zoo), ("LambertW(0.0, -1)", -oo),
    ("N(1 + erfinv(1.0))", oo), ("N(2*elliptic_k(1.0))", oo), ("N(li(1.0) + 1)", -oo),
    ("N(besselj(0, 1.0)*erfinv(1.0))", oo),
    ("Float(0.0)**-1", zoo), ("Float(0.0)**-2", zoo), ("Float(0.0)**Rational(-1, 2)", zoo),
    ("Float(0.0)**(-1.0)", zoo),
]


@pytest.mark.parametrize("src, expected", SINGULAR)
def test_float_singularities_match_sympy(src, expected):
    value = eval(src)
    assert value == expected, (src, value)


def test_float_singularities_through_subs():
    x = Symbol('x')
    assert erfinv(x).evalf(subs={x: 1.0}) == oo
    assert erfinv(x).subs(x, 1.0) == oo
    assert atanh(x).evalf(subs={x: 1.0}) == oo


@pytest.mark.parametrize("src, message", [
    ("loggamma(-3.0)", "gamma function pole"),
    ("N(zeta(1 + x), subs={x: 0.0})", r"zeta\(1\) pole"),
    ("zeta(1.0)", r"zeta\(1\) pole"),
    ("erfinv(2.7)", "defined only for"),
    ("N(erfinv(Rational(27, 10)))", "defined only for"),
])
def test_singular_errors_match_sympy(src, message):
    with pytest.raises(ValueError, match=message):
        eval(src)


# Intentional corrections: exact arguments next to a pole.  SymPy evaluates
# the argument at the working precision only, so it lands on the pole and
# raises or returns an infinity; symox refines the argument and returns the
# value, checked here against mpmath at 60 digits.
def _mp(fn):
    with mpmath.workdps(60):
        return fn(mpmath.mpf(10) ** -30)


@pytest.mark.parametrize("expr, dps, reference", [
    (lambda: zeta(1 + Rational(1, 10**30)), 15, lambda t: mpmath.zeta(1 + t)),
    (lambda: elliptic_k(1 - Rational(1, 10**30)), 15, lambda t: mpmath.ellipk(1 - t)),
    (lambda: erfinv(-1 + Rational(1, 10**30)), 15, lambda t: mpmath.erfinv(-1 + t)),
    (lambda: gamma(-3 + Rational(1, 10**30)), 30, lambda t: mpmath.gamma(-3 + t)),
    (lambda: loggamma(-3 + Rational(1, 10**30)), 15, lambda t: mpmath.loggamma(-3 + t)),
    # SymPy: oo at 15 digits, 6.60158062235514296... (wrong digits) at 30
    (lambda: erfinv(1 - Rational(1, 10**20)), 15, lambda t: mpmath.erfinv(1 - t*10**10)),
    (lambda: erfinv(1 - Rational(1, 10**20)), 30, lambda t: mpmath.erfinv(1 - t*10**10)),
])
def test_exact_near_pole_values(expr, dps, reference):
    value = N(expr(), dps)
    want = _mp(reference)
    with mpmath.workdps(dps + 20):
        got = mpmath.mpc(*_parts(value))
        assert abs(got - want) <= mpmath.mpf(10) ** (1 - dps) * abs(want)


def test_exact_near_pole_values_through_subs():
    x = Symbol('x')
    with mpmath.workdps(60):
        want = mpmath.erfinv(1 - mpmath.mpf(10) ** -30)
        want_gamma = mpmath.gamma(-3 + mpmath.mpf(10) ** -25)
    got = erfinv(x).evalf(subs={x: 1 - Rational(1, 10**30)})
    assert abs(mpmath.mpf(str(got)) - want) <= 1e-13 * abs(want)
    got = N(gamma(-3 + x), subs={x: Rational(1, 10**25)})
    assert abs(mpmath.mpf(str(got)) - want_gamma) <= 1e-13 * abs(want_gamma)


def test_meijerg_polar_branch_is_stable_across_precisions():
    # G(z) = gamma(1/3)*z**(1/3)*(1 + z)**(-1/3); on exp_polar(k*I*pi)*2 the
    # factor (1 + z) sits on its cut.  SymPy agrees for k = 1, 3, -3, -5
    # (principal side); for other k its side follows mpmath rounding noise
    # and changes with the precision (k = 7), symox keeps the principal side.
    for k, phase in [(1, 0), (3, Rational(2, 3)), (-3, Rational(2, 3)), (-5, 0), (5, Rational(4, 3))]:
        expr = meijerg([[1], []], [[Rational(1, 3)], []], exp_polar(k*I*pi)*2)
        with mpmath.workdps(80):
            want = mpmath.gamma(mpmath.mpf(1)/3) * mpmath.cbrt(2) * mpmath.expjpi(mpmath.mpf(phase.p)/phase.q if phase else 0)
        for dps in (15, 30, 50):
            with mpmath.workdps(dps + 20):
                got = mpmath.mpc(*_parts(N(expr, dps)))
                assert abs(got - want) <= mpmath.mpf(10) ** (1 - dps) * abs(want), (k, dps, got)


# The retry at higher argument precision must not swallow an interrupt:
# a non-Exception error, or an error raised by a signal handler during the
# first mpmath attempt, propagates at once and no retry is made.
import os  # noqa: E402
import signal  # noqa: E402

import symox._mpmath_special as _host  # noqa: E402


class _Interrupt(Exception):
    pass


def _pole_then(action, monkeypatch):
    calls = []

    def fake(*args):
        calls.append(args)
        if len(calls) == 1:
            action()
        raise ValueError("zeta(1) pole")

    monkeypatch.setitem(_host._GENERIC, "zeta", fake)
    return calls


def _raise_keyboard_interrupt():
    raise KeyboardInterrupt


@pytest.mark.parametrize("error", [KeyboardInterrupt, SystemExit, GeneratorExit])
def test_base_exception_during_first_attempt_propagates(error, monkeypatch):
    def action():
        raise error

    calls = _pole_then(action, monkeypatch)
    with pytest.raises(error):
        N(zeta(1 + Rational(1, 10**30)))
    assert len(calls) == 1
    # Later evaluations are not affected.
    monkeypatch.undo()
    assert abs(N(zeta(1 + Rational(1, 10**30))) - Float(10)**30) < 10**16


def _with_handler(handler, body):
    previous = signal.signal(signal.SIGUSR1, handler)
    try:
        body()
    finally:
        signal.signal(signal.SIGUSR1, previous)


def test_signal_handler_error_during_first_attempt_propagates(monkeypatch):
    def handler(signum, frame):
        raise _Interrupt("alarm")

    def body():
        calls = _pole_then(lambda: os.kill(os.getpid(), signal.SIGUSR1), monkeypatch)
        with pytest.raises(_Interrupt):
            N(zeta(1 + Rational(1, 10**30)))
        assert len(calls) == 1

    _with_handler(handler, body)


def test_ordinary_error_is_still_retried(monkeypatch):
    # An ordinary mpmath error at the rounded point is retried with more
    # argument bits, and the value is found.
    calls = []
    real = mpmath.zeta

    def fake(*args):
        calls.append(args)
        if len(calls) == 1:
            raise ValueError("zeta(1) pole")
        return real(*args)

    monkeypatch.setitem(_host._GENERIC, "zeta", fake)
    value = N(zeta(1 + Rational(1, 10**30)))
    assert len(calls) >= 2
    assert abs(value - Float(10)**30) < 10**16
