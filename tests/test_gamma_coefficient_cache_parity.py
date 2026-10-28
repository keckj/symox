"""Spouge coefficient reuse preserves precision, signs and host fallbacks."""
import mpmath as mp
import pytest

import symox as sx


def assert_accurate(value, reference, digits):
    with mp.workdps(digits + 70):
        real, imag = value.as_real_imag()
        actual = mp.mpc(mp.mpf(str(real)), mp.mpf(str(imag)))
        expected = reference()
        assert abs(actual - expected) <= max(abs(expected), mp.mpf(1)) * mp.mpf(10) ** (1 - digits)


REAL_CASES = [
    ('1/3', lambda: mp.mpf(1)/3),
    ('2/3', lambda: mp.mpf(2)/3),
    ('7/3', lambda: mp.mpf(7)/3),
    ('3001/3', lambda: mp.mpf(3001)/3),
    ('-1/3', lambda: -mp.mpf(1)/3),
    ('-7/3', lambda: -mp.mpf(7)/3),
    ('-3 + 1/10**24', lambda: -3 + mp.mpf(10)**-24),
    ('-3 - 1/10**24', lambda: -3 - mp.mpf(10)**-24),
    ('1/10**25', lambda: mp.mpf(10)**-25),
]


@pytest.mark.parametrize('source,argument', REAL_CASES)
@pytest.mark.parametrize('digits', [30, 90])
def test_real_gamma_matches_independent_high_precision(source, argument, digits):
    x = sx.sympify(source)
    value = sx.gamma(x, evaluate=False).evalf(digits)
    assert_accurate(value, lambda: mp.gamma(argument()), digits)


def test_repeated_arguments_and_precision_changes_keep_all_digits():
    # Separate expression allocations prevent the completed-value memo from
    # answering the calls intended to exercise coefficient reuse.
    for digits in [35, 80, 36, 35, 100, 80]:
        value = sx.gamma(sx.Rational(11, 7), evaluate=False).evalf(digits)
        assert_accurate(value, lambda: mp.gamma(mp.mpf(11)/7), digits)


@pytest.mark.parametrize('options', [{}, {'strict': True}, {'maxn': 150}])
def test_float_argument_preserves_precision_and_evalf_options(options):
    value = sx.gamma(sx.Float('0.3', 120), evaluate=False).evalf(75, **options)
    assert_accurate(value, lambda: mp.gamma(mp.mpf(3)/10), 75)


@pytest.mark.parametrize('sign', [-1, 1])
def test_complex_gamma_keeps_host_path_and_pole_sides(sign):
    sx.gamma(sx.Rational(2, 3)).evalf(50)
    x = -2 + sign*sx.I/sx.Integer(10)**20
    value = sx.gamma(x, evaluate=False).evalf(60)
    assert_accurate(value, lambda: mp.gamma(-2 + sign*mp.j*mp.mpf(10)**-20), 60)


def test_real_gamma_does_not_require_the_numeric_host(monkeypatch):
    import symox._mpmath_special as host

    def unavailable(*args):
        raise AssertionError('native real Gamma must not call the numeric host')

    monkeypatch.setitem(host._GENERIC, 'gamma', unavailable)
    for numerator in [5, 8, 11, 5]:
        value = sx.gamma(sx.Rational(numerator, 7), evaluate=False).evalf(45)
        assert_accurate(value, lambda: mp.gamma(mp.mpf(numerator)/7), 45)


@pytest.mark.parametrize('error', [KeyboardInterrupt, SystemExit])
def test_warm_real_gamma_does_not_suppress_complex_host_interrupt(monkeypatch, error):
    import symox._mpmath_special as host
    sx.gamma(sx.Rational(1, 3)).evalf(40)

    def interrupted(*args):
        raise error

    monkeypatch.setitem(host._GENERIC, 'gamma', interrupted)
    with pytest.raises(error):
        sx.gamma(2 + sx.I, evaluate=False).evalf(43)
